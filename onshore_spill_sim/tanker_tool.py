# -*- coding: utf-8 -*-
"""Movable tanker: click anywhere on the map, set the load, and run the spill model.

The tool reads 00_Config/project_config.json beside the project file to find the DEM and
the ground, river and receptor layers, so the Processing form opens already filled in.
"""
import json
import os
import re

from qgis.core import (QgsProject, QgsPointXY, QgsFeature, QgsGeometry, QgsCoordinateTransform,
                       QgsVectorLayer, QgsField)
from qgis.gui import QgsMapToolEmitPoint, QgsVertexMarker
from qgis.PyQt.QtCore import Qt, QVariant
from qgis.PyQt.QtGui import QColor
from qgis.PyQt.QtWidgets import (QDialog, QFormLayout, QLineEdit, QComboBox, QDoubleSpinBox,
                                 QCheckBox, QDialogButtonBox, QLabel, QVBoxLayout)

from .engine.oils import OILS
from .engine import climate

TANKER_LAYER = "Movable tanker"


def project_config():
    home = QgsProject.instance().homePath()
    p = os.path.join(home, "00_Config", "project_config.json")
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            cfg = json.load(f)
        cfg["_home"] = home
        return cfg
    return {"_home": home, "layers": {}}


def _layer_id(name):
    if not name:
        return None
    ls = QgsProject.instance().mapLayersByName(name)
    return ls[0].id() if ls else None


class TankerDialog(QDialog):
    def __init__(self, parent, pt, crs_authid, defaults):
        super().__init__(parent)
        self.setWindowTitle("Movable tanker - spill at the clicked point")
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel("Point: E %.1f  N %.1f  (%s)" % (pt.x(), pt.y(), crs_authid)))
        form = QFormLayout()
        d = defaults
        self.name = QLineEdit(d.get("name", "Tanker spill"))
        self.product = QComboBox()
        self.keys = list(OILS.keys())
        for k in self.keys:
            self.product.addItem(OILS[k]["label"], k)
        self.product.setCurrentIndex(self.keys.index(d.get("product", "DIESEL")))

        def spin(v, lo, hi, step=1.0, dec=1, suffix=""):
            s = QDoubleSpinBox()
            s.setRange(lo, hi)
            s.setDecimals(dec)
            s.setSingleStep(step)
            s.setValue(v)
            if suffix:
                s.setSuffix(suffix)
            return s
        self.volume = spin(d.get("volume_m3", 36.0), 0.1, 5000, 1, 1, " m3")
        self.release = spin(d.get("release_min", 20.0), 0.5, 1440, 5, 1, " min")
        self.month = QComboBox()
        self.month.addItems(climate.MONTHS)
        self.month.setCurrentIndex(int(d.get("month", 7)))
        self.start = QLineEdit(d.get("start", "2027-08-20T15:00"))
        self.storms = QLineEdit(d.get("storms", "0.5, 2, 40"))
        self.storms.setPlaceholderText("start h, duration h, depth mm; ...  (blank = dry)")
        self.snap = spin(d.get("snap_m", 0.0), 0, 500, 10, 0, " m")
        self.hours = spin(d.get("sim_hours", 24.0), 1, 240, 6, 0, " h")
        self.frame = spin(d.get("frame_min", 15.0), 1, 240, 5, 0, " min")
        self.heat = QCheckBox("Heat-transfer mode (temperature-dependent rheology)")
        self.heat.setChecked(bool(d.get("heat", False)))
        self.keep = QCheckBox("Add the point to the '%s' layer" % TANKER_LAYER)
        self.keep.setChecked(True)
        self.direct = QCheckBox("Run straight away (otherwise the filled-in form opens for checking)")
        self.direct.setChecked(False)
        for lab, w in (("Scenario name", self.name), ("Product", self.product), ("Volume lost", self.volume),
                       ("Time to empty", self.release), ("Month", self.month), ("Start (local time)", self.start),
                       ("Storms", self.storms), ("Snap onto stream within", self.snap),
                       ("Simulated period", self.hours), ("Frame interval", self.frame)):
            form.addRow(lab, w)
        lay.addLayout(form)
        lay.addWidget(self.heat)
        lay.addWidget(self.keep)
        lay.addWidget(self.direct)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setText("Run simulation")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)

    def values(self):
        return dict(name=self.name.text().strip() or "Tanker spill",
                    product=self.product.currentData(), volume_m3=self.volume.value(),
                    release_min=self.release.value(), month=self.month.currentIndex(),
                    start=self.start.text().strip(), storms=self.storms.text().strip(),
                    snap_m=self.snap.value(), sim_hours=self.hours.value(), frame_min=self.frame.value(),
                    heat=self.heat.isChecked(), keep=self.keep.isChecked(), direct=self.direct.isChecked())


class TankerTool(QgsMapToolEmitPoint):
    def __init__(self, iface):
        self.iface = iface
        super().__init__(iface.mapCanvas())
        self.marker = None
        self.setCursor(Qt.CrossCursor)

    def canvasReleaseEvent(self, e):
        canvas = self.iface.mapCanvas()
        pt = self.toMapCoordinates(e.pos())
        crs = canvas.mapSettings().destinationCrs()
        if self.marker is None:
            self.marker = QgsVertexMarker(canvas)
            self.marker.setIconType(QgsVertexMarker.ICON_BOX)
            self.marker.setColor(QColor(255, 140, 0))
            self.marker.setIconSize(14)
            self.marker.setPenWidth(3)
        self.marker.setCenter(pt)
        cfg = project_config()
        dlg = TankerDialog(self.iface.mainWindow(), pt, crs.authid(), cfg.get("tanker_defaults", {}))
        if dlg.exec_() != QDialog.Accepted:
            return
        v = dlg.values()
        if v["keep"]:
            self._store(pt, crs, v)
        self._run(pt, crs, v, cfg)

    def deactivate(self):
        if self.marker is not None:
            self.iface.mapCanvas().scene().removeItem(self.marker)
            self.marker = None
        super().deactivate()

    # ------------------------------------------------------------------
    def _store(self, pt, crs, v):
        ls = QgsProject.instance().mapLayersByName(TANKER_LAYER)
        lyr = ls[0] if ls else None
        if lyr is None:
            lyr = QgsVectorLayer("Point?crs=%s" % crs.authid(), TANKER_LAYER, "memory")
            QgsProject.instance().addMapLayer(lyr)
        names = [f.name() for f in lyr.fields()]
        need = [("name", QVariant.String), ("volume_m3", QVariant.Double), ("product", QVariant.String),
                ("release_min", QVariant.Double), ("month", QVariant.Int), ("start", QVariant.String),
                ("storms", QVariant.String), ("snap_m", QVariant.Double)]
        add = [QgsField(n, t) for n, t in need if n not in names]
        if add:
            lyr.dataProvider().addAttributes(add)
            lyr.updateFields()
        p = QgsPointXY(pt)
        if lyr.crs() != crs:
            p = QgsCoordinateTransform(crs, lyr.crs(), QgsProject.instance()).transform(p)
        f = QgsFeature(lyr.fields())
        f.setGeometry(QgsGeometry.fromPointXY(p))
        for k in ("name", "volume_m3", "product", "release_min", "start", "storms", "snap_m"):
            f[k] = v[k]
        f["month"] = v["month"] + 1
        lyr.dataProvider().addFeatures([f])
        lyr.triggerRepaint()

    def _run(self, pt, crs, v, cfg):
        import processing
        L = cfg.get("layers", {})
        home = cfg.get("_home", "")
        clean = re.sub(r"[^A-Za-z0-9_\-]+", "_", v["name"]).strip("_") or "Tanker"
        out = os.path.join(home, cfg.get("tanker_output", os.path.join("04_Outputs", "Tanker_Runs")), clean)
        ids = lambda key: [i for i in (_layer_id(n) for n in L.get(key, [])) if i]
        params = {
            "NAME": v["name"], "SOURCE": "%f,%f [%s]" % (pt.x(), pt.y(), crs.authid()),
            "SOURCE_LABEL": "Tanker", "KIND": 1, "PRODUCT": list(OILS.keys()).index(v["product"]),
            "VOLUME": v["volume_m3"], "RELEASE_MIN": v["release_min"], "CONTAINED": 0.0,
            "MONTH": v["month"], "START": v["start"], "STORMS": v["storms"], "SIM_HOURS": v["sim_hours"],
            "FRAME_MIN": v["frame_min"], "SNAP": v["snap_m"], "HEAT": v["heat"],
            "ROADS": ids("ROADS"), "WATER": ids("WATER"), "WOODLAND": ids("WOODLAND"),
            "RIVERS": ids("RIVERS"), "RECEPTORS": ids("RECEPTORS"),
            "ANIMATE": True, "HTML": True, "LOAD": True,
            "PREPARED_BY": cfg.get("prepared_by", ""), "OUT": out,
            "HALF_WIDTH": cfg.get("tanker_half_width_m", 12000.0)}
        dem = _layer_id(L.get("DEM"))
        if dem:
            params["DEM"] = dem
        for key in ("PADS", "SITEPLAN"):
            i = _layer_id(L.get(key))
            if i:
                params[key] = i
        if cfg.get("wind_csv"):
            params["WIND_CSV"] = os.path.join(home, cfg["wind_csv"])
        if not v["direct"]:
            processing.execAlgorithmDialog("onshorespill:simulate", params)
            return
        from qgis.core import QgsProcessingFeedback, Qgis
        from qgis.PyQt.QtWidgets import QProgressDialog, QApplication
        pd = QProgressDialog("Simulating the tanker spill...", "Cancel", 0, 100, self.iface.mainWindow())
        pd.setWindowModality(Qt.WindowModal)
        pd.setMinimumDuration(0)

        class FB(QgsProcessingFeedback):
            def __init__(self):
                super().__init__()
                self.progressChanged.connect(lambda p: (pd.setValue(int(p)), QApplication.processEvents()))

            def pushInfo(self, msg):
                pd.setLabelText(msg[:120])
                QApplication.processEvents()
                if pd.wasCanceled():
                    self.cancel()
        try:
            r = processing.run("onshorespill:simulate", params, feedback=FB())
            self.iface.messageBar().pushMessage("Tanker spill", "Finished - results in %s" % r.get("OUT", out),
                                                Qgis.Success, 10)
            if r.get("HTML"):
                from qgis.PyQt.QtGui import QDesktopServices
                from qgis.PyQt.QtCore import QUrl
                QDesktopServices.openUrl(QUrl.fromLocalFile(r["HTML"]))
        except Exception as ex:
            self.iface.messageBar().pushMessage("Tanker spill", str(ex), Qgis.Critical, 15)
        finally:
            pd.close()
