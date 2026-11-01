# -*- coding: utf-8 -*-
"""Processing algorithm: run the spill model for every point of a 'movable tanker' layer.

Each feature is one scenario.  Optional fields override the dialog values:
  name, volume_m3, product (key, e.g. DIESEL), release_min, month (1-12),
  start (YYYY-MM-DDTHH:MM), storms ('start h, dur h, mm; ...'), snap_m.
"""
import os
from qgis.core import (QgsProcessingParameterFeatureSource, QgsProcessingParameterBoolean,
                       QgsProcessing, QgsProcessingException, QgsCoordinateTransform,
                       QgsProject, QgsPointXY)

from .alg_simulate import SimulateSpillAlgorithm, _clean, PRODUCTS
from .engine.oils import OILS


class TankerLayerAlgorithm(SimulateSpillAlgorithm):

    def name(self):
        return "tankerlayer"

    def displayName(self):
        return "Simulate spills at movable tanker points"

    def createInstance(self):
        return TankerLayerAlgorithm()

    def shortHelpString(self):
        return ("Runs one tank / tanker release for each point of the chosen point layer (for example "
                "the 'Movable tanker' layer of the project: move or add points with the editing tools, "
                "then run). Field values override the form for that point: name, volume_m3, product "
                "(DIESEL, LIGHT_CRUDE, ...), release_min, month (1-12), start, storms, snap_m. "
                "Each run is written to its own sub-folder of the output folder.")

    def initAlgorithm(self, config=None):
        super().initAlgorithm(config)
        self.removeParameter("SOURCE")
        self.addParameter(QgsProcessingParameterFeatureSource(
            "TANKERS", "Movable tanker points", [QgsProcessing.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterBoolean("SELECTED", "Only the selected points", False))
        # tanker runs default to a tank release of diesel
        self.parameterDefinition("KIND").setDefaultValue(1)
        self.parameterDefinition("PRODUCT").setDefaultValue(PRODUCTS.index("DIESEL"))
        self.parameterDefinition("VOLUME").setDefaultValue(36.0)
        self.parameterDefinition("RELEASE_MIN").setDefaultValue(20.0)
        self.parameterDefinition("CONTAINED").setDefaultValue(0.0)
        self.parameterDefinition("SIM_HOURS").setDefaultValue(24.0)
        self.parameterDefinition("FRAME_MIN").setDefaultValue(15.0)
        self.parameterDefinition("NAME").setDefaultValue("Tanker")
        self.parameterDefinition("SOURCE_LABEL").setDefaultValue("Tanker")

    def source_point(self, parameters, context, crs):
        return self._pt

    def processAlgorithm(self, parameters, context, feedback):
        src = self.parameterAsSource(parameters, "TANKERS", context)
        if src is None:
            raise QgsProcessingException("Choose the tanker point layer.")
        dem = self.parameterAsRasterLayer(parameters, "DEM", context)
        ct = QgsCoordinateTransform(src.sourceCrs(), dem.crs(), QgsProject.instance())
        feats = list(src.getFeatures())
        if self.parameterAsBool(parameters, "SELECTED", context):
            lyr = self.parameterAsVectorLayer(parameters, "TANKERS", context)
            if lyr is not None:
                ids = set(lyr.selectedFeatureIds())
                feats = [f for f in feats if f.id() in ids]
        if not feats:
            raise QgsProcessingException("No tanker points to run.")
        base_out = self.parameterAsString(parameters, "OUT", context)
        names = [f.name() for f in src.fields()]
        results = {"OUT": base_out, "RUNS": []}
        for i, f in enumerate(feats):
            if feedback.isCanceled():
                break
            g = f.geometry()
            if g is None or g.isEmpty():
                continue
            p = g.asPoint() if not g.isMultipart() else g.asMultiPoint()[0]
            self._pt = ct.transform(QgsPointXY(p))
            pr = dict(parameters)

            def fv(key):
                if key in names:
                    v = f[key]
                    if v is not None and str(v) not in ("", "NULL"):
                        return v
                return None
            nm = fv("name") or "%s %d" % (self.parameterAsString(parameters, "NAME", context), i + 1)
            pr["NAME"] = str(nm)
            pr["SOURCE_LABEL"] = "Tanker spill"
            if fv("volume_m3") is not None:
                pr["VOLUME"] = float(fv("volume_m3"))
            if fv("release_min") is not None:
                pr["RELEASE_MIN"] = float(fv("release_min"))
            if fv("product") is not None and str(fv("product")).upper() in OILS:
                pr["PRODUCT"] = PRODUCTS.index(str(fv("product")).upper())
            if fv("month") is not None:
                pr["MONTH"] = int(fv("month")) - 1
            if fv("start") is not None:
                pr["START"] = str(fv("start"))
            if fv("storms") is not None:
                pr["STORMS"] = str(fv("storms"))
            if fv("snap_m") is not None:
                pr["SNAP"] = float(fv("snap_m"))
            pr["OUT"] = os.path.join(base_out, _clean(str(nm)))
            feedback.pushInfo("---- Tanker point %d of %d: %s" % (i + 1, len(feats), nm))
            cfg = self.prepare(pr, context, feedback)
            r = self.execute(cfg, pr, context, feedback)
            results["RUNS"].append(r.get("RUN_JSON"))
        return results
