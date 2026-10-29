# -*- coding: utf-8 -*-
import os
from qgis.core import QgsApplication
from qgis.PyQt.QtGui import QIcon, QDesktopServices
from qgis.PyQt.QtWidgets import QAction
from qgis.PyQt.QtCore import QUrl

from .provider import SpillProvider

MENU = "&Onshore Spill Simulator"


class OnshoreSpillPlugin:
    def __init__(self, iface):
        self.iface = iface
        self.provider = None
        self.actions = []
        self.tool = None

    def initProcessing(self):
        if self.provider is None:
            self.provider = SpillProvider()
            QgsApplication.processingRegistry().addProvider(self.provider)

    def initGui(self):
        self.initProcessing()
        here = os.path.dirname(__file__)
        icon = QIcon(os.path.join(here, "icon.png"))
        a = QAction(icon, "Simulate onshore spill...", self.iface.mainWindow())
        a.triggered.connect(lambda: self._open("onshorespill:simulate"))
        self.iface.addToolBarIcon(a)
        self.iface.addPluginToMenu(MENU, a)
        self.actions.append(a)
        b = QAction(icon, "Render animation from a finished run...", self.iface.mainWindow())
        b.triggered.connect(lambda: self._open("onshorespill:animate"))
        self.iface.addPluginToMenu(MENU, b)
        self.actions.append(b)
        t = QAction(QIcon(os.path.join(here, "icon_tanker.png")), "Place a movable tanker and simulate",
                    self.iface.mainWindow())
        t.setCheckable(True)
        t.triggered.connect(self._tanker)
        self.iface.addToolBarIcon(t)
        self.iface.addPluginToMenu(MENU, t)
        self.actions.append(t)
        self.tanker_action = t
        tl = QAction(QIcon(os.path.join(here, "icon_tanker.png")), "Simulate spills at movable tanker points...",
                     self.iface.mainWindow())
        tl.triggered.connect(lambda: self._open("onshorespill:tankerlayer", self._tanker_params()))
        self.iface.addPluginToMenu(MENU, tl)
        self.actions.append(tl)
        h = QAction(QIcon(os.path.join(here, "icon_h2s.png")), "Simulate H2S leakage at a well...",
                    self.iface.mainWindow())
        h.triggered.connect(lambda: self._open("onshorespill:h2s"))
        self.iface.addToolBarIcon(h)
        self.iface.addPluginToMenu(MENU, h)
        self.actions.append(h)
        d = QAction(QIcon(os.path.join(here, "icon_panel.png")), "Simulation panel (set up, run and play in QGIS)",
                    self.iface.mainWindow())
        d.setCheckable(True)
        d.triggered.connect(self._panel)
        self.iface.addToolBarIcon(d)
        self.iface.addPluginToMenu(MENU, d)
        self.actions.append(d)
        self.panel_action = d
        self.dock = None
        c = QAction("Method and user guide", self.iface.mainWindow())
        c.triggered.connect(lambda: QDesktopServices.openUrl(
            QUrl.fromLocalFile(os.path.join(here, "help", "index.html"))))
        self.iface.addPluginToMenu(MENU, c)
        self.actions.append(c)

    def _panel(self, checked=True):
        from qgis.PyQt.QtCore import Qt
        if self.dock is None:
            from .sim_dock import SimulationDock
            self.dock = SimulationDock(self.iface)
            self.iface.addDockWidget(Qt.RightDockWidgetArea, self.dock)
            self.dock.visibilityChanged.connect(self.panel_action.setChecked)
        self.dock.setVisible(checked)
        if checked:
            self.dock.raise_()

    def _open(self, alg_id, params=None):
        import processing
        processing.execAlgorithmDialog(alg_id, params or {})

    def _tanker(self, checked):
        from .tanker_tool import TankerTool
        canvas = self.iface.mapCanvas()
        if checked:
            if self.tool is None:
                self.tool = TankerTool(self.iface)
                self.tool.setAction(self.tanker_action)
            canvas.setMapTool(self.tool)
            self.iface.messageBar().pushInfo("Movable tanker", "Click the map where the tanker loses its load.")
        elif self.tool is not None:
            canvas.unsetMapTool(self.tool)

    def _tanker_params(self):
        from .tanker_tool import project_config, _layer_id, TANKER_LAYER
        cfg = project_config()
        L = cfg.get("layers", {})
        ids = lambda key: [i for i in (_layer_id(n) for n in L.get(key, [])) if i]
        p = {"ROADS": ids("ROADS"), "WATER": ids("WATER"), "WOODLAND": ids("WOODLAND"),
             "RIVERS": ids("RIVERS"), "RECEPTORS": ids("RECEPTORS"),
             "PREPARED_BY": cfg.get("prepared_by", "")}
        for key in ("DEM", "PADS", "SITEPLAN"):
            i = _layer_id(L.get(key))
            if i:
                p[key] = i
        t = _layer_id(TANKER_LAYER)
        if t:
            p["TANKERS"] = t
        if cfg.get("_home"):
            p["OUT"] = os.path.join(cfg["_home"], cfg.get("tanker_output", "04_Outputs/Tanker_Runs"))
        return p

    def unload(self):
        if getattr(self, "dock", None) is not None:
            try:
                self.dock.cleanup()
                self.iface.removeDockWidget(self.dock)
                self.dock.deleteLater()
            except Exception:
                pass
            self.dock = None
        if self.tool is not None:
            try:
                self.iface.mapCanvas().unsetMapTool(self.tool)
            except Exception:
                pass
            self.tool = None
        for a in self.actions:
            self.iface.removePluginMenu(MENU, a)
            self.iface.removeToolBarIcon(a)
        if self.provider is not None:
            QgsApplication.processingRegistry().removeProvider(self.provider)
            self.provider = None
