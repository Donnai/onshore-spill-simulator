# -*- coding: utf-8 -*-
"""
Simulation panel: set up, run and play back oil-spill and H2S scenarios inside QGIS.

  Set up       hazard, place (well, fuel store, stream crossing or any point clicked on the map),
               date and time of the release - the season, temperatures, rain, wind, stream flow,
               soil wetness and groundwater depth follow the date - storms and release figures.
  Play         steps through the results on the map canvas with a time slider, play / pause,
               and subtitles that explain what is happening in each frame.
  Settlements  towns and villages near the chosen site with their distance and direction, and how
               each is affected in the run being played; new settlements can be digitised.
  Log          messages from the model.

The panel uses the study set-up in <project folder>/00_Config/model_inputs.json, so a run made
here uses the same layers and rules as the reference scenarios.
"""
import logging
import csv
import datetime as _dt
import glob
import json
import math
import os
import queue
import re
import sys

from qgis.core import (QgsApplication, QgsProject, QgsTask, QgsRasterLayer, QgsPointXY,
                       QgsCoordinateReferenceSystem, QgsCoordinateTransform, QgsRectangle,
                       QgsVectorLayer, QgsMessageLog, Qgis)
from qgis.gui import QgsMapToolEmitPoint, QgsVertexMarker
from qgis.PyQt.QtCore import Qt, QTimer, QDateTime, QDate, QTime, QObject, QEvent, pyqtSignal, QUrl
from qgis.PyQt.QtGui import QColor, QFont, QDesktopServices
from qgis.PyQt.QtWidgets import (QDockWidget, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QGridLayout,
                                 QTabWidget, QGroupBox, QLabel, QComboBox, QPushButton, QToolButton,
                                 QDoubleSpinBox, QSpinBox, QDateTimeEdit, QLineEdit, QCheckBox, QSlider,
                                 QProgressBar, QPlainTextEdit, QScrollArea, QStackedWidget, QFrame,
                                 QTableWidget, QTableWidgetItem, QHeaderView, QFileDialog, QInputDialog,
                                 QMessageBox, QSizePolicy, QApplication, QProgressDialog, QAbstractItemView)

from .engine import climate, study, outputs
from .engine.oils import OILS

GREEN = "#0b5d3b"
GOLD = "#c9a227"
HAZARDS = [("BLOWOUT", "Oil - well oil spillage"),
           ("H2S", "H2S - gas leakage at the well"),
           ("TANK", "Oil - rig diesel tank failure"),
           ("TANKER", "Oil - road tanker spill"),
           ("CONTINUOUS", "Oil - continuous leak (flowline, hose)")]
STAB = {"A": "very unstable", "B": "unstable", "C": "slightly unstable", "D": "neutral", "E": "stable",
        "F": "very stable"}


def _clean(s):
    return re.sub(r"[^A-Za-z0-9_\-]+", "_", s).strip("_") or "run"


def _log(msg):
    QgsMessageLog.logMessage(msg, "Spill simulator", Qgis.Info)


# ---------------------------------------------------------------------------
class _Feedback(object):
    """Progress / cancel adapter between the engine and a QgsTask."""

    def __init__(self, task):
        self.task = task

    def setProgress(self, v):
        self.task.setProgress(max(0.0, min(99.0, float(v))))

    def isCanceled(self):
        return self.task.isCanceled()

    def pushInfo(self, s):
        self.task.q.put(s)


class RunTask(QgsTask):
    def __init__(self, desc, kind, cfg):
        QgsTask.__init__(self, desc, QgsTask.CanCancel)
        self.kind, self.cfg = kind, cfg
        self.q = queue.Queue()
        self.result = None
        self.error = None

    def run(self):
        fb = _Feedback(self)
        try:
            if self.kind == "H2S":
                from .engine.h2s import run_h2s
                r = run_h2s(self.cfg, log=self.q.put, feedback=fb)
                self.result = dict(files=r["files"], summary=r["summary"])
            else:
                from .engine.runner import run_scenario
                r = run_scenario(self.cfg, log=self.q.put, feedback=fb)
                self.result = dict(files=r.files, summary=r.files.get("summary"))
            return not self.isCanceled()
        except Exception as ex:
            import traceback
            self.error = "%s\n%s" % (ex, traceback.format_exc())
            return False


# ---------------------------------------------------------------------------
class CaptionOverlay(QFrame):
    """Subtitle panel drawn over the bottom of the map canvas."""

    def __init__(self, canvas):
        QFrame.__init__(self, canvas)
        self.canvas = canvas
        self.setObjectName("spillcap")
        self.setStyleSheet("#spillcap{background:rgba(16,28,22,215);border-radius:8px}"
                           "QLabel{color:white;background:transparent}")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 8, 12, 9)
        lay.setSpacing(2)
        top = QHBoxLayout()
        self.clock = QLabel()
        self.clock.setStyleSheet("color:#ffffff;font-weight:bold;font-size:15px")
        self.title = QLabel()
        self.title.setStyleSheet("color:#ffd43b;font-weight:bold;font-size:14px")
        self.wind = QLabel()
        self.wind.setStyleSheet("color:#c7d3cc;font-size:12px")
        self.title.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        top.addWidget(self.clock)
        top.addSpacing(10)
        top.addWidget(self.title, 1)
        top.addWidget(self.wind)
        lay.addLayout(top)
        self.text = QLabel()
        self.text.setWordWrap(True)
        self.text.setStyleSheet("font-size:13px")
        lay.addWidget(self.text)
        canvas.installEventFilter(self)
        self.hide()

    def eventFilter(self, obj, ev):
        if obj is self.canvas and ev.type() == QEvent.Resize:
            self.place()
        return False

    def place(self):
        w = min(820, max(320, self.canvas.width() - 40))
        self.setFixedWidth(w)
        h = self.layout().totalHeightForWidth(w) if self.layout().hasHeightForWidth() else self.sizeHint().height()
        self.setFixedHeight(max(h, 60) + 4)
        self.move(16, self.canvas.height() - self.height() - 34)

    def show_caption(self, clock, title, text, wind):
        self.clock.setText(clock or "")
        self.title.setText(title or "")
        self.wind.setText(wind or "")
        self.text.setText(text or "")
        self.show()
        self.raise_()
        self.place()


# ---------------------------------------------------------------------------
class SimulationDock(QDockWidget):
    def __init__(self, iface, parent=None):
        QDockWidget.__init__(self, "Spill and H2S Simulator", parent or iface.mainWindow())
        self.setObjectName("OnshoreSpillSimulationDock")
        self.iface = iface
        self.canvas = iface.mapCanvas()
        self.root = None
        self.M = None
        self.pick_xy = None
        self.marker = None
        self.task = None
        self.play_layers = []
        self.play_rec = None
        self.play_dir = None
        self.frame = 0
        self.timer = QTimer(self)
        self.timer.timeout.connect(self._tick)
        self.logtimer = QTimer(self)
        self.logtimer.timeout.connect(self._drain_log)
        self.overlay = CaptionOverlay(self.canvas)
        self._build()
        self.reload_project()
        self._rp = lambda *a: self.reload_project()
        QgsProject.instance().readProject.connect(self._rp)

    # ------------------------------------------------------------------ UI
    def _build(self):
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(6, 6, 6, 6)
        head = QLabel()
        head.setTextFormat(Qt.RichText)
        head.setStyleSheet("background:%s;color:white;padding:8px 10px;border-bottom:3px solid %s;"
                           "font-size:13px" % (GREEN, GOLD))
        self.head = head
        v.addWidget(head)
        self.tabs = QTabWidget()
        self.tabs.addTab(self._tab_setup(), "Set up and run")
        self.tabs.addTab(self._tab_play(), "Play")
        self.tabs.addTab(self._tab_settle(), "Settlements")
        self.logbox = QPlainTextEdit()
        self.logbox.setReadOnly(True)
        self.tabs.addTab(self.logbox, "Log")
        v.addWidget(self.tabs, 1)
        self.setWidget(w)
        self.setMinimumWidth(400)

    def _spin(self, lo, hi, val, dec=1, step=1.0, suffix=""):
        s = QDoubleSpinBox()
        s.setRange(lo, hi)
        s.setDecimals(dec)
        s.setSingleStep(step)
        s.setValue(val)
        if suffix:
            s.setSuffix(" " + suffix)
        return s

    def _tab_setup(self):
        page = QWidget()
        v = QVBoxLayout(page)
        v.setSpacing(8)
        # what
        g = QGroupBox("1  What happens")
        f = QFormLayout(g)
        self.hazard = QComboBox()
        for k, lab in HAZARDS:
            self.hazard.addItem(lab, k)
        f.addRow("Hazard", self.hazard)
        self.product = QComboBox()
        for k, o in OILS.items():
            self.product.addItem(o["label"], k)
        f.addRow("Product", self.product)
        v.addWidget(g)
        # where
        g = QGroupBox("2  Where")
        f = QFormLayout(g)
        self.place = QComboBox()
        self.place.currentIndexChanged.connect(self._place_changed)
        f.addRow("Location", self.place)
        h = QHBoxLayout()
        self.pickbtn = QPushButton("Pick a point on the map")
        self.pickbtn.setCheckable(True)
        self.pickbtn.toggled.connect(self._pick_toggled)
        self.xylab = QLabel("-")
        self.xylab.setStyleSheet("color:#5b6b63")
        h.addWidget(self.pickbtn)
        h.addWidget(self.xylab, 1)
        f.addRow(h)
        self.snap = QCheckBox("Move the point onto the nearest stream (within 90 m)")
        f.addRow(self.snap)
        v.addWidget(g)
        # when
        g = QGroupBox("3  When")
        f = QFormLayout(g)
        self.start = QDateTimeEdit(QDateTime(QDate(2027, 3, 15), QTime(6, 0)))
        self.start.setCalendarPopup(True)
        self.start.setDisplayFormat("ddd dd MMM yyyy  HH:mm")
        self.start.dateTimeChanged.connect(self._date_changed)
        f.addRow("Release starts", self.start)
        h = QHBoxLayout()
        self.simh = self._spin(1, 240, 96, 0, 6, "h")
        self.framem = self._spin(5, 180, 60, 0, 5, "min")
        h.addWidget(QLabel("Simulate"))
        h.addWidget(self.simh)
        h.addWidget(QLabel("frame every"))
        h.addWidget(self.framem)
        f.addRow(h)
        self.season = QLabel()
        self.season.setWordWrap(True)
        self.season.setTextFormat(Qt.RichText)
        self.season.setStyleSheet("background:#f3f6f4;border:1px solid #d9e0dc;border-left:4px solid %s;"
                                  "padding:7px;color:#14211b" % GOLD)
        f.addRow(self.season)
        h = QHBoxLayout()
        self.storms = QLineEdit()
        self.storms.setPlaceholderText("start h, duration h, depth mm; ...  (blank = no rain)")
        b1 = QPushButton("Typical for the month")
        b1.clicked.connect(self._typical_storms)
        h.addWidget(self.storms, 1)
        h.addWidget(b1)
        f.addRow("Storms", h)
        v.addWidget(g)
        # release
        g = QGroupBox("4  Release")
        gl = QVBoxLayout(g)
        self.stack = QStackedWidget()
        # well spillage / continuous
        p = QWidget()
        f = QFormLayout(p)
        self.rate = self._spin(1, 100000, 1500, 0, 100, "bbl/d")
        self.rdur = self._spin(0.5, 720, 72, 1, 6, "h")
        self.spray = self._spin(0, 100, 40, 0, 5, "% as spray")
        f.addRow("Flow rate", self.rate)
        f.addRow("Flows for", self.rdur)
        f.addRow("Spray", self.spray)
        self.stack.addWidget(p)
        # tank / tanker
        p = QWidget()
        f = QFormLayout(p)
        self.vol = self._spin(0.1, 10000, 60, 1, 1, "m3")
        self.rmin = self._spin(0.5, 600, 10, 1, 5, "min")
        f.addRow("Volume lost", self.vol)
        f.addRow("Escapes in", self.rmin)
        self.stack.addWidget(p)
        # H2S
        p = QWidget()
        f = QFormLayout(p)
        self.gas = self._spin(0.1, 500, 10, 1, 1, "MMscf/d")
        self.ppm = self._spin(1, 500000, 10000, 0, 1000, "ppm H2S")
        self.jet = QComboBox()
        self.jet.addItem("Impinged on the substructure / BOP (3 m, no plume rise)", False)
        self.jet.addItem("Free vertical jet (10 m, with plume rise)", True)
        self.gdur = self._spin(0.5, 720, 72, 1, 6, "h")
        f.addRow("Gas rate", self.gas)
        f.addRow("H2S content", self.ppm)
        f.addRow("Release", self.jet)
        f.addRow("Flows for", self.gdur)
        self.gbuild = self._spin(0, 48, 0, 1, 0.5, "h")
        self.gbuild.setSpecialValueText("full rate at once (worst case)")
        self.gbuild.setToolTip("Time for the leak to grow from 10 % to the full rate. 0 = full rate from the start.")
        f.addRow("Leak builds up over", self.gbuild)
        self.stack.addWidget(p)
        gl.addWidget(self.stack)
        hh = QHBoxLayout()
        self.heat = QCheckBox("Heat transfer (sun, air, wind, ground)")
        self.heat.setChecked(True)
        self.srct = self._spin(-1, 150, 60, 0, 5, "C at source")
        self.srct.setSpecialValueText("air temperature")
        hh.addWidget(self.heat)
        hh.addWidget(self.srct)
        self.heatrow = QWidget()
        self.heatrow.setLayout(hh)
        gl.addWidget(self.heatrow)
        v.addWidget(g)
        # outputs + run
        g = QGroupBox("5  Run")
        f = QFormLayout(g)
        self.name = QLineEdit()
        f.addRow("Name", self.name)
        h = QHBoxLayout()
        self.html = QCheckBox("HTML page")
        self.html.setChecked(True)
        self.anim = QCheckBox("GIF / MP4")
        self.load = QCheckBox("Add result layers")
        self.load.setChecked(True)
        for c in (self.html, self.anim, self.load):
            h.addWidget(c)
        f.addRow(h)
        h = QHBoxLayout()
        self.runbtn = QPushButton("Run simulation")
        self.runbtn.setStyleSheet("QPushButton{background:%s;color:white;font-weight:bold;padding:7px 14px;"
                                  "border-radius:4px}QPushButton:disabled{background:#9bb3a8}" % GREEN)
        self.runbtn.clicked.connect(self.run)
        self.cancelbtn = QPushButton("Cancel")
        self.cancelbtn.setEnabled(False)
        self.cancelbtn.clicked.connect(self._cancel)
        h.addWidget(self.runbtn, 1)
        h.addWidget(self.cancelbtn)
        f.addRow(h)
        self.bar = QProgressBar()
        self.bar.setRange(0, 100)
        self.bar.setValue(0)
        self.status = QLabel("")
        self.status.setWordWrap(True)
        f.addRow(self.bar)
        f.addRow(self.status)
        v.addWidget(g)
        v.addStretch(1)
        self.hazard.currentIndexChanged.connect(self._hazard_changed)
        sc = QScrollArea()
        sc.setWidgetResizable(True)
        sc.setFrameShape(QFrame.NoFrame)
        sc.setWidget(page)
        return sc

    def _tab_play(self):
        page = QWidget()
        v = QVBoxLayout(page)
        h = QHBoxLayout()
        self.runs = QComboBox()
        self.runs.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.runs.setMinimumContentsLength(20)
        self.runs.currentIndexChanged.connect(self._run_selected)
        rb = QToolButton()
        rb.setText("Refresh")
        rb.clicked.connect(self.refresh_runs)
        h.addWidget(self.runs, 1)
        h.addWidget(rb)
        v.addLayout(h)
        h = QHBoxLayout()
        h.addWidget(QLabel("Show"))
        self.field = QComboBox()
        self.field.currentIndexChanged.connect(self._field_changed)
        h.addWidget(self.field, 1)
        v.addLayout(h)
        self.clock = QLabel("T+00:00")
        self.clock.setStyleSheet("font-size:24px;font-weight:bold;color:%s" % GREEN)
        self.date = QLabel("")
        self.date.setStyleSheet("color:#5b6b63")
        v.addWidget(self.clock)
        v.addWidget(self.date)
        self.slider = QSlider(Qt.Horizontal)
        self.slider.valueChanged.connect(self.show_frame)
        v.addWidget(self.slider)
        h = QHBoxLayout()
        for txt, fn in (("|<", lambda: self.slider.setValue(0)), ("<", lambda: self.slider.setValue(self.frame - 1))):
            b = QPushButton(txt)
            b.setFixedWidth(34)
            b.clicked.connect(fn)
            h.addWidget(b)
        self.playbtn = QPushButton("Play")
        self.playbtn.setStyleSheet("QPushButton{background:%s;color:white;font-weight:bold;padding:4px 14px}" % GREEN)
        self.playbtn.clicked.connect(self.toggle_play)
        h.addWidget(self.playbtn)
        for txt, fn in ((">", lambda: self.slider.setValue(self.frame + 1)),
                        (">|", lambda: self.slider.setValue(self.slider.maximum()))):
            b = QPushButton(txt)
            b.setFixedWidth(34)
            b.clicked.connect(fn)
            h.addWidget(b)
        self.speed = QComboBox()
        for lab, ms in (("Slow", 1200), ("Normal", 600), ("Fast", 300), ("Very fast", 120)):
            self.speed.addItem(lab, ms)
        self.speed.setCurrentIndex(1)
        self.speed.currentIndexChanged.connect(lambda *_: self.timer.setInterval(self.speed.currentData()))
        h.addWidget(self.speed)
        self.loop = QCheckBox("Loop")
        h.addWidget(self.loop)
        v.addLayout(h)
        self.subs = QCheckBox("Show subtitles on the map")
        self.subs.setChecked(True)
        self.subs.toggled.connect(lambda on: self.show_frame(self.frame) if on else self.overlay.hide())
        v.addWidget(self.subs)
        cap = QFrame()
        cap.setStyleSheet("QFrame{background:#14211b;border-radius:6px}QLabel{color:white}")
        cl = QVBoxLayout(cap)
        self.captitle = QLabel()
        self.captitle.setStyleSheet("color:#ffd43b;font-weight:bold;font-size:13px")
        self.captext = QLabel()
        self.captext.setWordWrap(True)
        cl.addWidget(self.captitle)
        cl.addWidget(self.captext)
        v.addWidget(cap)
        self.stats = QLabel()
        self.stats.setTextFormat(Qt.RichText)
        self.stats.setWordWrap(True)
        v.addWidget(self.stats)
        h = QGridLayout()
        b = QPushButton("Zoom to the run")
        b.clicked.connect(self.zoom_run)
        h.addWidget(b, 0, 0)
        b = QPushButton("Open the HTML page")
        b.clicked.connect(self.open_html)
        h.addWidget(b, 0, 1)
        b = QPushButton("Make GIF / MP4 from this map")
        b.clicked.connect(self.make_video)
        h.addWidget(b, 1, 0)
        b = QPushButton("Open the run folder")
        b.clicked.connect(lambda: self.play_dir and QDesktopServices.openUrl(QUrl.fromLocalFile(self.play_dir)))
        h.addWidget(b, 1, 1)
        b = QPushButton("Stop playing and remove the animation layers")
        b.clicked.connect(self.clear_player)
        h.addWidget(b, 2, 0, 1, 2)
        v.addLayout(h)
        v.addStretch(1)
        return page

    def _tab_settle(self):
        page = QWidget()
        v = QVBoxLayout(page)
        h = QHBoxLayout()
        h.addWidget(QLabel("Around"))
        self.ssite = QComboBox()
        self.ssite.currentIndexChanged.connect(self.fill_settlements)
        h.addWidget(self.ssite, 1)
        h.addWidget(QLabel("within"))
        self.srad = QSpinBox()
        self.srad.setRange(1, 50)
        self.srad.setValue(10)
        self.srad.setSuffix(" km")
        self.srad.valueChanged.connect(self.fill_settlements)
        h.addWidget(self.srad)
        v.addLayout(h)
        self.stable = QTableWidget(0, 6)
        self.stable.setHorizontalHeaderLabels(["Settlement", "Distance", "Direction", "District",
                                               "In the run being played", "Source"])
        self.stable.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.stable.horizontalHeader().setStretchLastSection(True)
        self.stable.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.stable.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.stable.cellDoubleClicked.connect(self._zoom_settlement)
        self.stable.verticalHeader().setVisible(False)
        v.addWidget(self.stable, 1)
        self.ssum = QLabel()
        self.ssum.setWordWrap(True)
        v.addWidget(self.ssum)
        h = QHBoxLayout()
        self.addbtn = QPushButton("Add a settlement (click it on the map)")
        self.addbtn.setCheckable(True)
        self.addbtn.toggled.connect(self._add_toggled)
        h.addWidget(self.addbtn)
        self.rings = QCheckBox("Distance rings")
        self.rings.toggled.connect(self._rings)
        h.addWidget(self.rings)
        v.addLayout(h)
        note = QLabel("Double-click a row to zoom to it. Settlements come from the configured towns layers. "
                      "Points added here are saved in the "
                      "'settlements_digitised' layer and are used as receptors in every later run. Turn on the "
                      "Google Hybrid basemap to see the compounds before digitising.")
        note.setWordWrap(True)
        note.setStyleSheet("color:#5b6b63;font-size:11px")
        v.addWidget(note)
        return page

    # ------------------------------------------------------------- project
    def reload_project(self):
        fn = QgsProject.instance().fileName()
        root = study.find_root(fn) if fn else None
        self.root, self.M = root, None
        if root:
            try:
                self.M = study.load(root)
            except Exception as ex:
                self._msg("Could not read 00_Config/model_inputs.json: %s" % ex)
        if self.M:
            wc = self.M.get("wells_config", {})
            self.head.setText("<b>%s</b><br><span style='font-size:11px'>%s</span>" % (
                os.path.basename(fn), wc.get("prepared_by", self.M.get("prepared_by", ""))))
            self.runbtn.setEnabled(True)
        else:
            self.head.setText("<b>Open the study project</b><br><span style='font-size:11px'>The panel needs a "
                              "project saved in a folder with 00_Config/model_inputs.json</span>")
            self.runbtn.setEnabled(False)
        self._fill_places()
        self._hazard_changed()
        self._date_changed()
        self.refresh_runs()
        self.ssite.blockSignals(True)
        self.ssite.clear()
        for k, xy in self._wells().items():
            self.ssite.addItem(k, xy)
        self.ssite.blockSignals(False)
        self.fill_settlements()

    def _wells(self):
        if not self.M:
            return {}
        return {k: tuple(v["xy"]) for k, v in self.M.get("wells_config", {}).get("wells", {}).items()}

    def _fill_places(self):
        self.place.blockSignals(True)
        self.place.clear()
        wc = (self.M or {}).get("wells_config", {}).get("wells", {})
        for k, w in wc.items():
            self.place.addItem("%s - wellhead" % k, ("WELL", k, tuple(w["xy"])))
        for k, w in wc.items():
            if w.get("fuel_storage_xy"):
                self.place.addItem("%s - rig fuel storage" % k, ("FUEL", k, tuple(w["fuel_storage_xy"])))
        for k, w in wc.items():
            if w.get("stream_crossing_xy"):
                self.place.addItem("%s - %s" % (k, w.get("stream_crossing_desc", "stream crossing")),
                                   ("CROSS", k, tuple(w["stream_crossing_xy"])))
        self.place.addItem("Point picked on the map", ("PICK", None, None))
        self.place.blockSignals(False)
        self._place_changed()

    # --------------------------------------------------------------- set up
    def _basis(self):
        return (self.M or {}).get("wells_config", {}).get("release_basis", {})

    def _hazard_changed(self, *a):
        k = self.hazard.currentData()
        B = self._basis()
        self.product.setEnabled(k != "H2S")
        self.heatrow.setVisible(k != "H2S")
        if k in ("BLOWOUT", "CONTINUOUS"):
            self.stack.setCurrentIndex(0)
            self.spray.setEnabled(k == "BLOWOUT")
            self._set_product("LIGHT_CRUDE")
            if k == "BLOWOUT":
                self.rate.setValue(B.get("blowout_rate_bbl_d", 1500))
                self.rdur.setValue(B.get("blowout_duration_h", 72))
                self.spray.setValue(100 * B.get("spray_fraction", 0.4))
                self.srct.setValue(B.get("wellhead_temp_c", 60))
            else:
                self.rate.setValue(200)
                self.rdur.setValue(12)
                self.srct.setValue(-1)
            self._choose_place("WELL")
            self.simh.setValue(96 if k == "BLOWOUT" else 48)
            self.framem.setValue(60)
        elif k in ("TANK", "TANKER"):
            self.stack.setCurrentIndex(1)
            self._set_product("DIESEL")
            self.srct.setValue(-1)
            if k == "TANK":
                self.vol.setValue(B.get("diesel_tank_m3", 60))
                self.rmin.setValue(10)
                self._choose_place("FUEL")
                self.simh.setValue(48)
                self.framem.setValue(30)
            else:
                self.vol.setValue(B.get("road_tanker_m3", 36))
                self.rmin.setValue(20)
                self._choose_place("CROSS")
                self.simh.setValue(24)
                self.framem.setValue(15)
        else:
            self.stack.setCurrentIndex(2)
            self.gas.setValue(B.get("gas_rate_mmscfd", 10))
            self.ppm.setValue(B.get("h2s_ppm", 10000))
            self.gdur.setValue(B.get("h2s_duration_h", 72))
            self._choose_place("WELL")
            self.simh.setValue(B.get("h2s_duration_h", 72) + 2)
            self.framem.setValue(15)
        self._auto_name()

    def _set_product(self, key):
        i = self.product.findData(key)
        if i >= 0:
            self.product.setCurrentIndex(i)

    def _choose_place(self, kind):
        cur = self.place.currentData()
        if cur and cur[0] == kind:
            return
        well = cur[1] if cur and cur[1] else None
        for i in range(self.place.count()):
            d = self.place.itemData(i)
            if d[0] == kind and (well is None or d[1] == well):
                self.place.setCurrentIndex(i)
                return
        for i in range(self.place.count()):
            if self.place.itemData(i)[0] == kind:
                self.place.setCurrentIndex(i)
                return

    def _place_changed(self, *a):
        d = self.place.currentData()
        if not d:
            return
        self.snap.setChecked(d[0] == "CROSS")
        xy = d[2] if d[0] != "PICK" else self.pick_xy
        self.xylab.setText("E %.0f  N %.0f" % xy if xy else "click 'Pick a point on the map'")
        self._mark(xy)
        self._auto_name()

    def _mark(self, xy):
        if self.marker is None:
            self.marker = QgsVertexMarker(self.canvas)
            self.marker.setIconType(QgsVertexMarker.ICON_CIRCLE)
            self.marker.setColor(QColor(220, 30, 30))
            self.marker.setFillColor(QColor(255, 212, 0, 200))
            self.marker.setIconSize(14)
            self.marker.setPenWidth(3)
        if xy is None:
            self.marker.hide()
            return
        p = self._to_canvas(QgsPointXY(*xy))
        self.marker.setCenter(p)
        self.marker.show()

    def _model_crs(self):
        c = QgsCoordinateReferenceSystem("EPSG:32630")
        if self.M and self.M.get("dem"):
            l = QgsRasterLayer(self.M["dem"], "dem")
            if l.isValid():
                c = l.crs()
        return c

    def _to_canvas(self, p):
        src, dst = self._model_crs(), self.canvas.mapSettings().destinationCrs()
        if src == dst:
            return p
        return QgsCoordinateTransform(src, dst, QgsProject.instance()).transform(p)

    def _from_canvas(self, p):
        src, dst = self.canvas.mapSettings().destinationCrs(), self._model_crs()
        if src == dst:
            return p
        return QgsCoordinateTransform(src, dst, QgsProject.instance()).transform(p)

    def _pick_toggled(self, on):
        if on:
            self.tool = QgsMapToolEmitPoint(self.canvas)
            self.tool.canvasClicked.connect(self._picked)
            self.canvas.setMapTool(self.tool)
            self.iface.messageBar().pushInfo("Simulator", "Click the map where the release happens.")
        elif getattr(self, "tool", None) is not None:
            self.canvas.unsetMapTool(self.tool)

    def _picked(self, pt, button):
        p = self._from_canvas(pt)
        self.pick_xy = (p.x(), p.y())
        self.pickbtn.setChecked(False)
        for i in range(self.place.count()):
            if self.place.itemData(i)[0] == "PICK":
                self.place.setCurrentIndex(i)
        self._place_changed()

    def _date_changed(self, *a):
        d = self.start.dateTime().toPyDateTime()
        m = d.month - 1
        s = climate.season_summary(m, (self.M or {}).get("wind_csv"))
        hints = "".join("<li>%s</li>" % x for x in s["hints"])
        self.season.setText(
            "<b style='color:%s;font-size:13px'>%s - %s</b><br>%s<br>%s<br>%s<br>%s<br>%s<br>%s%s" % (
                GREEN, d.strftime("%B"), s["season"], s["temperature"], s["rain"], s["wind"], s["streams"],
                s["soil"], s["sky"], ("<ul style='margin:4px 0 0 -18px'>%s</ul>" % hints) if hints else ""))
        self._auto_name()

    def _typical_storms(self):
        d = self.start.dateTime().toPyDateTime()
        st = climate.typical_storms(d.month - 1, self.simh.value(), d.hour + d.minute / 60.0)
        self.storms.setText("; ".join("%g, %g, %g" % s for s in st))
        if not st:
            self.status.setText("No storm is expected in %.0f h in %s (Tamale normals)." % (
                self.simh.value(), d.strftime("%B")))

    def _auto_name(self):
        d = self.place.currentData()
        k = self.hazard.currentData()
        if not d:
            return
        site = d[1] or "Point"
        what = {"BLOWOUT": "well spillage", "H2S": "H2S", "TANK": "diesel tank", "TANKER": "tanker",
                "CONTINUOUS": "leak"}[k]
        dt = self.start.dateTime().toPyDateTime()
        self.name.setText("%s %s %s" % (site, what, dt.strftime("%d %b %Y %Hh")))

    # ------------------------------------------------------------------ run
    def _storms(self):
        out = []
        txt = self.storms.text().strip()
        for part in [p for p in txt.split(";") if p.strip()]:
            vals = [float(x) for x in re.split(r"[,\s]+", part.strip()) if x]
            if len(vals) != 3:
                raise ValueError("Storms must be written as 'start h, duration h, depth mm; ...'")
            out.append(tuple(vals))
        return out

    def build_cfg(self):
        d = self.place.currentData()
        xy = d[2] if d[0] != "PICK" else self.pick_xy
        if xy is None:
            raise ValueError("Pick the release point on the map first.")
        k = self.hazard.currentData()
        M = self.M
        dt = self.start.dateTime().toPyDateTime()
        name = self.name.text().strip() or "Panel run"
        prefix = _clean(name)
        out = os.path.join(self.root, "04_Outputs", "Panel_Runs", prefix)
        well = d[1] or "Point"
        W = M.get("wells_config", {}).get("wells", {}).get(well, {})
        label = {"WELL": "%s (wellhead)" % well, "FUEL": "Rig fuel storage", "CROSS": "Stream crossing",
                 "PICK": "Release point"}[d[0]]
        common = dict(name=name, prefix=prefix, title=name, well=well, source_xy=xy, source_label=label,
                      month=dt.month - 1, start_iso=dt.strftime("%Y-%m-%dT%H:%M:00"),
                      sim_hours=self.simh.value(), frame_minutes=self.framem.value(), out_dir=out,
                      prepared_by=M.get("prepared_by", ""), html=self.html.isChecked(),
                      notes=["Run from the simulation panel in QGIS on %s." % _dt.date.today().strftime("%d %b %Y")])
        if k == "H2S":
            jet = bool(self.jet.currentData())
            cfg = dict(common, dem=M["dem"], wind_csv=M.get("wind_csv"), gas_mmscfd=self.gas.value(),
                       h2s_ppm=self.ppm.value(), jet=jet, release_h_m=10.0 if jet else 3.0,
                       duration_h=self.gdur.value(), buildup_h=self.gbuild.value(), half_width_m=15000.0, cell_m=50.0,
                       receptor_specs=M["receptor_specs"], html_layers=M["html_layers"],
                       zone_m=W.get("h2s_zone_m"), zone_label=W.get("h2s_zone_label"),
                       subtitle="H2S dispersion from a gas leakage at the well")
            return "H2S", cfg
        cfg = dict(study.common_cfg(M), **common)
        cfg.update(product=self.product.currentData(), rain=self._storms(), heat_mode=self.heat.isChecked(),
                   discharge_temp_c=None if self.srct.value() < 0 else self.srct.value())
        if k in ("BLOWOUT", "CONTINUOUS"):
            cfg.update(kind=k, rate_bbl_d=self.rate.value(), duration_h=self.rdur.value(),
                       spray_fraction=self.spray.value() / 100.0 if k == "BLOWOUT" else 0.0)
        else:
            cfg.update(kind="TANK", volume_m3=self.vol.value(), release_minutes=self.rmin.value(), contained_m3=0.0)
        if self.snap.isChecked():
            cfg["snap_to_channel_m"] = 90.0
        if k == "TANKER":
            cfg["half_width_m"] = 12000.0
        return "OIL", cfg

    def run(self):
        if self.task is not None:
            return
        try:
            kind, cfg = self.build_cfg()
        except Exception as ex:
            QMessageBox.warning(self, "Simulator", str(ex))
            return
        self._last_cfg = (kind, cfg)
        self.logbox.appendPlainText("\n==== %s  (%s)" % (cfg["name"], _dt.datetime.now().strftime("%H:%M:%S")))
        self.task = RunTask("Simulating %s" % cfg["name"], kind, cfg)
        self.task.progressChanged.connect(lambda v: self.bar.setValue(int(v)))
        self.task.taskCompleted.connect(self._done)
        self.task.taskTerminated.connect(self._failed)
        self.runbtn.setEnabled(False)
        self.cancelbtn.setEnabled(True)
        self.status.setText("Running... the map stays usable while the model works.")
        self.bar.setValue(0)
        QgsApplication.taskManager().addTask(self.task)
        self.logtimer.start(300)

    def _cancel(self):
        if self.task is not None:
            self.task.cancel()

    def _drain_log(self):
        t = self.task
        if t is None:
            return
        try:
            while True:
                self.logbox.appendPlainText(t.q.get_nowait())
        except queue.Empty:
            logging.getLogger(__name__).debug("Optional operation failed in SimulationDock._drain_log; existing fallback retained.", exc_info=True)

    def _failed(self):
        self._drain_log()
        t = self.task
        self.task = None
        self.logtimer.stop()
        self.runbtn.setEnabled(True)
        self.cancelbtn.setEnabled(False)
        msg = (t.error or "Cancelled.") if t else "Stopped."
        self.status.setText("The run did not finish: %s" % msg.split("\n")[0])
        self.logbox.appendPlainText(msg)

    def _done(self):
        self._drain_log()
        t = self.task
        self.task = None
        self.logtimer.stop()
        self.runbtn.setEnabled(True)
        self.cancelbtn.setEnabled(False)
        self.bar.setValue(100)
        if t is None or t.result is None:
            self._failed()
            return
        kind, cfg = self._last_cfg
        files = t.result["files"]
        files["summary"] = t.result.get("summary") or files.get("summary")
        try:
            if self.load.isChecked():
                if kind == "H2S":
                    from .qgis_style import load_h2s_results
                    g = load_h2s_results(files, cfg["name"])
                else:
                    from .qgis_style import load_results
                    g = load_results(files, cfg["name"])
                if g is not None:
                    g.setExpanded(False)
                    g.setItemVisibilityChecked(False)   # the animation player shows the run; tick to see the maxima
        except Exception as ex:
            self.logbox.appendPlainText("Results not added to the project: %s" % ex)
        self.status.setText("Finished: %s" % self._one_line(kind, files["summary"]))
        self.refresh_runs(select=files["run_json"])
        self.tabs.setCurrentIndex(1)
        if self.anim.isChecked():
            self.make_video()

    def _one_line(self, kind, s):
        if not s:
            return ""
        if kind == "H2S":
            return "AEGL-2 reaches %.0f m, IDLH %.0f m, highest 10-min %.0f ppm." % (
                s["dist_AEGL2_m"], s["dist_IDLH_m"], s["max_ppm"])
        return "%.1f m3 released, %.2f ha oiled, %.0f m3 evaporated, %.0f m3 soaked in." % (
            s["released_m3"], s["footprint_ha"], s["evaporated_m3"], s["infiltrated_m3"])

    # ----------------------------------------------------------------- play
    def refresh_runs(self, select=None):
        cur = select or self.runs.currentData()
        self.runs.blockSignals(True)
        self.runs.clear()
        if self.root:
            items = []
            for rj in glob.glob(os.path.join(self.root, "04_Outputs", "**", "*_run.json"), recursive=True):
                try:
                    with open(rj) as f:
                        r = json.load(f)
                except Exception:
                    logging.getLogger(__name__).debug("Optional operation failed in SimulationDock.refresh_runs; existing fallback retained.", exc_info=True); continue
                grp = os.path.relpath(rj, os.path.join(self.root, "04_Outputs")).split(os.sep)[0]
                items.append(("%s  |  %s" % (r.get("title") or r.get("name"), grp.replace("_", " ")), rj,
                              os.path.getmtime(rj)))
            items.sort(key=lambda x: (not x[1].count("Panel_Runs"), x[0]))
            for lab, rj, _ in items:
                self.runs.addItem(lab, rj)
        self.runs.blockSignals(False)
        i = self.runs.findData(cur) if cur else -1
        if i >= 0:
            self.runs.setCurrentIndex(i)
            self._run_selected()
        elif select is None and self.play_rec is None:
            pass

    def _run_selected(self, *a):
        rj = self.runs.currentData()
        if not rj or not os.path.exists(rj):
            return
        self.timer.stop()
        self.playbtn.setText("Play")
        with open(rj) as f:
            rec = json.load(f)
        self.play_rec, self.play_dir = rec, os.path.dirname(rj)
        self.play_prefix = os.path.basename(rj).replace("_run.json", "")
        if not rec.get("captions"):
            self._make_captions(rj, rec)
        self.field.blockSignals(True)
        self.field.clear()
        h2s = rec.get("model") == "h2s"
        if h2s:
            self.field.addItem("H2S at ground level (ppm)", "h2s")
        else:
            self.field.addItem("Oil on the ground (mm) and spray", "oil")
            if os.path.exists(self._ts("oil_temperature_timeseries_C")):
                self.field.addItem("Oil temperature (C)", "temp")
            if os.path.exists(self._ts("yield_stress_timeseries_Pa")):
                self.field.addItem("Yield stress - gelled oil (Pa)", "tau")
        self.field.blockSignals(False)
        self._field_changed()
        self.fill_settlements()

    def _make_captions(self, rj, rec):
        """Older runs: write the captions now from the stored results."""
        try:
            from .engine.runner import add_captions
            recs = self._run_receptors()
            rec["captions"] = add_captions(rj, os.path.join(self.play_dir, rec["timeseries"]), recs)
        except Exception as ex:
            self.logbox.appendPlainText("Captions could not be made for this run: %s" % ex)

    def _run_receptors(self):
        """Receptor results of the run being played (from its CSV)."""
        h2s = self.play_rec.get("model") == "h2s"
        p = os.path.join(self.play_dir, self.play_prefix + ("_h2s_receptors.csv" if h2s else "_receptors.csv"))
        out = []
        if not os.path.exists(p):
            return out
        with open(p, newline="") as f:
            for r in csv.DictReader(f):
                d = dict(name=r["name"], type=r["type"], x=float(r["x"]), y=float(r["y"]))
                if h2s:
                    d.update(hazard=r.get("hazard"), peak_ppm=float(r.get("peak_ppm") or 0))
                else:
                    d.update(arrival_h=float(r["arrival_h"]) if r.get("arrival_h") not in (None, "") else None,
                             status=r.get("status"))
                out.append(d)
        return out

    def _ts(self, key):
        return os.path.join(self.play_dir, "%s_%s.tif" % (self.play_prefix, key))

    def _remove_play_layers(self):
        proj = QgsProject.instance()
        for lid in self.play_layers:
            if proj.mapLayer(lid) is not None:
                proj.removeMapLayer(lid)
        self.play_layers = []
        g = proj.layerTreeRoot().findGroup("Animation player")
        if g is not None and not g.children():
            proj.layerTreeRoot().removeChildNode(g)

    def _add_play_layer(self, path, title, ramp, fmt, linear=False):
        from .qgis_style import ramp_renderer, lin_renderer
        l = QgsRasterLayer(path, title)
        if not l.isValid():
            return None
        try:
            from qgis.core import QgsBilinearRasterResampler
            l.resampleFilter().setZoomedInResampler(QgsBilinearRasterResampler())
        except Exception:
            logging.getLogger(__name__).debug("Optional operation failed in SimulationDock._add_play_layer; existing fallback retained.", exc_info=True)
        if linear:
            lin_renderer(l, ramp, fmt)
        else:
            ramp_renderer(l, ramp, fmt)
        proj = QgsProject.instance()
        root = proj.layerTreeRoot()
        g = root.findGroup("Animation player")
        if g is None:
            g = root.insertGroup(0, "Animation player")
        proj.addMapLayer(l, False)
        g.addLayer(l)
        self.play_layers.append(l.id())
        return l

    def _field_changed(self, *a):
        if not self.play_rec:
            return
        self._remove_play_layers()
        rec = self.play_rec
        key = self.field.currentData()
        title = rec.get("title") or rec.get("name")
        if key == "h2s":
            from .engine.h2s import H2S_RAMP
            self._add_play_layer(os.path.join(self.play_dir, rec["timeseries"]), "H2S (ppm) - " + title,
                                 [tuple(r) for r in rec.get("ramp", H2S_RAMP)], "%g ppm")
        elif key == "temp":
            self._add_play_layer(self._ts("oil_temperature_timeseries_C"), "Oil temperature (C) - " + title,
                                 outputs.TEMP_RAMP, "%g C", linear=True)
        elif key == "tau":
            self._add_play_layer(self._ts("yield_stress_timeseries_Pa"), "Yield stress (Pa) - " + title,
                                 outputs.YIELD_RAMP, "%g Pa")
        else:
            self._add_play_layer(os.path.join(self.play_dir, rec["timeseries"]), "Oil on the ground (mm) - " + title,
                                 outputs.OIL_RAMP, "%g mm")
            if rec.get("fallout_timeseries"):
                fl = self._add_play_layer(os.path.join(self.play_dir, rec["fallout_timeseries"]),
                                          "Spray fallout (mm) - " + title, outputs.FALL_RAMP, "%g mm")
                if fl is not None:
                    # the spray film is thin and wide: keep it see-through so the ground and the oil show
                    fl.renderer().setOpacity(0.45)
        n = max(len(rec["frames"]), 1)
        self.slider.blockSignals(True)
        self.slider.setRange(0, n - 1)
        self.slider.blockSignals(False)
        self.show_frame(min(self.frame, n - 1))

    def show_frame(self, i):
        rec = self.play_rec
        if not rec:
            return
        n = len(rec["frames"])
        i = max(0, min(int(i), n - 1))
        self.frame = i
        if self.slider.value() != i:
            self.slider.blockSignals(True)
            self.slider.setValue(i)
            self.slider.blockSignals(False)
        proj = QgsProject.instance()
        for lid in self.play_layers:
            l = proj.mapLayer(lid)
            if l is None:
                continue
            b = min(i + 1, l.bandCount())
            r = l.renderer()
            if hasattr(r, "setBand"):
                r.setBand(b)
            elif hasattr(r, "setInputBand"):
                r.setInputBand(b)
            l.triggerRepaint()
        fr = rec["frames"][i]
        t = fr["t_h"]
        hh = int(t + 1e-6)
        mm = int(round((t - hh) * 60)) % 60
        clock = "T+%02d:%02d" % (hh, mm)
        self.clock.setText(clock)
        caps = rec.get("captions") or []
        cap = caps[i] if i < len(caps) else {}
        self.date.setText(cap.get("clock") or "")
        self.captitle.setText(cap.get("title", ""))
        self.captext.setText(cap.get("text", ""))
        w = fr.get("wind")
        wtxt = ("Wind from %s %.1f m/s" % (climate.compass((w[0] + 180) % 360), w[1])) if w else ""
        if rec.get("model") == "h2s":
            s = rec["stats"][min(i, len(rec["stats"]) - 1)]
            self.stats.setText("<b>%.1f ppm</b> highest at ground &nbsp; AEGL-2 reach <b>%.0f m</b> &nbsp; IDLH reach "
                               "<b>%.0f m</b><br>Stability %s (%s), mixing height %.0f m" % (
                                   s["max_ppm"], s.get("AEGL2", 0), s.get("IDLH", 0), s.get("stability", "-"),
                                   STAB.get(s.get("stability"), ""), s.get("mixing_m", 0)))
        else:
            r = rec["mass_balance"][min(i, len(rec["mass_balance"]) - 1)]
            extra = ""
            if r.get("oil_t_mean") is not None:
                extra = "<br>Oil %.0f C, air %.0f C%s" % (r["oil_t_mean"], r.get("t_air", 0),
                                                         (", %.0f m3 gelled" % r["gelled_m3"]) if r.get("gelled_m3", 0) >= 0.5 else "")
            self.stats.setText("Released <b>%.1f m3</b> &nbsp; oiled <b>%.2f ha</b> &nbsp; evaporated <b>%.1f m3</b> "
                               "&nbsp; soaked in <b>%.1f m3</b>%s" % (r["released"], r["area_ha"], r["evaporated"],
                                                                     r["infiltrated"], extra))
        if self.subs.isChecked():
            self.overlay.show_caption("%s   %s" % (clock, cap.get("clock") or ""), cap.get("title", ""),
                                      cap.get("text", ""), wtxt)

    def toggle_play(self):
        if self.timer.isActive():
            self.timer.stop()
            self.playbtn.setText("Play")
        elif self.play_rec:
            if self.frame >= self.slider.maximum():
                self.slider.setValue(0)
            self.timer.start(self.speed.currentData())
            self.playbtn.setText("Pause")

    def _tick(self):
        if self.frame >= self.slider.maximum():
            if self.loop.isChecked():
                self.slider.setValue(0)
            else:
                self.toggle_play()
            return
        self.slider.setValue(self.frame + 1)

    def clear_player(self):
        self.timer.stop()
        self.playbtn.setText("Play")
        self._remove_play_layers()
        self.overlay.hide()
        self.canvas.refresh()

    def zoom_run(self):
        if not self.play_layers:
            return
        l = QgsProject.instance().mapLayer(self.play_layers[0])
        rec = self.play_rec
        if l is None:
            return
        # the oiled / gassed cells over the whole run, plus a margin
        from osgeo import gdal
        import numpy as np
        ds = gdal.Open(l.source())
        gt = ds.GetGeoTransform()
        thr = 0.51 if rec.get("model") == "h2s" else 0.05
        acc = None
        for b in range(1, ds.RasterCount + 1, max(ds.RasterCount // 24, 1)):
            a = ds.GetRasterBand(b).ReadAsArray() >= thr
            acc = a if acc is None else (acc | a)
        sx, sy = rec["source_xy"]
        if acc is not None and acc.any():
            rr, cc = np.nonzero(acc)
            x0, x1 = gt[0] + cc.min() * gt[1], gt[0] + (cc.max() + 1) * gt[1]
            y1, y0 = gt[3] + rr.min() * gt[5], gt[3] + (rr.max() + 1) * gt[5]
        else:
            x0, x1, y0, y1 = sx - 300, sx + 300, sy - 300, sy + 300
        pad = max(x1 - x0, y1 - y0) * 0.25 + 300
        rect = QgsRectangle(x0 - pad, y0 - pad, x1 + pad, y1 + pad)
        src = l.crs()
        dst = self.canvas.mapSettings().destinationCrs()
        if src != dst:
            rect = QgsCoordinateTransform(src, dst, QgsProject.instance()).transformBoundingBox(rect)
        self.canvas.setExtent(rect)
        self.canvas.refresh()

    def open_html(self):
        if not self.play_dir:
            return
        h = glob.glob(os.path.join(self.play_dir, "*_animation.html"))
        if h:
            QDesktopServices.openUrl(QUrl.fromLocalFile(h[0]))
        else:
            QMessageBox.information(self, "Simulator", "This run has no HTML page (tick 'HTML page' before running).")

    def make_video(self):
        if not self.play_rec:
            return
        from .qgis_anim import render_animation
        root = QgsProject.instance().layerTreeRoot()
        skip = set(self.play_layers)
        res = root.findGroup("Results")
        if res is not None:
            skip |= set(n.layerId() for n in res.findLayers())
        vis = [n.layer() for n in root.findLayers()
               if n.isVisible() and n.layer() is not None and n.layerId() not in skip
               and not n.layer().name().startswith(("Spill - ", "H2S - "))]
        rj = os.path.join(self.play_dir, self.play_prefix + "_run.json")
        dlg = QProgressDialog("Drawing the animation frames from the visible layers...", "Cancel", 0, 100, self)
        dlg.setWindowModality(Qt.WindowModal)
        dlg.setMinimumDuration(0)

        class FB(object):
            def setProgress(self_, v):
                dlg.setValue(int(v))
                QApplication.processEvents()

            def isCanceled(self_):
                return dlg.wasCanceled()

        cfg = (self.M or {}).get("wells_config", {})
        try:
            out = render_animation(rj, vis, width=1600, height=1000,
                                   fps=6 if self.play_rec.get("model") == "h2s" else 5,
                                   prepared_by=cfg.get("prepared_by", ""), log=self.logbox.appendPlainText,
                                   feedback=FB(), extent=None)
        finally:
            dlg.close()
        msg = ", ".join(os.path.basename(v) for k, v in out.items() if k in ("gif", "mp4"))
        self.iface.messageBar().pushSuccess("Simulator", "Animation written: %s" % (msg or "frames only"))

    # ---------------------------------------------------------- settlements
    def _settlement_rows(self):
        rows = []
        if not self.M:
            return rows
        from osgeo import ogr
        ds = ogr.Open(self.M["gpkg"])
        if ds is None:
            return rows
        for lname, src in (("towns", None), ("settlements_digitised", "Digitised")):
            L = ds.GetLayerByName(lname)
            if L is None:
                continue
            fields = [L.GetLayerDefn().GetFieldDefn(k).GetName() for k in range(L.GetLayerDefn().GetFieldCount())]
            for f in L:
                g = f.GetGeometryRef()
                if g is None or ("kind" in fields and (f["kind"] or "").startswith("Dugout")):
                    continue
                if g.GetGeometryCount():
                    g = g.GetGeometryRef(0)
                rows.append(dict(name=f["name"] or "(no name)", x=g.GetX(), y=g.GetY(),
                                 district=f["district"] if "district" in fields else "",
                                 source=src or (f["source"] if "source" in fields else "")))
        return rows

    def fill_settlements(self, *a):
        xy = self.ssite.currentData()
        if not xy:
            self.stable.setRowCount(0)
            return
        rad = self.srad.value() * 1000.0
        res = {}
        if self.play_rec:
            for r in self._run_receptors():
                res[r["name"].lower()] = r
        rows = []
        for s in self._settlement_rows():
            d = math.hypot(s["x"] - xy[0], s["y"] - xy[1])
            if d <= rad:
                rows.append((d, s))
        rows.sort(key=lambda x: x[0])
        self.stable.setRowCount(len(rows))
        n1 = n2 = n5 = 0
        h2s = self.play_rec and self.play_rec.get("model") == "h2s"
        for i, (d, s) in enumerate(rows):
            brg = math.degrees(math.atan2(s["x"] - xy[0], s["y"] - xy[1])) % 360
            n1 += d < 1000
            n2 += d < 2000
            n5 += d < 5000
            r = res.get(s["name"].lower())
            if r is None:
                eff = "outside the run area" if self.play_rec else ""
            elif h2s:
                eff = "%s (peak %.1f ppm)" % (r.get("hazard"), r.get("peak_ppm", 0))
            else:
                eff = ("oil arrives T+%.1f h" % r["arrival_h"]) if r.get("arrival_h") is not None else (r.get("status") or "")
            vals = [s["name"], "%.2f km" % (d / 1000.0), "%s (%.0f deg)" % (climate.compass(brg), brg),
                    s.get("district") or "", eff, s.get("source") or ""]
            for c, v in enumerate(vals):
                it = QTableWidgetItem(v)
                if c == 0:
                    it.setData(Qt.UserRole, (s["x"], s["y"]))
                if d < 1000:
                    it.setBackground(QColor(255, 230, 230))
                elif d < 2000:
                    it.setBackground(QColor(255, 244, 224))
                if c == 4 and eff and eff not in ("-", "Not reached"):
                    it.setForeground(QColor(180, 20, 30))
                self.stable.setItem(i, c, it)
        self.ssum.setText("Settlements around %s - within 1 km: <b>%d</b>; 2 km: <b>%d</b>; 5 km: <b>%d</b>; "
                          "%d km: <b>%d</b>." % (self.ssite.currentText(), n1, n2, n5, self.srad.value(), len(rows)))

    def _zoom_settlement(self, row, col):
        it = self.stable.item(row, 0)
        if it is None:
            return
        x, y = it.data(Qt.UserRole)
        p = self._to_canvas(QgsPointXY(x, y))
        r = 600.0
        if self.canvas.mapSettings().destinationCrs().isGeographic():
            r = 0.006
        self.canvas.setExtent(QgsRectangle(p.x() - r, p.y() - r, p.x() + r, p.y() + r))
        self.canvas.refresh()
        self._flash(p)

    def _flash(self, p):
        m = QgsVertexMarker(self.canvas)
        m.setCenter(p)
        m.setIconType(QgsVertexMarker.ICON_BOX)
        m.setIconSize(22)
        m.setPenWidth(3)
        m.setColor(QColor(220, 30, 30))
        QTimer.singleShot(2500, lambda: self.canvas.scene().removeItem(m))

    def _add_toggled(self, on):
        if on:
            self.atool = QgsMapToolEmitPoint(self.canvas)
            self.atool.canvasClicked.connect(self._add_settlement)
            self.canvas.setMapTool(self.atool)
            self.iface.messageBar().pushInfo("Simulator", "Click the centre of the settlement on the map.")
        elif getattr(self, "atool", None) is not None:
            self.canvas.unsetMapTool(self.atool)

    def _add_settlement(self, pt, button):
        self.addbtn.setChecked(False)
        name, ok = QInputDialog.getText(self, "New settlement", "Name of the settlement (or a description):")
        if not ok or not name.strip():
            return
        kind, ok = QInputDialog.getItem(self, "New settlement", "What is it?",
                                        ["Village", "Hamlet / compound cluster", "Farmstead", "School",
                                         "Market", "Clinic", "Water point", "Other"], 0, False)
        if not ok:
            return
        p = self._from_canvas(pt)
        from osgeo import ogr
        ds = ogr.Open(self.M["gpkg"], 1)
        L = ds.GetLayerByName("settlements_digitised")
        f = ogr.Feature(L.GetLayerDefn())
        f["name"] = name.strip()
        f["kind"] = kind
        f["source"] = "Digitised in QGIS on %s" % _dt.date.today().strftime("%d %b %Y")
        f.SetGeometry(ogr.CreateGeometryFromWkt("POINT (%f %f)" % (p.x(), p.y())))
        L.CreateFeature(f)
        ds = None
        for l in QgsProject.instance().mapLayers().values():
            if isinstance(l, QgsVectorLayer) and "settlements_digitised" in l.source():
                l.reload()
                l.triggerRepaint()
        self.fill_settlements()
        self.iface.messageBar().pushSuccess("Simulator", "'%s' saved; it is used as a receptor in the next run." % name)

    def _rings(self, on):
        for l in QgsProject.instance().mapLayers().values():
            if "well_distance_rings" in l.source():
                n = QgsProject.instance().layerTreeRoot().findLayer(l.id())
                if n is not None:
                    n.setItemVisibilityChecked(on)

    def _msg(self, s):
        self.logbox.appendPlainText(s)

    # ---------------------------------------------------------------- close
    def cleanup(self):
        self.timer.stop()
        self.overlay.hide()
        self.overlay.deleteLater()
        if self.marker is not None:
            self.canvas.scene().removeItem(self.marker)
            self.marker = None
        self._remove_play_layers()
        try:
            QgsProject.instance().readProject.disconnect(self._rp)
        except Exception:
            logging.getLogger(__name__).debug("Optional operation failed in SimulationDock.cleanup; existing fallback retained.", exc_info=True)
