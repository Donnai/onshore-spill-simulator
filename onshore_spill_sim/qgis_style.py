# -*- coding: utf-8 -*-
"""Load and style the results of a run in the current QGIS project."""
import logging
import os

from qgis.core import (QgsProject, QgsRasterLayer, QgsVectorLayer, QgsColorRampShader,
                       QgsRasterShader, QgsSingleBandPseudoColorRenderer, QgsCategorizedSymbolRenderer,
                       QgsRendererCategory, QgsSymbol, QgsFillSymbol, QgsLineSymbol, QgsMarkerSymbol,
                       QgsRuleBasedRenderer, QgsWkbTypes)
from qgis.PyQt.QtGui import QColor

from .engine import outputs


def ramp_renderer(layer, ramp, label_fmt="%g mm"):
    items = []
    for s in ramp:
        items.append(QgsColorRampShader.ColorRampItem(s[0], QColor(s[1], s[2], s[3], max(s[4], 160)),
                                                      label_fmt % s[0]))
    sh = QgsColorRampShader(ramp[0][0], ramp[-1][0])
    sh.setColorRampType(QgsColorRampShader.Interpolated)
    sh.setColorRampItemList(items)
    sh.setClip(True)          # values below the first class (no oil) stay transparent
    rs = QgsRasterShader()
    rs.setRasterShaderFunction(sh)
    r = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, rs)
    r.setClassificationMin(ramp[0][0])
    r.setClassificationMax(ramp[-1][0])
    layer.setRenderer(r)
    try:
        layer.dataProvider().setUserNoDataValue(1, [])
    except Exception:
        logging.getLogger(__name__).debug("Optional operation failed in ramp_renderer; existing fallback retained.", exc_info=True)


def arrival_renderer(layer, hmax):
    stops = [(0, "#7a0177"), (hmax * 0.1, "#c51b8a"), (hmax * 0.25, "#f768a1"), (hmax * 0.5, "#fdae61"),
             (hmax * 0.75, "#fee08b"), (hmax, "#ffffbf")]
    items = [QgsColorRampShader.ColorRampItem(v, QColor(c), "%.0f h" % v) for v, c in stops]
    sh = QgsColorRampShader(0, hmax)
    sh.setColorRampType(QgsColorRampShader.Interpolated)
    sh.setColorRampItemList(items)
    rs = QgsRasterShader()
    rs.setRasterShaderFunction(sh)
    layer.setRenderer(QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, rs))


def _status_renderer(layer):
    cols = {"Oil reaches groundwater": "#6f2c00", "Oiled": "#d73a49", "Spray fallout only": "#c9a227",
            "Not reached": "#ffffff", "Outside model domain": "#bfbfbf"}
    cats = []
    for k, c in cols.items():
        sym = QgsMarkerSymbol.createSimple({"name": "circle", "color": c, "outline_color": "#222222",
                                            "size": "3.2"})
        cats.append(QgsRendererCategory(k, sym, k))
    layer.setRenderer(QgsCategorizedSymbolRenderer("status", cats))


def load_results(files, name, context=None):
    proj = QgsProject.instance()
    root = proj.layerTreeRoot()
    grp = root.insertGroup(0, "Spill - %s" % name)
    summ = files.get("summary") or {}

    def add_vec(layer_name, style):
        uri = "%s|layername=%s" % (files["gpkg"], layer_name)
        lyr = QgsVectorLayer(uri, layer_name.replace("_", " ").capitalize(), "ogr")
        if not lyr.isValid() or lyr.featureCount() == 0:
            return None
        style(lyr)
        proj.addMapLayer(lyr, False)
        grp.addLayer(lyr)
        return lyr

    add_vec("receptor_impacts", _status_renderer)
    add_vec("spill_source", lambda l: l.setRenderer(
        QgsCategorizedSymbolRenderer("kind", [QgsRendererCategory(
            k, QgsMarkerSymbol.createSimple({"name": "star", "color": "#ffd400", "outline_color": "#000000",
                                             "size": "6"}), k) for k in ("BLOWOUT", "TANK", "CONTINUOUS")])))
    add_vec("arrival_isochrones", lambda l: l.renderer().setSymbol(QgsFillSymbol.createSimple(
        {"color": "0,0,0,0", "outline_color": "#7a1fa2", "outline_width": "0.4", "outline_style": "dash"})))
    add_vec("oil_footprint", lambda l: l.renderer().setSymbol(QgsFillSymbol.createSimple(
        {"color": "0,0,0,0", "outline_color": "#b42318", "outline_width": "0.6"})))
    add_vec("groundwater_concern", lambda l: l.renderer().setSymbol(QgsFillSymbol.createSimple(
        {"color": "255,140,0,70", "outline_color": "#6f2c00", "outline_width": "0.5"})))
    add_vec("spray_fallout", lambda l: l.renderer().setSymbol(QgsFillSymbol.createSimple(
        {"color": "90,62,38,40", "outline_color": "#5a3e26", "outline_width": "0.3"})))
    add_vec("drainage_network_dem", lambda l: l.renderer().setSymbol(QgsLineSymbol.createSimple(
        {"line_color": "#1c7ed6", "line_width": "0.45"})))

    def add_ras(key, title, styler):
        p = files.get(key)
        if not p or not os.path.exists(p):
            return None
        lyr = QgsRasterLayer(p, title)
        if not lyr.isValid():
            return None
        styler(lyr)
        proj.addMapLayer(lyr, False)
        grp.addLayer(lyr)
        return lyr

    add_ras("fallout_mm", "Spray fallout (mm)", lambda l: ramp_renderer(l, outputs.FALL_RAMP))
    add_ras("arrival_time_h", "Oil arrival time (h)",
            lambda l: arrival_renderer(l, max(1.0, float(summ.get("sim_hours", 0) or 0) or 96.0)))
    add_ras("yield_stress_max_Pa", "Highest yield stress (Pa)",
            lambda l: ramp_renderer(l, outputs.YIELD_RAMP, "%g Pa"))
    add_ras("oil_temp_min_C", "Lowest oil temperature (C)", lambda l: lin_renderer(l, outputs.TEMP_RAMP, "%g C"))
    add_ras("max_thickness_mm", "Maximum oil thickness (mm)", lambda l: ramp_renderer(l, outputs.OIL_RAMP))
    for n in grp.findLayers():
        if n.name() in ("Oil arrival time (h)", "Spray fallout (mm)", "Arrival isochrones",
                        "Groundwater concern", "Spray fallout", "Highest yield stress (Pa)",
                        "Lowest oil temperature (C)"):
            n.setItemVisibilityChecked(False)
    return grp


def lin_renderer(layer, ramp, fmt):
    items = [QgsColorRampShader.ColorRampItem(s[0], QColor(s[1], s[2], s[3], max(s[4], 180)), fmt % s[0])
             for s in ramp]
    sh = QgsColorRampShader(ramp[0][0], ramp[-1][0])
    sh.setColorRampType(QgsColorRampShader.Interpolated)
    sh.setColorRampItemList(items)
    rs = QgsRasterShader()
    rs.setRasterShaderFunction(sh)
    layer.setRenderer(QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, rs))


H2S_COLOURS = {"500 ppm+": "#3b0a45", "IDLH": "#9c36b5", "AEGL-3": "#e03131", "AEGL-2": "#f08c00",
               "AEGL-1": "#2f9e44", "Odour only": "#74c0fc", "Not reached": "#ffffff"}


def load_h2s_results(files, name, context=None):
    from .engine.h2s import H2S_RAMP
    proj = QgsProject.instance()
    root = proj.layerTreeRoot()
    grp = root.insertGroup(0, "H2S - %s" % name)

    def add_vec(layer_name, title, style):
        lyr = QgsVectorLayer("%s|layername=%s" % (files["gpkg"], layer_name), title, "ogr")
        if not lyr.isValid() or lyr.featureCount() == 0:
            return None
        style(lyr)
        proj.addMapLayer(lyr, False)
        grp.addLayer(lyr)
        return lyr

    def rec_style(l):
        cats = [QgsRendererCategory(k, QgsMarkerSymbol.createSimple(
            {"name": "circle", "color": c, "outline_color": "#222222", "size": "3.2"}), k)
            for k, c in H2S_COLOURS.items()]
        l.setRenderer(QgsCategorizedSymbolRenderer("hazard", cats))

    def zone_style(l):
        root = QgsRuleBasedRenderer.Rule(None)
        for f in l.getFeatures():
            c = QColor(f["colour"])
            sym = QgsFillSymbol.createSimple({"outline_color": f["colour"], "outline_width": "0.6"})
            c.setAlpha(35)
            sym.setColor(c)
            sym.symbolLayer(0).setStrokeColor(QColor(f["colour"]))
            root.appendChild(QgsRuleBasedRenderer.Rule(sym, 0, 0, "\"threshold\" = '%s'" % f["threshold"],
                                                       f["threshold"]))
        l.setRenderer(QgsRuleBasedRenderer(root))

    add_vec("h2s_receptor_exposure", "H2S receptor exposure", rec_style)
    add_vec("h2s_source", "H2S leak point", lambda l: l.renderer().setSymbol(QgsMarkerSymbol.createSimple(
        {"name": "star", "color": "#ffd400", "outline_color": "#000000", "size": "6"})))
    add_vec("h2s_planned_zone", "Planned H2S zone", lambda l: l.renderer().setSymbol(QgsFillSymbol.createSimple(
        {"color": "0,0,0,0", "outline_color": "#000000", "outline_width": "0.6", "outline_style": "dash"})))
    add_vec("h2s_threshold_zones", "H2S threshold zones", zone_style)
    for key, title in (("h2s_max_60min_ppm", "H2S highest 60-min mean (ppm)"),
                       ("h2s_max_10min_ppm", "H2S highest 10-min mean (ppm)")):
        p = files.get(key)
        if p and os.path.exists(p):
            lyr = QgsRasterLayer(p, title)
            if lyr.isValid():
                ramp_renderer(lyr, H2S_RAMP, "%g ppm")
                proj.addMapLayer(lyr, False)
                grp.addLayer(lyr)
    for n in grp.findLayers():
        if n.name() == "H2S highest 10-min mean (ppm)":
            n.setItemVisibilityChecked(False)
    return grp
