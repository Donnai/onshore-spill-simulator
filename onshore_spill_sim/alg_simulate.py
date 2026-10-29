# -*- coding: utf-8 -*-
"""Processing algorithm: Simulate onshore spill."""
import logging
import os
import re
import tempfile

from qgis.core import (QgsProcessing, QgsProcessingAlgorithm, QgsProcessingException,
                       QgsProcessingParameterRasterLayer, QgsProcessingParameterPoint,
                       QgsProcessingParameterEnum, QgsProcessingParameterNumber,
                       QgsProcessingParameterString, QgsProcessingParameterFile,
                       QgsProcessingParameterMultipleLayers, QgsProcessingParameterVectorLayer,
                       QgsProcessingParameterBoolean, QgsProcessingParameterFolderDestination,
                       QgsProcessingParameterDefinition, QgsCoordinateTransform,
                       QgsCoordinateReferenceSystem, QgsVectorFileWriter, QgsProject,
                       QgsRasterLayer, QgsVectorLayer, QgsWkbTypes)
from qgis.PyQt.QtGui import QIcon

from .engine.oils import OILS
from .engine import climate
from .engine.runner import run_scenario

PLUGIN_DIR = os.path.dirname(__file__)
KINDS = [("BLOWOUT", "Well oil spillage (pool at the cellar + spray fallout)"),
         ("TANK", "Tank / bund failure, road tanker (fixed volume)"),
         ("CONTINUOUS", "Continuous leak - flowline, pipe or hose (rate x duration)")]
PRODUCTS = list(OILS.keys())
NAME_FIELDS = ["name", "Name", "NAME", "facility_name", "O_Name", "TOWN_NAME", "label", "code"]


def _clean(s):
    return re.sub(r"[^A-Za-z0-9_\-]+", "_", s).strip("_") or "run"


class SimulateSpillAlgorithm(QgsProcessingAlgorithm):

    def name(self):
        return "simulate"

    def displayName(self):
        return "Simulate onshore spill"

    def group(self):
        return ""

    def groupId(self):
        return ""

    def icon(self):
        return QIcon(os.path.join(PLUGIN_DIR, "icon.png"))

    def flags(self):
        # runs in the main thread so that the animation can draw the project layers
        return super().flags() | QgsProcessingAlgorithm.FlagNoThreading

    def createInstance(self):
        return SimulateSpillAlgorithm()

    def shortHelpString(self):
        return ("Simulates an onshore oil or fuel spill over the DEM and animates it.\n\n"
                "1  Choose the DEM and click the spill point on the map (or type E,N).\n"
                "2  Choose the release type and product, then the size of the release.\n"
                "3  Set the month (climate, wind, stream flow, soil wetness and groundwater depth "
                "follow the month) and any storms as 'start h, duration h, depth mm' separated by ';'.\n"
                "4  Add the layers that describe the ground: rig pads (compacted laterite), roads, bunds "
                "and lined pits, ponds and dugouts, forest reserves; a site-plan polygon layer can be "
                "classified from its feature names.\n"
                "5  Add mapped rivers/streams (they are burned into the DEM) and the receptor layers "
                "(towns, water points, health facilities...).\n"
                "6  Pick an output folder. The run writes GeoTIFF time series, a GeoPackage of "
                "footprints, isochrones and receptor impacts, CSV/Excel tables, an interactive HTML "
                "animation and, if ticked, a GIF/MP4 drawn from the layers currently visible in the "
                "Layers panel.\n\n"
                "Defaults are set for the northern Voltaian Basin (Tamale climate). The release "
                "figures are placeholders: replace them with the worst-case discharge and the actual "
                "tank inventories for the rig selected. See Plugins > Onshore Spill Simulator > "
                "Method and user guide.")

    # ------------------------------------------------------------------
    def initAlgorithm(self, config=None):
        P = self.addParameter
        P(QgsProcessingParameterString("NAME", "Scenario name", "Well spillage scenario"))
        P(QgsProcessingParameterRasterLayer("DEM", "Digital elevation model (projected CRS, metres)"))
        P(QgsProcessingParameterPoint("SOURCE", "Spill point (click on the map)"))
        P(QgsProcessingParameterString("SOURCE_LABEL", "Label for the spill point", "Well cellar"))
        P(QgsProcessingParameterEnum("KIND", "Release type", [k[1] for k in KINDS], defaultValue=0))
        P(QgsProcessingParameterEnum("PRODUCT", "Product", [OILS[k]["label"] for k in PRODUCTS],
                                     defaultValue=0))
        P(QgsProcessingParameterNumber("RATE", "Release rate (bbl/d) - well spillage or continuous leak",
                                       QgsProcessingParameterNumber.Double, 1500.0, minValue=0.0))
        P(QgsProcessingParameterNumber("DURATION", "Release duration (h) - well spillage or continuous leak",
                                       QgsProcessingParameterNumber.Double, 72.0, minValue=0.0))
        P(QgsProcessingParameterNumber("VOLUME", "Volume released (m3) - tank / tanker",
                                       QgsProcessingParameterNumber.Double, 60.0, minValue=0.0))
        P(QgsProcessingParameterNumber("RELEASE_MIN", "Time for the volume to escape (min) - tank / tanker",
                                       QgsProcessingParameterNumber.Double, 10.0, minValue=0.1))
        P(QgsProcessingParameterNumber("CONTAINED", "Volume held before overflow (m3) - cellar or intact bund",
                                       QgsProcessingParameterNumber.Double, 3.0, minValue=0.0))
        P(QgsProcessingParameterEnum("MONTH", "Month", climate.MONTHS, defaultValue=0))
        P(QgsProcessingParameterString("START", "Start date and time (YYYY-MM-DDTHH:MM, optional)",
                                       "2027-01-15T06:00", optional=True))
        P(QgsProcessingParameterString("STORMS", "Storms: start h, duration h, depth mm; ... (blank = none)",
                                       "", optional=True))
        P(QgsProcessingParameterNumber("SIM_HOURS", "Simulated period (h)",
                                       QgsProcessingParameterNumber.Double, 96.0, minValue=0.5))
        P(QgsProcessingParameterNumber("FRAME_MIN", "Animation frame interval (min)",
                                       QgsProcessingParameterNumber.Double, 60.0, minValue=1.0))
        P(QgsProcessingParameterVectorLayer("PADS", "Rig pads / compacted laterite areas (polygons)",
                                            [QgsProcessing.TypeVectorPolygon], optional=True))
        P(QgsProcessingParameterVectorLayer("SITEPLAN", "Site plan polygons, classified by feature name "
                                            "(pad, laterite, bund, pit, cellar...)",
                                            [QgsProcessing.TypeVectorPolygon], optional=True))
        P(QgsProcessingParameterMultipleLayers("ROADS", "Roads and tracks (lines)",
                                               QgsProcessing.TypeVectorLine, optional=True))
        P(QgsProcessingParameterMultipleLayers("WATER", "Ponds, dugouts, reservoirs (polygons)",
                                               QgsProcessing.TypeVectorPolygon, optional=True))
        P(QgsProcessingParameterMultipleLayers("WOODLAND", "Forest reserves / woodland (polygons)",
                                               QgsProcessing.TypeVectorPolygon, optional=True))
        P(QgsProcessingParameterMultipleLayers("RIVERS", "Mapped rivers and streams (lines, burned into the DEM)",
                                               QgsProcessing.TypeVectorLine, optional=True))
        P(QgsProcessingParameterMultipleLayers("RECEPTORS", "Receptor layers (towns, water points, health "
                                               "facilities, dugouts...)", QgsProcessing.TypeVectorAnyGeometry,
                                               optional=True))
        P(QgsProcessingParameterBoolean("ANIMATE", "Render GIF / MP4 from the visible project layers", True))
        P(QgsProcessingParameterBoolean("HTML", "Write the interactive HTML animation", True))
        P(QgsProcessingParameterBoolean("LOAD", "Add the results to the project", True))
        P(QgsProcessingParameterString("PREPARED_BY", "Prepared by", "", optional=True))
        P(QgsProcessingParameterBoolean("HEAT", "Heat-transfer mode: oil temperature from sun, air, wind and ground; "
                                        "temperature-dependent density, viscosity and yield stress", False))
        P(QgsProcessingParameterNumber("SRC_TEMP", "Oil temperature at the source (C, blank = air temperature; "
                                       "e.g. 60 for wellhead fluid)", QgsProcessingParameterNumber.Double, None,
                                       optional=True))
        P(QgsProcessingParameterMultipleLayers("ROUGH", "Roughness zones (polygons) - Manning n override",
                                               QgsProcessing.TypeVectorPolygon, optional=True))
        P(QgsProcessingParameterMultipleLayers("BARRIERS", "Barriers - berms, bund walls, sandbag lines "
                                               "(lines or polygons)", QgsProcessing.TypeVectorAnyGeometry,
                                               optional=True))
        P(QgsProcessingParameterMultipleLayers("SINKS", "Recovery points, sumps, free-outflow cells "
                                               "(oil removed there)", QgsProcessing.TypeVectorAnyGeometry,
                                               optional=True))
        P(QgsProcessingParameterFolderDestination("OUT", "Output folder"))
        # ---- advanced
        adv = [
            QgsProcessingParameterNumber("SPRAY", "Well spillage: share of the flow atomised as spray (0-1)",
                                         QgsProcessingParameterNumber.Double, 0.40, minValue=0, maxValue=1),
            QgsProcessingParameterNumber("HEIGHT", "Well spillage: effective release height (m)",
                                         QgsProcessingParameterNumber.Double, 40.0, minValue=1),
            QgsProcessingParameterNumber("D50", "Well spillage: median droplet diameter (micrometres)",
                                         QgsProcessingParameterNumber.Double, 250.0, minValue=10),
            QgsProcessingParameterFile("WIND_CSV", "Monthly wind table (month,month_name,u_mps,v_mps,mean_speed_mps)",
                                       extension="csv", optional=True,
                                       defaultValue=os.path.join(PLUGIN_DIR, "data",
                                                                 "Tamale_DGLE_monthly_wind_1973-2026.csv")),
            QgsProcessingParameterNumber("HALF_WIDTH", "Model area half-width around the spill (m)",
                                         QgsProcessingParameterNumber.Double, 9000.0, minValue=500),
            QgsProcessingParameterNumber("CELL", "Model cell size (m, 0 = DEM cell)",
                                         QgsProcessingParameterNumber.Double, 0.0, minValue=0),
            QgsProcessingParameterNumber("SMOOTH", "DEM smoothing (cells, 0 = none)",
                                         QgsProcessingParameterNumber.Double, 1.0, minValue=0),
            QgsProcessingParameterNumber("PIT_TOL", "Depressions shallower than this are DEM noise (m)",
                                         QgsProcessingParameterNumber.Double, 3.0, minValue=0),
            QgsProcessingParameterNumber("CHANNEL_KM2", "Catchment area that starts a stream (km2)",
                                         QgsProcessingParameterNumber.Double, 1.0, minValue=0.01),
            QgsProcessingParameterNumber("SNAP", "Snap the spill point onto a stream within (m, 0 = off)",
                                         QgsProcessingParameterNumber.Double, 0.0, minValue=0),
            QgsProcessingParameterNumber("TEMP", "Ground temperature (C, blank = monthly normal)",
                                         QgsProcessingParameterNumber.Double, None, optional=True),
            QgsProcessingParameterNumber("QSPEC", "Stream specific discharge (L/s/km2, blank = monthly default)",
                                         QgsProcessingParameterNumber.Double, None, optional=True),
            QgsProcessingParameterNumber("SOILSAT", "Soil pores already water-filled (0-1, blank = monthly default)",
                                         QgsProcessingParameterNumber.Double, None, optional=True),
            QgsProcessingParameterNumber("DTW", "Depth to groundwater (m, blank = monthly default)",
                                         QgsProcessingParameterNumber.Double, None, optional=True),
            QgsProcessingParameterNumber("CN", "SCS curve number for storm runoff",
                                         QgsProcessingParameterNumber.Double, 78.0, minValue=30, maxValue=98),
            QgsProcessingParameterNumber("RIVULET", "Width of oil rivulets on sloping ground (m)",
                                         QgsProcessingParameterNumber.Double, 6.0, minValue=0.5),
            QgsProcessingParameterNumber("STRAND", "Bank stranding coefficient (1/km)",
                                         QgsProcessingParameterNumber.Double, 0.15, minValue=0),
            QgsProcessingParameterString("ROUGH_FIELD", "Field holding Manning n in the roughness zones",
                                         "manning_n", optional=True),
            QgsProcessingParameterNumber("ROUGH_N", "Manning n used when the zone has no such field",
                                         QgsProcessingParameterNumber.Double, 0.10, minValue=0.005),
            QgsProcessingParameterNumber("BARRIER_H", "Barrier height (m)",
                                         QgsProcessingParameterNumber.Double, 0.6, minValue=0.05),
            QgsProcessingParameterNumber("CLOUD", "Cloud / haze fraction for heat transfer (0-1, blank = monthly)",
                                         QgsProcessingParameterNumber.Double, None, optional=True,
                                         minValue=0, maxValue=1),
            QgsProcessingParameterNumber("RH", "Relative humidity (%, blank = monthly normal)",
                                         QgsProcessingParameterNumber.Double, None, optional=True,
                                         minValue=1, maxValue=100),
            QgsProcessingParameterNumber("ANIM_WIDTH", "Animation width (pixels)",
                                         QgsProcessingParameterNumber.Integer, 1600, minValue=640),
            QgsProcessingParameterNumber("FPS", "Animation frames per second",
                                         QgsProcessingParameterNumber.Integer, 4, minValue=1),
            QgsProcessingParameterNumber("SEED", "Wind meander random seed",
                                         QgsProcessingParameterNumber.Integer, 7),
        ]
        for a in adv:
            a.setFlags(a.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
            self.addParameter(a)

    # ------------------------------------------------------------------
    def _export(self, layer, crs, tmpdir, tag, feedback):
        """Write a layer to a GeoPackage in the model CRS; returns (path, layer name)."""
        if layer is None:
            return None
        path = os.path.join(tmpdir, "%s.gpkg" % tag)
        opts = QgsVectorFileWriter.SaveVectorOptions()
        opts.driverName = "GPKG"
        opts.layerName = "lyr"
        opts.ct = QgsCoordinateTransform(layer.crs(), crs, QgsProject.instance())
        ctx = QgsProject.instance().transformContext()
        if hasattr(QgsVectorFileWriter, "writeAsVectorFormatV3"):
            res = QgsVectorFileWriter.writeAsVectorFormatV3(layer, path, ctx, opts)
        else:
            res = QgsVectorFileWriter.writeAsVectorFormatV2(layer, path, ctx, opts)
        if res[0] != QgsVectorFileWriter.NoError:
            feedback.reportError("Could not export %s: %s" % (layer.name(), res[1]))
            return None
        return path, "lyr"

    def _name_field(self, layer):
        names = [f.name() for f in layer.fields()]
        for n in NAME_FIELDS:
            if n in names:
                return n
        for f in layer.fields():
            if f.type() == 10:   # QVariant.String
                return f.name()
        return names[0] if names else None

    def processAlgorithm(self, parameters, context, feedback):
        cfg = self.prepare(parameters, context, feedback)
        return self.execute(cfg, parameters, context, feedback)

    def source_point(self, parameters, context, crs):
        return self.parameterAsPoint(parameters, "SOURCE", context, crs)

    def prepare(self, parameters, context, feedback):
        dem = self.parameterAsRasterLayer(parameters, "DEM", context)
        if dem is None:
            raise QgsProcessingException("Choose a DEM.")
        crs = dem.crs()
        if crs.isGeographic():
            raise QgsProcessingException("The DEM must be in a projected CRS with metre units "
                                         "(e.g. WGS 84 / UTM zone 30N, EPSG:32630).")
        pt = self.source_point(parameters, context, crs)
        name = self.parameterAsString(parameters, "NAME", context)
        out = self.parameterAsString(parameters, "OUT", context)
        os.makedirs(out, exist_ok=True)
        tmpdir = tempfile.mkdtemp(prefix="spill_inputs_")
        log = feedback.pushInfo

        def multi(key):
            return self.parameterAsLayerList(parameters, key, context) or []

        surf, rivers, recs, html_layers = [], [], [], []
        pads = self.parameterAsVectorLayer(parameters, "PADS", context)
        site = self.parameterAsVectorLayer(parameters, "SITEPLAN", context)
        k = 0
        # later rules overwrite earlier ones: woodland < roads < pads < site plan < open water
        groups = [(multi("WOODLAND"), 6, "woodland", None), (multi("ROADS"), 2, "roads", 4.0),
                  ([pads] if pads else [], 1, "pads", None), ([site] if site else [], -1, "siteplan", None),
                  (multi("WATER"), 4, "water", None)]
        for lys, cls, tag, buf in groups:
            for lyr in lys:
                e = self._export(lyr, crs, tmpdir, "%s_%d" % (tag, k), feedback)
                k += 1
                if not e:
                    continue
                if cls == -1:
                    d = dict(path=e[0], layer=e[1], name_rules=True,
                             name_field=self._name_field(lyr) or "name")
                else:
                    d = dict(path=e[0], layer=e[1], cls=cls)
                    if buf:
                        d["buffer_m"] = buf
                    if cls == 4:
                        d["all_touched"] = True
                surf.append(d)
                html_layers.append(self._html_style(lyr, e))
        for lyr in multi("RIVERS"):
            e = self._export(lyr, crs, tmpdir, "river_%d" % k, feedback)
            k += 1
            if e:
                rivers.append(dict(path=e[0], layer=e[1], min_area_km2=5.0))
                html_layers.append(self._html_style(lyr, e))
        for lyr in multi("RECEPTORS"):
            e = self._export(lyr, crs, tmpdir, "rec_%d" % k, feedback)
            k += 1
            if e:
                nf = self._name_field(lyr)
                recs.append(dict(path=e[0], layer=e[1], type=lyr.name(), name_field=nf))
                st = self._html_style(lyr, e)
                st["label"] = nf
                html_layers.append(st)

        storms = []
        txt = self.parameterAsString(parameters, "STORMS", context).strip()
        if txt:
            for part in txt.split(";"):
                vals = [float(v) for v in re.split(r"[,\s]+", part.strip()) if v]
                if len(vals) != 3:
                    raise QgsProcessingException("Storms must be written as 'start h, duration h, depth mm; ...'")
                storms.append(tuple(vals))

        def opt(key):
            v = parameters.get(key)
            if v is None or v == "":
                return None
            return self.parameterAsDouble(parameters, key, context)

        kind = KINDS[self.parameterAsEnum(parameters, "KIND", context)][0]
        start = self.parameterAsString(parameters, "START", context).strip() or None
        if start and len(start) == 16:
            start += ":00"
        cfg = dict(
            name=name, title=name, well=self.parameterAsString(parameters, "SOURCE_LABEL", context),
            prefix=_clean(name), kind=kind, source_xy=(pt.x(), pt.y()),
            source_label=self.parameterAsString(parameters, "SOURCE_LABEL", context),
            dem=dem.source(), product=PRODUCTS[self.parameterAsEnum(parameters, "PRODUCT", context)],
            rate_bbl_d=self.parameterAsDouble(parameters, "RATE", context),
            duration_h=self.parameterAsDouble(parameters, "DURATION", context),
            volume_m3=self.parameterAsDouble(parameters, "VOLUME", context),
            release_minutes=self.parameterAsDouble(parameters, "RELEASE_MIN", context),
            contained_m3=self.parameterAsDouble(parameters, "CONTAINED", context),
            spray_fraction=self.parameterAsDouble(parameters, "SPRAY", context),
            release_height_m=self.parameterAsDouble(parameters, "HEIGHT", context),
            d50_um=self.parameterAsDouble(parameters, "D50", context),
            month=self.parameterAsEnum(parameters, "MONTH", context), start_iso=start,
            sim_hours=self.parameterAsDouble(parameters, "SIM_HOURS", context),
            frame_minutes=self.parameterAsDouble(parameters, "FRAME_MIN", context),
            half_width_m=self.parameterAsDouble(parameters, "HALF_WIDTH", context),
            cell_m=self.parameterAsDouble(parameters, "CELL", context) or None,
            smooth_sigma=self.parameterAsDouble(parameters, "SMOOTH", context),
            pit_tolerance_m=self.parameterAsDouble(parameters, "PIT_TOL", context),
            channel_area_km2=self.parameterAsDouble(parameters, "CHANNEL_KM2", context),
            snap_to_channel_m=self.parameterAsDouble(parameters, "SNAP", context),
            temp_c=opt("TEMP"), q_spec=opt("QSPEC"), soil_sat=opt("SOILSAT"), dtw_m=opt("DTW"),
            cn=self.parameterAsDouble(parameters, "CN", context),
            rivulet_width_m=self.parameterAsDouble(parameters, "RIVULET", context),
            strand_per_m=self.parameterAsDouble(parameters, "STRAND", context) / 1000.0,
            wind_csv=self.parameterAsFile(parameters, "WIND_CSV", context) or None,
            wind_seed=self.parameterAsInt(parameters, "SEED", context),
            rain=storms, rivers=rivers, surface_rules=surf, receptor_specs=recs,
            html=self.parameterAsBool(parameters, "HTML", context), html_layers=html_layers,
            prepared_by=self.parameterAsString(parameters, "PREPARED_BY", context),
            out_dir=out)
        # heat transfer and control features
        cfg["heat_mode"] = self.parameterAsBool(parameters, "HEAT", context)
        cfg["discharge_temp_c"] = opt("SRC_TEMP")
        cfg["cloud"] = opt("CLOUD")
        cfg["rh_pct"] = opt("RH")
        rf = (self.parameterAsString(parameters, "ROUGH_FIELD", context) or "").strip()
        zones = []
        for lyr in multi("ROUGH"):
            e = self._export(lyr, crs, tmpdir, "rough_%d" % k, feedback)
            k += 1
            if e:
                z = dict(path=e[0], layer=e[1], n=self.parameterAsDouble(parameters, "ROUGH_N", context))
                if rf and rf in [f.name() for f in lyr.fields()]:
                    z["attribute"] = rf
                zones.append(z)
        cfg["roughness_zones"] = zones
        cfg["barriers"] = []
        for lyr in multi("BARRIERS"):
            e = self._export(lyr, crs, tmpdir, "barrier_%d" % k, feedback)
            k += 1
            if e:
                cfg["barriers"].append(dict(path=e[0], layer=e[1],
                                            height_m=self.parameterAsDouble(parameters, "BARRIER_H", context)))
                html_layers.append(self._html_style(lyr, e))
        cfg["sinks"] = []
        for lyr in multi("SINKS"):
            e = self._export(lyr, crs, tmpdir, "sink_%d" % k, feedback)
            k += 1
            if e:
                cfg["sinks"].append(dict(path=e[0], layer=e[1]))
                html_layers.append(self._html_style(lyr, e))
        cfg["_name"] = name
        return cfg

    def execute(self, cfg, parameters, context, feedback):
        log = feedback.pushInfo
        name = cfg.pop("_name", cfg["name"])
        out = cfg["out_dir"]
        res = run_scenario(cfg, log=log, feedback=feedback)
        files = res.files
        results = {"OUT": out, "RUN_JSON": files["run_json"], "GPKG": files["gpkg"]}
        if files.get("html"):
            results["HTML"] = files["html"]
            log("HTML animation: %s" % files["html"])
        if self.parameterAsBool(parameters, "ANIMATE", context):
            from .qgis_anim import render_animation
            vis = self._visible_layers()
            log("Drawing the animation over %d visible project layers..." % len(vis))
            a = render_animation(files["run_json"], vis, width=self.parameterAsInt(parameters, "ANIM_WIDTH", context),
                                 height=int(self.parameterAsInt(parameters, "ANIM_WIDTH", context) * 0.625),
                                 fps=self.parameterAsInt(parameters, "FPS", context), log=log, feedback=feedback,
                                 prepared_by=cfg["prepared_by"])
            results.update({k.upper(): v for k, v in a.items() if k != "frames"})
        if self.parameterAsBool(parameters, "LOAD", context):
            from .qgis_style import load_results
            load_results(files, name, context)
        return results

    def _visible_layers(self):
        root = QgsProject.instance().layerTreeRoot()
        return [n.layer() for n in root.findLayers() if n.isVisible() and n.layer() is not None]

    def _html_style(self, layer, exported):
        """Carry the layer's main colour into the HTML map."""
        st = dict(name=layer.name(), path=exported[0], layer=exported[1], visible=True)
        col = "#555555"
        try:
            sym = layer.renderer().symbol()
            col = sym.color().name()
        except Exception:
            logging.getLogger(__name__).debug("Optional operation failed in SimulateSpillAlgorithm._html_style; existing fallback retained.", exc_info=True)
        g = layer.geometryType()
        if g == QgsWkbTypes.PointGeometry:
            st.update(kind="point", style=dict(color="#222222", fill=col, radius=5), permanent=False)
        elif g == QgsWkbTypes.LineGeometry:
            st.update(style=dict(color=col, weight=2))
        else:
            st.update(style=dict(color=col, weight=1.2, fill=col, fillOpacity=0.25))
        st["label"] = self._name_field(layer) or "name"
        return st
