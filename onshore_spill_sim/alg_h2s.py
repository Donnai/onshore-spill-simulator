# -*- coding: utf-8 -*-
"""Processing algorithm: H2S dispersion from a gas leakage at a well (Gaussian puff model)."""
import os
import re
import tempfile

from qgis.core import (QgsProcessing, QgsProcessingAlgorithm, QgsProcessingException,
                       QgsProcessingParameterRasterLayer, QgsProcessingParameterPoint,
                       QgsProcessingParameterEnum, QgsProcessingParameterNumber,
                       QgsProcessingParameterString, QgsProcessingParameterFile,
                       QgsProcessingParameterMultipleLayers, QgsProcessingParameterBoolean,
                       QgsProcessingParameterFolderDestination, QgsProcessingParameterDefinition)
from qgis.PyQt.QtGui import QIcon

from .engine import climate
from .engine.h2s import run_h2s
from .alg_simulate import SimulateSpillAlgorithm, _clean, PLUGIN_DIR

MODES = ["Free vertical jet from the rig floor (momentum and buoyancy rise)",
         "Impinged leak - jet stopped by the substructure / BOP, no plume rise (conservative)"]


class H2SAlgorithm(QgsProcessingAlgorithm):

    def name(self):
        return "h2s"

    def displayName(self):
        return "Simulate H2S leakage at a well"

    def icon(self):
        return QIcon(os.path.join(PLUGIN_DIR, "icon_h2s.png"))

    def flags(self):
        return super().flags() | QgsProcessingAlgorithm.FlagNoThreading

    def createInstance(self):
        return H2SAlgorithm()

    def shortHelpString(self):
        return ("Ground-level H2S concentrations downwind of a sour-gas leakage, hour by hour, with the "
                "wind, stability and mixing height of the chosen month.\n\n"
                "Source: gas rate (MMscf/d) x H2S content (ppm by volume). The figures in the form are "
                "placeholders - use the drilling programme's worst-case gas rate and the H2S content "
                "expected for the well.\n\n"
                "Outputs: H2S time series and worst 10-min / 60-min concentration rasters, AEGL-1/2/3, "
                "IDLH and 500 ppm zones, receptor exposure table, interactive HTML animation and, if "
                "ticked, a GIF/MP4 drawn from the visible layers.")

    def initAlgorithm(self, config=None):
        P = self.addParameter
        P(QgsProcessingParameterString("NAME", "Scenario name", "H2S leakage"))
        P(QgsProcessingParameterRasterLayer("DEM", "Digital elevation model (projected CRS, metres)"))
        P(QgsProcessingParameterPoint("SOURCE", "Leak point - the well (click on the map)"))
        P(QgsProcessingParameterNumber("GAS", "Gas flow rate (MMscf/d)", QgsProcessingParameterNumber.Double,
                                       10.0, minValue=0.001))
        P(QgsProcessingParameterNumber("H2S_PPM", "H2S content of the gas (ppm by volume; 10,000 = 1%)",
                                       QgsProcessingParameterNumber.Double, 10000.0, minValue=0.1))
        P(QgsProcessingParameterEnum("MODE", "Leak mode", MODES, defaultValue=1))
        P(QgsProcessingParameterNumber("HEIGHT", "Leak height above ground (m)",
                                       QgsProcessingParameterNumber.Double, 3.0, minValue=0.5))
        P(QgsProcessingParameterNumber("DURATION", "Leak duration (h)", QgsProcessingParameterNumber.Double,
                                       24.0, minValue=0.1))
        P(QgsProcessingParameterNumber("BUILDUP", "Leak build-up time (h): grows from 10 % to the full rate; 0 = full rate "
                                       "from the start (worst case)", QgsProcessingParameterNumber.Double, 0.0, minValue=0.0))
        P(QgsProcessingParameterNumber("SIM_HOURS", "Simulated period (h)", QgsProcessingParameterNumber.Double,
                                       24.0, minValue=0.5))
        P(QgsProcessingParameterNumber("FRAME_MIN", "Frame interval (min)", QgsProcessingParameterNumber.Double,
                                       10.0, minValue=2.0))
        P(QgsProcessingParameterEnum("MONTH", "Month", climate.MONTHS, defaultValue=0))
        P(QgsProcessingParameterString("START", "Start date and time (YYYY-MM-DDTHH:MM)", "2027-01-15T18:00",
                                       optional=True))
        P(QgsProcessingParameterMultipleLayers("RECEPTORS", "Receptor layers (towns, camps, muster points...)",
                                               QgsProcessing.TypeVectorAnyGeometry, optional=True))
        P(QgsProcessingParameterNumber("ZONE", "Planned H2S zone radius to compare with (m, 0 = none)",
                                       QgsProcessingParameterNumber.Double, 310.0, minValue=0))
        P(QgsProcessingParameterBoolean("ANIMATE", "Render GIF / MP4 from the visible project layers", True))
        P(QgsProcessingParameterBoolean("HTML", "Write the interactive HTML animation", True))
        P(QgsProcessingParameterBoolean("LOAD", "Add the results to the project", True))
        P(QgsProcessingParameterString("PREPARED_BY", "Prepared by", "", optional=True))
        P(QgsProcessingParameterFolderDestination("OUT", "Output folder"))
        adv = [QgsProcessingParameterNumber("DIAM", "Exit diameter of the flow path (m)",
                                            QgsProcessingParameterNumber.Double, 0.20, minValue=0.01),
               QgsProcessingParameterNumber("GAS_T", "Gas temperature at the exit (C)",
                                            QgsProcessingParameterNumber.Double, 40.0),
               QgsProcessingParameterNumber("GAS_MW", "Molecular weight of the hydrocarbon gas (g/mol)",
                                            QgsProcessingParameterNumber.Double, 19.0, minValue=2),
               QgsProcessingParameterNumber("HALF_WIDTH", "Model area half-width (m)",
                                            QgsProcessingParameterNumber.Double, 15000.0, minValue=1000),
               QgsProcessingParameterNumber("CELL", "Cell size (m)", QgsProcessingParameterNumber.Double,
                                            50.0, minValue=10),
               QgsProcessingParameterNumber("PUFF", "Puff interval (s)", QgsProcessingParameterNumber.Double,
                                            90.0, minValue=10),
               QgsProcessingParameterFile("WIND_CSV", "Monthly wind table", extension="csv", optional=True,
                                          defaultValue=os.path.join(PLUGIN_DIR, "data",
                                                                    "Tamale_DGLE_monthly_wind_1973-2026.csv")),
               QgsProcessingParameterNumber("SEED", "Wind meander random seed",
                                            QgsProcessingParameterNumber.Integer, 7),
               QgsProcessingParameterNumber("ANIM_WIDTH", "Animation width (pixels)",
                                            QgsProcessingParameterNumber.Integer, 1600, minValue=640),
               QgsProcessingParameterNumber("FPS", "Animation frames per second",
                                            QgsProcessingParameterNumber.Integer, 6, minValue=1)]
        for a in adv:
            a.setFlags(a.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
            self.addParameter(a)

    def processAlgorithm(self, parameters, context, feedback):
        dem = self.parameterAsRasterLayer(parameters, "DEM", context)
        if dem is None or dem.crs().isGeographic():
            raise QgsProcessingException("Choose a DEM in a projected CRS (metres).")
        crs = dem.crs()
        pt = self.parameterAsPoint(parameters, "SOURCE", context, crs)
        out = self.parameterAsString(parameters, "OUT", context)
        os.makedirs(out, exist_ok=True)
        tmp = tempfile.mkdtemp(prefix="h2s_inputs_")
        helper = SimulateSpillAlgorithm()
        recs, html_layers = [], []
        for i, lyr in enumerate(self.parameterAsLayerList(parameters, "RECEPTORS", context) or []):
            e = helper._export(lyr, crs, tmp, "rec_%d" % i, feedback)
            if e:
                nf = helper._name_field(lyr)
                recs.append(dict(path=e[0], layer=e[1], type=lyr.name(), name_field=nf))
                st = helper._html_style(lyr, e)
                st["label"] = nf
                html_layers.append(st)
        start = (self.parameterAsString(parameters, "START", context) or "").strip() or None
        if start and len(start) == 16:
            start += ":00"
        name = self.parameterAsString(parameters, "NAME", context)
        mode = self.parameterAsEnum(parameters, "MODE", context)
        zone = self.parameterAsDouble(parameters, "ZONE", context)
        cfg = dict(name=name, title=name, well=name, dem=dem.source(), source_xy=(pt.x(), pt.y()),
                   gas_mmscfd=self.parameterAsDouble(parameters, "GAS", context),
                   h2s_ppm=self.parameterAsDouble(parameters, "H2S_PPM", context),
                   jet=(mode == 0), release_h_m=self.parameterAsDouble(parameters, "HEIGHT", context),
                   exit_d_m=self.parameterAsDouble(parameters, "DIAM", context),
                   gas_temp_c=self.parameterAsDouble(parameters, "GAS_T", context),
                   gas_mw=self.parameterAsDouble(parameters, "GAS_MW", context),
                   duration_h=self.parameterAsDouble(parameters, "DURATION", context),
                   buildup_h=self.parameterAsDouble(parameters, "BUILDUP", context),
                   sim_hours=self.parameterAsDouble(parameters, "SIM_HOURS", context),
                   frame_minutes=self.parameterAsDouble(parameters, "FRAME_MIN", context),
                   month=self.parameterAsEnum(parameters, "MONTH", context), start_iso=start,
                   wind_csv=self.parameterAsFile(parameters, "WIND_CSV", context) or None,
                   wind_seed=self.parameterAsInt(parameters, "SEED", context),
                   half_width_m=self.parameterAsDouble(parameters, "HALF_WIDTH", context),
                   cell_m=self.parameterAsDouble(parameters, "CELL", context),
                   puff_s=self.parameterAsDouble(parameters, "PUFF", context),
                   receptor_specs=recs, html_layers=html_layers,
                   zone_m=zone or None, zone_label=("Planned H2S zone %.0f m" % zone) if zone else None,
                   out_dir=out, prefix=_clean(name), html=self.parameterAsBool(parameters, "HTML", context),
                   prepared_by=self.parameterAsString(parameters, "PREPARED_BY", context),
                   subtitle="H2S dispersion from a gas leakage at the well")
        res = run_h2s(cfg, log=feedback.pushInfo, feedback=feedback)
        files = res["files"]
        results = {"OUT": out, "RUN_JSON": files["run_json"], "GPKG": files["gpkg"]}
        if files.get("html"):
            results["HTML"] = files["html"]
        if self.parameterAsBool(parameters, "ANIMATE", context):
            from .qgis_anim import render_animation
            vis = helper._visible_layers()
            w = self.parameterAsInt(parameters, "ANIM_WIDTH", context)
            a = render_animation(files["run_json"], vis, width=w, height=int(w * 0.625),
                                 fps=self.parameterAsInt(parameters, "FPS", context), log=feedback.pushInfo,
                                 feedback=feedback, prepared_by=cfg["prepared_by"])
            results.update({k.upper(): v for k, v in a.items() if k != "frames"})
        if self.parameterAsBool(parameters, "LOAD", context):
            from .qgis_style import load_h2s_results
            load_h2s_results(files, name, context)
        return results
