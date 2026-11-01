# -*- coding: utf-8 -*-
"""Processing algorithm: redraw the animation of a finished run with the current map."""
import os

from qgis.core import (QgsProcessingAlgorithm, QgsProcessingParameterFile, QgsProcessingParameterNumber,
                       QgsProcessingParameterExtent, QgsProcessingParameterString,
                       QgsProcessingParameterBoolean, QgsProcessingException, QgsProject)
from qgis.PyQt.QtGui import QIcon

PLUGIN_DIR = os.path.dirname(__file__)


class AnimateSpillAlgorithm(QgsProcessingAlgorithm):

    def name(self):
        return "animate"

    def displayName(self):
        return "Render animation from a finished run"

    def group(self):
        return ""

    def groupId(self):
        return ""

    def icon(self):
        return QIcon(os.path.join(PLUGIN_DIR, "icon.png"))

    def flags(self):
        return super().flags() | QgsProcessingAlgorithm.FlagNoThreading

    def createInstance(self):
        return AnimateSpillAlgorithm()

    def shortHelpString(self):
        return ("Draws the time steps of a finished run (the *_run.json file in its output folder) over "
                "the layers that are visible in the Layers panel, using their current symbology, and "
                "writes PNG frames, an animated GIF and an MP4 (when ffmpeg is installed).\n\n"
                "Style and order the map as you want it to appear, optionally set the map extent, "
                "then run. The project CRS is used for the drawing.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFile("RUN", "Run record (*_run.json)", extension="json"))
        self.addParameter(QgsProcessingParameterExtent("EXTENT", "Map extent (blank = fit the spill)",
                                                       optional=True))
        self.addParameter(QgsProcessingParameterNumber("WIDTH", "Width (pixels)",
                                                       QgsProcessingParameterNumber.Integer, 1600, minValue=640))
        self.addParameter(QgsProcessingParameterNumber("HEIGHT", "Height (pixels)",
                                                       QgsProcessingParameterNumber.Integer, 1000, minValue=400))
        self.addParameter(QgsProcessingParameterNumber("FPS", "Frames per second",
                                                       QgsProcessingParameterNumber.Integer, 4, minValue=1))
        self.addParameter(QgsProcessingParameterString("TITLE", "Title (blank = scenario title)", optional=True))
        self.addParameter(QgsProcessingParameterString("PREPARED_BY", "Prepared by", optional=True))
        self.addParameter(QgsProcessingParameterBoolean("GIF", "Write GIF", True))
        self.addParameter(QgsProcessingParameterBoolean("MP4", "Write MP4 (needs ffmpeg)", True))

    def processAlgorithm(self, parameters, context, feedback):
        from .qgis_anim import render_animation
        run = self.parameterAsFile(parameters, "RUN", context)
        if not run or not os.path.exists(run):
            raise QgsProcessingException("Choose the *_run.json file of a finished run.")
        crs = QgsProject.instance().crs()
        ext = None
        if parameters.get("EXTENT"):
            ext = self.parameterAsExtent(parameters, "EXTENT", context, crs)
            if ext.isEmpty():
                ext = None
        root = QgsProject.instance().layerTreeRoot()
        layers = [n.layer() for n in root.findLayers() if n.isVisible() and n.layer() is not None]
        out = render_animation(run, layers, width=self.parameterAsInt(parameters, "WIDTH", context),
                               height=self.parameterAsInt(parameters, "HEIGHT", context), extent=ext,
                               dest_crs=crs, fps=self.parameterAsInt(parameters, "FPS", context),
                               make_gif=self.parameterAsBool(parameters, "GIF", context),
                               make_mp4=self.parameterAsBool(parameters, "MP4", context),
                               log=feedback.pushInfo, feedback=feedback,
                               title=self.parameterAsString(parameters, "TITLE", context) or None,
                               prepared_by=self.parameterAsString(parameters, "PREPARED_BY", context) or None)
        return {k.upper(): v for k, v in out.items() if k != "frames"}
