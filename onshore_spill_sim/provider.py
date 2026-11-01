# -*- coding: utf-8 -*-
import os
from qgis.core import QgsProcessingProvider
from qgis.PyQt.QtGui import QIcon


class SpillProvider(QgsProcessingProvider):

    def loadAlgorithms(self):
        from .alg_simulate import SimulateSpillAlgorithm
        from .alg_animate import AnimateSpillAlgorithm
        from .alg_tanker import TankerLayerAlgorithm
        from .alg_h2s import H2SAlgorithm
        self.addAlgorithm(SimulateSpillAlgorithm())
        self.addAlgorithm(TankerLayerAlgorithm())
        self.addAlgorithm(H2SAlgorithm())
        self.addAlgorithm(AnimateSpillAlgorithm())

    def id(self):
        return "onshorespill"

    def name(self):
        return "Onshore Spill Simulator"

    def icon(self):
        return QIcon(os.path.join(os.path.dirname(__file__), "icon.png"))

    def longName(self):
        return self.name()
