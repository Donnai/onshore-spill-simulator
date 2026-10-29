# Onshore Spill Simulator

QGIS plugin by Donald Annor Appiah for onshore oil/fuel spill and H2S simulations, with Processing algorithms and map animations.

## Installation

Use the installable `onshore_spill_sim.zip` in QGIS: **Plugins > Manage and Install Plugins > Install from ZIP**. The ZIP must contain the single `onshore_spill_sim/` folder. Exclude all `__pycache__` folders and compiled `.pyc`/`.pyo` files when packaging it.

The four algorithms are **Simulate onshore spill**, **Simulate spills at movable tanker points**, **Simulate H2S leakage at a well**, and **Render animation from a finished run**. Read `onshore_spill_sim/help/index.html` for the method and user guide.

## Version 2.1.1

Maintenance release correcting ZIP packaging and metadata. Internal project references have been removed from public documentation and one interface label; model calculations are unchanged. Software smoke tests used QGIS 3.34.4 and synthetic data; they do not validate physical predictions or compatibility with every QGIS version.

## Dependencies

QGIS/PyQGIS, NumPy, and GDAL/OGR/OSR are required. Optional: openpyxl for Excel, Pillow for GIFs, and ffmpeg or imageio-ffmpeg for MP4s.

Source: https://github.com/Donnai/onshore-spill-simulator

Issues: https://github.com/Donnai/onshore-spill-simulator/issues

The author must select a GPL-compatible license before submission to the official QGIS Plugins Directory.
