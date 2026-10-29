# Onshore Spill Simulator

QGIS plugin by Donald Annor Appiah for onshore oil/fuel spill and H2S simulations,
with Processing algorithms and map animations.

## Installation and use

In QGIS use **Plugins > Manage and Install Plugins > Install from ZIP** with
`onshore_spill_sim.zip`. The archive contains the single `onshore_spill_sim/`
folder. Exclude Python cache folders and compiled bytecode when packaging.

The Processing Toolbox contains **Simulate onshore spill**, **Simulate spills at
movable tanker points**, **Simulate H2S leakage at a well**, and **Render animation
from a finished run**. Open **Plugins > Onshore Spill Simulator > Method and user
guide** for the full guide in `onshore_spill_sim/help/index.html`.

Use a projected DEM in metres, supply release and environmental inputs, select
an output folder and run. The simulation panel also reads a study's
`00_Config/model_inputs.json`; operational projects and data are not included.

## Dependencies

QGIS/PyQGIS, NumPy and GDAL/OGR/OSR. Optional exporters: openpyxl for Excel,
Pillow for GIFs, and ffmpeg or imageio-ffmpeg for MP4s.

## Version 2.1.2

Adds diagnostic logging and validated, shell-free video encoding. Also corrects ZIP packaging and metadata and removes internal project references from
public help and one interface label. Model calculations are unchanged. Software
checks passed in QGIS 3.34.4 with synthetic data; these do not validate physical
predictions or every QGIS version.

## License

This plugin is distributed under the GNU General Public License, version 3 or
any later version (GPL-3.0-or-later). See `LICENSE` and
`onshore_spill_sim/LICENSE`. Bundled Leaflet 1.9.4 retains its BSD 2-Clause
license in `onshore_spill_sim/engine/web/LEAFLET_LICENSE.txt`.

Source: https://github.com/Donnai/onshore-spill-simulator

Issues: https://github.com/Donnai/onshore-spill-simulator/issues
