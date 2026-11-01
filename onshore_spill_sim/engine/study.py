# -*- coding: utf-8 -*-
"""
Study set-up read from <project folder>/00_Config/model_inputs.json: the DEM, wind table,
ground-surface rules, rivers burned into the DEM, receptor layers and the layers drawn on the
HTML page.  Paths in the file are relative to the project folder, so the folder can be moved.
Used by the reference-scenario script and by the simulation panel in QGIS, so both run the
model on exactly the same inputs.
"""
import json
import os

KEYS_WITH_PATH = ("surface_rules", "receptor_specs", "rivers", "html_layers")


def find_root(start):
    """Walk up from a file or folder until a folder holding 00_Config/model_inputs.json is found."""
    d = os.path.abspath(start)
    if os.path.isfile(d):
        d = os.path.dirname(d)
    for _ in range(6):
        if os.path.exists(os.path.join(d, "00_Config", "model_inputs.json")):
            return d
        nd = os.path.dirname(d)
        if nd == d:
            break
        d = nd
    return None


def load(root):
    with open(os.path.join(root, "00_Config", "model_inputs.json")) as f:
        m = json.load(f)
    absp = lambda p: p if (p is None or os.path.isabs(p)) else os.path.normpath(os.path.join(root, p))
    for k in ("dem", "wind_csv", "gpkg"):
        if m.get(k):
            m[k] = absp(m[k])
    for k in KEYS_WITH_PATH:
        for d in m.get(k, []):
            d["path"] = absp(d.get("path") or m["gpkg"])
    wc = os.path.join(root, "00_Config", "wells_config.json")
    if os.path.exists(wc):
        with open(wc) as f:
            m["wells_config"] = json.load(f)
    m["root"] = root
    return m


def common_cfg(m):
    """Keys every oil run shares."""
    return dict(dem=m["dem"], wind_csv=m.get("wind_csv"), surface_rules=m["surface_rules"],
                receptor_specs=m["receptor_specs"], rivers=m["rivers"], html_layers=m["html_layers"],
                prepared_by=m.get("prepared_by", ""), half_width_m=m.get("half_width_m", 9000.0),
                margin_m=m.get("margin_m", 6000.0), subtitle=m.get("subtitle"))
