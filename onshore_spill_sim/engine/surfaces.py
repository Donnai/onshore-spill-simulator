# -*- coding: utf-8 -*-
"""Build the surface-class raster from vector layers (later rules win)."""
import numpy as np
from osgeo import ogr, gdal
from .terrain import rasterize_layer

# keyword rules applied to site-plan polygon names (case-insensitive)
SITE_RULES = [
    (("drilling pad", "rig pad", "rig area", "pad (", "camp envelope", "car park",
      "laydown", "helipad"), 1),
    (("laterite", "carriageway"), 2),
    (("bund", "fuel storage", "cellar", "pit", "mud system", "septic", "sump",
      "chemical", "hazardous"), 5),
]


def classify_site_name(name):
    n = (name or "").lower()
    cls = None
    for keys, c in SITE_RULES:
        if any(k in n for k in keys):
            cls = c
    return cls


def build_surface_classes(grid, rules, log=print):
    """rules: list of dicts {path, layer, where, cls, name_rules(bool), buffer_m}
    Returns int16 array of SURFACES keys (0 = natural savannah default)."""
    out = np.zeros((grid.nrows, grid.ncols), dtype=np.int16)
    for r in rules:
        ds = ogr.Open(r["path"])
        if ds is None:
            log("  surface rule skipped (cannot open %s)" % r["path"])
            continue
        lyr = ds.GetLayerByName(r["layer"]) if r.get("layer") else ds.GetLayer(0)
        if lyr is None:
            log("  surface rule skipped (no layer %s)" % r.get("layer"))
            continue
        if r.get("name_rules"):
            # one pass per class, using the polygon names
            mem = ogr.GetDriverByName("Memory").CreateDataSource("m")
            ml = mem.CreateLayer("x", lyr.GetSpatialRef(), ogr.wkbMultiPolygon)
            ml.CreateField(ogr.FieldDefn("cls", ogr.OFTInteger))
            fld = r.get("name_field", "name")
            for f in lyr:
                c = classify_site_name(f.GetField(fld))
                g = f.GetGeometryRef()
                if c is None or g is None:
                    continue
                nf = ogr.Feature(ml.GetLayerDefn())
                nf.SetField("cls", c)
                nf.SetGeometry(g.Clone())
                ml.CreateFeature(nf)
            for c in (1, 2, 5):
                ml.SetAttributeFilter("cls = %d" % c)
                if ml.GetFeatureCount() == 0:
                    continue
                arr = rasterize_layer(mem, "x", grid, all_touched=(c == 5))
                out[arr > 0] = c
            continue
        src = ds
        lname = r.get("layer")
        if r.get("buffer_m"):
            mem = ogr.GetDriverByName("Memory").CreateDataSource("b")
            ml = mem.CreateLayer("x", lyr.GetSpatialRef(), ogr.wkbMultiPolygon)
            if r.get("where"):
                lyr.SetAttributeFilter(r["where"])
            for f in lyr:
                g = f.GetGeometryRef()
                if g is None:
                    continue
                nf = ogr.Feature(ml.GetLayerDefn())
                nf.SetGeometry(g.Buffer(r["buffer_m"]))
                ml.CreateFeature(nf)
            src, lname = mem, "x"
            arr = rasterize_layer(src, lname, grid, all_touched=r.get("all_touched", False))
        else:
            arr = rasterize_layer(src, lname, grid, where=r.get("where"),
                                  all_touched=r.get("all_touched", False))
        out[arr > 0] = r["cls"]
    return out
