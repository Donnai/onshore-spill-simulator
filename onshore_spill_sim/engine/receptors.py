# -*- coding: utf-8 -*-
"""Collect receptor points (towns, water points, health facilities, dugouts,
forest reserves, ...) from vector layers into plain dicts."""
import logging
from osgeo import ogr, osr


def load_receptors(specs, target_wkt, extent=None, log=print):
    """specs: list of dicts {path, layer, type, name_field, where}
    Polygons are represented by their nearest-to-centre point on surface.
    extent: (xmin, ymin, xmax, ymax) in target CRS to keep."""
    tgt = osr.SpatialReference()
    tgt.ImportFromWkt(target_wkt)
    tgt.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    out = []
    for s in specs:
        ds = ogr.Open(s["path"])
        if ds is None:
            log("  receptors: cannot open %s" % s["path"])
            continue
        lyr = ds.GetLayerByName(s["layer"]) if s.get("layer") else ds.GetLayer(0)
        if lyr is None:
            continue
        if s.get("where"):
            lyr.SetAttributeFilter(s["where"])
        src = lyr.GetSpatialRef()
        tr = None
        if src is not None:
            src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            if not src.IsSame(tgt):
                tr = osr.CoordinateTransformation(src, tgt)
        nf = s.get("name_field", "name")
        idx = lyr.GetLayerDefn().GetFieldIndex(nf)
        n = 0
        for f in lyr:
            g = f.GetGeometryRef()
            if g is None:
                continue
            g = g.Clone()
            if tr:
                g.Transform(tr)
            if g.IsEmpty():
                continue
            gt = ogr.GT_Flatten(g.GetGeometryType())
            try:
                if gt in (ogr.wkbPolygon, ogr.wkbMultiPolygon):
                    p = g.PointOnSurface()
                elif gt in (ogr.wkbLineString, ogr.wkbMultiLineString):
                    # nearest-to-midpoint vertex keeps the point on the line
                    pts = []
                    for i in range(max(g.GetGeometryCount(), 1)):
                        part = g.GetGeometryRef(i) if g.GetGeometryCount() else g
                        pts += part.GetPoints() or []
                    if not pts:
                        continue
                    p = ogr.CreateGeometryFromWkt("POINT (%f %f)" % pts[len(pts) // 2][:2])
                else:
                    p = g.Centroid() if gt == ogr.wkbMultiPoint else g
                x, y = p.GetX(), p.GetY()
            except Exception:
                logging.getLogger(__name__).debug("Optional operation failed in load_receptors; existing fallback retained.", exc_info=True); continue
            if extent and not (extent[0] <= x <= extent[2] and extent[1] <= y <= extent[3]):
                continue
            nm = f.GetField(idx) if idx >= 0 else None
            if not nm:
                nm = "%s %d" % (s["type"], f.GetFID())
            out.append(dict(name=str(nm), type=s["type"], x=x, y=y))
            n += 1
        log("  receptors: %d x %s" % (n, s["type"]))
    return out
