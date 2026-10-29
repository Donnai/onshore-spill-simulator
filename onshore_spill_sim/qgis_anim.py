# -*- coding: utf-8 -*-
"""
Frame renderer: draws every time step of a finished run on top of the
project layers (same symbology as the QGIS map), adds a title block, clock,
legend, wind arrow, rain indicator, scale bar and north arrow, then writes
PNG frames, an animated GIF and (when ffmpeg is available) an MP4.
"""
import datetime as _dt
import json
import math
import os
import shutil
import subprocess
import tempfile

import numpy as np
from osgeo import gdal

from qgis.core import (QgsCoordinateReferenceSystem, QgsCoordinateTransform, QgsMapRendererParallelJob,
                       QgsMapSettings, QgsProject, QgsRasterLayer, QgsRectangle, QgsVectorLayer,
                       QgsMultiBandColorRenderer, QgsRuleBasedRenderer, QgsLineSymbol, QgsPointXY)
from qgis.PyQt.QtCore import QSize, Qt, QRectF, QPointF, QRect
from qgis.PyQt.QtGui import (QColor, QFont, QImage, QPainter, QPen, QBrush, QLinearGradient,
                             QPolygonF, QFontMetrics)

from .engine import outputs

DIRS = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]


STAB = {"A": "very unstable", "B": "unstable", "C": "slightly unstable", "D": "neutral",
        "E": "slightly stable", "F": "stable"}


def _render(settings):
    job = QgsMapRendererParallelJob(settings)
    job.start()
    job.waitForFinished()
    return job.renderedImage()


def _rgba_tif(arr_mm, ramp, gt, wkt, path, upsample=4):
    """Colour-map a thickness grid. The grid is first resampled (bilinear on the
    log of the thickness) to a finer cell so that the frames read as smooth
    blobs rather than 30 m squares; the values are not changed otherwise."""
    drv = gdal.GetDriverByName("MEM")
    a = np.asarray(arr_mm, dtype=np.float32)
    if upsample and upsample > 1 and a.size:
        la = np.log10(np.maximum(a, 1e-4)).astype(np.float32)
        src = drv.Create("", a.shape[1], a.shape[0], 1, gdal.GDT_Float32)
        src.SetGeoTransform(gt)
        src.SetProjection(wkt)
        src.GetRasterBand(1).WriteArray(la)
        dst = gdal.Warp("", src, format="MEM", xRes=abs(gt[1]) / upsample, yRes=abs(gt[5]) / upsample,
                        resampleAlg="bilinear")
        near = gdal.Warp("", src, format="MEM", xRes=abs(gt[1]) / upsample, yRes=abs(gt[5]) / upsample,
                         resampleAlg="near")
        # bilinear for smooth edges, nearest underneath so thin threads stay continuous
        a = 10 ** np.maximum(dst.GetRasterBand(1).ReadAsArray(), near.GetRasterBand(1).ReadAsArray() - 0.3)
        a[a < 1.5e-4] = 0
        gt = dst.GetGeoTransform()
    rgba = outputs.ramp_rgba(a, ramp)
    ds = gdal.GetDriverByName("GTiff").Create(path, a.shape[1], a.shape[0], 4, gdal.GDT_Byte)
    ds.SetGeoTransform(gt)
    ds.SetProjection(wkt)
    for b in range(4):
        ds.GetRasterBand(b + 1).WriteArray(rgba[..., b])
    ds.GetRasterBand(4).SetColorInterpretation(gdal.GCI_AlphaBand)
    ds = None


def _overlay_layer(path, name):
    lyr = QgsRasterLayer(path, name)
    r = QgsMultiBandColorRenderer(lyr.dataProvider(), 1, 2, 3)
    r.setAlphaBand(4)
    lyr.setRenderer(r)
    return lyr


def _nice(v):
    e = 10 ** math.floor(math.log10(v))
    for m in (1, 2, 2.5, 5, 10):
        if m * e >= v:
            return m * e
    return 10 * e


def _with_run_drainage(run_dir, prefix, layers):
    """Use the drainage network traced for this run (the one the oil follows) in
    place of any project-wide drainage layer, so that the map and the oil agree."""
    gp = os.path.join(run_dir, prefix + "_results.gpkg")
    if not os.path.exists(gp):
        return layers
    lyr = QgsVectorLayer("%s|layername=drainage_network_dem" % gp, "Drainage network (this run)", "ogr")
    if not lyr.isValid() or lyr.featureCount() == 0:
        return layers
    rend = QgsRuleBasedRenderer(QgsLineSymbol.createSimple({"line_color": "#3b8fd9", "line_width": "0.3"}))
    root = rend.rootRule()
    root.removeChildAt(0)
    for expr, w in (('"area_km2" >= 100', "0.9"), ('"area_km2" >= 10 AND "area_km2" < 100', "0.55"),
                    ('"area_km2" >= 3 AND "area_km2" < 10', "0.35"), ('"area_km2" < 3', "0.2")):
        root.appendChild(QgsRuleBasedRenderer.Rule(QgsLineSymbol.createSimple(
            {"line_color": "#3b8fd9", "line_width": w}), filterExp=expr))
    lyr.setRenderer(rend)
    _KEEP.append(lyr)
    out, placed = [], False
    for l in layers:
        if "drainage" in l.name().lower() and l.type() == 0:
            if not placed:
                out.append(lyr)
                placed = True
            continue
        if not placed and l.type() == 1:      # first raster: put the streams above the terrain
            out.append(lyr)
            placed = True
        out.append(l)
    if not placed:
        out.append(lyr)
    return out


_KEEP = []


def _bbox(run_dir, rec, key, thr):
    """Bounding box (raster CRS) of cells that ever exceed thr in a time-series GeoTIFF."""
    if not rec.get(key):
        return None
    ds = gdal.Open(os.path.join(run_dir, rec[key]))
    gt = ds.GetGeoTransform()
    m = None
    for b in range(ds.RasterCount):
        a = ds.GetRasterBand(b + 1).ReadAsArray() >= thr
        m = a if m is None else (m | a)
    if m is None or not m.any():
        return None
    rr, cc = np.nonzero(m)
    x0 = gt[0] + cc.min() * gt[1]
    x1 = gt[0] + (cc.max() + 1) * gt[1]
    y1 = gt[3] + rr.min() * gt[5]
    y0 = gt[3] + (rr.max() + 1) * gt[5]
    return QgsRectangle(x0, y0, x1, y1)


def _fit(rect, aspect, margin, min_half):
    w, h = rect.width() * (1 + 2 * margin), rect.height() * (1 + 2 * margin)
    w, h = max(w, 2 * min_half), max(h, 2 * min_half)
    if w / h < aspect:
        w = h * aspect
    else:
        h = w / aspect
    c = rect.center()
    return QgsRectangle(c.x() - w / 2, c.y() - h / 2, c.x() + w / 2, c.y() + h / 2)


def default_extents(run_dir, rec, dest_crs, aspect):
    """Main map extent (oil + spray fallout) and an inset extent (oil pool only)."""
    src = QgsCoordinateReferenceSystem()
    src.createFromWkt(rec["crs_wkt"])
    sx, sy = rec["source_xy"]
    srect = QgsRectangle(sx - 1, sy - 1, sx + 1, sy + 1)
    h2s = rec.get("model") == "h2s"
    oil = _bbox(run_dir, rec, "timeseries", 5.0 if h2s else 0.1) or srect
    oil.combineExtentWith(srect)
    allr = QgsRectangle(oil)
    fall = _bbox(run_dir, rec, "fallout_timeseries", 0.02)
    if fall:
        allr.combineExtentWith(fall)
    tr = None
    if dest_crs is not None and dest_crs.isValid() and dest_crs != src:
        tr = QgsCoordinateTransform(src, dest_crs, QgsProject.instance())
        allr = tr.transformBoundingBox(allr)
        oil = tr.transformBoundingBox(oil)
    main = _fit(allr, aspect, 0.15, 1500.0 if h2s else 800.0)
    inset = None
    if not h2s and oil.width() < 0.3 * main.width() and oil.height() < 0.3 * main.height():
        inset = _fit(oil, 1.0, 0.35, 250.0)
    return main, inset


def render_animation(run_json, layers, out_dir=None, width=1600, height=1000, extent=None,
                     dest_crs=None, fps=4, make_gif=True, make_mp4=True, log=print, feedback=None,
                     title=None, subtitle=None, prepared_by=None):
    """layers: list of QgsMapLayer, top first (as in the Layers panel)."""
    run_dir = os.path.dirname(os.path.abspath(run_json))
    with open(run_json) as f:
        rec = json.load(f)
    prefix = os.path.basename(run_json).replace("_run.json", "")
    out_dir = out_dir or os.path.join(run_dir, prefix + "_animation")
    os.makedirs(out_dir, exist_ok=True)
    if dest_crs is None or not dest_crs.isValid():
        dest_crs = QgsProject.instance().crs()
    if not dest_crs.isValid():
        dest_crs = QgsCoordinateReferenceSystem()
        dest_crs.createFromWkt(rec["crs_wkt"])
    layers = _with_run_drainage(run_dir, prefix, [l for l in layers if l is not None and l.isValid()])
    map_h = height - 190          # title block on top (90 px), info strip at the bottom (100 px)
    inset_ext = None
    if extent is None:
        extent, inset_ext = default_extents(run_dir, rec, dest_crs, width / float(map_h))
    # ---- base map, rendered once
    ms = QgsMapSettings()
    ms.setDestinationCrs(dest_crs)
    ms.setOutputSize(QSize(width, map_h))
    ms.setExtent(extent)
    ms.setOutputDpi(96)
    ms.setBackgroundColor(QColor(236, 240, 232))
    ms.setFlag(QgsMapSettings.Antialiasing, True)
    ms.setLayers([l for l in layers if l is not None and l.isValid()])
    try:
        ms.setTransformContext(QgsProject.instance().transformContext())
    except Exception:
        pass
    base = _render(ms)
    ext = ms.visibleExtent()
    ins = None
    isz = 0
    if inset_ext is not None:
        isz = int(min(map_h * 0.46, width * 0.28))
        ins = QgsMapSettings(ms)
        ins.setOutputSize(QSize(isz, isz))
        ins.setExtent(inset_ext)
        ins_base = _render(ins)
        inset_ext = ins.visibleExtent()
    # spill source in the map CRS
    _src_crs = QgsCoordinateReferenceSystem()
    _src_crs.createFromWkt(rec["crs_wkt"])
    src_pt = QgsPointXY(*rec["source_xy"])
    if _src_crs != dest_crs:
        src_pt = QgsCoordinateTransform(_src_crs, dest_crs, QgsProject.instance()).transform(src_pt)
    gdal.PushErrorHandler("CPLQuietErrorHandler")
    # ---- frame stacks
    ts = gdal.Open(os.path.join(run_dir, rec["timeseries"]))
    gt = ts.GetGeoTransform()
    wkt = ts.GetProjection()
    fts = None
    if rec.get("fallout_timeseries"):
        fts = gdal.Open(os.path.join(run_dir, rec["fallout_timeseries"]))
    tmp = tempfile.mkdtemp(prefix="spillframes_")
    start = None
    if rec.get("start_iso"):
        try:
            start = _dt.datetime.fromisoformat(rec["start_iso"])
        except Exception:
            start = None
    h2s = rec.get("model") == "h2s"
    ramp = [tuple(r) for r in rec["ramp"]] if h2s else outputs.OIL_RAMP
    mb = rec["stats"] if h2s else rec["mass_balance"]
    frames = rec["frames"]
    events = _events(rec)
    paths = []
    n = ts.RasterCount
    for i in range(n):
        if feedback is not None:
            if feedback.isCanceled():
                break
            feedback.setProgress(100.0 * i / max(n, 1))
        oil = ts.GetRasterBand(i + 1).ReadAsArray()
        lay = []
        p1 = os.path.join(tmp, "oil_%03d.tif" % i)
        _rgba_tif(oil, ramp, gt, wkt, p1)
        lay.append(_overlay_layer(p1, "oil"))
        if fts is not None:
            fall = fts.GetRasterBand(i + 1).ReadAsArray()
            p2 = os.path.join(tmp, "fall_%03d.tif" % i)
            _rgba_tif(fall, outputs.FALL_RAMP, gt, wkt, p2)
            lay.append(_overlay_layer(p2, "fall"))
        ov = QgsMapSettings(ms)
        ov.setLayers(lay)
        ov.setBackgroundColor(QColor(0, 0, 0, 0))
        top = _render(ov)
        if ins is not None:
            iov = QgsMapSettings(ins)
            iov.setLayers(lay)
            iov.setBackgroundColor(QColor(0, 0, 0, 0))
            itop = _render(iov)
        img = QImage(width, height, QImage.Format_ARGB32_Premultiplied)
        img.fill(QColor(255, 255, 255))
        p = QPainter(img)
        p.setRenderHint(QPainter.Antialiasing, True)
        p.setRenderHint(QPainter.TextAntialiasing, True)
        p.drawImage(0, 90, base)
        p.drawImage(0, 90, top)
        fr = frames[min(i, len(frames) - 1)]
        row = mb[min(i, len(mb) - 1)]
        if fr.get("rain_mmh", 0) > 0:
            _rain(p, width, map_h, fr["rain_mmh"])
        _source_marker(p, ms, 0, 90, src_pt, rec.get("source_label") or
                       ("Leak point" if h2s else "Spill source"), True)
        if ins is not None:
            _inset(p, width, isz, ins_base, itop, ext, inset_ext, map_h)
            _source_marker(p, ins, width - isz - 16, 104, src_pt, None, False)
        _decorate(p, width, height, map_h, rec, fr, row, start, ext, title, subtitle,
                  prepared_by, events)
        caps = rec.get("captions") or []
        if i < len(caps) and caps[i]:
            _caption(p, width, map_h, caps[i])
        p.end()
        fp = os.path.join(out_dir, "%s_frame_%03d.png" % (prefix, i))
        img.save(fp)
        paths.append(fp)
    gdal.PopErrorHandler()
    shutil.rmtree(tmp, ignore_errors=True)
    log("%d frames drawn to %s" % (len(paths), out_dir))
    out = dict(frames=paths)
    if make_gif and paths:
        g = _gif(paths, os.path.join(run_dir, prefix + "_animation.gif"), fps, log)
        if g:
            out["gif"] = g
    if make_mp4 and paths:
        m = _mp4(out_dir, prefix, os.path.join(run_dir, prefix + "_animation.mp4"), fps, log)
        if m:
            out["mp4"] = m
    return out


def _source_marker(p, settings, ox, oy, pt, label, big):
    q = settings.mapToPixel().transform(pt)
    x, y = ox + q.x(), oy + q.y()
    r = 9 if big else 6
    p.setPen(QPen(QColor(0, 0, 0), 2))
    p.setBrush(QBrush(QColor(255, 212, 0)))
    p.drawEllipse(QPointF(x, y), r, r)
    p.setBrush(QBrush(QColor(0, 0, 0)))
    p.drawEllipse(QPointF(x, y), 2.5, 2.5)
    if label:
        p.setFont(QFont("Segoe UI", 9, QFont.Bold))
        fm = QFontMetrics(p.font())
        tw = fm.horizontalAdvance(label) if hasattr(fm, "horizontalAdvance") else fm.width(label)
        rect = QRectF(x - tw - 22, y - 10, tw + 10, 20)
        p.setPen(Qt.NoPen)
        p.setBrush(QBrush(QColor(255, 255, 255, 225)))
        p.drawRoundedRect(rect, 4, 4)
        p.setPen(QColor(20, 33, 27))
        p.drawText(rect, Qt.AlignCenter, label)


def _inset(p, W, isz, base, top, ext, iext, map_h):
    x0, y0 = W - isz - 16, 90 + 14
    # footprint of the inset on the main map
    sx = W / ext.width()
    rx = (iext.xMinimum() - ext.xMinimum()) * sx
    ry = 90 + (ext.yMaximum() - iext.yMaximum()) * sx
    rw = iext.width() * sx
    p.setBrush(Qt.NoBrush)
    p.setPen(QPen(QColor(20, 33, 27), 1.4, Qt.DashLine))
    p.drawRect(QRectF(rx, ry, rw, rw))
    p.setPen(QPen(QColor(20, 33, 27, 120), 1, Qt.DotLine))
    p.drawLine(QPointF(rx + rw, ry), QPointF(x0, y0))
    p.drawLine(QPointF(rx + rw, ry + rw), QPointF(x0, y0 + isz))
    p.setPen(Qt.NoPen)
    p.setBrush(QBrush(QColor(0, 0, 0, 60)))
    p.drawRect(QRectF(x0 + 3, y0 + 3, isz, isz))
    p.drawImage(int(x0), int(y0), base)
    p.drawImage(int(x0), int(y0), top)
    p.setBrush(Qt.NoBrush)
    p.setPen(QPen(QColor(255, 255, 255), 3))
    p.drawRect(QRectF(x0, y0, isz, isz))
    p.setPen(QPen(QColor(20, 33, 27), 1))
    p.drawRect(QRectF(x0, y0, isz, isz))
    # label and scale
    p.fillRect(QRectF(x0, y0 + isz - 22, isz, 22), QColor(255, 255, 255, 215))
    p.setPen(QColor(20, 33, 27))
    p.setFont(QFont("Segoe UI", 8))
    target = _nice(iext.width() * 0.3)
    lab = ("%g km" % (target / 1000.0)) if target >= 1000 else ("%g m" % target)
    p.drawText(QRectF(x0 + 8, y0 + isz - 22, isz - 16, 22), Qt.AlignLeft | Qt.AlignVCenter,
               "Zoom on the spill source")
    bw = target / iext.width() * isz
    bx = x0 + isz - 12 - bw
    p.fillRect(QRectF(bx, y0 + isz - 12, bw, 4), QColor(20, 33, 27))
    p.drawText(QRectF(bx - 60, y0 + isz - 24, bw + 60, 12), Qt.AlignRight | Qt.AlignVCenter, lab)


def _events(rec):
    ev = []
    if rec["kind"] in ("BLOWOUT", "H2S") and rec.get("duration_h") and rec["duration_h"] < rec.get("sim_hours", 1e9):
        ev.append((rec["duration_h"], "Well under control - spillage stops"))
    for st, du, dep in rec.get("storms", []):
        ev.append((st, "Storm: %.0f mm in %.0f h" % (dep, du)))
    s = rec.get("summary", {})
    fc = s.get("first_channel")
    if fc:
        ev.append((fc["t_h"], "Oil enters the drainage network"))
    if s.get("first_export_h") is not None:
        ev.append((s["first_export_h"], "Oil leaves the model area downstream"))
    for r in (s.get("receptors_reached") or [])[:6]:
        ev.append((r["arrival_h"], "Oil reaches %s" % r["name"]))
    return sorted(ev)


def _rain(p, w, h, mmh):
    a = int(min(150, 40 + mmh * 4))
    pen = QPen(QColor(90, 130, 230, a))
    pen.setWidthF(1.2)
    p.setPen(pen)
    rng = np.random.RandomState(int(mmh * 100) % 9973)
    nlines = int(min(900, 150 + mmh * 25))
    for _ in range(nlines):
        x = rng.uniform(0, w)
        y = rng.uniform(90, 90 + h)
        p.drawLine(QPointF(x, y), QPointF(x - 6, y + 16))
    p.fillRect(QRectF(0, 90, w, h), QColor(60, 90, 160, int(min(60, 10 + mmh * 2))))


def _decorate(p, W, H, map_h, rec, fr, row, start, ext, title, subtitle, prepared_by, events):
    green = QColor(11, 93, 59)
    gold = QColor(201, 162, 39)
    ink = QColor(20, 33, 27)
    # title block
    p.fillRect(QRectF(0, 0, W, 90), green)
    p.fillRect(QRectF(0, 86, W, 4), gold)
    p.setPen(QColor(255, 255, 255))
    f = QFont("Segoe UI", 20)
    f.setBold(True)
    p.setFont(f)
    p.drawText(QRectF(24, 10, W * 0.66, 40), Qt.AlignLeft | Qt.AlignVCenter, title or rec.get("title") or rec["name"])
    p.setFont(QFont("Segoe UI", 11))
    p.drawText(QRectF(24, 50, W * 0.66, 28), Qt.AlignLeft | Qt.AlignVCenter,
               subtitle or "  |  ".join(x for x in (rec.get("subtitle", ""), rec.get("product", "")) if x))
    # clock
    t = fr["t_h"]
    hh = int(t + 1e-6)
    mm = int(round((t - hh) * 60)) % 60
    f = QFont("Segoe UI", 26)
    f.setBold(True)
    p.setFont(f)
    p.drawText(QRectF(W - 360, 6, 336, 46), Qt.AlignRight | Qt.AlignVCenter, "T+%02d:%02d" % (hh, mm))
    p.setFont(QFont("Segoe UI", 11))
    if start:
        d = start + _dt.timedelta(hours=t)
        p.drawText(QRectF(W - 420, 52, 396, 26), Qt.AlignRight | Qt.AlignVCenter, d.strftime("%a %d %b %Y, %H:%M"))
    y0 = 90 + map_h
    # bottom strip
    p.fillRect(QRectF(0, y0, W, H - y0), QColor(248, 250, 248))
    p.setPen(QColor(217, 224, 220))
    p.drawLine(0, y0, W, y0)
    p.setPen(ink)
    h2s = rec.get("model") == "h2s"
    if h2s:
        kp = [("Highest at ground now", "%.1f ppm" % row["max_ppm"]),
              ("AEGL-2 (27 ppm) reach", "%.0f m" % row.get("AEGL2", 0)),
              ("AEGL-1 (0.51 ppm) reach", "%.0f m" % row.get("AEGL1", 0)),
              ("IDLH (100 ppm) reach", "%.0f m" % row.get("IDLH", 0)),
              ("Stability class", row.get("stability", "-")),
              ("Mixing height", "%.0f m" % row.get("mixing_m", 0))]
    else:
        kp = [("Released", "%.1f m3" % row["released"]), ("Oiled ground now", "%.2f ha" % row["area_ha"]),
              ("Evaporated", "%.1f m3" % row["evaporated"]), ("Soaked into soil", "%.1f m3" % row["infiltrated"]),
              ("Stranded on banks", "%.1f m3" % row["stranded"]),
              ("Left model area", "%.1f m3" % (row["exported"] + row["spray_out"]))]
    x = 24
    for k, v in kp:
        p.setFont(QFont("Segoe UI", 9))
        p.setPen(QColor(91, 107, 99))
        p.drawText(QRectF(x, y0 + 10, 170, 18), Qt.AlignLeft, k)
        f = QFont("Segoe UI", 14)
        f.setBold(True)
        p.setFont(f)
        p.setPen(ink)
        p.drawText(QRectF(x, y0 + 28, 170, 26), Qt.AlignLeft, v)
        x += 165
    # legend
    if h2s:
        _ramp_legend(p, W - 290, y0 + 10, 250, "H2S at ground level, 10-min mean (ppm)",
                     [tuple(r) for r in rec["ramp"]])
    elif rec["kind"] == "BLOWOUT":
        _ramp_legend(p, W - 590, y0 + 10, 240, "Surface oil (mm)", outputs.OIL_RAMP)
        _ramp_legend(p, W - 300, y0 + 10, 250, "Spray fallout, cumulative (mm)", outputs.FALL_RAMP)
    else:
        _ramp_legend(p, W - 290, y0 + 10, 250, "Surface oil (mm)", outputs.OIL_RAMP)
    # footer
    p.setFont(QFont("Segoe UI", 8))
    p.setPen(QColor(91, 107, 99))
    foot = "Prepared by %s" % (prepared_by or rec.get("prepared_by") or "-")
    p.drawText(QRectF(24, H - 22, W - 48, 18), Qt.AlignLeft | Qt.AlignVCenter, foot)
    p.drawText(QRectF(24, H - 22, W - 48, 18), Qt.AlignRight | Qt.AlignVCenter,
               "Screening simulation for emergency planning - not a forecast")
    # wind + rain box (top-left of the map)
    bx, by = 16, 104
    p.setBrush(QBrush(QColor(255, 255, 255, 225)))
    p.setPen(QColor(200, 208, 204))
    p.drawRoundedRect(QRectF(bx, by, 250, 92), 8, 8)
    w = fr.get("wind")
    cx, cy = bx + 46, by + 46
    p.setPen(QPen(QColor(150, 160, 155), 1))
    p.setBrush(QBrush(QColor(243, 246, 244)))
    p.drawEllipse(QPointF(cx, cy), 34, 34)
    p.setFont(QFont("Segoe UI", 7))
    for k, lab in enumerate("NESW"):
        a = math.radians(k * 90)
        p.drawText(QRectF(cx + 27 * math.sin(a) - 6, cy - 27 * math.cos(a) - 6, 12, 12), Qt.AlignCenter, lab)
    txt = "Calm"
    if w:
        to = math.radians(w[0])
        dx, dy = math.sin(to), -math.cos(to)
        p.setPen(QPen(green, 4, Qt.SolidLine, Qt.RoundCap))
        p.drawLine(QPointF(cx - 18 * dx, cy - 18 * dy), QPointF(cx + 16 * dx, cy + 16 * dy))
        tip = QPointF(cx + 26 * dx, cy + 26 * dy)
        lft = QPointF(cx + 12 * dx - 8 * dy, cy + 12 * dy + 8 * dx)
        rgt = QPointF(cx + 12 * dx + 8 * dy, cy + 12 * dy - 8 * dx)
        p.setBrush(QBrush(green))
        p.setPen(Qt.NoPen)
        p.drawPolygon(QPolygonF([tip, lft, rgt]))
        frm = (w[0] + 180) % 360
        txt = "Wind from %s, %.1f m/s" % (DIRS[int(round(frm / 22.5)) % 16], w[1])
    p.setPen(ink)
    p.setFont(QFont("Segoe UI", 9))
    p.drawText(QRectF(bx + 90, by + 10, 156, 22), Qt.AlignLeft | Qt.AlignVCenter, txt)
    rain = fr.get("rain_mmh", 0) or 0
    if h2s:
        stab = row.get("stability", "D")
        p.drawText(QRectF(bx + 90, by + 34, 156, 22), Qt.AlignLeft | Qt.AlignVCenter,
                   "Stability %s - %s" % (stab, STAB.get(stab, "")))
        p.drawText(QRectF(bx + 90, by + 58, 156, 22), Qt.AlignLeft | Qt.AlignVCenter,
                   "Plume height %.0f m" % row.get("plume_h_m", 0))
    else:
        p.drawText(QRectF(bx + 90, by + 34, 156, 22), Qt.AlignLeft | Qt.AlignVCenter,
                   ("Rain %.1f mm/h" % rain) if rain > 0 else "No rain")
        p.drawText(QRectF(bx + 90, by + 58, 156, 22), Qt.AlignLeft | Qt.AlignVCenter,
                   "Streams %s" % rec.get("stream_state", ""))
    # event banner
    fh = rec["frame_minutes"] / 60.0
    ev = [e for e in events if e[0] <= fr["t_h"] + 1e-6 and fr["t_h"] - e[0] < max(3 * fh, 1.0)]
    if ev:
        p.setFont(QFont("Segoe UI", 12, QFont.Bold))
        msg = ev[-1][1]
        fm = QFontMetrics(p.font())
        tw = fm.horizontalAdvance(msg) if hasattr(fm, "horizontalAdvance") else fm.width(msg)
        r = QRectF((W - tw) / 2 - 16, 104, tw + 32, 32)
        p.setBrush(QBrush(QColor(20, 33, 27, 215)))
        p.setPen(Qt.NoPen)
        p.drawRoundedRect(r, 16, 16)
        p.setPen(QColor(255, 255, 255))
        p.drawText(r, Qt.AlignCenter, msg)
    # scale bar and north arrow (bottom-left of map)
    mupp = ext.width() / float(W)
    target = _nice(mupp * 160)
    px = target / mupp
    sx, sy = 70, 90 + map_h - 30
    p.setBrush(QBrush(QColor(255, 255, 255, 220)))
    p.setPen(Qt.NoPen)
    p.drawRoundedRect(QRectF(sx - 12, sy - 26, px + 24, 42), 6, 6)
    p.setPen(QPen(ink, 1))
    p.setBrush(QBrush(ink))
    p.drawRect(QRectF(sx, sy, px / 2, 6))
    p.setBrush(QBrush(QColor(255, 255, 255)))
    p.drawRect(QRectF(sx + px / 2, sy, px / 2, 6))
    p.setFont(QFont("Segoe UI", 8))
    lab = ("%g km" % (target / 1000.0)) if target >= 1000 else ("%g m" % target)
    p.drawText(QRectF(sx, sy - 20, px, 16), Qt.AlignCenter, lab)
    nx, ny = 34, 90 + map_h - 34
    p.setBrush(QBrush(QColor(255, 255, 255, 220)))
    p.setPen(Qt.NoPen)
    p.drawRoundedRect(QRectF(nx - 16, ny - 28, 32, 54), 6, 6)
    p.setBrush(QBrush(ink))
    p.setPen(Qt.NoPen)
    p.drawPolygon(QPolygonF([QPointF(nx, ny - 22), QPointF(nx - 9, ny + 8), QPointF(nx, ny + 2),
                             QPointF(nx + 9, ny + 8)]))
    p.setPen(ink)
    p.setFont(QFont("Segoe UI", 9, QFont.Bold))
    p.drawText(QRectF(nx - 10, ny + 8, 20, 16), Qt.AlignCenter, "N")


def _caption(p, W, map_h, cap):
    """Subtitle box at the bottom of the map: headline and what is happening in this frame."""
    x0 = 345
    w = min(W - x0 - 24, 980)
    pad = 12
    ft = QFont("Segoe UI", 12, QFont.Bold)
    fb = QFont("Segoe UI", 11)
    fmb = QFontMetrics(fb)
    txt = cap.get("text") or ""
    br = fmb.boundingRect(QRect(0, 0, int(w - 2 * pad), 400), int(Qt.TextWordWrap), txt)
    h = 22 + br.height() + 2 * pad
    y0 = 90 + map_h - h - 14
    p.setPen(Qt.NoPen)
    p.setBrush(QBrush(QColor(16, 28, 22, 215)))
    p.drawRoundedRect(QRectF(x0, y0, w, h), 8, 8)
    p.setFont(ft)
    p.setPen(QColor(255, 212, 59))
    p.drawText(QRectF(x0 + pad, y0 + pad - 2, w - 2 * pad, 22), Qt.AlignLeft | Qt.AlignVCenter, cap.get("title") or "")
    p.setFont(fb)
    p.setPen(QColor(255, 255, 255))
    p.drawText(QRectF(x0 + pad, y0 + pad + 22, w - 2 * pad, br.height() + 4), int(Qt.TextWordWrap | Qt.AlignLeft), txt)


def _ramp_legend(p, x, y, w, label, ramp):
    p.setPen(QColor(91, 107, 99))
    p.setFont(QFont("Segoe UI", 9))
    p.drawText(QRectF(x, y, w, 16), Qt.AlignLeft, label)
    g = QLinearGradient(x, 0, x + w, 0)
    n = len(ramp)
    for i, s in enumerate(ramp):
        g.setColorAt(i / float(n - 1), QColor(s[1], s[2], s[3], max(s[4], 150)))
    p.setPen(Qt.NoPen)
    p.setBrush(QBrush(g))
    p.drawRoundedRect(QRectF(x, y + 20, w, 12), 3, 3)
    p.setPen(QColor(91, 107, 99))
    p.setFont(QFont("Segoe UI", 7))
    for i, s in enumerate(ramp):
        if i % 2 and i != n - 1:
            continue
        xx = x + w * i / float(n - 1)
        p.drawText(QRectF(xx - 20, y + 34, 40, 12), Qt.AlignCenter, "%g" % s[0])


def _gif(paths, out, fps, log):
    try:
        from PIL import Image
    except ImportError:
        log("Pillow not available - GIF skipped (PNG frames kept).")
        return None
    ims = []
    step = max(1, int(math.ceil(len(paths) / 150.0)))   # keep GIFs to about 150 frames
    if step > 1:
        paths = paths[::step] + ([paths[-1]] if (len(paths) - 1) % step else [])
        fps = max(1, int(round(fps / float(step) * 2)))
    for pth in paths:
        im = Image.open(pth).convert("RGB")
        if im.width > 1280:
            im = im.resize((1280, int(im.height * 1280.0 / im.width)), Image.LANCZOS)
        ims.append(im.convert("P", palette=Image.ADAPTIVE, colors=255))
    d = int(1000 / max(fps, 1))
    durations = [d] * len(ims)
    durations[-1] = d * 6
    ims[0].save(out, save_all=True, append_images=ims[1:], duration=durations, loop=0, optimize=True)
    log("GIF written: %s" % out)
    return out


def _ffmpeg():
    exe = shutil.which("ffmpeg")
    if exe:
        return exe
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except Exception:
        return None


def _mp4(frame_dir, prefix, out, fps, log):
    exe = _ffmpeg()
    if not exe:
        log("ffmpeg not found - MP4 skipped (GIF and PNG frames written). Install ffmpeg or "
            "imageio-ffmpeg to get MP4 files.")
        return None
    cmd = [exe, "-y", "-loglevel", "error", "-framerate", str(fps), "-i",
           os.path.join(frame_dir, prefix + "_frame_%03d.png"), "-c:v", "libx264", "-pix_fmt", "yuv420p",
           "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-movflags", "+faststart", out]
    try:
        kw = {}
        if os.name == "nt":
            kw["creationflags"] = 0x08000000
        subprocess.check_call(cmd, **kw)
        log("MP4 written: %s" % out)
        return out
    except Exception as e:
        log("MP4 encoding failed: %s" % e)
        return None
