# -*- coding: utf-8 -*-
"""
Plain-language captions (subtitles) for every animation frame.

The captions are written from the model results only: the clock and date, the release state,
the wind, rain and stability of the frame, the oil temperature and gelled volume (heat mode),
the extent of the oil or gas, and the settlements and other receptors nearest to it.
They are stored in the run record (<prefix>_run.json, key "captions") and shown in the QGIS
simulation panel, the GIF / MP4 animations and the HTML page.
"""
import datetime as _dt
import math

import numpy as np
from osgeo import gdal

from .climate import compass

SETTLEMENT_WORDS = ("town", "village", "settlement", "communit", "compound")
STAB_TEXT = {"A": "very unstable air (strong sun) - the gas mixes upwards quickly",
             "B": "unstable air (sunny) - the gas mixes upwards quickly",
             "C": "slightly unstable air - good mixing",
             "D": "neutral air (cloud or stronger wind) - moderate mixing",
             "E": "stable night air - the gas stays near the ground",
             "F": "very stable night air with light wind - the gas hugs the ground and travels furthest"}


def is_settlement(rtype):
    t = (rtype or "").lower()
    return any(w in t for w in SETTLEMENT_WORDS)


def part_of_day(h):
    if h < 5 or h >= 20:
        return "Night"
    if h < 9:
        return "Early morning"
    if h < 12:
        return "Morning"
    if h < 15:
        return "Midday"
    if h < 18:
        return "Afternoon"
    return "Evening"


def _clock(start, t):
    if start is None:
        return None, None
    d = start + _dt.timedelta(hours=t)
    return d, d.strftime("%a %d %b %Y %H:%M")


def _dist_text(m):
    return ("%.1f km" % (m / 1000.0)) if m >= 1000 else ("%.0f m" % (round(m / 10.0) * 10))


def _cells_xy(mask, gt):
    r, c = np.nonzero(mask)
    x = gt[0] + (c + 0.5) * gt[1]
    y = gt[3] + (r + 0.5) * gt[5]
    return x, y


def _nearest(px, py, cx, cy):
    """Distance from a point to the nearest of the cells (brute force on a subsample)."""
    if cx.size == 0:
        return None
    if cx.size > 20000:
        k = np.linspace(0, cx.size - 1, 20000).astype(int)
        cx, cy = cx[k], cy[k]
    return float(np.sqrt(((cx - px) ** 2 + (cy - py) ** 2).min()))


def build(rec, ts_path, receptors, threshold=None, log=print):
    """rec: run record dict (oil or H2S); ts_path: GeoTIFF time series (one band per frame);
    receptors: list of dicts with name, type, x, y (+ arrival_h for oil, series for H2S).
    Returns one caption dict per frame: {t_h, clock, title, text}."""
    h2s = rec.get("model") == "h2s"
    start = None
    if rec.get("start_iso"):
        try:
            start = _dt.datetime.fromisoformat(rec["start_iso"])
        except Exception:
            start = None
    ds = gdal.Open(ts_path)
    if ds is None:
        return []
    gt = ds.GetGeoTransform()
    sx, sy = rec["source_xy"]
    thr = threshold if threshold is not None else (0.51 if h2s else 0.1)
    rows = rec["stats"] if h2s else rec["mass_balance"]
    frames = rec["frames"]
    settle = [r for r in receptors if is_settlement(r.get("type"))]
    others = [r for r in receptors if not is_settlement(r.get("type"))]
    storms = rec.get("storms") or []
    dur = rec.get("duration_h")
    kind = rec.get("kind")
    out = []
    reached_before = set()
    for i in range(ds.RasterCount):
        fr = frames[min(i, len(frames) - 1)]
        row = rows[min(i, len(rows) - 1)]
        t = fr["t_h"]
        d, clock = _clock(start, t)
        pod = part_of_day(d.hour + d.minute / 60.0) if d else None
        a = ds.GetRasterBand(i + 1).ReadAsArray()
        mask = a >= thr
        cx, cy = _cells_xy(mask, gt)
        parts = []
        title = None
        # ---- release state
        if h2s:
            rf = row.get("releasing", 1.0)
            if 0 < rf < 0.999:
                parts.append("H2S is leaking from the well; the leak is building up (%.0f%% of the full rate)." % (rf * 100))
            elif rf > 0:
                parts.append("H2S is leaking from the well at the full rate.")
            else:
                parts.append("The leak has stopped; the remaining gas is clearing.")
        elif kind == "BLOWOUT":
            if dur is not None and t < dur:
                parts.append("The well is flowing: %.0f m3 of oil released so far." % row["released"]
                             if row["released"] >= 0.5 else "The well has just started to flow.")
            else:
                parts.append("The well was brought under control at T+%.0f h; no further oil is released." % dur)
        else:
            parts.append("%.1f m3 of %s has escaped." % (row["released"], (rec.get("product") or "product").lower()))
        # ---- weather
        w = fr.get("wind")
        if w:
            frm = (w[0] + 180.0) % 360.0
            wtxt = "Wind from the %s at %.1f m/s" % (compass(frm), w[1])
        else:
            wtxt = "Calm"
        rain = fr.get("rain_mmh") or 0.0
        storm_now = [s for s in storms if s[0] <= t < s[0] + s[1] + 1e-6]
        if h2s:
            st = row.get("stability") or fr.get("stab") or "D"
            parts.append("%s; %s." % (wtxt, STAB_TEXT.get(st, "")))
        elif rain > 0:
            parts.append("%s. Storm in progress (%.0f mm/h): run-off is washing oil downslope." % (wtxt, rain))
        else:
            spray = ""
            if kind == "BLOWOUT" and dur and t < dur and w:
                spray = "; the oil spray drifts to the %s" % compass(w[0])
            parts.append("%s%s." % (wtxt, spray))
        # ---- heat (oil)
        if not h2s and row.get("oil_t_mean") is not None and row.get("surface", 0) > 0.05:
            ht = "Oil on the ground averages %.0f C (air %.0f C" % (row["oil_t_mean"], row.get("t_air", 0))
            if row.get("solar_wm2", 0) > 20:
                ht += ", sunshine %.0f W/m2" % row["solar_wm2"]
            ht += ")."
            if row.get("gelled_m3", 0) >= 0.5:
                ht += " %.0f m3 has cooled below its pour point and gelled - it stops moving." % row["gelled_m3"]
            parts.append(ht)
        # ---- extent
        if cx.size:
            dd = np.hypot(cx - sx, cy - sy)
            k = int(np.argmax(dd))
            reach = float(dd[k])
            brg = math.degrees(math.atan2(cx[k] - sx, cy[k] - sy)) % 360
            if h2s:
                z = []
                if row.get("AEGL1", 0) > 0:
                    z.append("AEGL-1 %s" % _dist_text(row["AEGL1"]))
                if row.get("AEGL2", 0) > 0:
                    z.append("AEGL-2 %s" % _dist_text(row["AEGL2"]))
                if row.get("IDLH", 0) > 0:
                    z.append("IDLH (100 ppm, 10-minute) %s" % _dist_text(row["IDLH"]))
                parts.append("Gas at 0.51 ppm or more (10-minute mean) reaches %s to the %s; %s." % (
                    _dist_text(reach), compass(brg),
                    ("hazard zones now (AEGL on the 1-hour mean): " + ", ".join(z)) if z else
                    "it stays below the 1-hour AEGL levels and IDLH"))
            else:
                parts.append("Oiled ground %.2f ha; the furthest oil is %s from the source to the %s." % (
                    row.get("area_ha", 0.0), _dist_text(reach), compass(brg)))
        else:
            parts.append("No gas above the odour level at ground level." if h2s else "No oil on the ground surface.")
        # ---- settlements and receptors
        if h2s:
            hit = []
            for r in receptors:
                ser = r.get("series") or []
                v = ser[i] if i < len(ser) else 0.0
                if v >= 0.51:
                    lvl = ("IDLH level" if v >= 100 else "AEGL-2 level" if v >= 27 else "AEGL-1 level")
                    hit.append((v, "%s %.0f ppm (%s)" % (r["name"], v, lvl), is_settlement(r.get("type")), r["name"]))
            hit.sort(reverse=True)
            if hit:
                parts.append("In the gas now (10-minute mean): " + "; ".join(h[1] for h in hit[:5]) + ".")
                if any(h[2] for h in hit):
                    title = "Gas over %s" % ", ".join(h[3] for h in hit if h[2])[:80]
        else:
            new = [r for r in receptors if r.get("arrival_h") is not None and r["arrival_h"] <= t + 1e-6
                   and r["name"] not in reached_before]
            for r in new:
                reached_before.add(r["name"])
            if new:
                parts.append("Oil has now reached " + ", ".join("%s (%s)" % (r["name"], r["type"].lower())
                                                                for r in new[:4]) + ".")
                title = "Oil reaches %s" % ", ".join(r["name"] for r in new[:3])
        if settle and cx.size:
            best = None
            for r in settle:
                dm = _nearest(r["x"], r["y"], cx, cy)
                if dm is not None and (best is None or dm < best[0]):
                    best = (dm, r)
            if best is not None:
                dm, r = best
                brg = math.degrees(math.atan2(r["x"] - sx, r["y"] - sy)) % 360
                what = "the gas cloud (0.51 ppm)" if h2s else "the oil"
                if dm <= 30:
                    parts.append("Settlement %s (%s of the source) is inside %s." % (r["name"], compass(brg), what))
                else:
                    parts.append("Nearest settlement: %s, %s from the edge of %s." % (r["name"], _dist_text(dm), what))
        elif settle:
            dists = sorted((math.hypot(r["x"] - sx, r["y"] - sy), r["name"]) for r in settle)
            parts.append("Nearest settlement to the source: %s, %s away." % (dists[0][1], _dist_text(dists[0][0])))
        # ---- headline
        if title is None:
            if storm_now and not h2s:
                title = "Storm run-off"
            elif not h2s and kind == "BLOWOUT" and dur is not None and abs(t - dur) < (rec["frame_minutes"] / 60.0) / 2 + 1e-6:
                title = "Well brought under control"
            elif h2s:
                st = row.get("stability") or "D"
                title = "%s - %s" % (pod or "T+%.1f h" % t, {"A": "strong mixing", "B": "strong mixing",
                                                               "C": "good mixing", "D": "moderate mixing",
                                                               "E": "gas stays low", "F": "gas hugs the ground"}[st])
            elif row.get("gelled_m3", 0) >= 0.5:
                title = "%s - oil cooling and gelling" % (pod or "")
            elif row.get("solar_wm2", 0) > 600:
                title = "%s - sun heating the oil" % (pod or "")
            else:
                title = pod or ("T+%.1f h" % t)
        out.append(dict(t_h=round(t, 3), clock=clock, title=title.strip(" -"), text=" ".join(parts)))
    ds = None
    return out

