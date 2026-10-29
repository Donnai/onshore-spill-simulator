# -*- coding: utf-8 -*-
"""
One-call scenario runner used by the QGIS algorithm and by batch scripts.

    from engine.runner import run_scenario
    out = run_scenario(cfg, log=print)

cfg is a plain dict (see DEFAULTS); every key can be overridden.
"""
import logging
import datetime as _dt
import json
import math
import os

import numpy as np

from .terrain import Terrain, rasterize_layer
from .surfaces import build_surface_classes
from .model import Simulator, Scenario
from .oils import get_oil, OILS, BBL_M3
from .receptors import load_receptors
from . import climate, fallout, outputs, html_export

DEFAULTS = dict(
    name="Scenario",
    well="WELL",
    kind="BLOWOUT",                 # BLOWOUT | TANK | CONTINUOUS
    source_xy=None,
    source_label="Spill source",
    dem=None,
    product="LIGHT_CRUDE",
    product_overrides=None,
    rate_bbl_d=1500.0,
    duration_h=72.0,
    volume_m3=60.0,
    release_minutes=10.0,
    contained_m3=3.0,
    spray_fraction=0.40,
    release_height_m=40.0,
    d50_um=250.0,
    month=0,                        # 0 = January
    start_iso=None,                 # e.g. "2027-01-15T06:00:00"
    sim_hours=96.0,
    frame_minutes=60.0,
    half_width_m=9000.0,
    margin_m=6000.0,
    cell_m=None,
    smooth_sigma=1.0,
    pit_tolerance_m=3.0,
    channel_area_km2=1.0,
    rivers=None,                    # dict(path, layer, where, min_area_km2)
    surface_rules=[],
    receptor_specs=[],
    rain=None,                      # list of (start_h, duration_h, depth_mm); None = none
    cn=78.0,
    wind_csv=None,
    wind_seed=7,
    temp_c=None, q_spec=None, soil_sat=None, dtw_m=None,
    stream_kv=0.70, strand_per_m=1.5e-4, rain_carry=0.5, rain_flush=0.5, s_or=0.10,
    thickness_threshold_mm=0.1,
    rivulet_width_m=6.0,
    snap_to_channel_m=0.0,
    burn_depth_m=3.0,
    receptor_radius_m=60.0,
    out_dir=".",
    prefix=None,
    prepared_by="",
    html=True,
    html_layers=[],                 # list of dicts for html_export
    notes=[],
    # heat transfer / rheology (OilFlow2D-type)
    heat_mode=False,
    discharge_temp_c=None,          # oil temperature at the source; None = air temperature
    cloud=None,                     # cloud / haze fraction; None = monthly default
    rh_pct=None,
    t_max_c=None,
    ground_temp_c=None,
    h_ground=10.0,
    # control features
    roughness_zones=[],             # list of dict(path, layer, n, where) - Manning n override
    barriers=[],                    # list of dict(path, layer, height_m, where) - berms / bund walls
    sinks=[],                       # list of dict(path, layer, where) - sumps, recovery pits, outflow
    extra_sources=[],               # list of dict(x, y, volume_m3, start_h, release_min, temp_c)
)


def _merge(cfg):
    c = dict(DEFAULTS)
    c.update({k: v for k, v in cfg.items() if v is not None or k in ("rain",)})
    return c


def build_scenario(c, log=print):
    clim = climate.month_climate(c["month"])
    if c.get("wind_csv") and os.path.exists(c["wind_csv"]):
        u, v, s = climate.load_wind_csv(c["wind_csv"])
        clim["u"], clim["v"], clim["speed"] = u[c["month"]], v[c["month"]], s[c["month"]]
    start_hour = 0
    doy = 15 + 30 * int(c["month"])
    if c.get("start_iso"):
        d0 = _dt.datetime.fromisoformat(c["start_iso"])
        start_hour = d0.hour
        doy = d0.timetuple().tm_yday
    hours = int(c["sim_hours"]) + 2
    wind = fallout.wind_series(clim["u"], clim["v"], clim["speed"], hours + start_hour,
                               seed=c["wind_seed"])[start_hour:]
    oil = get_oil(c["product"], c.get("product_overrides"))
    sc = Scenario(
        name=c["name"], kind=c["kind"], source_xy=tuple(c["source_xy"]), oil=oil,
        rate_bbl_d=c["rate_bbl_d"], duration_h=c["duration_h"], volume_m3=c["volume_m3"],
        release_minutes=c["release_minutes"], contained_m3=c["contained_m3"],
        spray_fraction=c["spray_fraction"], release_height_m=c["release_height_m"],
        d50_um=c["d50_um"], sim_hours=c["sim_hours"], frame_minutes=c["frame_minutes"],
        month=c["month"], month_name=clim["month"],
        temp_c=c["temp_c"] if c["temp_c"] is not None else clim["temp_c"],
        wind=wind, rain=list(c["rain"] or []), cn=c["cn"],
        q_spec=c["q_spec"] if c["q_spec"] is not None else clim["q_spec"],
        soil_sat=c["soil_sat"] if c["soil_sat"] is not None else clim["soil_sat"],
        dtw_m=c["dtw_m"] if c["dtw_m"] is not None else clim["dtw"],
        stream_kv=c["stream_kv"], strand_per_m=c["strand_per_m"], rain_carry=c["rain_carry"],
        rain_flush=c["rain_flush"], s_or=c["s_or"], rivulet_width_m=c["rivulet_width_m"],
        thickness_threshold_m=c["thickness_threshold_mm"] / 1000.0,
        receptor_radius_m=c["receptor_radius_m"], seed=c["wind_seed"],
        heat_mode=bool(c.get("heat_mode")), discharge_temp_c=c.get("discharge_temp_c"),
        t_max_c=c["t_max_c"] if c.get("t_max_c") is not None else clim["tmax_c"],
        rh_pct=c["rh_pct"] if c.get("rh_pct") is not None else clim["rh"],
        cloud=c["cloud"] if c.get("cloud") is not None else clim["cloud"],
        ground_temp_c=c.get("ground_temp_c"), h_ground=c.get("h_ground", 10.0),
        doy=doy, start_hour=float(start_hour),
        extra_sources=list(c.get("extra_sources") or []))
    if c.get("_latlon"):
        sc.lat, sc.lon = c["_latlon"]
    sc.prepared_by = c.get("prepared_by", "")
    sc.climate = clim
    return sc


def run_scenario(cfg, log=print, feedback=None, terrain=None):
    c = _merge(cfg)
    if not c["dem"] or not os.path.exists(c["dem"]):
        raise IOError("DEM not found: %s" % c["dem"])
    os.makedirs(c["out_dir"], exist_ok=True)
    prefix = c["prefix"] or c["name"].replace(" ", "_")
    sx, sy = c["source_xy"]
    if terrain is None:
        terrain = Terrain(c["dem"], sx, sy, c["half_width_m"], cell_m=c["cell_m"],
                          pit_tolerance_m=c["pit_tolerance_m"],
                          channel_area_km2=c["channel_area_km2"], margin_m=c["margin_m"],
                          smooth_sigma_cells=c["smooth_sigma"],
                          burn=[dict(rv, depth_m=c.get("burn_depth_m", 3.0)) for rv in (c["rivers"] or [])],
                          log=log)
        for rv in (c["rivers"] or []):
            n = terrain.burn_rivers(rv["path"], rv.get("layer"), rv.get("min_area_km2", 300.0),
                                    rv.get("where"))
            log("Mapped watercourses burned into the drainage network: %d cells (%s)." %
                (n, rv.get("layer") or os.path.basename(rv["path"])))
    if c.get("snap_to_channel_m"):
        # move the source onto the nearest mapped/derived watercourse cell
        g = terrain.grid
        r, cc = g.xy_to_rc(sx, sy)
        rad = int(np.ceil(c["snap_to_channel_m"] / g.dx))
        r0, r1 = max(r - rad, 0), min(r + rad + 1, g.nrows)
        c0, c1 = max(cc - rad, 0), min(cc + rad + 1, g.ncols)
        ch = terrain.channel[r0:r1, c0:c1]
        if ch.any():
            rr, kk = np.nonzero(ch)
            # the largest watercourse within reach (the stream the crossing spans)
            i = int(np.argmax(terrain.area_km2[rr + r0, kk + c0]))
            nx, ny = g.rc_to_xy(rr[i] + r0, kk[i] + c0)
            log("Source snapped %.0f m onto the watercourse (catchment %.1f km2)." %
                (np.hypot(nx - sx, ny - sy), terrain.area_km2[rr[i] + r0, kk[i] + c0]))
            c["source_xy"] = (nx, ny)
            sx, sy = nx, ny
    cls = build_surface_classes(terrain.grid, c["surface_rules"], log=log)
    try:
        tr = html_export._to_ll(terrain.grid.wkt)
        lo, la, _ = tr.TransformPoint(sx, sy)
        c["_latlon"] = (la, lo)
    except Exception:
        logging.getLogger(__name__).debug("Optional operation failed in run_scenario; existing fallback retained.", exc_info=True)
    sc = build_scenario(c, log)
    controls = build_controls(terrain.grid, c, log)
    ext = terrain.grid.extent()
    recs = load_receptors(c["receptor_specs"], terrain.grid.wkt, ext, log=log)
    sim = Simulator(terrain, cls, sc, receptors=recs, log=log, feedback=feedback,
                    controls=controls)
    res = sim.run()
    res.terrain = terrain
    # ---- outputs -----------------------------------------------------------
    folder = c["out_dir"]
    files, window = outputs.write_rasters(res, folder, prefix, log=log)
    streams = terrain.stream_lines()
    gpkg = outputs.write_vectors(res, os.path.join(folder, prefix + "_results.gpkg"), prefix,
                                 res.receptors, streams=streams, log=log)
    files["gpkg"] = gpkg
    files.update(outputs.write_tables(res, folder, prefix, res.receptors, log=log))
    fdir = os.path.join(folder, prefix + "_frames")
    pngs, bounds, fpngs = outputs.frame_pngs_3857(res, window, fdir, prefix, log=log)
    hill_png = os.path.join(fdir, prefix + "_hillshade.png")
    hb = outputs.hillshade_png(terrain, window, hill_png)
    files.update(frames=pngs, fall_frames=fpngs, frame_bounds=bounds, hillshade=hill_png,
                 hill_bounds=hb, window=window)
    if sc.heat_mode:
        files["temp_frames"] = outputs.field_pngs_3857(res, window, fdir, prefix, "temp", outputs.TEMP_RAMP,
                                                       linear=True)
        files["tau_frames"] = outputs.field_pngs_3857(res, window, fdir, prefix, "tau", outputs.YIELD_RAMP)
    summary = summarise(res, c, files)
    files["run_json"] = write_run_record(res, c, files, summary, prefix)
    try:
        add_captions(files["run_json"], files["timeseries"], res.receptors, log=log)
    except Exception as ex:          # captions are a presentation aid; never stop a run
        log("Captions not written: %s" % ex)
    with open(os.path.join(folder, prefix + "_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=float)
    files["summary"] = summary
    if c["html"]:
        files["html"] = write_html(res, c, files, summary, log=log)
    res.files = files
    res.config = c
    return res


def build_controls(grid, c, log=print):
    """Rasterise roughness zones, barriers and sinks onto the simulation grid."""
    ctl = {}
    shape = (grid.nrows, grid.ncols)
    if c.get("roughness_zones"):
        n = np.full(shape, np.nan)
        for z in c["roughness_zones"]:
            if z.get("attribute"):
                a = rasterize_layer(z["path"], z.get("layer"), grid, where=z.get("where"),
                                    attribute=z["attribute"], all_touched=False)
            else:
                a = rasterize_layer(z["path"], z.get("layer"), grid, where=z.get("where"),
                                    burn=float(z["n"]), all_touched=False)
            n = np.where(a > 0, a, n)
        ctl["n"] = n
        log("Roughness zones applied: %d cells." % int(np.isfinite(n).sum()))
    if c.get("barriers"):
        b = np.zeros(shape)
        for z in c["barriers"]:
            a = rasterize_layer(z["path"], z.get("layer"), grid, where=z.get("where"),
                                burn=float(z.get("height_m", 1.0)), all_touched=True)
            b = np.maximum(b, a)
        ctl["barrier_m"] = b
        log("Barriers (berms / bunds) applied: %d cells, up to %.1f m high." % (int((b > 0).sum()), b.max()))
    if c.get("sinks"):
        m = np.zeros(shape, bool)
        for z in c["sinks"]:
            m |= rasterize_layer(z["path"], z.get("layer"), grid, where=z.get("where"),
                                 all_touched=True) > 0
        ctl["sink"] = m
        log("Recovery / outflow cells: %d." % int(m.sum()))
    return ctl


def summarise(res, c, files):
    sc = res.scenario
    mb = res.mass_balance[-1]
    reached = [r for r in res.receptors if r["arrival_h"] is not None]
    reached.sort(key=lambda r: r["arrival_h"])
    return dict(
        name=sc.name, kind=sc.kind, product=sc.oil["label"], month=sc.month_name, sim_hours=sc.sim_hours,
        temp_c=sc.temp_c, released_m3=mb["released"], evaporated_m3=mb["evaporated"],
        infiltrated_m3=mb["infiltrated"], surface_m3=mb["surface"], stranded_m3=mb["stranded"],
        on_water_m3=mb["on_water"], left_domain_m3=mb["exported"] + mb["spray_out"],
        recovered_m3=mb.get("recovered", 0.0), heat_mode=bool(sc.heat_mode),
        oil_temp_min_c=float(np.nanmin(res.temp_min)) if np.isfinite(res.temp_min).any() else None,
        oil_temp_max_c=float(np.nanmax(res.temp_max)) if np.isfinite(res.temp_max).any() else None,
        yield_max_pa=float(res.yield_max.max()),
        footprint_ha=mb["footprint_ha"], max_thickness_mm=float(res.hmax.max() * 1000),
        fallout_area_ha=float((res.fallout >= 1e-5).sum() * res.grid.cell_area / 1e4),
        fallout_max_mm=float(res.fallout.max() * 1000),
        gw_area_ha=float((np.nan_to_num(res.penetration, posinf=0) >= sc.dtw_m).sum() *
                         res.grid.cell_area / 1e4),
        first_channel=res.first_channel,
        first_export_h=res.first_export_h,
        receptors_reached=[dict(name=r["name"], type=r["type"], arrival_h=r["arrival_h"],
                                dist_km=r["dist_km"]) for r in reached],
        receptors_fallout=[dict(name=r["name"], type=r["type"], fallout_mm=r["fallout_mm"])
                           for r in res.receptors if r["fallout_mm"] >= 0.001 and r["arrival_h"] is None],
        steps=res.steps)


def write_html(res, c, files, summary, log=print):
    sc = res.scenario
    from osgeo import osr
    tr = html_export._to_ll(res.grid.wkt)
    lon, lat, _ = tr.TransformPoint(*sc.source_xy)
    b = files["frame_bounds"]
    hb = files["hill_bounds"]
    frames = []
    caps = []
    try:
        with open(files["run_json"]) as f_:
            caps = json.load(f_).get("captions") or []
    except Exception:
        logging.getLogger(__name__).debug("Optional operation failed in write_html; existing fallback retained.", exc_info=True)
    for i, fr in enumerate(res.frames):
        frames.append(dict(cap=caps[i] if i < len(caps) else None,
                           t=round(fr["t_h"], 3), oil=html_export._uri(files["frames"][i]),
                           fall=html_export._uri(files["fall_frames"][i]),
                           rain=round(fr["rain_mmh"], 2),
                           wind=[round(fr["wind"][0], 1), round(fr["wind"][1], 2)] if fr["wind"] else None))
    recs = []
    for r in res.receptors:
        rl, rb, _ = tr.TransformPoint(r["x"], r["y"])
        if not (hb[0] - 0.05 <= rb <= hb[2] + 0.05 and hb[1] - 0.05 <= rl <= hb[3] + 0.05):
            continue
        recs.append(dict(name=r["name"], type=r["type"], lat=round(rb, 6), lon=round(rl, 6),
                         dist=round(r["dist_km"], 3), arrival=r["arrival_h"],
                         max_mm=round(r["max_mm"], 3), fall_mm=round(r["fallout_mm"], 4),
                         status=outputs.receptor_status(r), gw=bool(r["gw"])))
    layers = []
    clip = [hb[0] - 0.01, hb[1] - 0.01, hb[2] + 0.01, hb[3] + 0.01]
    # model-derived layers first
    gp = files["gpkg"]
    model_layers = [
        dict(name="Drainage network (from DEM)", path=gp, layer="drainage_network_dem",
             style=dict(color="#1c7ed6", weight=1.6, opacity=0.9), visible=True, label="order"),
        dict(name="Arrival isochrones", path=gp, layer="arrival_isochrones",
             style=dict(color="#7a1fa2", weight=1, fillOpacity=0.0, dash="4 3"), visible=False,
             label="label"),
        dict(name="Final oiled footprint", path=gp, layer="oil_footprint",
             style=dict(color="#b42318", weight=1.2, fillOpacity=0.0), visible=False, label="class"),
        dict(name="Groundwater concern zone", path=gp, layer="groundwater_concern",
             style=dict(color="#6f2c00", weight=1.2, fill="#ff8c00", fillOpacity=0.25), visible=False,
             label="class"),
    ]
    for L in model_layers + list(c["html_layers"]):
        gj = html_export.layer_geojson(L["path"], L.get("layer"), clip_ll=clip,
                                       fields=L.get("fields"), simplify_m=L.get("simplify", 0.0),
                                       where=L.get("where"))
        if not gj or not gj["features"]:
            continue
        st = dict(L.get("style", {}))
        st["label"] = L.get("label", "name")
        if L.get("popup"):
            st["popup"] = L["popup"]
        if L.get("permanent"):
            st["permanent"] = True
        layers.append(dict(name=L["name"], kind=L.get("kind", "vector"), geojson=gj, style=st,
                           visible=L.get("visible", True)))
    fields = []
    if files.get("temp_frames"):
        fields.append(dict(key="temp", name="Oil temperature", unit="C", ramp=outputs.TEMP_RAMP, linear=True,
                           t=[round(fr["t_h"], 3) for fr in res.frames], bounds=b,
                           frames=[html_export._uri(p) for p in files["temp_frames"]]))
        fields.append(dict(key="tau", name="Yield stress (gelled oil)", unit="Pa", ramp=outputs.YIELD_RAMP,
                           linear=False, t=[round(fr["t_h"], 3) for fr in res.frames], bounds=b,
                           frames=[html_export._uri(p) for p in files["tau_frames"]]))
    if c.get("h2s_field"):
        fields.append(c["h2s_field"])
    h2s_rec = c.get("h2s_receptors") or {}
    for rr_ in recs:
        if rr_["name"] in h2s_rec:
            rr_["h2s"] = h2s_rec[rr_["name"]]
    sandbox = None
    try:
        sandbox = html_export.sandbox_payload(res.terrain, sc.source_xy, c.get("sandbox_half_m", 6000.0),
                                              sc.q_spec)
    except Exception as ex:
        log("Sandbox tracer not embedded: %s" % ex)
    if sc.kind == "TANK":
        release = "%.0f m3 %s, released over %.0f min" % (sc.volume_m3, sc.oil["label"].lower(),
                                                         sc.release_minutes)
    else:
        release = "%.0f bbl/d (%.1f m3/h) for %.0f h, %.0f%% as spray" % (
            sc.rate_bbl_d, sc.rate_bbl_d * BBL_M3 / 24.0, sc.duration_h, sc.spray_fraction * 100)
    flowing = sc.q_spec > 0
    stream_state = ("flowing (wet season, %.0f L/s per km2)" % sc.q_spec) if flowing else \
        "dry beds - oil collects in the channels (dry season)"
    events = []
    if sc.kind == "BLOWOUT" and sc.duration_h < sc.sim_hours:
        events.append(dict(t=sc.duration_h, text="Well brought under control - spillage stops at T+%.0f h" % sc.duration_h))
    for st_, du, dep in sc.rain:
        events.append(dict(t=st_, text="Storm begins: %.0f mm over %.0f h" % (dep, du)))
    fc = res.first_channel
    if fc:
        events.append(dict(t=fc["t_h"], text="Oil enters the drainage network at T+%.1f h" % fc["t_h"]))
    if res.first_export_h is not None:
        events.append(dict(t=res.first_export_h, text="Oil leaves the model area downstream at T+%.1f h" %
                           res.first_export_h))
    for r in summary["receptors_reached"][:6]:
        events.append(dict(t=r["arrival_h"], text="Oil reaches %s (%s)" % (r["name"], r["type"])))
    events.sort(key=lambda e: e["t"])
    facts = [["Release", release], ["Product", sc.oil["label"]],
             ["Heat transfer", ("on - source oil %s, pour point %.0f &deg;C, oil %.0f-%.0f &deg;C over the run" % (
                 ("%.0f &deg;C" % sc.discharge_temp_c) if sc.discharge_temp_c is not None else "at air temperature",
                 sc.oil["pour_point_c"], summary.get("oil_temp_min_c") or 0, summary.get("oil_temp_max_c") or 0))
              if sc.heat_mode else "off (constant properties)"],
             ["Month / ground temperature", "%s / %.1f &deg;C" % (sc.month_name, sc.temp_c)],
             ["Storms", "; ".join("T+%.0f h, %.0f mm in %.0f h" % (a, d, b) for a, b, d in sc.rain) or "none"],
             ["Depth to groundwater (assumed)", "%.1f m" % sc.dtw_m],
             ["Final footprint", "%.2f ha" % summary["footprint_ha"]],
             ["Spray fallout area (>= 0.01 mm)", "%.1f ha" % summary["fallout_area_ha"]],
             ["Evaporated / soaked in", "%.0f / %.0f m3" % (summary["evaporated_m3"], summary["infiltrated_m3"])],
             ["Recovered at sumps / recovery points", "%.1f m3" % summary.get("recovered_m3", 0.0)],
             ["First entry to drainage", ("T+%.1f h" % fc["t_h"]) if fc else "none in %.0f h" % sc.sim_hours],
             ["Prepared by", c.get("prepared_by") or "-"]]
    notes = list(c.get("notes") or []) + [
        "Terrain: %s (cell %.0f m), smoothed and with depressions shallower than %.1f m filled; "
        "deeper depressions are kept as ponds." % (os.path.basename(c["dem"]), res.grid.dx,
                                                    c["pit_tolerance_m"]),
        "Overland flow: free-surface multiple-flow-direction routing; speed is the lesser of the "
        "viscous thin-film and Manning values on the oil depth above the surface retention depth.",
        "Evaporation: Fingas (2004) generic equations; viscosity rises with evaporation "
        "(Mackay). Infiltration: Green-Ampt with conductivity scaled for oil.",
        "Spray: droplets released at an effective height of %.0f m, carried downwind and settled at "
        "their terminal velocity (Atlas et al. 1973)." % sc.release_height_m,
        ("Heat transfer: oil temperature from solar radiation, long-wave exchange, wind convection and "
         "ground conduction; density, viscosity and yield stress follow the temperature, and gelled oil "
         "below its pour point moves only where gravity overcomes the yield stress (Bingham film)."
         if sc.heat_mode else
         "Constant-property mode: oil properties fixed at the mean air temperature (heat transfer off)."),
        "This is a screening simulation for emergency planning. It does not model fire; H2S is modelled "
        "separately (see the H2S runs). Berms, sumps and recovery points are included only where given."]
    # initial map view: the oiled area (thickness >= threshold) plus 400 m
    m = res.hmax >= sc.thickness_threshold_m
    if not m.any():
        m = res.fallout >= 1e-5
    if m.any():
        rr, cc = np.nonzero(m)
        x0, y0 = res.grid.rc_to_xy(rr.max(), cc.min())
        x1, y1 = res.grid.rc_to_xy(rr.min(), cc.max())
    else:
        x0, y0 = x1, y1 = sc.source_xy
    pad = 400.0
    a_lon, a_lat, _ = tr.TransformPoint(x0 - pad, y0 - pad)
    b_lon, b_lat, _ = tr.TransformPoint(x1 + pad, y1 + pad)
    focus = [a_lat, a_lon, b_lat, b_lon]
    payload = dict(
        focus=focus,
        meta=dict(title=c.get("title") or ("%s - %s" % (c["well"], sc.name)),
                  subtitle=c.get("subtitle") or "Onshore oil spill simulation",
                  kind=sc.kind, source=[lat, lon], source_label=c["source_label"],
                  start_iso=c.get("start_iso"), facts=facts, notes=notes,
                  stream_state=stream_state, events=events,
                  frame_h=sc.frame_minutes / 60.0),
        frames=frames, bounds=b, hill=dict(uri=html_export._uri(files["hillshade"]), bounds=hb),
        mb=[{k: (round(v, 3) if isinstance(v, float) and math.isfinite(v) else (None if isinstance(v, float) else v))
             for k, v in r.items()} for r in res.mass_balance],
        receptors=recs, layers=layers, fields=fields, sandbox=sandbox, links=c.get("links") or [],
        ramps=dict(oil=outputs.OIL_RAMP, fall=outputs.FALL_RAMP))
    p = os.path.join(c["out_dir"], (c["prefix"] or c["name"].replace(" ", "_")) + "_animation.html")
    html_export.build_html(p, payload)
    log("Interactive HTML written: %s (%.1f MB)" % (p, os.path.getsize(p) / 1e6))
    return p


def add_captions(run_json, ts_path, receptors, log=print):
    """Write the per-frame captions into the run record and return them."""
    from . import narrative
    with open(run_json) as f:
        rec = json.load(f)
    caps = narrative.build(rec, ts_path, receptors, log=log)
    rec["captions"] = caps
    with open(run_json, "w") as f:
        json.dump(rec, f, indent=1, default=float)
    return caps


def write_run_record(res, c, files, summary, prefix):
    """Everything the animation renderer needs to redraw a finished run."""
    sc = res.scenario
    rec = dict(
        version=1, name=sc.name, well=c.get("well"), title=c.get("title") or sc.name,
        subtitle=c.get("subtitle") or "Onshore oil spill simulation",
        kind=sc.kind, product=sc.oil["label"], month=sc.month_name, start_iso=c.get("start_iso"),
        source_xy=list(sc.source_xy), source_label=c.get("source_label"),
        crs_wkt=res.grid.wkt, frame_minutes=sc.frame_minutes,
        timeseries=os.path.basename(files["timeseries"]),
        fallout_timeseries=os.path.basename(files["fallout_timeseries"]) if files.get("fallout_timeseries") else None,
        frames=[dict(t_h=fr["t_h"], rain_mmh=fr["rain_mmh"],
                     wind=list(fr["wind"]) if fr["wind"] else None) for fr in res.frames],
        mass_balance=res.mass_balance,
        stream_state=("flowing" if sc.q_spec > 0 else "dry"),
        storms=[list(r) for r in sc.rain], duration_h=sc.duration_h, sim_hours=sc.sim_hours,
        summary=summary, prepared_by=c.get("prepared_by", ""))
    p = os.path.join(c["out_dir"], prefix + "_run.json")
    with open(p, "w") as f:
        json.dump(rec, f, indent=1, default=float)
    return p
