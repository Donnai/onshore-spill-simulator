# -*- coding: utf-8 -*-
"""Synthetic example environmental inputs for software demonstrations.

All twelve months use the same invented values. These values are not station
observations, climate normals, site measurements or a regional calibration.
Supply reviewed wind, temperature, runoff, soil and groundwater inputs for any
assessment. The example profile is not suitable for operational predictions.
"""
import logging
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

RAIN_MM = [0.0] * 12
RAIN_DAYS = [0.0] * 12
TEMP_MEAN_C = [20.0] * 12
TEMP_MAX_C = [25.0] * 12
RH_PCT = [50] * 12

WIND_U = [3.0] * 12
WIND_V = [0.0] * 12
WIND_SPEED = [3.0] * 12

# Synthetic cloud fraction, constant throughout the example year.
CLOUD = [0.5] * 12

# Synthetic stream specific discharge, L/s per km2.
STREAM_Q_SPEC = [0.0] * 12

# soil moisture state used for oil infiltration (fraction of pore space filled by water)
SOIL_WATER_SAT = [0.5] * 12

# Synthetic depth to groundwater (m), constant throughout the example year.
DTW_M = [5.0] * 12


def month_climate(month_index):
    m = int(month_index)
    return dict(month=MONTHS[m], rain_mm=RAIN_MM[m], rain_days=RAIN_DAYS[m],
                temp_c=TEMP_MEAN_C[m], tmax_c=TEMP_MAX_C[m], rh=RH_PCT[m],
                u=WIND_U[m], v=WIND_V[m], speed=WIND_SPEED[m],
                cloud=CLOUD[m], q_spec=STREAM_Q_SPEC[m], soil_sat=SOIL_WATER_SAT[m], dtw=DTW_M[m])


def load_wind_csv(path):
    """Read a monthly wind table (month,month_name,u_mps,v_mps,mean_speed_mps)."""
    import csv
    u, v, s = list(WIND_U), list(WIND_V), list(WIND_SPEED)
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            i = int(row["month"]) - 1
            u[i] = float(row["u_mps"])
            v[i] = float(row["v_mps"])
            s[i] = float(row["mean_speed_mps"])
    return u, v, s


# ---------------------------------------------------------------------------
# Season descriptions (used by the simulation panel and the captions)
# ---------------------------------------------------------------------------
_DIR16 = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]

SEASON_NAME = ["Synthetic example - supply reviewed environmental inputs"] * 12


def compass(deg):
    return _DIR16[int(round((deg % 360) / 22.5)) % 16]


def wind_from_deg(u, v):
    """Direction the wind blows FROM (degrees from north) for air moving with (u, v)."""
    import math
    return (math.degrees(math.atan2(-u, -v)) + 360.0) % 360.0


def typical_storms(month_index, sim_hours, start_hour=6.0):
    """Example storms calculated from the configured monthly rain-day inputs,
    each with the mean rain-day depth, falling as late-afternoon convective storms (16:00, 2 h).
    Returns a list of (start_h, duration_h, depth_mm)."""
    m = int(month_index)
    days = sim_hours / 24.0
    n = int(round(RAIN_DAYS[m] / 30.4 * days))
    if n == 0:
        return []
    depth = round(RAIN_MM[m] / max(RAIN_DAYS[m], 0.1), 0)
    out = []
    step = days / n
    for k in range(n):
        day = int(k * step + step / 2.0)
        t = day * 24.0 + (16.0 - start_hour) % 24.0
        if t + 2.0 <= sim_hours:
            out.append((round(t, 1), 2.0, depth))
    return out


def season_summary(month_index, wind_csv=None):
    """Plain-language summary of the conditions a release in this month meets."""
    c = month_climate(month_index)
    m = int(month_index)
    if wind_csv:
        try:
            u, v, s = load_wind_csv(wind_csv)
            c["u"], c["v"], c["speed"] = u[m], v[m], s[m]
        except Exception:
            logging.getLogger(__name__).debug("Optional operation failed in season_summary; existing fallback retained.", exc_info=True)
    frm = wind_from_deg(c["u"], c["v"])
    to = (frm + 180.0) % 360.0
    steady = min((c["u"] ** 2 + c["v"] ** 2) ** 0.5 / max(c["speed"], 0.1), 1.0)
    storms_4d = RAIN_DAYS[m] / 30.4 * 4.0
    dry = c["q_spec"] <= 0.0
    lines = dict(
        season=SEASON_NAME[m],
        temperature="Mean %.1f C, afternoon maximum about %.0f C; relative humidity about %d%%." % (
            c["temp_c"], c["tmax_c"], c["rh"]),
        rain="%.0f mm in the month on about %.1f rain days (%.0f mm per rain day); "
             "about %.1f storms expected in a 4-day period." % (
                 c["rain_mm"], c["rain_days"], c["rain_mm"] / max(c["rain_days"], 0.1), storms_4d),
        wind="Resultant wind from the %s (%.0f deg), mean speed %.1f m/s, %s (steadiness %.0f%%) - gas and "
             "spray drift mostly towards the %s." % (
                 compass(frm), frm, c["speed"], "steady" if steady >= 0.6 else "fairly steady" if steady >= 0.35
                 else "variable in direction", steady * 100, compass(to)),
        streams=("Streams are dry; spilled oil collects in the channel beds." if dry else
                 "Streams are flowing (about %.1f L/s per km2 of catchment); oil reaching a channel is "
                 "carried downstream." % c["q_spec"]),
        soil="Soil pores about %.0f%% water-filled (%s); groundwater about %.1f m below ground." % (
            c["soil_sat"] * 100, "dry - oil soaks in readily" if c["soil_sat"] < 0.35 else
            "wet - little oil soaks in, more runs off", c["dtw"]),
        sky="Example cloud / haze about %.0f%%." % (c["cloud"] * 100),
    )
    hints = ["Synthetic example inputs are not site observations. Review all environmental inputs before assessment."]
    if c["tmax_c"] >= 36:
        hints.append("Hot afternoons heat oil on the ground well above air temperature; waxy oil re-mobilises by day.")
    if c["temp_c"] <= 27:
        hints.append("Clear, still nights cool the oil quickly; waxy crude can gel overnight.")
    if c["rain_days"] >= 7:
        hints.append("Frequent storms wash oil off the pad into the drainage lines.")
    elif c["rain_days"] >= 2:
        hints.append("Isolated storms are possible; add one to test run-off.")
    lines["hints"] = hints
    lines["month"] = c["month"]
    lines["wind_from_deg"] = frm
    lines["climate"] = c
    return lines
