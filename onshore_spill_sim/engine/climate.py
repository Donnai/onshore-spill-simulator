# -*- coding: utf-8 -*-
"""
Climate defaults for the northern Voltaian Basin (Tamale area).

Rainfall and temperature: Tamale 1991-2020 normals (NOAA, as published on the
Tamale climate table).  Winds: Tamale (DGLE) synoptic station monthly
vector-mean winds 1973-2026 compiled for the rig-pad orientation study
(u = eastward, v = northward component of air motion, m/s).
Stream specific discharge: indicative seasonal values for small Voltaian
tributaries (dry Nov-Apr, peak Aug-Oct); to be replaced by Ghana Hydrological
Authority gauge data for the Nasia / Daka sub-basins when available.
"""
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

RAIN_MM = [3.2, 10.4, 34.9, 82.1, 114.4, 139.7, 178.1, 189.3, 206.1, 106.2, 6.0, 3.4]
RAIN_DAYS = [0.2, 0.9, 2.4, 4.8, 7.2, 9.5, 10.5, 11.7, 14.3, 8.8, 0.6, 0.2]
TEMP_MEAN_C = [27.9, 30.1, 30.9, 29.8, 28.6, 26.8, 26.0, 25.8, 25.9, 27.3, 28.2, 27.5]
TEMP_MAX_C = [36.0, 38.2, 38.9, 37.0, 34.6, 32.2, 30.6, 30.2, 31.2, 33.6, 36.5, 36.3]
RH_PCT = [27, 32, 45, 61, 70, 77, 80, 79, 80, 74, 59, 40]

WIND_U = [-1.039038, -0.799912, 0.084095, 0.508494, 0.545411, 0.459206,
          0.892590, 0.818345, 0.422849, 0.183994, -0.358540, -0.717207]
WIND_V = [-1.071848, -0.661258, 0.577042, 1.154736, 1.177098, 1.089735,
          1.138358, 1.015346, 0.604763, 0.749095, 0.214328, -0.573651]
WIND_SPEED = [1.978711, 2.029398, 2.024875, 2.184169, 1.910125, 1.722347,
              1.942316, 1.800111, 1.415169, 1.468280, 1.524691, 1.765500]

# cloud / haze fraction used for solar and long-wave radiation (indicative:
# harmattan haze Dec-Feb, heavy convective cloud Jul-Sep)
CLOUD = [0.35, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.80, 0.75, 0.60, 0.40, 0.35]

# L/s per km2 - indicative baseflow-plus-quickflow for small tributaries
STREAM_Q_SPEC = [0.0, 0.0, 0.0, 0.1, 0.5, 2.0, 5.0, 12.0, 20.0, 10.0, 1.0, 0.0]

# soil moisture state used for oil infiltration (fraction of pore space filled by water)
SOIL_WATER_SAT = [0.15, 0.12, 0.15, 0.30, 0.45, 0.60, 0.70, 0.80, 0.85, 0.70, 0.40, 0.20]

# indicative depth to groundwater (m) - shallow regolith aquifer, seasonal
DTW_M = [8.0, 9.0, 9.5, 9.0, 7.0, 5.0, 3.5, 2.5, 2.0, 3.0, 5.0, 7.0]


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

SEASON_NAME = ["Dry season (Harmattan)", "Dry season (Harmattan)", "Late dry season - change to the rains",
               "Early rainy season", "Early rainy season", "Main rainy season", "Main rainy season",
               "Main rainy season (peak)", "Main rainy season (peak)", "End of the rainy season",
               "Start of the dry season", "Dry season (Harmattan)"]


def compass(deg):
    return _DIR16[int(round((deg % 360) / 22.5)) % 16]


def wind_from_deg(u, v):
    """Direction the wind blows FROM (degrees from north) for air moving with (u, v)."""
    import math
    return (math.degrees(math.atan2(-u, -v)) + 360.0) % 360.0


def typical_storms(month_index, sim_hours, start_hour=6.0):
    """Indicative storms for a run: the mean number of rain days in the period (Tamale normals),
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
            pass
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
        sky="Cloud / haze about %.0f%%%s." % (c["cloud"] * 100, " (Harmattan dust haze)" if m in (0, 1, 11) else ""),
    )
    hints = []
    if c["tmax_c"] >= 36:
        hints.append("Hot afternoons heat oil on the ground well above air temperature; waxy oil re-mobilises by day.")
    if c["temp_c"] <= 27 or m in (0, 1, 11):
        hints.append("Clear, still nights cool the oil quickly; waxy crude can gel overnight.")
    if m in (0, 1, 10, 11):
        hints.append("Stable Harmattan nights with light wind give the longest gas reach.")
    if c["rain_days"] >= 7:
        hints.append("Frequent storms wash oil off the pad into the drainage lines.")
    elif c["rain_days"] >= 2:
        hints.append("Isolated storms are possible; add one to test run-off.")
    lines["hints"] = hints
    lines["month"] = c["month"]
    lines["wind_from_deg"] = frm
    lines["climate"] = c
    return lines
