# -*- coding: utf-8 -*-
"""
Well spray ("oil rain") fallout.

A gas-driven surface spillage at the well breaks part of the liquid stream into droplets.
Each droplet class is released at an effective height H, is carried downwind
at the wind speed at that height and settles at its terminal velocity, so it
lands at x = u(H) * H / vt.  Scatter about the mean landing point comes from
wind-direction meander (sigma_theta) and from the spread of release heights.

Terminal velocity of a drop in air: Atlas, Srivastava & Sekhon (1973)
  vt = 9.65 - 10.3 exp(-0.6 D)   (m/s, D in mm, water drops, D > 0.1 mm)
scaled by sqrt(rho_oil / rho_water); Stokes law below 0.1 mm.
Wind profile: power law u(z) = u10 (z / 10)^p, p = 0.14 (open savannah).
"""
import math
import numpy as np

RHO_AIR = 1.16      # kg/m3 at ~30 C
MU_AIR = 1.86e-5    # Pa.s


def terminal_velocity(d_m, rho_oil):
    d_mm = d_m * 1000.0
    if d_mm < 0.1:
        return (rho_oil - RHO_AIR) * 9.81 * d_m ** 2 / (18.0 * MU_AIR)
    vt_w = max(9.65 - 10.3 * math.exp(-0.6 * d_mm), 0.27)
    return vt_w * math.sqrt(rho_oil / 1000.0)


def droplet_classes(d50_um=600.0, gsd=2.2, n=12):
    """Log-normal mass distribution of droplet diameters.
    Returns list of (diameter_m, mass_fraction)."""
    lo = math.log(d50_um) - 2.5 * math.log(gsd)
    hi = math.log(d50_um) + 2.5 * math.log(gsd)
    edges = np.linspace(lo, hi, n + 1)
    mids = 0.5 * (edges[:-1] + edges[1:])
    s = math.log(gsd)
    cdf = lambda x: 0.5 * (1.0 + math.erf((x - math.log(d50_um)) / (s * math.sqrt(2))))
    out = []
    tot = 0.0
    for i in range(n):
        f = cdf(edges[i + 1]) - cdf(edges[i])
        out.append((math.exp(mids[i]) * 1e-6, f))
        tot += f
    # tails folded into the end classes
    out[0] = (out[0][0], out[0][1] + cdf(edges[0]))
    out[-1] = (out[-1][0], out[-1][1] + (1.0 - cdf(edges[-1])))
    return out


def fallout_kernel(grid, sx, sy, wind_to_deg, u10, sigma_theta_deg, rho_oil,
                   release_height_m=15.0, d50_um=600.0, gsd=2.2, exponent=0.14):
    """Fraction of sprayed oil landing in each grid cell for one wind state.

    wind_to_deg : direction the wind blows TOWARDS, degrees clockwise from grid north
    Returns (kernel array summing to <= 1, fraction leaving the grid).
    """
    nr, nc = grid.nrows, grid.ncols
    xs = grid.gt[0] + (np.arange(nc) + 0.5) * grid.gt[1]
    ys = grid.gt[3] + (np.arange(nr) + 0.5) * grid.gt[5]
    X, Y = np.meshgrid(xs - sx, ys - sy)
    th = math.radians(wind_to_deg)
    ex, ey = math.sin(th), math.cos(th)            # downwind unit vector
    along = X * ex + Y * ey
    cross = -X * ey + Y * ex
    k = np.zeros((nr, nc))
    uH = max(u10, 0.3) * (release_height_m / 10.0) ** exponent
    sig_t = math.radians(max(sigma_theta_deg, 3.0))
    cell = grid.cell_area
    for d, f in droplet_classes(d50_um, gsd):
        vt = terminal_velocity(d, rho_oil)
        x0 = uH * release_height_m / vt
        sx_ = 0.35 * x0 + 0.7 * grid.dx          # height spread + plume width
        sy_ = math.tan(min(sig_t, 1.2)) * x0 + 0.7 * grid.dx
        g = np.exp(-0.5 * ((along - x0) / sx_) ** 2 - 0.5 * (cross / sy_) ** 2)
        g *= cell / (2.0 * math.pi * sx_ * sy_)
        k += f * g
    # drop the far Gaussian tails (keeps the deposit compact and the run fast)
    k[k < 1e-7 * max(float(k.max()), 1e-30)] = 0.0
    inside = float(k.sum())
    if inside > 1.0:
        k /= inside
        inside = 1.0
    return k, max(0.0, 1.0 - inside)


def wind_series(month_u, month_v, month_speed, hours, seed=7, diurnal_amp=0.35):
    """Hourly wind for the simulation, built from the monthly vector-mean wind.

    Direction meanders about the resultant direction with a spread set by the
    steadiness (|vector mean| / scalar mean); speed follows a diurnal cycle
    (stronger by day, calm at night) about the scalar mean.
    Returns list of (wind_to_deg, speed_mps, sigma_theta_deg) for each hour.
    """
    rng = np.random.RandomState(seed)
    vec = math.hypot(month_u, month_v)
    steadiness = min(vec / max(month_speed, 1e-6), 0.99)
    sigma = math.degrees(math.sqrt(max(2.0 * (1.0 - steadiness), 0.0)))
    sigma = min(max(sigma, 10.0), 75.0)
    mean_to = math.degrees(math.atan2(month_u, month_v)) % 360.0
    out = []
    dev = 0.0
    for h in range(int(hours) + 1):
        dev = 0.8 * dev + rng.normal(0.0, sigma * 0.6)
        dev = max(-2.0 * sigma, min(2.0 * sigma, dev))
        hod = h % 24
        diurnal = 1.0 + diurnal_amp * math.sin(2 * math.pi * (hod - 9) / 24.0)
        spd = max(0.4, month_speed * diurnal * (1.0 + rng.normal(0, 0.12)))
        out.append(((mean_to + dev) % 360.0, spd, 0.5 * sigma + 8.0))
    return out
