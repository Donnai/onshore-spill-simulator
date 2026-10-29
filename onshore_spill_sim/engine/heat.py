# -*- coding: utf-8 -*-
"""
Heat transfer between spilled oil and its surroundings, and the
temperature-dependent rheology that follows from it.

Energy balance of a thin oil layer of depth h (W/m2):

  rho cp h dT/dt = (1 - a) S                        absorbed solar radiation
                 + e sigma (e_a Ta^4 - T^4)          net long-wave radiation
                 - h_c (T - Ta)                      convection, h_c = 5.7 + 3.8 U
                 - h_g (T - Tg)                      conduction to the ground

S       clear-sky global radiation (Haurwitz, 1945) reduced for cloud and
        harmattan haze (Kasten and Czeplak, 1980)
e_a     clear-sky emissivity of the air (Brutsaert, 1975) raised for cloud
h_c     wind convection coefficient (McAdams / Watmuff et al., 1977)

The balance is linearised about the current temperature and integrated with
an exponential relaxation, which is stable for any time step and any depth.

Rheology (as in temperature-dependent overland oil models such as OilFlow2D):
  density    rho(T) = rho_ref + k_rho (T - T_ref)
  viscosity  mu(T, Fe) = mu_ref exp(B (1/T - 1/T_ref)) exp(C Fe)   (Andrade x Mackay)
  yield      tau_y(T) = tau_pp 10^((T_pp - T) / dT10)              below about T_pp + 2 dT10
"""
import math
import numpy as np

SIGMA = 5.670374e-8
EMISS_OIL = 0.95
ALBEDO_OIL = 0.05


def solar_elevation_sin(lat_deg, doy, hod, lon_deg=0.0):
    """sin(solar elevation) for local clock hour hod (UTC+0, longitude correction applied)."""
    decl = math.radians(23.45) * math.sin(math.radians(360.0 / 365.0 * (284 + doy)))
    solar_hour = hod + lon_deg / 15.0
    w = math.radians(15.0 * (solar_hour - 12.0))
    phi = math.radians(lat_deg)
    return math.sin(phi) * math.sin(decl) + math.cos(phi) * math.cos(decl) * math.cos(w)


def global_radiation(lat_deg, doy, hod, cloud, lon_deg=0.0):
    """Global horizontal irradiance (W/m2): Haurwitz clear sky x Kasten-Czeplak cloud factor."""
    s = solar_elevation_sin(lat_deg, doy, hod, lon_deg)
    if s <= 0.01:
        return 0.0
    clear = 1098.0 * s * math.exp(-0.057 / s)
    return clear * (1.0 - 0.75 * cloud ** 3.4)


def air_temperature(t_mean, t_max, hod):
    """Diurnal air temperature: minimum near 03:00, maximum near 15:00."""
    return t_mean + (t_max - t_mean) * math.sin(2.0 * math.pi * (hod - 9.0) / 24.0)


def air_emissivity(ta_c, rh_pct, cloud):
    tk = ta_c + 273.15
    es = 6.112 * math.exp(17.67 * ta_c / (ta_c + 243.5))
    e = max(rh_pct, 1.0) / 100.0 * es
    clear = 1.24 * (e / tk) ** (1.0 / 7.0)
    return min(1.0, clear * (1.0 + 0.22 * cloud ** 2))


def relax_temperature(T, h, dt, Ta, Tg, S, wind, eps_a, rho, cp, hg=10.0):
    """Advance oil temperature arrays T (C) of depth h (m) by dt seconds."""
    Tk = T + 273.15
    Tak = Ta + 273.15
    hc = 5.7 + 3.8 * max(wind, 0.0)
    q0 = ((1.0 - ALBEDO_OIL) * S + EMISS_OIL * SIGMA * (eps_a * Tak ** 4 - Tk ** 4)
          - hc * (T - Ta) - hg * (T - Tg))
    b = 4.0 * EMISS_OIL * SIGMA * Tk ** 3 + hc + hg
    t_eq = T + q0 / b
    cap = rho * cp * np.maximum(h, 2e-4)
    return t_eq + (T - t_eq) * np.exp(-b * dt / cap)


def viscosity(oil, T, Fe):
    """Dynamic viscosity (Pa.s) of the liquid phase."""
    mu_ref = oil["viscosity_mpas"] * 1e-3
    B = oil.get("visc_act_k", 3000.0)
    tref = oil.get("visc_ref_c", 25.0) + 273.15
    Tk = np.maximum(T + 273.15, 200.0)
    mu = mu_ref * np.exp(B * (1.0 / Tk - 1.0 / tref)) * np.exp(oil["visc_evap_c"] * np.clip(Fe, 0, 1))
    return np.clip(mu, 1e-4, 1e4)


def yield_stress(oil, T):
    """Yield stress (Pa); zero well above the pour point."""
    pp = oil["pour_point_c"]
    d10 = max(oil.get("yield_decade_c", 5.0), 0.5)
    ty = oil.get("yield_ref_pa", 1.0) * 10.0 ** ((pp - T) / d10)
    ty = np.where(T > pp + 2.0 * d10, 0.0, ty)
    return np.clip(ty, 0.0, 5000.0)


def density(oil, T):
    return oil["density_kgm3"] + oil.get("thermal_exp", -0.65) * (T - oil.get("visc_ref_c", 25.0))
