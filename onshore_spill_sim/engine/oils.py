# -*- coding: utf-8 -*-
"""
Product library for the onshore spill model.

Each product carries the properties the engine needs:
  density_kgm3      density at the reference temperature
  thermal_exp       density change per deg C (kg/m3/C, negative)
  viscosity_mpas    dynamic viscosity at visc_ref_c (mPa.s = cP) when fluid
  visc_ref_c        reference temperature of the viscosity value
  visc_act_k        Andrade activation term B (K): mu = mu_ref exp(B (1/T - 1/T_ref))
  pour_point_c      pour point; below it the oil gels and develops a yield stress
  yield_ref_pa      yield stress at the pour point (Pa)
  yield_decade_c    temperature drop that multiplies the yield stress by ten
  evap_form         'ln' or 'sqrt' (Fingas 2004 generic equations)
  pct_distilled_180 mass % distilled at 180 C (input to the Fingas generic equations)
  evap_max_pct      ceiling on evaporative loss (non-volatile residue limit)
  visc_evap_c       Mackay constant C in mu = mu0 * exp(C * Fe)
  cp_jkgk           specific heat (J/kg/K)
  label             text used in reports

Values are typical published values for the product class and are meant
to be replaced with the assay of the actual product when available
(for the Voltaian prospects the hydrocarbon phase is still uncertain).
"""

_COMMON = dict(thermal_exp=-0.65, yield_ref_pa=1.0, yield_decade_c=5.0, cp_jkgk=2000.0,
               visc_ref_c=25.0)

OILS = {
    "LIGHT_CRUDE": dict(_COMMON,
        label="Light crude oil (API ~35, generic)",
        density_kgm3=850.0, viscosity_mpas=8.0, visc_act_k=3000.0, pour_point_c=-9.0,
        evap_form="ln", pct_distilled_180=25.0, evap_max_pct=40.0, visc_evap_c=10.0),
    "WAXY_CRUDE": dict(_COMMON,
        label="Waxy crude oil (API ~32, pour point 27 C)",
        density_kgm3=860.0, viscosity_mpas=12.0, visc_ref_c=40.0, visc_act_k=3500.0,
        pour_point_c=27.0, yield_ref_pa=2.0, yield_decade_c=4.0,
        evap_form="ln", pct_distilled_180=18.0, evap_max_pct=30.0, visc_evap_c=10.0),
    "MEDIUM_CRUDE": dict(_COMMON,
        label="Medium crude oil (API ~25, generic)",
        density_kgm3=900.0, viscosity_mpas=60.0, visc_act_k=4000.0, pour_point_c=-3.0,
        evap_form="ln", pct_distilled_180=15.0, evap_max_pct=28.0, visc_evap_c=10.0),
    "HEAVY_WAXY_CRUDE": dict(_COMMON,
        label="Heavy / waxy crude (API ~18, bitumen-type shows in SGST wells)",
        density_kgm3=950.0, viscosity_mpas=1500.0, visc_act_k=6500.0, pour_point_c=24.0,
        yield_ref_pa=5.0, evap_form="ln", pct_distilled_180=4.0, evap_max_pct=10.0,
        visc_evap_c=10.0),
    "CONDENSATE": dict(_COMMON,
        label="Gas condensate (API ~50, generic)",
        density_kgm3=780.0, viscosity_mpas=1.0, visc_act_k=1500.0, pour_point_c=-30.0,
        evap_form="ln", pct_distilled_180=55.0, evap_max_pct=80.0, visc_evap_c=1.0),
    "DIESEL": dict(_COMMON,
        label="Automotive gas oil / rig diesel",
        density_kgm3=840.0, viscosity_mpas=3.5, visc_act_k=2200.0, pour_point_c=-15.0,
        evap_form="sqrt", pct_distilled_180=10.0, evap_max_pct=60.0, visc_evap_c=1.0),
    "BASE_OIL_OBM": dict(_COMMON,
        label="Low-toxicity mineral base oil (oil-based mud)",
        density_kgm3=800.0, viscosity_mpas=2.5, visc_act_k=2500.0, pour_point_c=-20.0,
        evap_form="sqrt", pct_distilled_180=2.0, evap_max_pct=20.0, visc_evap_c=1.0),
}

BBL_M3 = 0.158987294928  # 1 US oil barrel in cubic metres


def get_oil(key, overrides=None):
    o = dict(OILS[key])
    o["key"] = key
    if overrides:
        o.update({k: v for k, v in overrides.items() if v is not None})
    return o


def fingas_coefficient(oil, temp_c):
    """Coefficient of the Fingas (2004) generic evaporation equation.

    ln form   : %Ev = [0.165 (%D) + 0.045 (T - 15)] ln(t)
    sqrt form : %Ev = [0.0254 (%D) + 0.01 (T - 15)] sqrt(t)
    t in minutes, T in deg C, %D = % distilled at 180 C.
    """
    d = oil["pct_distilled_180"]
    if oil["evap_form"] == "ln":
        a = 0.165 * d + 0.045 * (temp_c - 15.0)
    else:
        a = 0.0254 * d + 0.01 * (temp_c - 15.0)
    return max(a, 1e-4)
