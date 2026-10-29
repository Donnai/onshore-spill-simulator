# -*- coding: utf-8 -*-
"""
Onshore oil-spill simulation engine.

State per grid cell
  h   surface oil thickness (m = m3 of oil per m2 of ground)
  O   'original-equivalent' thickness: h / (1 - E), so that E = 1 - h/O is the
      fraction of the original oil already evaporated (mass-weighted when
      parcels of different age mix)
  B   oil stranded on stream banks / vegetation (immobile part of h)
  F   cumulative infiltrated oil (m)

Processes, in order, every time step
  1  source terms: pool at the wellhead/tank (after the cellar or bund is full)
     and well spray fallout from the hourly wind state
  2  evaporation (Fingas 2004 generic equations, pseudo-time method)
  3  infiltration (Green-Ampt with oil-scaled conductivity and suction)
  4  transport:
       overland  - multiple-flow-direction on the oil free surface (DEM + h),
                   speed = min(viscous thin-film, Manning) on the mobile depth
                   above the surface retention depth; ponds fill and spill
       rain      - during runoff the mobile oil is carried at a fraction of
                   the kinematic sheet-flow speed
       channels  - flowing streams carry the oil down the D8 path at the
                   stream surface speed with exponential bank stranding;
                   dry beds behave as overland surfaces
  5  bookkeeping: frames, mass balance, receptor arrival

Heat-transfer mode (heat_mode=True) adds an oil temperature state carried with
the oil (heat content Hc = h T).  Each step the temperature relaxes towards the
radiative-convective balance with the air and the ground (see heat.py), and
density, viscosity, yield stress and the evaporation rate are evaluated from the
local temperature.  Oil flows as a Bingham film: no movement while the driving
shear stress rho g h S is below the yield stress; above it the laminar speed is
reduced by (1 - 1.5 r + 0.5 r^3), r = tau_y / (rho g h S).  A warm wellhead
release therefore spreads freely, gels as it cools at night and moves again
in the afternoon heat.

Control features (optional arrays supplied by the caller):
  n          spatial Manning roughness override (NaN = keep surface class value)
  barrier_m  height of berms / bunds / walls added to the flow surface (m)
  sink       cells where oil is removed (sump, recovery pit, skimming point,
             free-outflow boundary); volume reported as 'recovered'
Extra sources: list of dict(x, y, volume_m3, start_h, release_min, temp_c).
"""
import math
import time as _time
import numpy as np

from .terrain import D8, D8_DIST
from .oils import fingas_coefficient, BBL_M3
from . import fallout as fo
from . import heat as ht

G = 9.81
MU_WATER = 0.8e-3        # Pa.s at ~30 C
SIGMA_RATIO = 0.4        # oil/air vs water/air surface tension (capillary scaling)

# surface classes ------------------------------------------------------------
SURFACES = {
    0: dict(name="Savannah / farmland", n=0.25, hret=0.006, k=2.0e-6, psi=0.11, por=0.40),
    1: dict(name="Rig pad (compacted laterite)", n=0.03, hret=0.002, k=1.0e-7, psi=0.21, por=0.30),
    2: dict(name="Laterite road", n=0.03, hret=0.002, k=2.0e-7, psi=0.21, por=0.30),
    3: dict(name="Stream bed (dry)", n=0.05, hret=0.003, k=1.0e-5, psi=0.05, por=0.38),
    4: dict(name="Open water (dugout / pond)", n=0.02, hret=0.0, k=0.0, psi=0.0, por=0.0),
    5: dict(name="Lined / bunded area", n=0.02, hret=0.001, k=0.0, psi=0.0, por=0.0),
    6: dict(name="Woodland / forest reserve", n=0.40, hret=0.010, k=5.0e-6, psi=0.09, por=0.45),
}


class Scenario(object):
    """Container for all scenario inputs (plain attributes, sensible defaults)."""

    def __init__(self, **kw):
        self.name = "Scenario"
        self.kind = "BLOWOUT"            # BLOWOUT | TANK | CONTINUOUS
        self.source_xy = (0.0, 0.0)
        self.oil = None                  # dict from oils.get_oil
        self.rate_bbl_d = 1500.0         # well spillage / continuous rate
        self.duration_h = 72.0           # release duration
        self.volume_m3 = 60.0            # tank inventory
        self.release_minutes = 10.0      # time for the tank inventory to escape
        self.contained_m3 = 0.0          # volume held by bund / cellar before overflow
        self.spray_fraction = 0.40       # well spillage share leaving as droplets
        self.release_height_m = 40.0
        self.d50_um = 250.0
        self.gsd = 2.2
        self.sim_hours = 72.0
        self.frame_minutes = 30.0
        self.month = 0
        self.temp_c = 28.0
        self.wind = None                 # list of hourly (to_deg, speed, sigma)
        self.rain = []                   # list of (start_h, duration_h, depth_mm)
        self.cn = 78.0                   # SCS curve number for storm runoff
        self.q_spec = 0.0                # stream specific discharge (L/s/km2)
        self.stream_kv = 0.70            # v = kv * Q^0.1 (hydraulic geometry)
        self.q_min_flow = 0.003          # m3/s below which a channel is treated as dry
        self.strand_per_m = 1.5e-4       # bank stranding coefficient (1/m)
        self.rain_carry = 0.5            # fraction of sheet-flow speed given to oil
        self.rain_flush = 0.5            # share of the retention depth flushed by runoff
        self.soil_sat = 0.3              # pore fraction already water-filled
        self.s_or = 0.10                 # residual oil saturation (CONCAWE 3/79: ~40 L/m3 in silty soils)
        self.dtw_m = 6.0                 # depth to groundwater
        self.cfl = 0.7
        self.dt_max = 60.0
        self.thickness_threshold_m = 1e-4   # 0.1 mm = visible sheen / oiled ground
        self.rivulet_width_m = 6.0       # width oil actually wets when running down a slope
        self.rivulet_min_slope = 0.005   # below this the oil spreads as a sheet over the cell
        self.receptor_radius_m = 60.0
        self.seed = 7
        # heat transfer (OilFlow2D-type) inputs
        self.heat_mode = False
        self.discharge_temp_c = None     # oil temperature at the source (None = air temperature)
        self.t_max_c = None              # mean daily maximum air temperature (None = temp_c + 6)
        self.rh_pct = 50.0
        self.cloud = 0.3                 # cloud / haze fraction 0-1
        self.lat = 9.3
        self.lon = -0.5
        self.doy = 15
        self.start_hour = 0.0            # local clock hour at T+0
        self.ground_temp_c = None        # None = mean air temperature + 2 C
        self.h_ground = 10.0             # oil-ground conduction coefficient (W/m2/K)
        self.extra_sources = []
        for k, v in kw.items():
            setattr(self, k, v)


class Result(object):
    pass


class Simulator(object):

    def __init__(self, terrain, surface_class, scenario, receptors=None, log=print,
                 feedback=None, controls=None):
        self.T = terrain
        self.g = terrain.grid
        self.sc = scenario
        self.log = log
        self.fb = feedback
        self.cls = surface_class.astype(np.int16).copy()
        # natural ground on the drainage network takes stream-bed properties
        self.cls[(terrain.channel) & (self.cls == 0)] = 3
        nr, nc = self.g.nrows, self.g.ncols
        self.nr, self.nc = nr, nc
        self.A = self.g.cell_area
        self.dx = self.g.dx
        self.receptors = receptors or []
        # static property maps
        lut = lambda key: np.vectorize(lambda c: SURFACES[int(c)][key])(np.arange(7))
        self.n_map = lut("n")[self.cls]
        self.hret_map = lut("hret")[self.cls]
        self.k_map = lut("k")[self.cls]
        self.psi_map = lut("psi")[self.cls] * SIGMA_RATIO
        self.por_map = lut("por")[self.cls]
        self.water = self.cls == 4
        ctl = controls or {}
        if ctl.get("n") is not None:
            nn = ctl["n"]
            self.n_map = np.where(np.isfinite(nn) & (nn > 0), nn, self.n_map)
        self.zf = terrain.z_flow.astype(np.float64).copy()
        bar = ctl.get("barrier_m")
        self.barrier = np.zeros((nr, nc), bool)
        if bar is not None and np.any(bar > 0):
            self.zf += np.where(bar > 0, bar, 0.0)
            self.barrier = bar > 0
        self.sink = ctl.get("sink")
        if self.sink is not None and not self.sink.any():
            self.sink = None
        # sub-grid rivulets: on sloping natural ground the oil runs in narrow
        # threads, so only part of each cell holds oil and soaks it up
        f = min(1.0, max(scenario.rivulet_width_m / self.dx, 0.02))
        natural = np.isin(self.cls, (0, 3, 6))
        self.fwet = np.where(natural & (terrain.slope >= scenario.rivulet_min_slope), f, 1.0)
        self.hret_map = self.hret_map * self.fwet
        # channel hydraulics
        self.chan = terrain.channel & ~self.water & ~self.barrier
        Q = scenario.q_spec * terrain.area_km2 / 1000.0
        self.Qbase = Q
        # no oil infiltration through a flowing (saturated) stream bed
        self.k_map = np.where(self.chan & (Q >= scenario.q_min_flow), 0.0, self.k_map)
        self.pen_div = np.where(self.por_map > 0, scenario.s_or * self.por_map, np.inf)
        self.rec = terrain.rec.reshape(nr, nc)
        # state
        self.h = np.zeros((nr, nc))
        self.O = np.zeros((nr, nc))
        self.B = np.zeros((nr, nc))
        self.F = np.zeros((nr, nc))
        self.Hc = np.zeros((nr, nc))     # heat content h * T (m C)
        self.Tmax = np.full((nr, nc), np.nan)
        self.Tmin = np.full((nr, nc), np.nan)
        self.mumax = np.zeros((nr, nc))
        self.taumax = np.zeros((nr, nc))
        self.hmax = np.zeros((nr, nc))
        self.fall_cum = np.zeros((nr, nc))
        self.t_arr = np.full((nr, nc), np.nan)
        self.cellar = 0.0
        self.src_rc = self.g.xy_to_rc(*scenario.source_xy)
        if not self.g.inside(*self.src_rc):
            raise ValueError("Source outside the simulation domain.")
        oil = scenario.oil
        self.rho = oil["density_kgm3"]
        self.mu0 = oil["viscosity_mpas"] * 1e-3
        self.Cv = oil["visc_evap_c"]
        self.emax = oil["evap_max_pct"] / 100.0
        self.evap_a = fingas_coefficient(oil, scenario.temp_c)
        self.evap_form = oil["evap_form"]
        self.oil = oil
        self.cp = oil.get("cp_jkgk", 2000.0)
        sc_ = scenario
        self.t_mean = sc_.temp_c
        self.t_max = sc_.t_max_c if sc_.t_max_c is not None else sc_.temp_c + 6.0
        self.t_ground = sc_.ground_temp_c if sc_.ground_temp_c is not None else sc_.temp_c + 2.0
        self.t_src = sc_.discharge_temp_c
        self._met = None
        self.extra = []
        for es in (sc_.extra_sources or []):
            rc = self.g.xy_to_rc(es["x"], es["y"])
            if self.g.inside(*rc):
                e = dict(es)
                e["rc"] = rc
                self.extra.append(e)
        # bookkeeping (m3)
        self.mb = dict(released=0.0, spray=0.0, spray_out=0.0, contained=0.0, surface=0.0,
                       stranded=0.0, infiltrated=0.0, evaporated=0.0, exported=0.0,
                       on_water=0.0, recovered=0.0)
        self.kernel_cache = {}
        self.frames = []
        self.mb_rows = []
        self.first_channel = None
        self.first_export_h = None
        self.rain_cum_mm = 0.0
        self.runoff_cum_mm = 0.0
        self._rain_on = 0.0
        self._ue = None

    # ------------------------------------------------------------------
    def _source_volume(self, t, dt):
        """Volume (m3) released between t and t + dt."""
        sc = self.sc
        if sc.kind in ("BLOWOUT", "CONTINUOUS"):
            dur = sc.duration_h * 3600.0
            q = sc.rate_bbl_d * BBL_M3 / 86400.0
        else:  # TANK: inventory released over release_minutes
            dur = max(sc.release_minutes * 60.0, 1.0)
            q = sc.volume_m3 / dur
        return q * max(0.0, min(t + dt, dur) - t)

    def source_end_s(self, t=None):
        sc = self.sc
        if sc.kind in ("BLOWOUT", "CONTINUOUS"):
            ends = [sc.duration_h * 3600.0]
        else:
            ends = [max(sc.release_minutes * 60.0, 1.0)]
        for e in self.extra:
            s0 = e.get("start_h", 0.0) * 3600.0
            ends += [s0, s0 + max(e.get("release_min", 10.0) * 60.0, 1.0)]
        if t is None:
            return max(ends)
        nxt = [x for x in ends if x > t + 1e-6]
        return min(nxt) if nxt else max(ends)

    # ------------------------------------------------------------------
    def _meteo(self, t):
        """Air temperature, solar radiation, sky emissivity and wind at time t."""
        sc = self.sc
        hod = (sc.start_hour + t / 3600.0) % 24.0
        day = sc.doy + int((sc.start_hour + t / 3600.0) // 24)
        Ta = ht.air_temperature(self.t_mean, self.t_max, hod)
        S = ht.global_radiation(sc.lat, day, hod, sc.cloud, sc.lon)
        ea = ht.air_emissivity(Ta, sc.rh_pct, sc.cloud)
        w = sc.wind[min(int(t // 3600), len(sc.wind) - 1)] if sc.wind else (0, 2.0, 0)
        return dict(hod=hod, Ta=Ta, S=S, ea=ea, U=w[1])

    def _rain_intensity(self, t):
        """Rain intensity (mm/h) at time t (s) from the storm list."""
        th = t / 3600.0
        for st, du, depth in self.sc.rain:
            if st <= th < st + du and du > 0:
                # triangular hyetograph peaking at one third of the duration
                tp = du / 3.0
                x = th - st
                peak = 2.0 * depth / du
                return peak * (x / tp if x < tp else (du - x) / (du - tp))
        return 0.0

    def _runoff_rate(self, i_mmh, dt):
        """Incremental SCS-CN runoff (mm/h) for the current storm."""
        S = 25400.0 / self.sc.cn - 254.0
        Ia = 0.2 * S
        P0 = self.rain_cum_mm
        P1 = P0 + i_mmh * dt / 3600.0
        q = lambda P: (P - Ia) ** 2 / (P + 0.8 * S) if P > Ia else 0.0
        dq = q(P1) - q(P0)
        self.rain_cum_mm = P1
        self.runoff_cum_mm += dq
        return dq / (dt / 3600.0)

    def _kernel(self, hour):
        w = self.sc.wind[min(hour, len(self.sc.wind) - 1)]
        key = (round(w[0] / 5.0) * 5, round(w[1], 1), round(w[2]))
        if key not in self.kernel_cache:
            sx, sy = self.sc.source_xy
            k, out = fo.fallout_kernel(self.g, sx, sy, key[0], key[1], key[2], self.rho,
                                       self.sc.release_height_m, self.sc.d50_um, self.sc.gsd)
            self.kernel_cache[key] = (k, out)
        return self.kernel_cache[key]

    # ------------------------------------------------------------------
    def _add_oil(self, rr, cc, dv_m3, temp_c=None):
        dh = dv_m3 / self.A
        self.h[rr, cc] += dh
        self.O[rr, cc] += dh
        self.Hc[rr, cc] += dh * (self._met["Ta"] if temp_c is None else temp_c)

    def _sources(self, t, dt):
        sc = self.sc
        v = self._source_volume(t, dt)
        r, c = self.src_rc
        if v <= 0:
            vp = 0.0
            self._extra_only(t, dt)
            return
        self.mb["released"] += v
        if sc.kind == "BLOWOUT":
            vs = v * sc.spray_fraction
            vp = v - vs
            k, out = self._kernel(int(t // 3600))
            self.mb["spray"] += vs
            self.mb["spray_out"] += vs * out
            add = k * vs / self.A
            self.h += add
            self.O += add
            # fine droplets reach air temperature before landing
            self.Hc += add * self._met["Ta"]
            self.fall_cum += add
        else:
            vp = v
        # cellar / bund first
        room = max(sc.contained_m3 - self.cellar, 0.0)
        take = min(room, vp)
        self.cellar += take
        self.mb["contained"] = self.cellar
        vp -= take
        if vp > 0:
            self._add_oil(r, c, vp, self.t_src)
        # additional point sources (second tank, tanker, flowline leak ...)
        for e in self.extra:
            s0 = e.get("start_h", 0.0) * 3600.0
            dur = max(e.get("release_min", 10.0) * 60.0, 1.0)
            q = e["volume_m3"] / dur
            ve = q * max(0.0, min(t + dt, s0 + dur) - max(t, s0))
            if ve > 0:
                self.mb["released"] += ve
                self._add_oil(e["rc"][0], e["rc"][1], ve, e.get("temp_c"))

    def _extra_only(self, t, dt):
        for e in self.extra:
            s0 = e.get("start_h", 0.0) * 3600.0
            dur = max(e.get("release_min", 10.0) * 60.0, 1.0)
            q = e["volume_m3"] / dur
            ve = q * max(0.0, min(t + dt, s0 + dur) - max(t, s0))
            if ve > 0:
                self.mb["released"] += ve
                self._add_oil(e["rc"][0], e["rc"][1], ve, e.get("temp_c"))

    def _temperature(self, sl):
        """Oil temperature (C) in the window; ground temperature where there is no oil."""
        h = self.h[sl]
        if not self.sc.heat_mode:
            return np.full(h.shape, self.sc.temp_c)
        return np.where(h > 1e-9, self.Hc[sl] / np.maximum(h, 1e-12), self._met["Ta"])

    def _heat(self, dt, sl):
        if not self.sc.heat_mode:
            return
        h = self.h[sl]
        on = h > 1e-9
        if not on.any():
            return
        T = self._temperature(sl)
        m = self._met
        rho = ht.density(self.oil, T)
        Tn = ht.relax_temperature(T, h, dt, m["Ta"], self.t_ground, m["S"], m["U"], m["ea"],
                                  rho, self.cp, self.sc.h_ground)
        self.Hc[sl] = np.where(on, h * Tn, 0.0)

    def _rheology(self, sl):
        """Viscosity (Pa.s), yield stress (Pa) and density (kg/m3) in the window."""
        h = self.h[sl]
        O = self.O[sl]
        E = np.where(h > 1e-9, 1.0 - h / np.maximum(O, 1e-12), 0.0)
        T = self._temperature(sl)
        if self.sc.heat_mode:
            mu = ht.viscosity(self.oil, T, E)
            rho = ht.density(self.oil, T)
        else:
            mu = self.mu0 * np.exp(self.Cv * np.clip(E, 0, 1))
            rho = np.full(h.shape, self.rho)
        tau = ht.yield_stress(self.oil, T)
        return mu, tau, rho

    def _evaporate(self, dt, sl):
        h = self.h[sl]
        O = self.O[sl]
        on = h > 1e-9
        if not on.any():
            return
        E = np.zeros_like(h)
        E[on] = 1.0 - h[on] / O[on]
        E = np.clip(E, 0.0, self.emax)
        if self.sc.heat_mode:
            T = self._temperature(sl)
            d = self.oil["pct_distilled_180"]
            if self.evap_form == "ln":
                a = np.maximum(0.165 * d + 0.045 * (T - 15.0), 1e-4)
            else:
                a = np.maximum(0.0254 * d + 0.01 * (T - 15.0), 1e-4)
        else:
            a = self.evap_a
        Ep = E * 100.0
        dtm = dt / 60.0
        if self.evap_form == "ln":
            teq = np.exp(np.minimum(Ep / a, 50.0))
            teq = np.maximum(teq, 1.0)
            Enew = a * np.log(teq + dtm)
        else:
            teq = (Ep / a) ** 2
            Enew = a * np.sqrt(teq + dtm)
        Enew = np.minimum(Enew / 100.0, self.emax)
        Enew = np.where(on, Enew, 0.0)
        hnew = np.where(on, O * (1.0 - Enew), 0.0)
        hnew = np.minimum(hnew, h)
        loss = h - hnew
        self.mb["evaporated"] += float(loss.sum()) * self.A
        # stranded oil shrinks in proportion
        Bs = self.B[sl]
        ratio = np.where(on, hnew / np.maximum(h, 1e-12), 0.0)
        self.B[sl] = Bs * ratio
        self.Hc[sl] = self.Hc[sl] * ratio
        self.h[sl] = hnew


    def _infiltrate(self, dt, sl, mu):
        h = self.h[sl]
        on = h > 1e-9
        if not on.any():
            return
        # wetted share of the cell: rivulet width, growing to the full cell as oil ponds
        fw = self.fwet[sl]
        hr = np.maximum(self.hret_map[sl] / np.maximum(fw, 1e-6), 1e-6)
        fw = np.clip(h / hr, fw, 1.0)
        k = self.k_map[sl] * fw
        kr = (1.0 - self.sc.soil_sat) ** 3          # water already in the pores
        ko = k * (self.rho / 1000.0) * (MU_WATER / mu) * kr
        dtheta = self.por_map[sl] * (1.0 - self.sc.soil_sat)
        F = self.F[sl]
        f = ko * (1.0 + self.psi_map[sl] * dtheta / np.maximum(F, 1e-3))
        dI = np.where(on, np.minimum(f * dt, h), 0.0)
        if self._rain_on > 0:
            dI *= 0.5
        frac = np.where(on, 1.0 - dI / np.maximum(h, 1e-12), 1.0)
        self.h[sl] = h - dI
        self.O[sl] *= frac
        self.B[sl] *= frac
        self.Hc[sl] *= frac
        self.F[sl] = F + dI
        self.mb["infiltrated"] += float(dI.sum()) * self.A

    def _runoff_speed(self, re_mmh, sl):
        """Kinematic sheet-flow speed (m/s) at each cell for runoff excess re."""
        re = re_mmh / 3600.0 / 1000.0
        acc = self.T.area_km2[sl] * 1e6
        q = re * acc / self.dx
        S = np.maximum(self.T.d8_slope[sl], 1e-3)
        n = self.n_map[sl]
        d = (q * n / np.sqrt(S)) ** 0.6
        return np.where(d > 1e-6, q / np.maximum(d, 1e-6), 0.0)

    def _transport(self, dt, sl, rheo, re_mmh):
        """Move mobile oil; returns max speed for the next time step."""
        r0, r1, c0, c1 = sl[0].start, sl[0].stop, sl[1].start, sl[1].stop
        h = self.h[sl]
        B = self.B[sl]
        hret = self.hret_map[sl]
        if re_mmh > 0:
            hret = hret * (1.0 - self.sc.rain_flush)
        hret = np.maximum(hret, B)
        # channels carrying water: the oil floats, only bank-stranded oil stays behind
        chan = self.chan[sl]
        Q = self.Qbase[sl].copy()
        if re_mmh > 0:
            Q += re_mmh / 3.6e6 * self.T.area_km2[sl] * 1e6
        flowing = chan & (Q >= self.sc.q_min_flow)
        hret = np.where(flowing, B, hret)
        hm = np.maximum(h - hret, 0.0)
        if not (hm > 1e-7).any():
            return 0.0
        mu, tau, rho = rheo
        nr, nc = h.shape
        eta = self.zf[sl] + h
        pad = self._eta_padded(r0, r1, c0, c1)
        # cells on the outer border of the full grid may drain out
        dh = []
        smax = np.zeros_like(h)
        for k, (dr, dc) in enumerate(D8):
            en = pad[1 + dr:1 + dr + nr, 1 + dc:1 + dc + nc]
            d = np.maximum(eta - en, 0.0)
            dh.append(d)
            smax = np.maximum(smax, d / (self.dx * D8_DIST[k]))
        # --- overland speed
        u_lam = rho * G * hm ** 2 * smax / (3.0 * mu)
        u_man = hm ** (2.0 / 3.0) * np.sqrt(smax) / self.n_map[sl]
        u = np.minimum(u_lam, u_man)
        # Bingham plastic: the film moves only where the bed shear stress exceeds
        # the yield stress (gelled oil below its pour point stops)
        if np.any(tau > 0):
            tb = rho * G * hm * smax
            with np.errstate(divide="ignore", invalid="ignore"):
                r = np.where(tb > 0, tau / tb, np.inf)
            rr_ = np.minimum(r, 1.0)
            fb = np.where(r < 1.0, 1.0 - 1.5 * rr_ + 0.5 * rr_ ** 3, 0.0)
            u = u * np.clip(fb, 0.0, 1.0)
        if re_mmh > 0:
            u = np.maximum(u, self.sc.rain_carry * self._runoff_speed(re_mmh, sl) *
                           (hm > 0))
        # --- channels
        u_ch = np.where(flowing, 1.1 * self.sc.stream_kv * np.maximum(Q, 1e-6) ** 0.1, 0.0)
        # oil on open water spreads but does not 'flow' downhill
        wat = self.water[sl]
        u = np.where(wat, np.minimum(u, 0.02), u)
        frac = np.minimum(u * dt / self.dx, 0.9)
        frac_ch = np.minimum(u_ch * dt / self.dx, 0.95)
        # overland outflow volumes (thickness units), MFD weights
        out_tot = np.where(flowing, 0.0, hm * frac)
        wsum = np.zeros_like(h)
        ws = []
        for k in range(8):
            w = (dh[k] / (self.dx * D8_DIST[k])) ** 1.1
            ws.append(w)
            wsum += w
        has_out = wsum > 0
        out_tot = np.where(has_out, out_tot, 0.0)
        wsum = np.where(has_out, wsum, 1.0)
        # limit by half the head difference to the lowest neighbour
        lim = np.full_like(h, np.inf)
        for k in range(8):
            wk = ws[k] / wsum
            with np.errstate(divide="ignore", invalid="ignore"):
                cap = np.where(wk > 0, 0.5 * dh[k] / wk, np.inf)
            lim = np.minimum(lim, cap)
        out_tot = np.minimum(out_tot, lim)
        moved_frac = np.where(h > 0, out_tot / np.maximum(h, 1e-12), 0.0)
        dO = self.O[sl] * moved_frac
        dHc = self.Hc[sl] * moved_frac
        exported = 0.0
        H = self.h
        Ofull = self.O
        for k, (dr, dc) in enumerate(D8):
            wk = ws[k] / wsum
            mv = out_tot * wk
            mo = dO * wk
            mt = dHc * wk
            if not (mv > 0).any():
                continue
            rr0, rr1 = r0 + dr, r1 + dr
            cc0, cc1 = c0 + dc, c1 + dc
            # destination slice clipped to the grid
            sr0, sr1 = max(rr0, 0), min(rr1, self.nr)
            scc0, scc1 = max(cc0, 0), min(cc1, self.nc)
            ar0, ar1 = sr0 - rr0, nr - (rr1 - sr1)
            ac0, ac1 = scc0 - cc0, nc - (cc1 - scc1)
            exported += float(mv.sum() - mv[ar0:ar1, ac0:ac1].sum())
            H[sr0:sr1, scc0:scc1] += mv[ar0:ar1, ac0:ac1]
            Ofull[sr0:sr1, scc0:scc1] += mo[ar0:ar1, ac0:ac1]
            self.Hc[sr0:sr1, scc0:scc1] += mt[ar0:ar1, ac0:ac1]
        # subtract outflow from sources (after adding inflows elsewhere)
        H[sl] -= out_tot
        Ofull[sl] -= dO
        self.Hc[sl] -= dHc
        # --- flowing channel cells: D8 downstream with bank stranding
        if flowing.any():
            idx_r, idx_c = np.nonzero(flowing & (hm > 1e-9))
            if idx_r.size:
                hm_c = hm[idx_r, idx_c]
                fr = frac_ch[idx_r, idx_c]
                uu = u_ch[idx_r, idx_c]
                strand = hm_c * (1.0 - np.exp(-self.sc.strand_per_m * uu * dt))
                mv = (hm_c - strand) * fr
                gr = idx_r + r0
                gc = idx_c + c0
                hcell = self.h[gr, gc]
                ofr = np.where(hcell > 0, mv / np.maximum(hcell, 1e-12), 0.0)
                mo = self.O[gr, gc] * ofr
                mt = self.Hc[gr, gc] * ofr
                self.B[gr, gc] += strand
                self.mb["stranded"] += float(strand.sum()) * self.A
                self.h[gr, gc] -= mv
                self.O[gr, gc] -= mo
                self.Hc[gr, gc] -= mt
                rec = self.rec[gr, gc]
                ok = rec >= 0
                exported += float(mv[~ok].sum())
                rr, cc = np.divmod(rec[ok], self.nc)
                np.add.at(self.h, (rr, cc), mv[ok])
                np.add.at(self.O, (rr, cc), mo[ok])
                np.add.at(self.Hc, (rr, cc), mt[ok])
        self.mb["exported"] += exported * self.A
        np.maximum(self.h, 0.0, out=self.h)
        np.maximum(self.O, self.h, out=self.O)
        self.Hc[self.h <= 0] = 0.0
        um = float(np.max(np.where(hm > 1e-7, u, 0.0)))
        if flowing.any():
            um = max(um, float(np.max(np.where(flowing & (hm > 1e-7), u_ch, 0.0))))
        return um

    def _eta_padded(self, r0, r1, c0, c1):
        """Free surface (DEM + oil) of the window plus a one-cell ring taken
        from the neighbouring cells; at the outer edge of the domain the
        ring is extrapolated (odd reflection) so oil can leave the domain."""
        R0, R1 = r0 - 1, r1 + 1
        C0, C1 = c0 - 1, c1 + 1
        a0, a1 = max(R0, 0), min(R1, self.nr)
        b0, b1 = max(C0, 0), min(C1, self.nc)
        e = self.zf[a0:a1, b0:b1] + self.h[a0:a1, b0:b1]
        pr = (a0 - R0, R1 - a1)
        pc = (b0 - C0, C1 - b1)
        if any(pr) or any(pc):
            e = np.pad(e, (pr, pc), mode="reflect", reflect_type="odd")
        return e

    def _sinks(self):
        if self.sink is None:
            return
        v = float(self.h[self.sink].sum()) * self.A
        if v > 0:
            self.mb["recovered"] += v
            for a in (self.h, self.O, self.B, self.Hc):
                a[self.sink] = 0.0

    # ------------------------------------------------------------------
    def _active_slice(self, mobile=False):
        if mobile:
            # cells that can hold mobile oil (rain may flush part of the retention)
            on = self.h > self.hret_map * (1.0 - self.sc.rain_flush) + 1e-7
            on |= self.chan & (self.h > 1e-7)
        else:
            on = self.h > 1e-9
        r, c = self.src_rc
        rows = np.nonzero(on.any(axis=1))[0]
        cols = np.nonzero(on.any(axis=0))[0]
        if rows.size == 0:
            r0, r1, c0, c1 = r, r + 1, c, c + 1
        else:
            r0, r1, c0, c1 = rows[0], rows[-1] + 1, cols[0], cols[-1] + 1
            r0, r1 = min(r0, r), max(r1, r + 1)
            c0, c1 = min(c0, c), max(c1, c + 1)
        m = 3
        return (slice(max(r0 - m, 0), min(r1 + m, self.nr)),
                slice(max(c0 - m, 0), min(c1 + m, self.nc)))

    def _snapshot(self, t):
        thr = self.sc.thickness_threshold_m
        h = self.h
        newly = np.isnan(self.t_arr) & (h >= thr)
        self.t_arr[newly] = t / 3600.0
        np.maximum(self.hmax, h, out=self.hmax)
        sel = h >= 1e-6
        rr, cc = np.nonzero(sel)
        fr, fc = np.nonzero(self.fall_cum > 1e-7)
        full = (slice(0, self.nr), slice(0, self.nc))
        Tw = self._temperature(full)[sel]
        mu, tau, _ = self._rheology(full)
        mu_s, tau_s = mu[sel], tau[sel]
        if rr.size:
            tm = self.Tmax[rr, cc]
            self.Tmax[rr, cc] = np.where(np.isnan(tm), Tw, np.maximum(tm, Tw))
            tn = self.Tmin[rr, cc]
            self.Tmin[rr, cc] = np.where(np.isnan(tn), Tw, np.minimum(tn, Tw))
            self.mumax[rr, cc] = np.maximum(self.mumax[rr, cc], mu_s)
            self.taumax[rr, cc] = np.maximum(self.taumax[rr, cc], tau_s)
        met = self._met or {}
        self.frames.append(dict(
            t_h=t / 3600.0,
            rows=rr.astype(np.int32), cols=cc.astype(np.int32),
            h=h[sel].astype(np.float32),
            frows=fr.astype(np.int32), fcols=fc.astype(np.int32),
            fall=self.fall_cum[fr, fc].astype(np.float32),
            temp=Tw.astype(np.float32), visc=mu_s.astype(np.float32),
            tau=tau_s.astype(np.float32),
            t_air=met.get("Ta"), solar=met.get("S"), hod=met.get("hod"),
            rain_mmh=self._rain_on,
            wind=self.sc.wind[min(int(t // 3600), len(self.sc.wind) - 1)] if self.sc.wind else None))
        A = self.A
        on_w = float(h[self.water].sum()) * A
        surf = float(h.sum()) * A
        strd = float(self.B.sum()) * A
        row = dict(t_h=t / 3600.0, released=self.mb["released"],
                   contained=self.mb["contained"],
                   surface=surf - on_w - strd, on_water=on_w,
                   stranded=strd,
                   infiltrated=self.mb["infiltrated"], evaporated=self.mb["evaporated"],
                   spray_out=self.mb["spray_out"], exported=self.mb["exported"],
                   recovered=self.mb["recovered"],
                   oil_t_mean=float(np.average(Tw, weights=h[sel])) if rr.size else float("nan"),
                   gelled_m3=float(h[sel][Tw < self.oil["pour_point_c"]].sum()) * A if rr.size else 0.0,
                   t_air=met.get("Ta", float("nan")), solar_wm2=met.get("S", 0.0),
                   area_ha=float((h >= thr).sum()) * A / 1e4,
                   footprint_ha=float((~np.isnan(self.t_arr)).sum()) * A / 1e4,
                   rain_mmh=self._rain_on, rain_cum_mm=self.rain_cum_mm)
        row["balance_err"] = row["released"] - (row["contained"] + surf + row["infiltrated"] +
                                                row["evaporated"] + row["spray_out"] +
                                                row["exported"] + row["recovered"])
        self.mb_rows.append(row)
        if self.first_export_h is None and self.mb["exported"] > max(0.05, 0.001 * self.mb["released"]):
            self.first_export_h = t / 3600.0
        # first channel entry
        if self.first_channel is None:
            hit = self.chan & (h >= thr)
            if hit.any():
                rr, cc = np.nonzero(hit)
                i = int(np.argmin(self.T.area_km2[rr, cc] * -1))
                x, y = self.g.rc_to_xy(rr[i], cc[i])
                self.first_channel = dict(t_h=t / 3600.0, x=x, y=y,
                                          area_km2=float(self.T.area_km2[rr[i], cc[i]]))

    def _receptor_report(self):
        g = self.g
        thr = self.sc.thickness_threshold_m
        rad = self.sc.receptor_radius_m
        pen = self.F / self.pen_div                         # oil penetration depth (m)
        out = []
        nr_ = int(math.ceil(rad / g.dx))
        for rc in self.receptors:
            r, c = g.xy_to_rc(rc["x"], rc["y"])
            d = math.hypot(rc["x"] - self.sc.source_xy[0], rc["y"] - self.sc.source_xy[1])
            item = dict(rc)
            item.update(dist_km=d / 1000.0, arrival_h=None, max_mm=0.0, fallout_mm=0.0,
                        infil_mm=0.0, pen_m=0.0, gw=False, inside=g.inside(r, c))
            if g.inside(r, c):
                r0, r1 = max(r - nr_, 0), min(r + nr_ + 1, g.nrows)
                c0, c1 = max(c - nr_, 0), min(c + nr_ + 1, g.ncols)
                ta = self.t_arr[r0:r1, c0:c1]
                if np.isfinite(ta).any():
                    item["arrival_h"] = float(np.nanmin(ta))
                item["max_mm"] = float(self.hmax[r0:r1, c0:c1].max()) * 1000.0
                item["fallout_mm"] = float(self.fall_cum[r0:r1, c0:c1].max()) * 1000.0
                item["infil_mm"] = float(self.F[r0:r1, c0:c1].max()) * 1000.0
                item["pen_m"] = float(pen[r0:r1, c0:c1].max())
                item["gw"] = item["pen_m"] >= self.sc.dtw_m
            out.append(item)
        return out

    # ------------------------------------------------------------------
    def run(self):
        sc = self.sc
        t = 0.0
        t_end = sc.sim_hours * 3600.0
        next_frame = 0.0
        dt = 5.0
        umax = 0.0
        steps = 0
        tic = _time.time()
        self.log("Running '%s': %s, %.0f h, frame every %.0f min." %
                 (sc.name, sc.oil["label"], sc.sim_hours, sc.frame_minutes))
        if sc.heat_mode:
            self.log("Heat-transfer mode: source oil at %s, air %.1f-%.1f C, RH %.0f%%, cloud %.2f; "
                     "pour point %.1f C." % (
                         "%.1f C" % self.t_src if self.t_src is not None else "air temperature",
                         2 * self.t_mean - self.t_max, self.t_max, sc.rh_pct, sc.cloud,
                         self.oil["pour_point_c"]))
        self._met = self._meteo(0.0)
        while t <= t_end + 1e-6:
            if t >= next_frame - 1e-6:
                self._snapshot(t)
                next_frame += sc.frame_minutes * 60.0
                if self.fb is not None:
                    self.fb.setProgress(100.0 * t / max(t_end, 1.0))
                    if self.fb.isCanceled():
                        break
            # time step: Courant limit on the previous max speed
            dt = sc.dt_max if umax <= 0 else min(sc.dt_max, sc.cfl * self.dx / umax)
            dt = max(dt, 1.0)
            # do not step over a frame time or the end of the source
            dt = min(dt, max(next_frame - t, 1.0))
            se = self.source_end_s(t)
            if t < se:
                dt = min(dt, max(se - t, 1.0))
            self._met = self._meteo(t)
            i = self._rain_intensity(t)
            self._rain_on = i
            re = self._runoff_rate(i, dt) if i > 0 else 0.0
            self._sources(t, dt)
            sl = self._active_slice()
            self._heat(dt, sl)
            self._evaporate(dt, sl)
            mu, _, _ = self._rheology(sl)
            self._infiltrate(dt, sl, mu)
            slm = self._active_slice(mobile=True)
            umax = self._transport(dt, slm, self._rheology(slm), re)
            self._sinks()
            t += dt
            steps += 1
        self.log("Finished in %.1f s (%d steps)." % (_time.time() - tic, steps))
        res = Result()
        res.frames = self.frames
        res.mass_balance = self.mb_rows
        res.hmax = self.hmax
        res.t_arrival = self.t_arr
        res.infiltrated = self.F
        res.penetration = self.F / self.pen_div
        res.fallout = self.fall_cum
        res.stranded = self.B
        res.final_h = self.h
        res.temp_max = self.Tmax
        res.temp_min = self.Tmin
        res.visc_max = self.mumax
        res.yield_max = self.taumax
        res.receptors = self._receptor_report()
        res.first_channel = self.first_channel
        res.first_export_h = self.first_export_h
        res.grid = self.g
        res.scenario = sc
        res.steps = steps
        return res
