# -*- coding: utf-8 -*-
"""
Terrain preparation: DEM window, depression handling, D8 flow directions,
flow accumulation and the drainage (channel) network.

Only numpy and GDAL are used so that the code runs inside any QGIS 3 install.
"""
import logging
import heapq
import math
import numpy as np
from osgeo import gdal, ogr, osr

# D8 neighbour offsets (row, col) and distance factors
D8 = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]
D8_DIST = [1.0, math.sqrt(2), 1.0, math.sqrt(2), 1.0, math.sqrt(2), 1.0, math.sqrt(2)]


class Grid(object):
    """Raster window with its georeference."""

    def __init__(self, z, gt, wkt, nodata_mask):
        self.z = z
        self.gt = gt
        self.wkt = wkt
        self.mask = nodata_mask          # True where DEM is missing
        self.nrows, self.ncols = z.shape
        self.dx = abs(gt[1])
        self.dy = abs(gt[5])
        self.cell_area = self.dx * self.dy

    # ---- coordinates -------------------------------------------------
    def rc_to_xy(self, r, c):
        x = self.gt[0] + (c + 0.5) * self.gt[1]
        y = self.gt[3] + (r + 0.5) * self.gt[5]
        return x, y

    def xy_to_rc(self, x, y):
        c = int(math.floor((x - self.gt[0]) / self.gt[1]))
        r = int(math.floor((y - self.gt[3]) / self.gt[5]))
        return r, c

    def inside(self, r, c):
        return 0 <= r < self.nrows and 0 <= c < self.ncols

    def extent(self):
        x0 = self.gt[0]
        y1 = self.gt[3]
        x1 = x0 + self.ncols * self.gt[1]
        y0 = y1 + self.nrows * self.gt[5]
        return x0, y0, x1, y1

    def sub(self, r0, r1, c0, c1):
        gt = list(self.gt)
        gt[0] = self.gt[0] + c0 * self.gt[1]
        gt[3] = self.gt[3] + r0 * self.gt[5]
        return Grid(self.z[r0:r1, c0:c1].copy(), tuple(gt), self.wkt,
                    self.mask[r0:r1, c0:c1].copy())


def read_dem_window(path, cx, cy, half_width_m, target_cell=None):
    """Read a square window of the DEM centred on (cx, cy).

    If target_cell is given (metres) and is coarser than the native cell,
    the window is resampled (average) to that cell size.
    """
    ds = gdal.Open(path)
    if ds is None:
        raise IOError("Cannot open DEM: %s" % path)
    gt = ds.GetGeoTransform()
    band = ds.GetRasterBand(1)
    nd = band.GetNoDataValue()
    px = abs(gt[1])
    c0 = int(math.floor((cx - half_width_m - gt[0]) / gt[1]))
    c1 = int(math.ceil((cx + half_width_m - gt[0]) / gt[1]))
    r0 = int(math.floor((cy + half_width_m - gt[3]) / gt[5]))
    r1 = int(math.ceil((cy - half_width_m - gt[3]) / gt[5]))
    c0 = max(c0, 0)
    r0 = max(r0, 0)
    c1 = min(c1, ds.RasterXSize)
    r1 = min(r1, ds.RasterYSize)
    if c1 <= c0 or r1 <= r0:
        raise ValueError("The source point lies outside the DEM.")
    ngt = (gt[0] + c0 * gt[1], gt[1], 0.0, gt[3] + r0 * gt[5], 0.0, gt[5])
    if target_cell and target_cell > px * 1.01:
        f = target_cell / px
        ox = int(round((c1 - c0) / f))
        oy = int(round((r1 - r0) / f))
        arr = band.ReadAsArray(c0, r0, c1 - c0, r1 - r0, buf_xsize=ox, buf_ysize=oy,
                               resample_alg=gdal.GRIORA_Average).astype(np.float64)
        ngt = (ngt[0], gt[1] * (c1 - c0) / ox, 0.0, ngt[3], 0.0, gt[5] * (r1 - r0) / oy)
    else:
        arr = band.ReadAsArray(c0, r0, c1 - c0, r1 - r0).astype(np.float64)
    mask = ~np.isfinite(arr)
    if nd is not None:
        mask |= np.isclose(arr, nd) | (arr > 1e30) | (arr < -1e4)
    wkt = ds.GetProjection()
    ds = None
    if mask.any():
        # fill holes with the mean of valid neighbours (few cells only)
        arr[mask] = np.nan
        valid_mean = np.nanmean(arr)
        arr[mask] = valid_mean
    return Grid(arr, ngt, wkt, mask)


def gaussian_smooth(z, sigma_cells):
    """Separable Gaussian smoothing (numpy only) to suppress DEM speckle
    (canopy / radar noise) that would otherwise trap oil in false pits."""
    if not sigma_cells or sigma_cells <= 0:
        return z
    rad = int(math.ceil(3 * sigma_cells))
    x = np.arange(-rad, rad + 1)
    k = np.exp(-0.5 * (x / float(sigma_cells)) ** 2)
    k /= k.sum()
    p = np.pad(z, rad, mode="reflect")
    tmp = np.zeros_like(p)
    for i, w in enumerate(k):
        tmp += w * np.roll(p, i - rad, axis=0)
    out = np.zeros_like(p)
    for i, w in enumerate(k):
        out += w * np.roll(tmp, i - rad, axis=1)
    return out[rad:-rad, rad:-rad]


def priority_flood(z, epsilon=1e-4):
    """Barnes et al. (2014) priority-flood with epsilon: returns a
    depression-free surface in which every cell drains to the edge."""
    nr, nc = z.shape
    filled = z.copy()
    closed = np.zeros((nr, nc), dtype=bool)
    heap = []
    for r in range(nr):
        for c in (0, nc - 1):
            if not closed[r, c]:
                closed[r, c] = True
                heap.append((filled[r, c], r, c))
    for c in range(1, nc - 1):
        for r in (0, nr - 1):
            if not closed[r, c]:
                closed[r, c] = True
                heap.append((filled[r, c], r, c))
    heapq.heapify(heap)
    fz = filled.tolist()
    cl = closed.tolist()
    push = heapq.heappush
    pop = heapq.heappop
    while heap:
        h, r, c = pop(heap)
        for dr, dc in D8:
            rr = r + dr
            cc = c + dc
            if rr < 0 or rr >= nr or cc < 0 or cc >= nc:
                continue
            row = cl[rr]
            if row[cc]:
                continue
            row[cc] = True
            v = fz[rr][cc]
            if v <= h:
                v = h + epsilon
                fz[rr][cc] = v
            push(heap, (v, rr, cc))
    return np.array(fz, dtype=np.float64)


def d8_receivers(zf, dx):
    """Index of the steepest-descent receiver of every cell (-1 = outlet)."""
    nr, nc = zf.shape
    best = np.zeros((nr, nc))
    rec = -np.ones((nr, nc), dtype=np.int64)
    idx = np.arange(nr * nc).reshape(nr, nc)
    pad = np.pad(zf, 1, mode="constant", constant_values=np.inf)
    ipad = np.pad(idx, 1, mode="constant", constant_values=-1)
    for k, (dr, dc) in enumerate(D8):
        zn = pad[1 + dr:1 + dr + nr, 1 + dc:1 + dc + nc]
        inn = ipad[1 + dr:1 + dr + nr, 1 + dc:1 + dc + nc]
        s = (zf - zn) / (dx * D8_DIST[k])
        better = (s > best) & (inn >= 0)
        best = np.where(better, s, best)
        rec = np.where(better, inn, rec)
    # edge cells with no lower neighbour drain out of the window
    return rec.ravel(), best


def flow_accumulation(zf, rec):
    """Number of cells draining through each cell (including itself)."""
    n = zf.size
    order = np.argsort(-zf.ravel(), kind="stable").tolist()
    acc = [1.0] * n
    r = rec.tolist()
    for i in order:
        j = r[i]
        if j >= 0:
            acc[j] += acc[i]
    return np.array(acc).reshape(zf.shape)


def slope_mfd_surface(z, dx):
    """Max downhill gradient (m/m) on a surface, used for display/QA."""
    nr, nc = z.shape
    pad = np.pad(z, 1, mode="edge")
    smax = np.zeros_like(z)
    for k, (dr, dc) in enumerate(D8):
        zn = pad[1 + dr:1 + dr + nr, 1 + dc:1 + dc + nc]
        s = (z - zn) / (dx * D8_DIST[k])
        smax = np.maximum(smax, s)
    return smax


def hillshade(z, dx, azimuth=315.0, altitude=45.0, zfactor=1.0):
    gy, gx = np.gradient(z * zfactor, dx)
    slope = np.pi / 2.0 - np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(-gx, gy)
    az = math.radians(360.0 - azimuth + 90.0)
    alt = math.radians(altitude)
    hs = (math.sin(alt) * np.sin(slope) +
          math.cos(alt) * np.cos(slope) * np.cos(az - aspect))
    return np.clip(hs, 0, 1)


def relief_shade(z, dx, zfactor=None, smooth_m=45.0):
    """Multidirectional hillshade (azimuths 225, 270, 315, 360 deg, weighted as in
    GDAL -multidirectional) with vertical exaggeration, contrast-stretched to 0-1.
    Used for display only. On the flat Voltaian plain (about 180 m of relief across
    the block) a z-factor of about 5 is needed for the drainage pattern to show."""
    z = np.asarray(z, float)
    ok = np.isfinite(z) & (z > -1000)
    z = np.where(ok, z, np.mean(z[ok]) if ok.any() else 0.0)
    if zfactor is None:
        rng = np.percentile(z[ok], 99) - np.percentile(z[ok], 1) if ok.any() else 1.0
        zfactor = float(np.clip(600.0 / max(rng, 1.0), 1.0, 6.0))
    if smooth_m:
        try:
            from scipy.ndimage import gaussian_filter
            z = gaussian_filter(z, smooth_m / float(dx))   # suppress canopy / speckle noise
        except ImportError:
            logging.getLogger(__name__).debug("Optional operation failed in relief_shade; existing fallback retained.", exc_info=True)
    # rows run north to south, columns west to east
    dzdy, dzdx = np.gradient(z * zfactor, dx)
    dzdy = -dzdy                                   # positive = rising to the north
    slope = np.arctan(np.hypot(dzdx, dzdy))
    aspect = np.arctan2(-dzdx, -dzdy)              # direction the slope faces, from north, clockwise
    alt = math.radians(45.0)
    out = np.zeros_like(z)
    wsum = 0.0
    for az in (225.0, 270.0, 315.0, 360.0):
        a = math.radians(az)
        h = math.sin(alt) * np.cos(slope) + math.cos(alt) * np.sin(slope) * np.cos(a - aspect)
        w = np.sin(aspect - math.radians(az - 90.0)) ** 2 + 0.05
        out += np.clip(h, 0, 1) * w
        wsum = wsum + w
    out = out / wsum
    lo, hi = np.percentile(out, 1), np.percentile(out, 99.5)
    return np.clip((out - lo) / max(hi - lo, 1e-6), 0, 1)


def rasterize_layer(path_or_ds, layer_name, grid, where=None, burn=1, all_touched=True,
                    attribute=None, dtype=gdal.GDT_Float32):
    """Burn a vector layer onto the grid. Returns a numpy array (0 elsewhere)."""
    drv = gdal.GetDriverByName("MEM")
    out = drv.Create("", grid.ncols, grid.nrows, 1, dtype)
    out.SetGeoTransform(grid.gt)
    out.SetProjection(grid.wkt)
    out.GetRasterBand(1).Fill(0)
    src = path_or_ds if not isinstance(path_or_ds, str) else ogr.Open(path_or_ds)
    if src is None:
        return np.zeros((grid.nrows, grid.ncols))
    lyr = src.GetLayerByName(layer_name) if layer_name else src.GetLayer(0)
    if lyr is None:
        return np.zeros((grid.nrows, grid.ncols))
    if where:
        lyr.SetAttributeFilter(where)
    opts = []
    if all_touched:
        opts.append("ALL_TOUCHED=TRUE")
    if attribute:
        opts.append("ATTRIBUTE=%s" % attribute)
        gdal.RasterizeLayer(out, [1], lyr, options=opts)
    else:
        gdal.RasterizeLayer(out, [1], lyr, burn_values=[burn], options=opts)
    if where:
        lyr.SetAttributeFilter(None)
    return out.GetRasterBand(1).ReadAsArray()


class Terrain(object):
    """All terrain-derived arrays for one simulation domain."""

    def __init__(self, dem_path, cx, cy, half_width_m, cell_m=None,
                 pit_tolerance_m=3.0, channel_area_km2=1.0, margin_m=6000.0,
                 smooth_sigma_cells=1.0, burn=None, log=print):
        self.log = log
        log("Reading DEM window (%.1f km wide plus %.1f km catchment margin)..." %
            (2 * half_width_m / 1000.0, margin_m / 1000.0))
        big = read_dem_window(dem_path, cx, cy, half_width_m + margin_m, cell_m)
        self.z_original = None
        if smooth_sigma_cells:
            big.z = gaussian_smooth(big.z, smooth_sigma_cells)
        # stream burning: mapped watercourses are cut into the DEM so that the
        # flow directions follow them (Hellweger 1997 'AGREE'-style burn-in)
        self.burned = None
        z_unburned = big.z.copy()
        if burn:
            mask = np.zeros(big.z.shape, dtype=bool)
            for b in burn:
                arr = rasterize_layer(b["path"], b.get("layer"), big, where=b.get("where"))
                mask |= arr > 0
            if mask.any():
                big.z = big.z - np.where(mask, float(burn[0].get("depth_m", 3.0)), 0.0)
                log("Burned %d mapped-watercourse cells into the DEM." % int(mask.sum()))
            self.burned = mask
        # simulation window inside the big one
        dx = big.dx
        m = int(round(margin_m / dx))
        r0 = min(m, big.nrows // 4)
        c0 = min(m, big.ncols // 4)
        r1 = big.nrows - r0
        c1 = big.ncols - c0
        log("Filling depressions and tracing flow paths (%d x %d cells)..." %
            (big.nrows, big.ncols))
        zf_big = priority_flood(big.z)
        rec_big, _ = d8_receivers(zf_big, dx)
        acc_big = flow_accumulation(zf_big, rec_big)
        self.grid = big.sub(r0, r1, c0, c1)
        g = self.grid
        self.dx = g.dx
        self.z_raw = g.z
        zf = zf_big[r0:r1, c0:c1]
        self.acc_cells = acc_big[r0:r1, c0:c1]
        self.area_km2 = self.acc_cells * big.cell_area / 1e6
        # depressions deeper than the tolerance are kept as real ponds; they are
        # judged on the un-burned DEM so that burning does not create false pits
        if self.burned is not None and self.burned.any():
            zu = z_unburned[r0:r1, c0:c1]
            depth = priority_flood(z_unburned)[r0:r1, c0:c1] - zu
        else:
            depth = zf - g.z
        pond = depth > pit_tolerance_m
        self.z_flow = np.where(pond, g.z, zf)
        self.pond_depth = np.where(pond, depth, 0.0)
        self.z_filled = zf
        # D8 receivers inside the sim window
        rec, s = d8_receivers(zf, g.dx)
        self.rec = rec
        self.d8_slope = np.maximum(s, 1e-4)
        self.channel = self.area_km2 >= channel_area_km2
        self.channel_area_km2 = channel_area_km2
        self.slope = slope_mfd_surface(g.z, g.dx)
        self.hillshade = relief_shade(self.z_raw if getattr(self, 'z_raw', None) is not None else g.z, g.dx)
        log("Drainage network: %d channel cells (threshold %.2f km2)." %
            (int(self.channel.sum()), channel_area_km2))

    def burn_rivers(self, path, layer=None, min_area_km2=300.0, where=None):
        """Mark mapped rivers as channels; named major rivers get a large
        contributing area because their catchments extend beyond the DEM."""
        arr = rasterize_layer(path, layer, self.grid, where=where)
        hit = arr > 0
        self.channel |= hit
        self.area_km2 = np.where(hit, np.maximum(self.area_km2, min_area_km2), self.area_km2)
        return int(hit.sum())

    def stream_lines(self, min_area_km2=None):
        """Vectorise the channel network as polylines (list of coordinate lists)
        with the contributing area of the downstream vertex."""
        g = self.grid
        thr = min_area_km2 if min_area_km2 is not None else self.channel_area_km2
        ch = (self.area_km2 >= thr).ravel()
        rec = self.rec
        nr, nc = g.nrows, g.ncols
        n_up = np.zeros(nr * nc, dtype=np.int32)
        idx = np.nonzero(ch)[0]
        for i in idx:
            j = rec[i]
            if j >= 0 and ch[j]:
                n_up[j] += 1
        lines = []
        starts = [i for i in idx if n_up[i] != 1]
        seen = set()
        area = self.area_km2.ravel()
        for s in starts:
            # start a line from heads and junctions
            i = s
            pts = []
            while True:
                r, c = divmod(int(i), nc)
                pts.append(g.rc_to_xy(r, c))
                j = rec[i]
                if j < 0 or not ch[j]:
                    break
                if n_up[j] != 1:
                    r2, c2 = divmod(int(j), nc)
                    pts.append(g.rc_to_xy(r2, c2))
                    break
                if j in seen:
                    break
                seen.add(j)
                i = j
            if len(pts) > 1:
                lines.append((pts, float(area[i])))
        return lines
