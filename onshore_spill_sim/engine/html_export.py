# -*- coding: utf-8 -*-
"""
Self-contained interactive HTML player for a simulation run.

Leaflet 1.9.4 is embedded, so the page opens without internet access
(satellite / street basemaps need a connection; the terrain hillshade,
the spill frames and every vector layer are inside the file).
"""
import logging
import base64
import json
import math
import os
from osgeo import ogr, osr

HERE = os.path.dirname(os.path.abspath(__file__))


def _uri(path):
    with open(path, "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode("ascii")


def _to_ll(wkt):
    s = osr.SpatialReference()
    s.ImportFromWkt(wkt)
    s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    t = osr.SpatialReference()
    t.ImportFromEPSG(4326)
    t.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return osr.CoordinateTransformation(s, t)


def layer_geojson(path, layer, clip_ll=None, fields=None, simplify_m=0.0, where=None):
    """Read a vector layer, reproject to WGS84 and return a GeoJSON dict."""
    ds = ogr.Open(path)
    if ds is None:
        return None
    lyr = ds.GetLayerByName(layer) if layer else ds.GetLayer(0)
    if lyr is None:
        return None
    if where:
        lyr.SetAttributeFilter(where)
    src = lyr.GetSpatialRef()
    t = osr.SpatialReference()
    t.ImportFromEPSG(4326)
    t.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = None
    if src is not None:
        src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        tr = osr.CoordinateTransformation(src, t)
    clip = None
    if clip_ll:
        s, w, n, e = clip_ll
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for x, y in ((w, s), (e, s), (e, n), (w, n), (w, s)):
            ring.AddPoint_2D(x, y)
        clip = ogr.Geometry(ogr.wkbPolygon)
        clip.AddGeometry(ring)
    feats = []
    defn = lyr.GetLayerDefn()
    names = [defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())]
    use = [n for n in names if (fields is None or n in fields)]
    for f in lyr:
        g = f.GetGeometryRef()
        if g is None:
            continue
        g = g.Clone()
        g.FlattenTo2D()
        if simplify_m > 0:
            g = g.SimplifyPreserveTopology(simplify_m)
        if tr:
            g.Transform(tr)
        if clip is not None:
            try:
                if not g.Intersects(clip):
                    continue
            except Exception:
                logging.getLogger(__name__).debug("Optional operation failed in layer_geojson; existing fallback retained.", exc_info=True); continue
            if ogr.GT_Flatten(g.GetGeometryType()) not in (ogr.wkbPoint, ogr.wkbMultiPoint):
                try:
                    gi = g.Intersection(clip)
                except Exception:
                    try:
                        gi = g.MakeValid().Intersection(clip)
                    except Exception:
                        gi = g
                g = gi
                if g is None or g.IsEmpty():
                    continue
        props = {}
        for n in use:
            v = f.GetField(n)
            if isinstance(v, float):
                v = round(v, 3)
            props[n] = v
        gj = json.loads(g.ExportToJson())
        _round_coords(gj)
        feats.append({"type": "Feature", "properties": props, "geometry": gj})
    return {"type": "FeatureCollection", "features": feats}


def _round_coords(gj, nd=6):
    def r(c):
        if isinstance(c[0], (int, float)):
            return [round(c[0], nd), round(c[1], nd)]
        return [r(x) for x in c]
    if gj.get("coordinates") is not None:
        gj["coordinates"] = r(gj["coordinates"])
    if gj.get("geometries"):
        for g in gj["geometries"]:
            _round_coords(g, nd)


def sandbox_payload(terrain, source_xy, half_m, q_spec):
    """Compact D8 grid for the in-browser 'sandbox tanker' path tracer."""
    import numpy as np
    from .terrain import D8
    g = terrain.grid
    r, c = g.xy_to_rc(*source_xy)
    n = int(half_m / g.dx)
    r0, r1 = max(r - n, 0), min(r + n + 1, g.nrows)
    c0, c1 = max(c - n, 0), min(c + n + 1, g.ncols)
    nr, nc = r1 - r0, c1 - c0
    rec = terrain.rec.reshape(g.nrows, g.ncols)[r0:r1, c0:c1]
    rr, cc = np.divmod(np.where(rec >= 0, rec, 0), g.ncols)
    rows = np.arange(r0, r1)[:, None]
    cols = np.arange(c0, c1)[None, :]
    dr, dc = rr - rows, cc - cols
    code = np.full((nr, nc), 255, np.uint8)
    for k, (a, b) in enumerate(D8):
        code[(rec >= 0) & (dr == a) & (dc == b)] = k
    ch = terrain.channel[r0:r1, c0:c1]
    code = np.where((code != 255) & ch, code + 8, code).astype(np.uint8)
    z = terrain.z_raw[r0:r1, c0:c1]
    zmin = float(np.floor(z.min()))
    zq = np.clip(np.round((z - zmin) * 10), 0, 65535).astype("<u2")
    la = np.log10(np.maximum(terrain.area_km2[r0:r1, c0:c1], 1e-4))
    aq = np.clip(np.round((la + 4.0) * 30.0), 0, 255).astype(np.uint8)
    # affine map between grid metres and WGS84 (exact to < 1 m over a few km)
    tr = _to_ll(g.wkt)
    x0 = g.gt[0] + c0 * g.dx
    y1 = g.gt[3] - r0 * g.dx
    xs = np.linspace(x0, x0 + nc * g.dx, 5)
    ys = np.linspace(y1 - nr * g.dx, y1, 5)
    P, LL = [], []
    for x in xs:
        for y in ys:
            lo, la_, _ = tr.TransformPoint(float(x), float(y))
            P.append((1.0, x, y))
            LL.append((la_, lo))
    P, LL = np.array(P), np.array(LL)
    fwd = np.linalg.lstsq(P, LL, rcond=None)[0]           # [1,x,y] -> [lat, lon]
    Q = np.column_stack([np.ones(len(LL)), LL[:, 0], LL[:, 1]])
    inv = np.linalg.lstsq(Q, P[:, 1:], rcond=None)[0]     # [1,lat,lon] -> [x, y]
    return dict(nr=nr, nc=nc, x0=x0, y1=y1, dx=g.dx, zmin=zmin, q_spec=q_spec,
                code=base64.b64encode(code.tobytes()).decode("ascii"),
                z=base64.b64encode(zq.tobytes()).decode("ascii"),
                a=base64.b64encode(aq.tobytes()).decode("ascii"),
                fwd=[[round(v, 12) for v in row] for row in fwd.tolist()],
                inv=[[round(v, 6) for v in row] for row in inv.tolist()])


def build_html(path, payload):
    with open(os.path.join(HERE, "web", "leaflet.js"), encoding="utf-8") as f:
        ljs = f.read()
    with open(os.path.join(HERE, "web", "leaflet.css"), encoding="utf-8") as f:
        lcss = f.read()
    html = TEMPLATE.replace("/*LEAFLET_CSS*/", lcss)
    html = html.replace("/*LEAFLET_JS*/", ljs)
    html = html.replace("__TITLE__", payload["meta"]["title"])
    html = html.replace("/*PAYLOAD*/", "const DATA = " + json.dumps(payload, separators=(",", ":")) + ";\n" + BASE_JS)
    html = html.replace("</style>", CAPTION_CSS + "</style>", 1)
    with open(path, "w", encoding="utf-8") as f:
        f.write(html)
    return path


BASE_JS = r"""function makeBase(hillUri, hillB){
 const esri=(svc,o)=>L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/'+svc+'/MapServer/tile/{z}/{y}/{x}',Object.assign({maxZoom:19,attribution:'Tiles &copy; Esri'},o||{}));
 const hb=[[hillB[0],hillB[1]],[hillB[2],hillB[3]]];
 const b={
  'Satellite imagery (Esri)': esri('World_Imagery'),
  'Satellite hybrid - imagery with roads and places': L.layerGroup([esri('World_Imagery'), esri('Reference/World_Transportation',{opacity:0.9}), esri('Reference/World_Boundaries_and_Places')]),
  'OpenStreetMap (CARTO Voyager)': L.tileLayer('https://{s}.basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}{r}.png',{maxZoom:20,subdomains:'abcd',attribution:'&copy; OpenStreetMap contributors &copy; CARTO'}),
  'Street map (Esri)': esri('World_Street_Map'),
  'Topographic map (Esri)': esri('World_Topo_Map'),
  'Terrain - world hillshade (Esri)': L.layerGroup([esri('Elevation/World_Hillshade'), esri('Reference/World_Boundaries_and_Places',{opacity:0.85})]),
  'Terrain - relief from the model DEM (offline)': L.layerGroup([L.imageOverlay(hillUri,hb,{opacity:1})])
 };
 return b;
}
function reliefOverlay(map, hillUri, hillB){
 if(!map.getPane('relief')){ const p=map.createPane('relief'); p.style.zIndex=250; p.style.mixBlendMode='multiply'; p.style.filter='grayscale(1) contrast(1.3) brightness(1.15)'; p.style.pointerEvents='none'; }
 return L.imageOverlay(hillUri,[[hillB[0],hillB[1]],[hillB[2],hillB[3]]],{opacity:0.6,pane:'relief',interactive:false});
}
function baseFallback(map, base, first, offline){
 let te=0,sw=false; const watch=l=>{ (l.eachLayer?(f=>l.eachLayer(f)):(f=>f(l)))(t=>{ if(t.on) t.on('tileerror',()=>{te++; if(te>6&&!sw){sw=true; map.removeLayer(base[first]); base[offline].addTo(map);}}); }); };
 watch(base[first]); base[first].addTo(map);
}
function subtitleBox(map){
 const C=L.Control.extend({options:{position:'bottomleft'},onAdd:function(){const d=L.DomUtil.create('div','captionbox');
  d.innerHTML='<div class="cap-h"><b id="cap_t"></b><span id="cap_c"></span><a href="#" id="cap_x" title="Hide subtitles">&times;</a></div><div id="cap_b"></div>';
  L.DomEvent.disableClickPropagation(d); return d;}});
 new C().addTo(map);
 document.getElementById('cap_x').onclick=e=>{e.preventDefault(); document.querySelector('.captionbox').classList.toggle('min');};
}
function setCaption(c){ const t=document.getElementById('cap_t'); if(!t) return;
 if(!c){ document.querySelector('.captionbox').style.display='none'; return; }
 document.querySelector('.captionbox').style.display='block';
 t.textContent=c.title||''; document.getElementById('cap_c').textContent=c.clock?('  '+c.clock):''; document.getElementById('cap_b').textContent=c.text||''; }
"""
CAPTION_CSS = r""".captionbox{background:rgba(16,28,22,.86);color:#fff;border-radius:8px;padding:8px 12px 9px;max-width:min(760px,calc(100vw - 120px));
font-size:13.5px;line-height:1.45;box-shadow:0 2px 10px rgba(0,0,0,.35);margin-bottom:26px!important}
.captionbox .cap-h{display:flex;gap:8px;align-items:baseline;margin-bottom:2px}
.captionbox b{color:#ffd43b;font-size:14px}
.captionbox #cap_c{color:#c7d3cc;font-size:12px;flex:1}
.captionbox a{color:#fff;text-decoration:none;font-size:18px;line-height:1;padding:0 2px}
.captionbox.min #cap_b{display:none}
.relief-mult{mix-blend-mode:multiply}
"""

TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
/*LEAFLET_CSS*/
:root{--ink:#14211b;--muted:#5b6b63;--line:#d9e0dc;--panel:#ffffff;--bg:#eef2ef;--brand:#0b5d3b;
--brand2:#c9a227;--warn:#b42318;--ok:#1a7f4b;--chip:#f3f6f4;}
*{box-sizing:border-box}
html,body{margin:0;height:100%;font-family:"Segoe UI",Calibri,Arial,sans-serif;color:var(--ink);background:var(--bg)}
#app{display:flex;height:100vh;width:100vw}
#map{flex:1;min-width:0;height:100%}
#side{width:390px;max-width:100%;height:100%;overflow-y:auto;background:var(--panel);border-left:1px solid var(--line);
padding:14px 16px 24px}
h1{font-size:17px;margin:0 0 2px;line-height:1.25;color:var(--brand)}
.sub{font-size:12px;color:var(--muted);margin-bottom:10px}
.card{border:1px solid var(--line);border-radius:8px;padding:10px 12px;margin:10px 0;background:#fff}
.card h2{font-size:12px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);margin:0 0 8px}
.clock{font-size:26px;font-weight:700;font-variant-numeric:tabular-nums;color:var(--ink)}
.clock small{font-size:13px;font-weight:500;color:var(--muted);display:block}
.controls{display:flex;gap:6px;align-items:center;margin-top:8px;flex-wrap:wrap}
button{border:1px solid var(--line);background:var(--chip);border-radius:6px;padding:6px 10px;font-size:13px;cursor:pointer;color:var(--ink)}
button.primary{background:var(--brand);color:#fff;border-color:var(--brand)}
select{border:1px solid var(--line);border-radius:6px;padding:5px;font-size:13px}
input[type=range]{width:100%;margin-top:8px;accent-color:var(--brand)}
.kpis{display:grid;grid-template-columns:1fr 1fr;gap:8px}
.kpi{background:var(--chip);border-radius:6px;padding:7px 9px}
.kpi b{display:block;font-size:17px;font-variant-numeric:tabular-nums}
.kpi span{font-size:11px;color:var(--muted)}
.row{display:flex;justify-content:space-between;font-size:12.5px;padding:3px 0;border-bottom:1px dashed #edf1ee}
.row:last-child{border-bottom:0}
.row span:first-child{color:var(--muted)}
.wx{display:flex;gap:12px;align-items:center}
#compass{width:70px;height:70px;flex:none}
.rain{font-size:12.5px}
.legend-bar{height:12px;border-radius:3px;margin:4px 0 2px}
.legend-scale{display:flex;justify-content:space-between;font-size:10.5px;color:var(--muted)}
.rlist{max-height:260px;overflow:auto;font-size:12.5px}
.ritem{display:flex;justify-content:space-between;gap:8px;padding:4px 2px;border-bottom:1px solid #f0f3f1;cursor:pointer}
.ritem.hit{color:var(--warn);font-weight:600}
.ritem .t{font-variant-numeric:tabular-nums;white-space:nowrap}
.tag{font-size:10.5px;color:var(--muted);font-weight:400}
.note{font-size:11.5px;color:var(--muted);line-height:1.4}
#chart{width:100%;height:150px}
.leaflet-container{background:#dfe6e1}
.maplabel{background:rgba(255,255,255,.85);border:0;box-shadow:none;font-size:11px;padding:1px 4px;color:#222}
.raincurtain{position:absolute;inset:0;pointer-events:none;z-index:450;opacity:0;transition:opacity .4s;
background:repeating-linear-gradient(105deg,rgba(120,160,255,.0) 0 7px,rgba(120,160,255,.22) 7px 8px);}
.badge{position:absolute;top:10px;left:50%;transform:translateX(-50%);z-index:800;background:rgba(20,33,27,.85);
color:#fff;padding:6px 12px;border-radius:16px;font-size:13px;font-weight:600;pointer-events:none}
.hz{display:inline-block;font-size:10px;font-weight:700;padding:1px 5px;border-radius:8px;color:#fff;margin-left:4px;white-space:nowrap}
#links{width:100%;margin:0 0 8px}
.sbx{display:none}
.sbx.on{display:block}
.tankicon{background:#ff8c00;border:2px solid #222;border-radius:4px;width:18px!important;height:18px!important;margin:-9px 0 0 -9px!important;box-shadow:0 1px 4px rgba(0,0,0,.5)}
@media (max-width:820px){#app{flex-direction:column}#map{height:58vh;flex:none}#side{width:100%;height:auto;border-left:0;border-top:1px solid var(--line)}}
</style>
</head>
<body>
<div id="app">
 <div id="map"><div class="raincurtain" id="curtain"></div><div class="badge" id="badge"></div></div>
 <aside id="side">
  <select id="links" style="display:none"></select>
  <h1 id="title"></h1>
  <div class="sub" id="subtitle"></div>
  <div class="card">
   <div class="clock" id="clock">T+00:00<small id="date"></small></div>
   <input type="range" id="slider" min="0" value="0" step="1">
   <div class="controls">
    <button id="back" title="Previous frame">&#9664;&#9664;</button>
    <button id="play" class="primary">Play</button>
    <button id="fwd" title="Next frame">&#9654;&#9654;</button>
    <select id="speed" title="Playback speed">
     <option value="900">Slow</option><option value="450" selected>Normal</option><option value="180">Fast</option><option value="70">Very fast</option>
    </select>
    <label style="font-size:12px"><input type="checkbox" id="loop" checked> Loop</label>
   </div>
  </div>
  <div class="card">
   <h2>Now</h2>
   <div class="kpis">
    <div class="kpi"><b id="k_rel">0</b><span>released (m&sup3;)</span></div>
    <div class="kpi"><b id="k_area">0</b><span>oiled ground now (ha)</span></div>
    <div class="kpi"><b id="k_evap">0</b><span>evaporated (m&sup3;)</span></div>
    <div class="kpi"><b id="k_inf">0</b><span>soaked into soil (m&sup3;)</span></div>
   </div>
   <div class="wx" style="margin-top:10px">
    <svg id="compass" viewBox="-50 -50 100 100"></svg>
    <div class="rain"><div id="wind_txt"></div><div id="rain_txt"></div><div id="stream_txt"></div></div>
   </div>
  </div>
  <div class="card">
   <h2>Where the oil goes (m&sup3;)</h2>
   <svg id="chart" viewBox="0 0 360 150" preserveAspectRatio="none"></svg>
   <div id="chart_legend" class="note"></div>
  </div>
  <div class="card">
   <h2>Receptors</h2>
   <div class="note" style="margin-bottom:6px">Map markers: <span style="color:#d73a49">&#9679;</span> oil has arrived &nbsp; <span style="color:#c9a227">&#9679;</span> spray fallout only &nbsp; <span style="color:#888">&#9675;</span> not reached. Click a name to zoom to it.</div>
   <div class="rlist" id="rlist"></div>
  </div>
  <div class="card sbx" id="sbx">
   <h2>Sandbox tanker (quick path trace)</h2>
   <div class="note" style="margin-bottom:6px">Drag the orange marker (or click the map) to see where a spill at that point would run, following the steepest descent of the DEM. Indicative only - run the full model in QGIS (movable tanker tool) for volumes, spreading, soaking-in and timing.</div>
   <div id="sbx_out" style="font-size:12.5px"></div>
  </div>
  <div class="card">
   <h2>Legend</h2>
   <div id="fieldsel" style="display:none;margin-bottom:8px"><span style="font-size:12px">Extra animated layer </span><select id="field"></select>
    <div style="font-size:12px;margin-top:6px" id="lg_f_t"></div><div class="legend-bar" id="lg_f"></div><div class="legend-scale" id="lg_f_s"></div></div>
   <div style="font-size:12px">Surface oil thickness (mm)</div>
   <div class="legend-bar" id="lg_oil"></div><div class="legend-scale" id="lg_oil_s"></div>
   <div style="font-size:12px;margin-top:8px">Well spray fallout, cumulative (mm)</div>
   <div class="legend-bar" id="lg_fall"></div><div class="legend-scale" id="lg_fall_s"></div>
  </div>
  <div class="card">
   <h2>Scenario</h2>
   <div id="facts"></div>
  </div>
  <div class="card">
   <h2>Basis and limitations</h2>
   <div class="note" id="notes"></div>
  </div>
 </aside>
</div>
<script>
/*LEAFLET_JS*/
</script>
<script>
/*PAYLOAD*/
(function(){
const M = DATA.meta, F = DATA.frames, MB = DATA.mb;
document.getElementById('title').textContent = M.title;
document.getElementById('subtitle').textContent = M.subtitle;
if((DATA.links||[]).length){ const ls=document.getElementById('links'); ls.style.display='block';
  ls.innerHTML='<option value="">Other scenarios...</option>'+DATA.links.map(l=>'<option value="'+l.href+'">'+l.name+'</option>').join('');
  ls.onchange=()=>{ if(ls.value) window.location.href=ls.value; }; }
const HZ={'500 ppm+':'#3b0a45','IDLH':'#9c36b5','AEGL-3':'#e03131','AEGL-2':'#f08c00','AEGL-1':'#2f9e44','Odour only':'#74c0fc','Not reached':'#adb5bd'};
const map = L.map('map',{zoomControl:true, preferCanvas:true});
const base = makeBase(DATA.hill.uri, DATA.hill.bounds);
baseFallback(map, base, 'Satellite imagery (Esri)', 'Terrain - relief from the model DEM (offline)');
subtitleBox(map);
const B = [[DATA.bounds[0],DATA.bounds[1]],[DATA.bounds[2],DATA.bounds[3]]];
const FOCUS = DATA.focus ? [[DATA.focus[0],DATA.focus[1]],[DATA.focus[2],DATA.focus[3]]] : B;
map.fitBounds(FOCUS,{padding:[20,20]});
const ZoomCtl = L.Control.extend({options:{position:'topleft'}, onAdd:function(){ const d=L.DomUtil.create('div','leaflet-bar');
  d.innerHTML='<a href="#" title="Zoom to the spill" style="width:auto;padding:0 6px;font-size:12px">Spill</a><a href="#" title="Whole model area" style="width:auto;padding:0 6px;font-size:12px">Area</a>';
  const a=d.querySelectorAll('a'); a[0].onclick=e=>{e.preventDefault();map.fitBounds(FOCUS,{padding:[20,20]});};
  a[1].onclick=e=>{e.preventDefault();map.fitBounds([[DATA.hill.bounds[0],DATA.hill.bounds[1]],[DATA.hill.bounds[2],DATA.hill.bounds[3]]]);};
  L.DomEvent.disableClickPropagation(d); return d; }});
new ZoomCtl().addTo(map);
map.createPane('vec'); map.getPane('vec').style.zIndex = 410;
map.createPane('fall'); map.getPane('fall').style.zIndex = 415;
map.createPane('oil'); map.getPane('oil').style.zIndex = 420;
map.createPane('pts'); map.getPane('pts').style.zIndex = 440;
// preload frames
const oilImgs = F.map(f=>{const i=new Image(); i.src=f.oil; return i;});
const fallImgs = F.map(f=>{const i=new Image(); i.src=f.fall; return i;});
const oilLayer = L.imageOverlay(F[0].oil, B, {pane:'oil', opacity:0.95, interactive:false});
const fallLayer = L.imageOverlay(F[0].fall, B, {pane:'fall', opacity:0.9, interactive:false});
const overlays = {'Terrain relief shading (DEM, over any basemap)': reliefOverlay(map, DATA.hill.uri, DATA.hill.bounds)};
map.createPane('fld'); map.getPane('fld').style.zIndex = 425;
const FL = DATA.fields||[];
let fieldLayer=null, fieldIdx=-1, fieldImgs=[];
if (FL.length){
  document.getElementById('fieldsel').style.display='block';
  const fs=document.getElementById('field');
  fs.innerHTML='<option value="-1">none</option>'+FL.map((f,i)=>'<option value="'+i+'">'+f.name+'</option>').join('');
  fieldImgs = FL.map(f=>f.frames.map(u=>{const i=new Image(); i.src=u; return i;}));
  fs.onchange=()=>{ fieldIdx=+fs.value; if(fieldLayer){ map.removeLayer(fieldLayer); fieldLayer=null; }
    const lt=document.getElementById('lg_f_t'), lb=document.getElementById('lg_f'), lsc=document.getElementById('lg_f_s');
    if(fieldIdx>=0){ const f=FL[fieldIdx]; fieldLayer=L.imageOverlay(f.frames[0],[[f.bounds[0],f.bounds[1]],[f.bounds[2],f.bounds[3]]],{pane:'fld',opacity:0.9,interactive:false}).addTo(map);
      lt.textContent=f.name+' ('+f.unit+')'; lb.style.background=grad(f.ramp); lsc.innerHTML=f.ramp.map(s=>'<span>'+s[0]+'</span>').join(''); }
    else { lt.textContent=''; lb.style.background='none'; lsc.innerHTML=''; }
    show(idx); };
}
function fieldFrame(f, t){ let k=0,best=1e9; f.t.forEach((v,i)=>{ const d=Math.abs(v-t); if(d<best){best=d;k=i;} }); return k; }
overlays['Oil on the ground (animated)'] = oilLayer.addTo(map);
overlays['Well spray fallout (animated)'] = fallLayer;
if (M.kind === 'BLOWOUT') fallLayer.addTo(map);
// vector layers
function styleFn(s){ return function(){ return {color:s.color||'#333', weight:s.weight||1.5, opacity:s.opacity==null?1:s.opacity,
  fillColor:s.fill||s.color||'#333', fillOpacity:s.fillOpacity==null?0.15:s.fillOpacity, dashArray:s.dash||null}; }; }
(DATA.layers||[]).forEach(Ly=>{
  const s = Ly.style||{};
  let lay;
  if (Ly.kind === 'point'){
    lay = L.geoJSON(Ly.geojson,{pane:'pts', pointToLayer:(f,ll)=>L.circleMarker(ll,{radius:s.radius||4,color:s.color||'#222',weight:1,
      fillColor:s.fill||s.color||'#fff',fillOpacity:s.fillOpacity==null?0.9:s.fillOpacity}),
      onEachFeature:(f,l)=>{ const p=f.properties||{}; const nm=p[s.label]||''; if(nm){ l.bindTooltip(String(nm),{permanent:!!s.permanent,direction:'right',className:'maplabel',offset:[4,0]}); }
        l.bindPopup('<b>'+(nm||Ly.name)+'</b><br>'+Ly.name+(s.popup?('<br>'+s.popup.map(k=>k+': '+(p[k]==null?'-':p[k])).join('<br>')):'')); }});
  } else {
    lay = L.geoJSON(Ly.geojson,{pane:'vec', style:styleFn(s), onEachFeature:(f,l)=>{ const p=f.properties||{}; const nm=p[s.label]||Ly.name; l.bindTooltip(String(nm),{sticky:true}); }});
  }
  overlays[Ly.name] = lay;
  if (Ly.visible) lay.addTo(map);
});
// source marker
const src = L.circleMarker([M.source[0],M.source[1]],{pane:'pts',radius:7,color:'#000',weight:2,fillColor:'#ffd400',fillOpacity:1}).bindTooltip(M.source_label,{permanent:true,direction:'left',className:'maplabel',offset:[-6,0]});
overlays['Spill source'] = src.addTo(map);
// receptors
const recLayer = L.layerGroup();
const recMarkers = [];
DATA.receptors.forEach((r,i)=>{
  const m = L.circleMarker([r.lat,r.lon],{pane:'pts',radius:5,color:'#333',weight:1,fillColor:'#ffffff',fillOpacity:1});
  m.bindPopup('<b>'+r.name+'</b><br>'+r.type+'<br>'+r.dist.toFixed(2)+' km from source<br>Status: '+r.status+
   (r.arrival!=null?('<br>Oil arrives: T+'+fmtH(r.arrival)):'')+'<br>Max thickness: '+r.max_mm.toFixed(2)+' mm'+
   (r.fall_mm>0?('<br>Spray fallout: '+r.fall_mm.toFixed(3)+' mm'):'')+(r.gw?'<br><b>Oil front reaches the water table</b>':'')+(r.h2s?'<br>H2S (leaking gas): <span class="hz" style="background:'+(HZ[r.h2s]||'#888')+'">'+r.h2s+'</span>':''));
  m.addTo(recLayer); recMarkers.push(m);
});
overlays['Receptor check points'] = recLayer.addTo(map);
L.control.layers(base, overlays, {collapsed: window.innerWidth<1600}).addTo(map);
L.control.scale({imperial:false}).addTo(map);
// legends
function grad(r){ return 'linear-gradient(90deg,'+r.map(s=>'rgba('+s[1]+','+s[2]+','+s[3]+','+(Math.max(s[4],140)/255).toFixed(2)+')').join(',')+')'; }
document.getElementById('lg_oil').style.background = grad(DATA.ramps.oil);
document.getElementById('lg_fall').style.background = grad(DATA.ramps.fall);
document.getElementById('lg_oil_s').innerHTML = DATA.ramps.oil.filter((s,i)=>i%2==0||i==DATA.ramps.oil.length-1).map(s=>'<span>'+s[0]+'</span>').join('');
document.getElementById('lg_fall_s').innerHTML = DATA.ramps.fall.map(s=>'<span>'+s[0]+'</span>').join('');
document.getElementById('facts').innerHTML = M.facts.map(f=>'<div class="row"><span>'+f[0]+'</span><span>'+f[1]+'</span></div>').join('');
document.getElementById('notes').innerHTML = M.notes.map(n=>'<p style="margin:0 0 6px">'+n+'</p>').join('');
// chart
const series = [
 {k:'evaporated',c:'#9aa5ad',n:'Evaporated'},
 {k:'infiltrated',c:'#8c5a2b',n:'Soaked into soil'},
 {k:'surface',c:'#1f6feb',n:'On the ground'},
 {k:'stranded',c:'#2aa198',n:'Stranded on banks'},
 {k:'on_water',c:'#6f42c1',n:'On open water'},
 {k:'contained',c:'#c9a227',n:'Held in cellar / bund'},
 {k:'left',c:'#d73a49',n:'Left the model area'}];
MB.forEach(r=>{ r.left = (r.exported||0)+(r.spray_out||0); });
const svg = document.getElementById('chart');
const W=360,H=150,PL=34,PB=16,PT=6,PR=6;
const tmax = MB[MB.length-1].t_h, vmax = Math.max(1, ...MB.map(r=>r.released))*1.05;
const X = t=>PL+(W-PL-PR)*t/Math.max(tmax,1e-6), Y = v=>H-PB-(H-PB-PT)*v/vmax;
let paths=''; let lower = MB.map(()=>0);
series.forEach(s=>{
  const up = MB.map((r,i)=>lower[i]+Math.max(r[s.k]||0,0));
  let d='M'+X(MB[0].t_h)+','+Y(lower[0]);
  MB.forEach((r,i)=>{ d+=' L'+X(r.t_h).toFixed(1)+','+Y(up[i]).toFixed(1); });
  for(let i=MB.length-1;i>=0;i--){ d+=' L'+X(MB[i].t_h).toFixed(1)+','+Y(lower[i]).toFixed(1); }
  paths+='<path d="'+d+' Z" fill="'+s.c+'" fill-opacity="0.85" stroke="none"/>';
  lower = up;
});
let rel='M'; MB.forEach((r,i)=>{ rel+=(i?' L':'')+X(r.t_h).toFixed(1)+','+Y(r.released).toFixed(1); });
const ticks=[0,0.25,0.5,0.75,1].map(f=>f*tmax);
let ax='<line x1="'+PL+'" y1="'+(H-PB)+'" x2="'+(W-PR)+'" y2="'+(H-PB)+'" stroke="#999" stroke-width="0.6"/>';
ticks.forEach(t=>{ ax+='<text x="'+X(t)+'" y="'+(H-3)+'" font-size="9" text-anchor="middle" fill="#666">'+Math.round(t)+' h</text>'; });
[0,0.5,1].forEach(f=>{ const v=f*vmax/1.05; ax+='<text x="'+(PL-3)+'" y="'+(Y(v)+3)+'" font-size="9" text-anchor="end" fill="#666">'+Math.round(v)+'</text><line x1="'+PL+'" x2="'+(W-PR)+'" y1="'+Y(v)+'" y2="'+Y(v)+'" stroke="#e5e9e6" stroke-width="0.6"/>'; });
svg.innerHTML = ax+paths+'<path d="'+rel+'" fill="none" stroke="#111" stroke-width="1.4" stroke-dasharray="4 2"/><line id="cursor" x1="0" x2="0" y1="'+PT+'" y2="'+(H-PB)+'" stroke="#111" stroke-width="1"/>';
document.getElementById('chart_legend').innerHTML = series.map(s=>'<span style="white-space:nowrap;margin-right:8px"><span style="display:inline-block;width:9px;height:9px;background:'+s.c+';border-radius:2px;margin-right:3px"></span>'+s.n+'</span>').join('')+'<span style="white-space:nowrap">- - released</span>';
// receptors list
const rl = document.getElementById('rlist');
const order = DATA.receptors.map((r,i)=>i).sort((a,b)=>{const A=DATA.receptors[a],Bb=DATA.receptors[b];
  if((A.arrival==null)!==(Bb.arrival==null)) return A.arrival==null?1:-1; if(A.arrival!=null) return A.arrival-Bb.arrival; return A.dist-Bb.dist;});
rl.innerHTML = order.slice(0,60).map(i=>{const r=DATA.receptors[i]; return '<div class="ritem" data-i="'+i+'"><span>'+r.name+' <span class="tag">'+r.type+', '+r.dist.toFixed(1)+' km</span>'+(r.h2s&&r.h2s!='Not reached'?'<span class="hz" title="H2S exposure in the leaking gas run" style="background:'+(HZ[r.h2s]||'#888')+'">H2S '+r.h2s+'</span>':'')+'</span><span class="t" id="rt'+i+'">'+(r.arrival!=null?('T+'+fmtH(r.arrival)):(r.fall_mm>0?'spray':'-'))+'</span></div>';}).join('');
rl.querySelectorAll('.ritem').forEach(el=>el.addEventListener('click',()=>{const i=+el.dataset.i; const r=DATA.receptors[i]; map.setView([r.lat,r.lon], Math.max(map.getZoom(),15)); recMarkers[i].openPopup();}));
// compass
function compass(toDeg, spd){
  const fromDeg = (toDeg+180)%360;
  const c = document.getElementById('compass');
  let s='<circle r="46" fill="#f3f6f4" stroke="#cfd8d2"/>';
  ['N','E','S','W'].forEach((d,i)=>{ const a=i*Math.PI/2; s+='<text x="'+(38*Math.sin(a))+'" y="'+(-38*Math.cos(a)+4)+'" font-size="10" text-anchor="middle" fill="#5b6b63">'+d+'</text>'; });
  s+='<g transform="rotate('+toDeg+')"><line x1="0" y1="26" x2="0" y2="-22" stroke="#0b5d3b" stroke-width="5" stroke-linecap="round"/><path d="M0,-32 L-10,-16 L10,-16 Z" fill="#0b5d3b"/></g>';
  c.innerHTML = s;
  const dirs=['N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW'];
  document.getElementById('wind_txt').innerHTML = '<b>Wind</b> from '+dirs[Math.round(fromDeg/22.5)%16]+' ('+Math.round(fromDeg)+'&deg;), '+spd.toFixed(1)+' m/s';
}
function fmtH(h){ const hh=Math.floor(h+1e-6), mm=Math.round((h-hh)*60); return String(hh).padStart(2,'0')+':'+String(mm%60).padStart(2,'0'); }
const start = M.start_iso ? new Date(M.start_iso) : null;
const slider = document.getElementById('slider'); slider.max = F.length-1;
let idx=0, timer=null;
function show(i){
  idx = Math.max(0, Math.min(F.length-1, i)); slider.value = idx;
  const f = F[idx], r = MB[Math.min(idx, MB.length-1)];
  oilLayer.setUrl(oilImgs[idx].src); fallLayer.setUrl(fallImgs[idx].src); setCaption(f.cap);
  if (fieldLayer && fieldIdx>=0){ const fd=FL[fieldIdx]; const k=fieldFrame(fd, f.t); fieldLayer.setUrl(fieldImgs[fieldIdx][k].src); }
  document.getElementById('clock').firstChild.nodeValue = 'T+'+fmtH(f.t);
  if (start){ const d=new Date(start.getTime()+f.t*3600e3); document.getElementById('date').textContent = d.toLocaleString('en-GB',{weekday:'short',day:'2-digit',month:'short',year:'numeric',hour:'2-digit',minute:'2-digit',timeZone:'UTC'})+' (local time)'; }
  document.getElementById('k_rel').textContent = r.released.toFixed(r.released<100?1:0);
  document.getElementById('k_area').textContent = r.area_ha.toFixed(2);
  document.getElementById('k_evap').textContent = r.evaporated.toFixed(r.evaporated<100?1:0);
  document.getElementById('k_inf').textContent = r.infiltrated.toFixed(r.infiltrated<100?1:0);
  if (f.wind) compass(f.wind[0], f.wind[1]);
  const rain = f.rain||0;
  document.getElementById('rain_txt').innerHTML = '<b>Rain</b> '+(rain>0?(rain.toFixed(1)+' mm/h (storm total so far '+r.rain_cum_mm.toFixed(0)+' mm)'):'none');
  document.getElementById('stream_txt').innerHTML = '<b>Streams</b> '+M.stream_state;
  document.getElementById('curtain').style.opacity = rain>0 ? Math.min(1, 0.35+rain/40) : 0;
  const badge = document.getElementById('badge');
  const ev = (M.events||[]).filter(e=>Math.abs(e.t-f.t)<=M.frame_h*1.01 && e.t<=f.t+1e-6).pop();
  badge.style.display = ev ? 'block':'none'; if(ev) badge.textContent = ev.text;
  const cx = X(f.t); const cur=document.getElementById('cursor'); cur.setAttribute('x1',cx); cur.setAttribute('x2',cx);
  DATA.receptors.forEach((rc,i)=>{ const hit = (rc.arrival!=null && rc.arrival<=f.t+1e-6);
    recMarkers[i].setStyle({fillColor: hit?'#d73a49':(rc.fall_mm>0 && f.t>0?'#c9a227':'#ffffff'), radius: hit?7:5});
    const el=document.querySelector('.ritem[data-i="'+i+'"]'); if(el) el.classList.toggle('hit',hit); });
}
function play(){ if(timer){ clearInterval(timer); timer=null; document.getElementById('play').textContent='Play'; return; }
  document.getElementById('play').textContent='Pause';
  timer = setInterval(()=>{ if(idx>=F.length-1){ if(document.getElementById('loop').checked) show(0); else { play(); return; } } else show(idx+1); }, +document.getElementById('speed').value); }
document.getElementById('play').onclick = play;
document.getElementById('back').onclick = ()=>show(idx-1);
document.getElementById('fwd').onclick = ()=>show(idx+1);
document.getElementById('speed').onchange = ()=>{ if(timer){ play(); play(); } };
slider.oninput = ()=>show(+slider.value);
document.addEventListener('keydown',e=>{ if(e.key===' '){e.preventDefault();play();} if(e.key==='ArrowRight')show(idx+1); if(e.key==='ArrowLeft')show(idx-1); });
// ---------------- sandbox tanker: D8 path tracer ----------------
const SB = DATA.sandbox;
if (SB){
  const dec=b=>{const s=atob(b), a=new Uint8Array(s.length); for(let i=0;i<s.length;i++) a[i]=s.charCodeAt(i); return a;};
  const code=dec(SB.code), zb=dec(SB.z), ab=dec(SB.a);
  const zv=new Uint16Array(zb.buffer), D8R=[-1,-1,0,1,1,1,0,-1], D8C=[0,1,1,1,0,-1,-1,-1];
  const toXY=ll=>[SB.inv[0][0]+SB.inv[1][0]*ll.lat+SB.inv[2][0]*ll.lng, SB.inv[0][1]+SB.inv[1][1]*ll.lat+SB.inv[2][1]*ll.lng];
  const toLL=(x,y)=>[SB.fwd[0][0]+SB.fwd[1][0]*x+SB.fwd[2][0]*y, SB.fwd[0][1]+SB.fwd[1][1]*x+SB.fwd[2][1]*y];
  const cellXY=(r,c)=>[SB.x0+(c+0.5)*SB.dx, SB.y1-(r+0.5)*SB.dx];
  let sbOn=false, sbMarker=null, sbLine=null, sbChan=null;
  const tankIcon=L.divIcon({className:'tankicon'});
  function trace(ll){
    const xy=toXY(ll); let c=Math.floor((xy[0]-SB.x0)/SB.dx), r=Math.floor((SB.y1-xy[1])/SB.dx);
    const out=document.getElementById('sbx_out');
    if(r<0||c<0||r>=SB.nr||c>=SB.nc){ out.innerHTML='Outside the traced area (within '+Math.round(SB.nr*SB.dx/2000)+' km of the source).'; return; }
    const pts=[], seen=new Set(); let land=0, strm=0, tmin=0, first=null, n=0, lastz=null, edge=false;
    while(n<6000){ const i=r*SB.nc+c; if(seen.has(i)) break; seen.add(i);
      const xyc=cellXY(r,c); pts.push(toLL(xyc[0],xyc[1]));
      const k=code[i]; if(k===255){ break; }
      const isCh=k>=8, d=k%8, rn=r+D8R[d], cn=c+D8C[d], step=SB.dx*((d%2)?1.41421:1);
      const z0=SB.zmin+zv[i]/10; let S=0.002;
      if(rn>=0&&cn>=0&&rn<SB.nr&&cn<SB.nc){ S=Math.max((z0-(SB.zmin+zv[rn*SB.nc+cn]/10))/step,0.002); }
      if(isCh){ if(first===null) first={i:pts.length-1, d:land, t:tmin}; strm+=step;
        const A=Math.pow(10, ab[i]/30-4), Q=SB.q_spec*A/1000; const u=(Q>0.003)?0.77*Math.pow(Q,0.1):Math.pow(0.005,2/3)*Math.sqrt(S)/0.05; tmin+=step/u/60; }
      else { land+=step; tmin+=step/(Math.pow(0.005,2/3)*Math.sqrt(S)/0.25)/60; }
      if(rn<0||cn<0||rn>=SB.nr||cn>=SB.nc){ edge=true; break; }
      r=rn; c=cn; n++; }
    if(sbLine) map.removeLayer(sbLine); if(sbChan) map.removeLayer(sbChan);
    sbLine=L.polyline(first?pts.slice(0,first.i+1):pts,{pane:'pts',color:'#ff8c00',weight:4,opacity:0.95,dashArray:'6 4'}).addTo(map);
    if(first) sbChan=L.polyline(pts.slice(first.i),{pane:'pts',color:'#e8590c',weight:5,opacity:0.95}).addTo(map);
    // receptors near the path
    const near=[];
    DATA.receptors.forEach(rc=>{ let best=1e9, at=0; pts.forEach((p,j)=>{ const dy=(p[0]-rc.lat)*111320, dx_=(p[1]-rc.lon)*111320*Math.cos(rc.lat*Math.PI/180); const dd=Math.hypot(dx_,dy); if(dd<best){best=dd;at=j;} });
      if(best<=120) near.push([rc.name, rc.type, at*SB.dx/1000]); });
    near.sort((a,b)=>a[2]-b[2]);
    const fmt=m=>m<60?Math.round(m)+' min':(m/60).toFixed(1)+' h';
    out.innerHTML='<div class="row"><span>Path traced</span><span>'+((land+strm)/1000).toFixed(2)+' km'+(edge?' (leaves the traced area)':'')+'</span></div>'+
      '<div class="row"><span>Overland before a watercourse</span><span>'+(first?Math.round(first.d)+' m':'no watercourse reached')+'</span></div>'+
      (first?'<div class="row"><span>Indicative time to the watercourse</span><span>'+fmt(first.t)+'</span></div>':'')+
      '<div class="row"><span>Along watercourses</span><span>'+(strm/1000).toFixed(2)+' km ('+(SB.q_spec>0?'flowing':'dry beds')+')</span></div>'+
      '<div class="row"><span>Indicative travel time, whole path</span><span>'+fmt(tmin)+'</span></div>'+
      '<div style="margin-top:6px"><b>Receptors within 120 m of the path</b><br>'+(near.length?near.slice(0,12).map(x=>x[0]+' <span class="tag">'+x[1]+', '+x[2].toFixed(2)+' km along</span>').join('<br>'):'none')+'</div>';
  }
  const SbCtl=L.Control.extend({options:{position:'topleft'}, onAdd:function(){ const d=L.DomUtil.create('div','leaflet-bar');
    d.innerHTML='<a href="#" title="Sandbox tanker: drag a tanker and trace the spill path" style="width:auto;padding:0 6px;font-size:12px;background:#fff3e0">Sandbox tanker</a>';
    d.querySelector('a').onclick=e=>{ e.preventDefault(); sbOn=!sbOn; document.getElementById('sbx').classList.toggle('on',sbOn);
      e.target.style.background=sbOn?'#ff8c00':'#fff3e0';
      if(sbOn){ const c0=map.getCenter(); sbMarker=L.marker(c0,{draggable:true,icon:tankIcon,pane:'pts'}).addTo(map); sbMarker.on('dragend',()=>trace(sbMarker.getLatLng())); trace(c0); }
      else { [sbMarker,sbLine,sbChan].forEach(l=>{ if(l) map.removeLayer(l); }); sbMarker=sbLine=sbChan=null; } };
    L.DomEvent.disableClickPropagation(d); return d; }});
  new SbCtl().addTo(map);
  map.on('click',e=>{ if(sbOn && sbMarker){ sbMarker.setLatLng(e.latlng); trace(e.latlng); } });
}
show(0);
})();
</script>
</body>
</html>
"""


def build_h2s_html(path, payload):
    with open(os.path.join(HERE, "web", "leaflet.js"), encoding="utf-8") as f:
        ljs = f.read()
    with open(os.path.join(HERE, "web", "leaflet.css"), encoding="utf-8") as f:
        lcss = f.read()
    css = TEMPLATE.split("<style>", 1)[1].split("</style>", 1)[0]
    html = H2S_TEMPLATE.replace("/*BASE_CSS*/", css)
    html = html.replace("/*LEAFLET_CSS*/", lcss)
    html = html.replace("/*LEAFLET_JS*/", ljs)
    html = html.replace("__TITLE__", payload["meta"]["title"])
    html = html.replace("/*PAYLOAD*/", "const DATA = " + json.dumps(payload, separators=(",", ":")) + ";\n" + BASE_JS)
    html = html.replace("</style>", CAPTION_CSS + "</style>", 1)
    with open(path, "w", encoding="utf-8") as f:
        f.write(html)
    return path


H2S_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
/*BASE_CSS*/
.hz{display:inline-block;font-size:10.5px;font-weight:700;padding:1px 6px;border-radius:9px;color:#fff;white-space:nowrap}
.thr{display:flex;align-items:center;gap:6px;font-size:12px;padding:2px 0}
.thr i{display:inline-block;width:22px;height:0;border-top:3px solid}
#links{width:100%;margin:4px 0 8px}
</style>
</head>
<body>
<div id="app">
 <div id="map"><div class="badge" id="badge"></div></div>
 <aside id="side">
  <select id="links" style="display:none"></select>
  <h1 id="title"></h1>
  <div class="sub" id="subtitle"></div>
  <div class="card">
   <div class="clock" id="clock">T+00:00<small id="date"></small></div>
   <input type="range" id="slider" min="0" value="0" step="1">
   <div class="controls">
    <button id="back">&#9664;&#9664;</button><button id="play" class="primary">Play</button><button id="fwd">&#9654;&#9654;</button>
    <select id="speed"><option value="700">Slow</option><option value="300" selected>Normal</option><option value="120">Fast</option></select>
    <label style="font-size:12px"><input type="checkbox" id="loop" checked> Loop</label>
   </div>
  </div>
  <div class="card">
   <h2>Now</h2>
   <div class="kpis">
    <div class="kpi"><b id="k_q">0</b><span>H2S leak rate (g/s)</span></div>
    <div class="kpi"><b id="k_max">0</b><span>highest ground level now (ppm)</span></div>
    <div class="kpi"><b id="k_a2">0</b><span>AEGL-2 reach now (m)</span></div>
    <div class="kpi"><b id="k_idlh">0</b><span>IDLH reach now (m)</span></div>
   </div>
   <div class="wx" style="margin-top:10px">
    <svg id="compass" viewBox="-50 -50 100 100"></svg>
    <div class="rain"><div id="wind_txt"></div><div id="stab_txt"></div><div id="mix_txt"></div></div>
   </div>
  </div>
  <div class="card">
   <h2>Downwind reach of each threshold (m)</h2>
   <svg id="chart" viewBox="0 0 360 150" preserveAspectRatio="none"></svg>
   <div id="chart_legend" class="note"></div>
  </div>
  <div class="card">
   <h2>Receptors</h2>
   <div class="note" style="margin-bottom:6px">Badge = worst exposure over the run. Marker colour = concentration at the displayed time. Click a name to zoom.</div>
   <div class="rlist" id="rlist"></div>
  </div>
  <div class="card">
   <h2>Legend</h2>
   <div style="font-size:12px">Ground-level H2S, 10-minute mean (ppm)</div>
   <div class="legend-bar" id="lg"></div><div class="legend-scale" id="lg_s"></div>
   <div id="thr" style="margin-top:8px"></div>
  </div>
  <div class="card"><h2>Scenario</h2><div id="facts"></div></div>
  <div class="card"><h2>Basis and limitations</h2><div class="note" id="notes"></div></div>
 </aside>
</div>
<script>
/*LEAFLET_JS*/
</script>
<script>
/*PAYLOAD*/
(function(){
const M=DATA.meta, F=DATA.frames, S=DATA.stats, T=DATA.thresholds;
document.getElementById('title').textContent=M.title;
document.getElementById('subtitle').textContent=M.subtitle;
if((M.links||[]).length){ const sel=document.getElementById('links'); sel.style.display='block';
  sel.innerHTML='<option value="">Other scenarios...</option>'+M.links.map(l=>'<option value="'+l.href+'">'+l.name+'</option>').join('');
  sel.onchange=()=>{ if(sel.value) window.location.href=sel.value; }; }
const map=L.map('map',{preferCanvas:true});
const base=makeBase(DATA.hill.uri, DATA.hill.bounds);
baseFallback(map, base, 'Satellite imagery (Esri)', 'Terrain - relief from the model DEM (offline)');
subtitleBox(map);
const B=[[DATA.bounds[0],DATA.bounds[1]],[DATA.bounds[2],DATA.bounds[3]]];
map.createPane('vec').style.zIndex=410; map.createPane('gas').style.zIndex=420; map.createPane('pts').style.zIndex=440;
const imgs=F.map(f=>{const i=new Image(); i.src=f.img; return i;});
const plume=L.imageOverlay(F[0].img,B,{pane:'gas',opacity:0.9,interactive:false}).addTo(map);
const overlays={'H2S at ground level (animated)':plume,'Terrain relief shading (DEM, over any basemap)':reliefOverlay(map, DATA.hill.uri, DATA.hill.bounds)};
let zoneBounds=null;
(DATA.layers||[]).forEach(Ly=>{ const s=Ly.style||{}; let lay;
 if(Ly.kind==='zones'){
  lay=L.geoJSON(Ly.geojson,{pane:'vec',style:f=>({color:f.properties.colour,weight:2,fillColor:f.properties.colour,fillOpacity:0.08,dashArray:'5 3'}),
   onEachFeature:(f,l)=>l.bindTooltip(f.properties.threshold+' - reaches '+Math.round(f.properties.max_dist_m)+' m ('+f.properties.basis+')',{sticky:true})});
  try{ lay.eachLayer(l=>{ const t=l.feature.properties.threshold||''; if(t.indexOf('AEGL-2')===0) zoneBounds=l.getBounds(); }); }catch(e){}
 } else if(Ly.kind==='point'){
  lay=L.geoJSON(Ly.geojson,{pane:'pts',pointToLayer:(f,ll)=>L.circleMarker(ll,{radius:s.radius||4,color:s.color||'#222',weight:1,fillColor:s.fill||'#fff',fillOpacity:0.9}),
   onEachFeature:(f,l)=>{const nm=(f.properties||{})[s.label]||''; if(nm) l.bindTooltip(String(nm),{permanent:!!s.permanent,direction:'right',className:'maplabel',offset:[4,0]});}});
 } else {
  lay=L.geoJSON(Ly.geojson,{pane:'vec',style:()=>({color:s.color||'#333',weight:s.weight||1.5,fillColor:s.fill||s.color||'#333',fillOpacity:s.fillOpacity==null?0.1:s.fillOpacity,dashArray:s.dash||null}),
   onEachFeature:(f,l)=>{const nm=(f.properties||{})[s.label]||Ly.name; l.bindTooltip(String(nm),{sticky:true});}});
 }
 overlays[Ly.name]=lay; if(Ly.visible) lay.addTo(map); });
const rings=L.layerGroup(); [500,1000,2000,5000].forEach(r=>L.circle([M.source[0],M.source[1]],{radius:r,pane:'vec',color:'#555',weight:0.8,fill:false,dashArray:'2 4'}).bindTooltip((r/1000)+' km',{sticky:true}).addTo(rings));
overlays['Distance rings']=rings.addTo(map);
const src=L.circleMarker(M.source,{pane:'pts',radius:7,color:'#000',weight:2,fillColor:'#ffd400',fillOpacity:1}).bindTooltip(M.source_label,{permanent:true,direction:'left',className:'maplabel',offset:[-6,0]});
overlays["Leak point"]=src.addTo(map);
const HZ={'500 ppm+':'#3b0a45','IDLH':'#9c36b5','AEGL-3':'#e03131','AEGL-2':'#f08c00','AEGL-1':'#2f9e44','Odour only':'#74c0fc','Not reached':'#adb5bd'};
function ppmCol(v){ if(v>=500)return HZ['500 ppm+']; if(v>=100)return HZ['IDLH']; if(v>=50)return HZ['AEGL-3']; if(v>=27)return HZ['AEGL-2']; if(v>=0.51)return HZ['AEGL-1']; if(v>=0.05)return HZ['Odour only']; return '#ffffff'; }
const recL=L.layerGroup(), mk=[];
DATA.receptors.forEach((r,i)=>{ const m=L.circleMarker([r.lat,r.lon],{pane:'pts',radius:5,color:'#333',weight:1,fillColor:'#fff',fillOpacity:1});
 m.bindPopup('<b>'+r.name+'</b><br>'+r.type+'<br>'+r.dist.toFixed(2)+' km from release<br>Worst exposure: <span class="hz" style="background:'+HZ[r.hazard]+'">'+r.hazard+'</span><br>Peak 10-min: '+r.peak.toFixed(2)+' ppm<br>Peak 60-min: '+r.m1.toFixed(2)+' ppm'+(r.a2!=null?'<br>AEGL-2 first exceeded: T+'+fmtH(r.a2):'')+(r.idlh!=null?'<br>IDLH first exceeded: T+'+fmtH(r.idlh):''));
 m.addTo(recL); mk.push(m); });
overlays['Receptors']=recL.addTo(map);
L.control.layers(base,overlays,{collapsed:window.innerWidth<1600}).addTo(map);
L.control.scale({imperial:false}).addTo(map);
const near=L.latLng(M.source).toBounds(4000);
if(zoneBounds && zoneBounds.isValid()){ zoneBounds.extend(near); map.fitBounds(zoneBounds.pad(0.1)); } else map.fitBounds(near);
const grad='linear-gradient(90deg,'+DATA.ramp.map(s=>'rgba('+s[1]+','+s[2]+','+s[3]+','+(Math.max(s[4],140)/255).toFixed(2)+')').join(',')+')';
document.getElementById('lg').style.background=grad;
document.getElementById('lg_s').innerHTML=DATA.ramp.map(s=>'<span>'+s[0]+'</span>').join('');
document.getElementById('thr').innerHTML=T.map(t=>'<div class="thr"><i style="border-color:'+t.color+'"></i>'+t.label+' <span class="tag">('+(t.basis=='1h'?'60-min mean':'10-min mean')+')</span></div>').join('')+(M.zone_m?'<div class="thr"><i style="border-color:#000;border-top-style:dashed"></i>Planned H2S zone, '+M.zone_m+' m</div>':'');
document.getElementById('facts').innerHTML=M.facts.map(f=>'<div class="row"><span>'+f[0]+'</span><span>'+f[1]+'</span></div>').join('');
document.getElementById('notes').innerHTML=M.notes.map(n=>'<p style="margin:0 0 6px">'+n+'</p>').join('');
// chart
const svg=document.getElementById('chart'); const W=360,H=150,PL=40,PB=16,PT=6,PR=6;
const tmax=S[S.length-1].t_h; const keys=['AEGL1','AEGL2','AEGL3','IDLH'];
const vmax=Math.max(100,...S.map(s=>Math.max(...keys.map(k=>s[k]||0))))*1.08;
const X=t=>PL+(W-PL-PR)*t/Math.max(tmax,1e-6), Y=v=>H-PB-(H-PB-PT)*v/vmax;
let g='<line x1="'+PL+'" y1="'+(H-PB)+'" x2="'+(W-PR)+'" y2="'+(H-PB)+'" stroke="#999" stroke-width="0.6"/>';
[0,0.25,0.5,0.75,1].forEach(f=>{const t=f*tmax; g+='<text x="'+X(t)+'" y="'+(H-3)+'" font-size="9" text-anchor="middle" fill="#666">'+Math.round(t)+' h</text>';});
[0,0.5,1].forEach(f=>{const v=f*vmax/1.08; g+='<text x="'+(PL-3)+'" y="'+(Y(v)+3)+'" font-size="9" text-anchor="end" fill="#666">'+Math.round(v)+'</text><line x1="'+PL+'" x2="'+(W-PR)+'" y1="'+Y(v)+'" y2="'+Y(v)+'" stroke="#e5e9e6" stroke-width="0.6"/>';});
if(M.zone_m) g+='<line x1="'+PL+'" x2="'+(W-PR)+'" y1="'+Y(M.zone_m)+'" y2="'+Y(M.zone_m)+'" stroke="#000" stroke-dasharray="4 3" stroke-width="0.8"/>';
keys.forEach(k=>{const t=T.find(x=>x.key==k); let d=''; S.forEach((s,i)=>{d+=(i?' L':'M')+X(s.t_h).toFixed(1)+','+Y(s[k]||0).toFixed(1);}); g+='<path d="'+d+'" fill="none" stroke="'+t.color+'" stroke-width="1.5"/>';});
svg.innerHTML=g+'<line id="cursor" x1="0" x2="0" y1="'+PT+'" y2="'+(H-PB)+'" stroke="#111" stroke-width="1"/>';
document.getElementById('chart_legend').innerHTML=keys.map(k=>{const t=T.find(x=>x.key==k); return '<span style="white-space:nowrap;margin-right:8px"><span style="display:inline-block;width:10px;height:3px;background:'+t.color+';margin-right:3px;vertical-align:middle"></span>'+t.label.split(' (')[0].split(' ')[0]+'</span>';}).join('')+(M.zone_m?'<span>- - planned zone</span>':'');
// receptors
const ord=DATA.receptors.map((r,i)=>i).sort((a,b)=>DATA.receptors[b].peak-DATA.receptors[a].peak||DATA.receptors[a].dist-DATA.receptors[b].dist);
const rl=document.getElementById('rlist');
rl.innerHTML=ord.slice(0,60).map(i=>{const r=DATA.receptors[i]; return '<div class="ritem" data-i="'+i+'"><span>'+r.name+' <span class="tag">'+r.type+', '+r.dist.toFixed(1)+' km</span></span><span class="t"><span class="hz" style="background:'+HZ[r.hazard]+'">'+r.hazard+'</span> <span id="rv'+i+'"></span></span></div>';}).join('');
rl.querySelectorAll('.ritem').forEach(el=>el.addEventListener('click',()=>{const i=+el.dataset.i,r=DATA.receptors[i]; map.setView([r.lat,r.lon],Math.max(map.getZoom(),14)); mk[i].openPopup();}));
function compass(toDeg,spd){ const fromDeg=(toDeg+180)%360; let s='<circle r="46" fill="#f3f6f4" stroke="#cfd8d2"/>';
 ['N','E','S','W'].forEach((d,i)=>{const a=i*Math.PI/2; s+='<text x="'+(38*Math.sin(a))+'" y="'+(-38*Math.cos(a)+4)+'" font-size="10" text-anchor="middle" fill="#5b6b63">'+d+'</text>';});
 s+='<g transform="rotate('+toDeg+')"><line x1="0" y1="26" x2="0" y2="-22" stroke="#0b5d3b" stroke-width="5" stroke-linecap="round"/><path d="M0,-32 L-10,-16 L10,-16 Z" fill="#0b5d3b"/></g>';
 document.getElementById('compass').innerHTML=s; const dirs=['N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW'];
 document.getElementById('wind_txt').innerHTML='<b>Wind</b> from '+dirs[Math.round(fromDeg/22.5)%16]+' ('+Math.round(fromDeg)+'&deg;), '+spd.toFixed(1)+' m/s'; }
const SD={A:'very unstable',B:'unstable',C:'slightly unstable',D:'neutral',E:'slightly stable',F:'stable (night, light wind)'};
function fmtH(h){const hh=Math.floor(h+1e-6),mm=Math.round((h-hh)*60); return String(hh).padStart(2,'0')+':'+String(mm%60).padStart(2,'0');}
const start=M.start_iso?new Date(M.start_iso):null; const sl=document.getElementById('slider'); sl.max=F.length-1;
let idx=0,timer=null;
function show(i){ idx=Math.max(0,Math.min(F.length-1,i)); sl.value=idx; const f=F[idx], s=S[idx];
 plume.setUrl(imgs[idx].src); setCaption(f.cap);
 document.getElementById('clock').firstChild.nodeValue='T+'+fmtH(f.t);
 if(start){const d=new Date(start.getTime()+f.t*3600e3); document.getElementById('date').textContent=d.toLocaleString('en-GB',{weekday:'short',day:'2-digit',month:'short',year:'numeric',hour:'2-digit',minute:'2-digit',timeZone:'UTC'})+' (local time)';}
 document.getElementById('k_q').textContent=((s.releasing||0)*M.q_gs).toFixed(2);
 document.getElementById('k_max').textContent=s.max_ppm>=10?Math.round(s.max_ppm):s.max_ppm.toFixed(2);
 document.getElementById('k_a2').textContent=Math.round(s.AEGL2||0);
 document.getElementById('k_idlh').textContent=Math.round(s.IDLH||0);
 compass(f.wind[0],Math.max(f.wind[1],1));
 document.getElementById('stab_txt').innerHTML='<b>Stability</b> class '+s.stability+' - '+SD[s.stability];
 document.getElementById('mix_txt').innerHTML='<b>Mixing height</b> '+Math.round(s.mixing_m)+' m; plume height '+Math.round(s.plume_h_m)+' m';
 const cx=X(f.t); const cu=document.getElementById('cursor'); cu.setAttribute('x1',cx); cu.setAttribute('x2',cx);
 const badge=document.getElementById('badge'); const stop=M.duration_h; const near=Math.abs(f.t-stop)<=M.frame_h*1.01&&f.t>=stop;
 badge.style.display=near?'block':'none'; if(near) badge.textContent='Release stops at T+'+fmtH(stop);
 DATA.receptors.forEach((r,i)=>{ const v=r.series?r.series[idx]||0:0; mk[i].setStyle({fillColor:ppmCol(v),radius:v>=27?8:(v>=0.51?6:5)});
   const el=document.getElementById('rv'+i); if(el) el.textContent=v>=0.05?(v>=10?Math.round(v):v.toFixed(1))+' ppm':''; });
}
function play(){ if(timer){clearInterval(timer);timer=null;document.getElementById('play').textContent='Play';return;}
 document.getElementById('play').textContent='Pause';
 timer=setInterval(()=>{ if(idx>=F.length-1){ if(document.getElementById('loop').checked) show(0); else {play();return;} } else show(idx+1); }, +document.getElementById('speed').value); }
document.getElementById('play').onclick=play; document.getElementById('back').onclick=()=>show(idx-1); document.getElementById('fwd').onclick=()=>show(idx+1);
document.getElementById('speed').onchange=()=>{if(timer){play();play();}}; sl.oninput=()=>show(+sl.value);
document.addEventListener('keydown',e=>{ if(e.key===' '){e.preventDefault();play();} if(e.key==='ArrowRight')show(idx+1); if(e.key==='ArrowLeft')show(idx-1); });
show(0);
})();
</script>
</body>
</html>
"""
