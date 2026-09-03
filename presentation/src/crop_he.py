import json, os
from PIL import Image, ImageDraw
Image.MAX_IMAGE_PIXELS = None
BASE = "/Users/youyun/Documents/HMS/PhD/beroukhimlab/dfci_mount/youyun/plgg/data/Xenium_annotations"
OUT = "/Users/youyun/plgg_deck/assets"; os.makedirs(OUT, exist_ok=True)
S = 1/7
# prettier, still high-contrast on pink/purple H&E: compact = teal-green, loose = blue
RED=(34,160,122); BLUE=(46,107,196)

def name_of(ft):
    p = ft.get("properties",{}) or {}; cls = p.get("classification") or {}
    return p.get("name") or (cls.get("name") if isinstance(cls,dict) else None)
def rings(geom):
    t=geom.get("type"); c=geom.get("coordinates",[]); out=[]
    if t=="Polygon" and c: out.append(c[0])
    elif t=="MultiPolygon":
        for poly in c:
            if poly: out.append(poly[0])
    return out
def ringinfo(ring):  # in JPG px
    xs=[p[0]*S for p in ring]; ys=[p[1]*S for p in ring]
    return dict(pts=list(zip(xs,ys)), cx=sum(xs)/len(xs), cy=sum(ys)/len(ys),
                bb=(min(xs),min(ys),max(xs),max(ys)),
                area=(max(xs)-min(xs))*(max(ys)-min(ys)))

img = Image.open(f"{BASE}/images/230918_Xenium_CytAssist_LGG1.jpg").convert("RGB")
W,H = img.size
gj = json.load(open(f"{BASE}/geojsons/230918_Xenium_CytAssist_LGG1.geojson"))
feats = gj["features"] if isinstance(gj,dict) else gj
comp=[]; loose=[]
for ft in feats:
    n=name_of(ft); g=ft.get("geometry",{}) or {}
    for r in rings(g):
        info=ringinfo(r)
        if info["area"]<2500: continue        # drop tiny slivers
        if n=="Compact. fibrillary component": comp.append(info)
        elif n=="Loose, myxoid component": loose.append(info)
print("compact",len(comp),"loose",len(loose))

# prefer LARGE adjacent compact+loose pairs (fuller, more intact tissue), not tiny slivers
pairs=[]
for a in comp:
    for b in loose:
        d=((a["cx"]-b["cx"])**2+(a["cy"]-b["cy"])**2)**0.5
        touch = d < 0.8*((a["area"]**0.5)+(b["area"]**0.5))   # roughly adjacent
        if not touch: continue
        pairs.append((a["area"]+b["area"], d, a, b))
pairs.sort(key=lambda t:-t[0])                                 # biggest combined area first

def sq_window(a,b,pad=1.18,lo=760,hi=1050):
    x0=min(a["bb"][0],b["bb"][0]); y0=min(a["bb"][1],b["bb"][1])
    x1=max(a["bb"][2],b["bb"][2]); y1=max(a["bb"][3],b["bb"][3])
    cx=(x0+x1)/2; cy=(y0+y1)/2
    side=max(x1-x0,y1-y0)*pad; side=max(lo,min(hi,side))
    X0=int(max(0,cx-side/2)); Y0=int(max(0,cy-side/2))
    X1=int(min(W,X0+side));  Y1=int(min(H,Y0+side))
    return X0,Y0,X1,Y1

allrings=[(r,RED) for r in comp]+[(r,BLUE) for r in loose]
made=0
for k,(sc,d,a,b) in enumerate(pairs[:4]):
    X0,Y0,X1,Y1=sq_window(a,b)
    crop=img.crop((X0,Y0,X1,Y1))
    crop.save(f"{OUT}/pa_sq{k+1}_clean.jpg",quality=93)
    ann=crop.copy(); dr=ImageDraw.Draw(ann)
    for r,color in allrings:
        bb=r["bb"]
        if bb[2]<X0 or bb[0]>X1 or bb[3]<Y0 or bb[1]>Y1: continue
        pts=[(px-X0,py-Y0) for px,py in r["pts"]]
        dr.line(pts+[pts[0]], fill=color, width=3)     # thin
    ann.save(f"{OUT}/pa_sq{k+1}_annot.jpg",quality=93)
    print(f"sq{k+1}: area={sc:.0f} dist={d:.0f}px window={X1-X0}x{Y1-Y0} at ({X0},{Y0})")
    made+=1
print("DONE",made,"square candidates ->",OUT)
