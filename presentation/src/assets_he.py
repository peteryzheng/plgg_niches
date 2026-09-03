import json, os
from PIL import Image, ImageDraw
Image.MAX_IMAGE_PIXELS = None
BASE = "/Users/youyun/Documents/HMS/PhD/beroukhimlab/dfci_mount/youyun/plgg/data/Xenium_annotations"
OUT = "/Users/youyun/plgg_deck/assets"; os.makedirs(OUT, exist_ok=True); S=1/7
GREEN=(34,160,122); BLUE=(46,107,196)

def name_of(ft):
    p=ft.get("properties",{}) or {}; cls=p.get("classification") or {}
    return p.get("name") or (cls.get("name") if isinstance(cls,dict) else None)
def rings(g):
    t=g.get("type"); c=g.get("coordinates",[]); out=[]
    if t=="Polygon" and c: out.append(c[0])
    elif t=="MultiPolygon":
        for poly in c:
            if poly: out.append(poly[0])
    return out

img1 = Image.open(f"{BASE}/images/230918_Xenium_CytAssist_LGG1.jpg").convert("RGB")
img2 = Image.open(f"{BASE}/images/230918_Xenium_CytAssist_LGG2.jpg").convert("RGB")

# 1) motivation biphasic crop, THICK outlines (width 4), box on LGG1
X0,Y0,X1,Y1 = 10650,3250,11150,3750
crop = img1.crop((X0,Y0,X1,Y1)); crop.save(f"{OUT}/pa_biphasic_slide_clean.jpg", quality=95)
gj1 = json.load(open(f"{BASE}/geojsons/230918_Xenium_CytAssist_LGG1.geojson"))
feats1 = gj1["features"] if isinstance(gj1,dict) else gj1
ann = crop.copy(); dr = ImageDraw.Draw(ann)
for ft in feats1:
    n=name_of(ft); col=GREEN if n=="Compact. fibrillary component" else BLUE if n=="Loose, myxoid component" else None
    if not col: continue
    for r in rings(ft.get("geometry",{}) or {}):
        pts=[(px*S-X0,py*S-Y0) for px,py in r]
        xs=[p[0] for p in pts]; ys=[p[1] for p in pts]
        if max(xs)<0 or min(xs)>500 or max(ys)<0 or min(ys)>500: continue
        dr.line(pts+[pts[0]], fill=col, width=4)   # 2x thicker
ann.save(f"{OUT}/pa_biphasic_slide_annot.jpg", quality=95)

# 2) four entity thumbnails (clean squares), by diagnostic region location
def sq(img, cx, cy, side, path):
    x0=int(cx-side/2); y0=int(cy-side/2)
    img.crop((x0,y0,x0+side,y0+side)).save(path, quality=92)
sq(img1, (X0+X1)/2, (Y0+Y1)/2, 600, f"{OUT}/thumb_PA.jpg")     # PA: biphasic (LGG1)
sq(img1, 1513, 1029, 620, f"{OUT}/thumb_GG.jpg")               # GG: dysplastic neurons/ganglion cells (LGG1)
sq(img2, 5712, 4807, 520, f"{OUT}/thumb_CN.jpg")               # CN: neuropil zone (LGG2) [verify]
sq(img2, 2740, 3690, 640, f"{OUT}/thumb_CLN.jpg")              # CLN: lipidized cells (LGG2), shifted off artifact
print("DONE: biphasic slide crop (thick) + 4 thumbnails ->", OUT)
