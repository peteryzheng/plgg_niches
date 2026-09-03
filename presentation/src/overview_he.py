import json, os
from PIL import Image, ImageDraw, ImageFont
Image.MAX_IMAGE_PIXELS = None
BASE = "/Users/youyun/Documents/HMS/PhD/beroukhimlab/dfci_mount/youyun/plgg/data/Xenium_annotations"
OUT = "/Users/youyun/plgg_deck/assets"; os.makedirs(OUT, exist_ok=True)
S = 1/7
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

img = Image.open(f"{BASE}/images/230918_Xenium_CytAssist_LGG1.jpg").convert("RGB")
W,H = img.size
OW = 2600; sd = OW/W; OH = int(H*sd)
ov = img.resize((OW,OH)); d = ImageDraw.Draw(ov)

gj = json.load(open(f"{BASE}/geojsons/230918_Xenium_CytAssist_LGG1.geojson"))
feats = gj["features"] if isinstance(gj,dict) else gj
for ft in feats:
    n=name_of(ft); col=GREEN if n=="Compact. fibrillary component" else BLUE if n=="Loose, myxoid component" else None
    if not col: continue
    for r in rings(ft.get("geometry",{}) or {}):
        pts=[(px*S*sd, py*S*sd) for px,py in r]
        d.line(pts+[pts[0]], fill=col, width=2)

try: font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 20)
except Exception:
    try: font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 20)
    except Exception: font = ImageFont.load_default()

# grid every 1000 JPG px, labeled in JPG-pixel coords (the crop coordinate system)
for gx in range(0, W, 1000):
    x=gx*sd; d.line([(x,0),(x,OH)], fill=(90,90,90), width=1)
    d.text((x+3,3), str(gx), fill=(0,0,0), font=font)
for gy in range(0, H, 1000):
    y=gy*sd; d.line([(0,y),(OW,y)], fill=(90,90,90), width=1)
    d.text((3,y+2), str(gy), fill=(0,0,0), font=font)

ov.save(f"{OUT}/overview_lgg1_grid.jpg", quality=90)
print("overview:", (OW,OH), "JPG full:", (W,H), "-> assets/overview_lgg1_grid.jpg")
print("green outlines = compact fibrillary ; blue = loose myxoid ; grid labels are JPG px")
