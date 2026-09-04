#!/usr/bin/env python
"""
Prepare square headshots for the Acknowledgements slide (slide 18) bubble grid.

Reads source portraits from `presentation/src/ack_src/<slug>.<ext>`, face-biased
center-crops each to a square, resizes to 512 px, and writes
`presentation/assets/ack_<slug>.jpg`. build_deck.js then renders each as a circle
via addImage({rounding:true}).

Source provenance (downloaded by curl into ack_src/, or user-pasted):
  beroukhim      user-supplied direct URL (googleusercontent)
  bandopadhayay  ispno2026.org  -> squarespace-cdn "Mimi pic"
  sievers        ghst.de/interview-philipp-sievers -> Philipp_Sievers_Hertie.jpg
  sahm           user-supplied portrait (~/Downloads/felix.jpg)
  savran         user-pasted (chat Image #1)
  collins        user-pasted (chat Image #2, MIT)

Run: /Users/youyun/miniforge3/envs/spatial/bin/python presentation/src/prep_ack_headshots.py
"""
import glob
from pathlib import Path
from PIL import Image, ImageDraw, ImageOps, ImageFont

HERE = Path(__file__).resolve().parent            # presentation/src
SRC = HERE / "ack_src"
ASSETS = HERE.parent / "assets"                    # presentation/assets
OUT = 512
# order matches the slide grid: row1 (DKFZ/German), row2
SLUGS = ["savran", "sievers", "sahm", "collins", "bandopadhayay", "beroukhim"]

# faces sit slightly above the vertical middle -> bias the square crop upward
FACE_CENTER_FRAC = 0.42

# per-slug pre-crop boxes (left,top,right,bottom) for sources that aren't clean
# portraits. Currently none — all sources are proper headshots.
PRECROP = {}


def square_crop(im, slug=None):
    im = ImageOps.exif_transpose(im).convert("RGB")   # honor EXIF rotation
    if slug in PRECROP:
        im = im.crop(PRECROP[slug])
    w, h = im.size
    side = min(w, h)
    if h >= w:                                         # portrait: bias toward the face
        top = int(h * FACE_CENTER_FRAC - side / 2)
        top = max(0, min(top, h - side))
        box = (0, top, side, top + side)
    else:                                              # landscape: center horizontally
        left = (w - side) // 2
        box = (left, 0, left + side, side)
    return im.crop(box).resize((OUT, OUT), Image.LANCZOS)


crops = {}
for slug in SLUGS:
    hits = [p for p in glob.glob(str(SRC / f"{slug}.*")) if not p.endswith(".jpg_out")]
    if not hits:
        print(f"MISSING source for {slug}")
        continue
    im = square_crop(Image.open(hits[0]), slug)
    outp = ASSETS / f"ack_{slug}.jpg"
    im.save(outp, quality=90)
    crops[slug] = im
    print(f"wrote {outp.name}  (from {Path(hits[0]).name})")

# verification montage: 3 cols x 2 rows with slug labels
cell, pad = OUT // 2, 16
cols, rows = 3, 2
W = cols * cell + (cols + 1) * pad
H = rows * (cell + 28) + pad
sheet = Image.new("RGB", (W, H), "white")
draw = ImageDraw.Draw(sheet)
for i, slug in enumerate(SLUGS):
    if slug not in crops:
        continue
    r, c = divmod(i, cols)
    x = pad + c * (cell + pad)
    y = pad + r * (cell + 28)
    sheet.paste(crops[slug].resize((cell, cell)), (x, y))
    draw.text((x, y + cell + 6), slug, fill="black")
montage = Path("/tmp/ack_montage.png")
sheet.save(montage)
print(f"wrote {montage}")
