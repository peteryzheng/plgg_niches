#!/usr/bin/env python3
# Measurement primitives for the nucleus-guided STAlign rerun, VENDORED (copied verbatim, no edits) so
# this pipeline does not depend on untracked files in another checkout. Sources, in the CL lipid-map
# session's sc_integration/ENVI/2_spatial_imputation/CL/segmentation_qc/ (untracked in the main checkout):
#   he_image.py                 md5 07d0432e570abb860aad47279e9dc984 : ROOT, SVS_PATH, MPP0, JPG_TO_L0, SVS
#   registration_nuclei_check.py md5 ea934ef128d2c4c794aeb8c58e036539 : DATA, GROUP, folder, DAPI_UM,
#                                  window constants, DAPI, peaks, raster, best_shift
#   lipid_detector_v2.py        md5 7d4da3868961d4b8bf38bfaadec10e01 : F = 2, UM = MPP0 * F
# Equivalence with the imported originals was checked by re-measuring a sample's windows (README).
#
# Frames: H&E level-0 px (SVS full resolution, MPP0 um/px; JPG px = level-0 / 7, the STAlign_x/_y
# frame); Xenium um (x_centroid, y_centroid); DAPI morphology_focus px = Xenium um / DAPI_UM.
from functools import lru_cache
from pathlib import Path
import numpy as np, pandas as pd, tifffile
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage import feature, filters

# ---- from he_image.py
ROOT = Path("/mnt/storage/dept/medonc/beroukhim/youyun/plgg/data/Xenium_annotations")
SVS_PATH = ROOT / "images/230918_Xenium_CytAssist_LGG2.svs"
MPP0 = 0.2522          # um per level-0 px (SVS header)

JPG_TO_L0 = 7          # JPG px * 7 = level-0 px

class SVS:
    """Region reader for one pyramid level of the SVS (default level 0); tiles decoded on demand
    and cached. Coordinates are in that level's pixels; `down` = level-0 px per level px."""
    def __init__(self, path=SVS_PATH, page=0):
        self.tif = tifffile.TiffFile(path)
        self.page = self.tif.pages[page]
        self.down = self.tif.pages[0].shape[1] / self.page.shape[1]
        self.H, self.W = self.page.shape[:2]
        self.th, self.tw = self.page.tilelength, self.page.tilewidth
        self.ntx = -(-self.W // self.tw)
        # fetch the decoder, tables and offsets ONCE: touching these page attributes can do file
        # I/O, which would move the shared file pointer between our seek and read
        self._dec, self._jt = self.page.decode, self.page.jpegtables
        self._off, self._bc = list(self.page.dataoffsets), list(self.page.databytecounts)
        self._tile = lru_cache(maxsize=4096)(self._decode)

    def _decode(self, ty, tx):
        i = ty * self.ntx + tx
        fh = self.tif.filehandle
        fh.seek(self._off[i]); data = fh.read(self._bc[i])
        return self._dec(data, i, jpegtables=self._jt)[0][0]

    def read(self, r0, c0, h, w):
        """Level-0 RGB uint8 block [r0:r0+h, c0:c0+w]; out-of-slide area filled white."""
        out = np.full((h, w, 3), 255, np.uint8)
        ra, rb = max(r0, 0), min(r0 + h, self.H); ca, cb = max(c0, 0), min(c0 + w, self.W)
        if ra >= rb or ca >= cb:
            return out
        for ty in range(ra // self.th, (rb - 1) // self.th + 1):
            for tx in range(ca // self.tw, (cb - 1) // self.tw + 1):
                t = self._tile(ty, tx)
                y0, x0 = ty * self.th, tx * self.tw
                ya, yb = max(ra, y0), min(rb, y0 + self.th); xa, xb = max(ca, x0), min(cb, x0 + self.tw)
                out[ya - r0:yb - r0, xa - c0:xb - c0] = t[ya - y0:yb - y0, xa - x0:xb - x0]
        return out

    def read_scaled(self, r0, c0, h, w, f):
        """Level-0 block downsampled by an integer factor f (block mean); h, w in level-0 px."""
        a = self.read(r0, c0, h - h % f, w - w % f).astype(np.float32)
        return a.reshape(a.shape[0] // f, f, a.shape[1] // f, f, 3).mean((1, 3))

# ---- from lipid_detector_v2.py
F = 2                      # level-0 px per output px -> 0.504 um/px
UM = MPP0 * F

# ---- from registration_nuclei_check.py
DATA = Path("/mnt/storage/dept/medonc/beroukhim/youyun/plgg/data")

# all 11 STAlign-annotated samples: H&E slide group (annotations/STalign/pathology_*_proseg.qmd) and
# Xenium folder (data/metadata/Xenium_PS.xlsx, the same lookup the STAlign notebooks use)
GROUP = {"CL42112": "LGG2", "CL67068": "LGG2", "CL87352": "LGG2", "CN51420": "LGG2", "CN53398": "LGG2",
         "GG240468": "LGG1", "GG241670": "LGG1", "PA188570": "LGG1", "PA258482": "LGG1", "PA320784": "LGG1",
         "PA328368": "LGG1"}

def folder(s):
    m = pd.read_excel(DATA / "metadata/Xenium_PS.xlsx")
    return m.loc[m["id"].astype(str) == "".join(c for c in s if c.isdigit()), "file"].iloc[0]

DAPI_UM = 0.2125

# SHIFT_UM: the first test found STAlign offsets of 75-120 um, so the search range is wide
N_WIN, WIN_UM, SHIFT_UM, MATCH_UM, NUC_MIN_DIST_UM = 30, 150.0, 150.0, 2.5, 4.0

class DAPI:
    """Region reader for the tiled JPEG2000 morphology_focus image (level 0)."""
    def __init__(self, path):
        self.tif = tifffile.TiffFile(path); p = self.tif.pages[0]
        self.H, self.W, self.th, self.tw = p.shape[0], p.shape[1], p.tilelength, p.tilewidth
        self.ntx = -(-self.W // self.tw)
        self._dec, self._off, self._bc = p.decode, list(p.dataoffsets), list(p.databytecounts)
        self._tile = lru_cache(maxsize=256)(self._decode)
    def _decode(self, ty, tx):
        i = ty * self.ntx + tx; fh = self.tif.filehandle
        fh.seek(self._off[i]); data = fh.read(self._bc[i])
        return np.squeeze(self._dec(data, i)[0])
    def read(self, r0, c0, h, w):
        out = np.zeros((h, w), np.float32)
        ra, rb = max(r0, 0), min(r0 + h, self.H); ca, cb = max(c0, 0), min(c0 + w, self.W)
        for ty in range(ra // self.th, (rb - 1) // self.th + 1):
            for tx in range(ca // self.tw, (cb - 1) // self.tw + 1):
                t = self._tile(ty, tx); y0, x0 = ty * self.th, tx * self.tw
                ya, yb = max(ra, y0), min(rb, y0 + self.th); xa, xb = max(ca, x0), min(cb, x0 + self.tw)
                out[ya - r0:yb - r0, xa - c0:xb - c0] = t[ya - y0:yb - y0, xa - x0:xb - x0]
        return out

def peaks(img, sigma_px, min_dist_px, mask=None):
    """Nucleus centres = local maxima of the smoothed image above its Otsu threshold."""
    sm = ndi.gaussian_filter(img, sigma_px)
    thr = filters.threshold_otsu(sm[mask] if mask is not None and mask.any() else sm)
    return feature.peak_local_max(sm, min_distance=int(min_dist_px), threshold_abs=thr, exclude_border=False)

def raster(p, lo, size, sigma=1.5):
    """Point set (um) -> blurred 1 um image of side `size` whose origin is `lo` um."""
    a = np.zeros((size, size)); q = np.round(p - lo).astype(int)
    q = q[(q >= 0).all(1) & (q < size).all(1)]
    a[q[:, 0], q[:, 1]] = 1
    return ndi.gaussian_filter(a, sigma)

def best_shift(a, b):
    """a: H&E nuclei in the window (um, 0..WIN_UM); b: mapped DAPI nuclei (um, window frame, padded by
    SHIFT_UM). Coarse: FFT cross-correlation of blurred point images over +-SHIFT_UM (1 um), with a
    z-score of the peak against the whole correlation surface (confidence that the shift is real, not
    a chance alignment). Fine: +-1.5 um at 0.25 um maximising the share of H&E nuclei with a DAPI partner
    within MATCH_UM. Returns (dy, dx, z, match@0, match@best, median residual of matched pairs)."""
    from scipy.signal import fftconvolve
    A = raster(a, 0, int(WIN_UM)); B = raster(b, -SHIFT_UM, int(WIN_UM + 2 * SHIFT_UM))
    cc = fftconvolve(B, A[::-1, ::-1], mode="valid")
    iy, ix = np.unravel_index(cc.argmax(), cc.shape)
    z = (cc.max() - cc.mean()) / cc.std()
    tree = cKDTree(b)
    frac = lambda d: (tree.query(a + d, distance_upper_bound=MATCH_UM)[0] < np.inf).mean()
    d0 = np.array([iy - SHIFT_UM, ix - SHIFT_UM], float)
    fine = np.arange(-1.5, 1.51, 0.25)
    sf = np.array([[frac(d0 + [fy, fx]) for fx in fine] for fy in fine])
    jy, jx = np.unravel_index(sf.argmax(), sf.shape); d = d0 + [fine[jy], fine[jx]]
    dist = tree.query(a + d, distance_upper_bound=MATCH_UM)[0]
    res = float(np.median(dist[np.isfinite(dist)])) if np.isfinite(dist).any() else np.nan
    return d[0], d[1], z, frac(np.zeros(2)), sf.max(), res
