# PA Imputation Gene-Universe Audit (022626)

## Context

This audit was run to resolve a discrepancy between earlier PA overlap results showing a `718`-gene shared imputation universe versus newly observed `18077`-gene full overlap.

- Data source directory:
  - `/Users/youyun/Documents/HMS/PhD/beroukhimlab/dfci_mount/youyun/plgg/data/sc_integration/ENVI_results_022626`
- Scope:
  - PA spatial ENVI objects only (4 samples).

## Files Audited

- `/Users/youyun/Documents/HMS/PhD/beroukhimlab/dfci_mount/youyun/plgg/data/sc_integration/ENVI_results_022626/188570_pilocytic_astrocytoma_spatial_envi.h5ad`
- `/Users/youyun/Documents/HMS/PhD/beroukhimlab/dfci_mount/youyun/plgg/data/sc_integration/ENVI_results_022626/258482_pilocytic_astrocytoma_spatial_envi.h5ad`
- `/Users/youyun/Documents/HMS/PhD/beroukhimlab/dfci_mount/youyun/plgg/data/sc_integration/ENVI_results_022626/320784_pilocytic_astrocytoma_spatial_envi.h5ad`
- `/Users/youyun/Documents/HMS/PhD/beroukhimlab/dfci_mount/youyun/plgg/data/sc_integration/ENVI_results_022626/328368_pilocytic_astrocytoma_spatial_envi.h5ad`

## Methods

Gene extraction logic:

- `var` genes:
  - `h5ad["var"]["_index"]` (fallback to `h5ad["var"]["index"]` if needed)
- Imputation genes:
  - `h5ad["obsm"]["imputation"].attrs["column-order"]`
  - fallback to `h5ad["obsm"]["imputation"].keys()` if no `column-order`
  - final fallback to `var` genes if `imputation` is not stored as a group

Normalization and overlap metrics:

- Gene symbols converted to uppercase
- Deduplication preserving first-seen order
- Per-sample gene counts recorded
- Global overlap:
  - intersection size
  - union size
  - Jaccard index (`intersection / union`)

## Results

### `var` gene universe

| Sample | n_genes |
|---|---:|
| 188570_pilocytic_astrocytoma | 266 |
| 258482_pilocytic_astrocytoma | 266 |
| 320784_pilocytic_astrocytoma | 266 |
| 328368_pilocytic_astrocytoma | 266 |

- intersection = `266`
- union = `266`
- jaccard = `1.000000`

### Imputation gene universe

| Sample | n_genes |
|---|---:|
| 188570_pilocytic_astrocytoma | 18077 |
| 258482_pilocytic_astrocytoma | 18077 |
| 320784_pilocytic_astrocytoma | 18077 |
| 328368_pilocytic_astrocytoma | 18077 |

- imputation_intersection = `18077`
- imputation_union = `18077`
- imputation_jaccard = `1.000000`

## Interpretation

The current redownloaded PA spatial ENVI objects are fully aligned for imputation genes. The shared imputation universe is now the full `18077` genes, not `718`.

Caution: the earlier mismatch came from mixed filename/object sets (space-named and underscore-named variants) in the data directory. The current directory state for these four PA files is now resolved and consistent.

## Downstream Implications for V2

- Expected effect:
  - `myeloid_pathway_v2_prepare.py` should now report `shared_universe_n: 18077` for PA.
- Practical next step:
  - rerun V2 Python prepare and R limma/fgsea stage to regenerate outputs under the updated shared universe.

## Reproducibility Snippet

Run from:

- `/Users/youyun/Documents/HMS/PhD/beroukhimlab/dfci_mount/youyun/plgg/data/sc_integration/ENVI_results_022626`

```python
import glob, os, h5py
from itertools import combinations

files = sorted(glob.glob("*astrocytoma_spatial_envi.h5ad"))

def to_str_list(arr):
    out = []
    for x in arr:
        if isinstance(x, bytes):
            out.append(x.decode("utf-8"))
        else:
            out.append(str(x))
    return out

def read_var(h):
    var = h["var"]
    key = "_index" if "_index" in var else "index"
    return list(dict.fromkeys([g.upper() for g in to_str_list(var[key][()])]))

def read_imputation(h):
    imp = h["obsm"]["imputation"]
    if isinstance(imp, h5py.Group):
        if "column-order" in imp.attrs:
            genes = to_str_list(imp.attrs["column-order"])
        else:
            genes = list(imp.keys())
    else:
        genes = read_var(h)
    return list(dict.fromkeys([g.upper() for g in genes]))

var_by = {}
imp_by = {}
for f in files:
    sid = os.path.basename(f).replace("_spatial_envi.h5ad", "")
    with h5py.File(f, "r") as h:
        var_by[sid] = read_var(h)
        imp_by[sid] = read_imputation(h)

def summarize(label, d):
    sets = {k: set(v) for k, v in d.items()}
    ids = sorted(sets)
    print(f"\n=== {label} ===")
    for s in ids:
        print(f"{s}\tn={len(d[s])}\tunique={len(set(d[s]))}")
    inter = set.intersection(*[sets[s] for s in ids])
    union = set.union(*[sets[s] for s in ids])
    print(f"{label}_intersection\t{len(inter)}")
    print(f"{label}_union\t{len(union)}")
    print(f"{label}_jaccard\t{len(inter)/len(union):.6f}")

summarize("var", var_by)
summarize("imputation", imp_by)
```
