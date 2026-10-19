#!/usr/bin/env python3
# Imputation QC: how much regional signal does ENVI erase?
# Compares the measured-panel fit against the imputed-genome fit on the shared genes,
# using identical statistics on both sides so the difference isolates imputation.
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
B = HERE.parent
COH = {"CL": ("lipid", "nonlipid"), "PA": ("loose", "compact")}
DISP = {"astro_glial":"Astrocytes","opc":"OPC","myeloid":"Myeloid","neuronal":"Neuronal",
        "radial_glia":"Radial glia","oligo_mature":"Oligodendrocytes"}

summary = []
for coh, (hi, lo) in COH.items():
    meas = pd.read_csv(HERE / f"measured_{coh}_limma_de.csv")
    imp  = pd.read_csv(B / coh / "morphology_state_gsea/morphology_state_limma_de.csv")
    imp = imp.rename(columns={f"log2FC_{hi}_vs_{lo}": "log2FC"})
    m = meas.merge(imp, on=["cell_type", "gene"], suffixes=("_meas", "_imp"))
    print("=" * 96)
    print(f"{coh}: {m.gene.nunique()} shared genes x {m.cell_type.nunique()} lineages "
          f"({len(m)} paired tests)")
    print("=" * 96)
    rows = []
    for ct, g in m.groupby("cell_type"):
        # effect attenuation on the same genes
        att = g.log2FC_meas.abs().median() / max(g.log2FC_imp.abs().median(), 1e-12)
        # do they even agree on direction / ranking?
        conc = (np.sign(g.log2FC_meas) == np.sign(g.log2FC_imp)).mean()
        rho = stats.spearmanr(g.log2FC_meas, g.log2FC_imp).correlation
        rows.append({"lineage": DISP.get(ct, ct),
                     "median_abs_log2FC_measured": g.log2FC_meas.abs().median(),
                     "median_abs_log2FC_imputed": g.log2FC_imp.abs().median(),
                     "attenuation_x": att,
                     "direction_concordance": conc,
                     "spearman_rho": rho,
                     "n_sig_measured": int((g.fdr_meas < 0.25).sum()),
                     "n_sig_imputed": int((g.fdr_imp < 0.25).sum())})
    t = pd.DataFrame(rows)
    print(t.round(4).to_string(index=False))
    print(f"\n  cohort median attenuation: {t.attenuation_x.median():.1f}x   "
          f"mean direction concordance: {t.direction_concordance.mean():.2f}   "
          f"mean Spearman rho: {t.spearman_rho.mean():.3f}")
    # recovery: what happens to measured-significant genes in the imputed arm?
    sig = m[m.fdr_meas < 0.25]
    if len(sig):
        print(f"\n  Of {len(sig)} measured-significant genes, in the IMPUTED arm:")
        print(f"    still FDR<0.25 : {int((sig.fdr_imp < 0.25).sum())}")
        print(f"    raw p<0.05     : {int((sig.pvalue_imp < 0.05).sum())}")
        print(f"    same direction : {int((np.sign(sig.log2FC_meas)==np.sign(sig.log2FC_imp)).sum())} / {len(sig)}")
        print(f"    median |log2FC| measured {sig.log2FC_meas.abs().median():.4f} "
              f"-> imputed {sig.log2FC_imp.abs().median():.4f}")
    t.insert(0, "cohort", coh); summary.append(t)
    print()

pd.concat(summary).to_csv(HERE / "measured_vs_imputed_attenuation.tsv", sep="\t", index=False)
print("wrote measured_vs_imputed_attenuation.tsv")
