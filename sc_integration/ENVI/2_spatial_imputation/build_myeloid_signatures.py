# Builds myeloid_signatures.csv (Signature,Gene) from three literature sources
# (raw supplementary tables kept in myeloid_signature_sources/ for provenance)
# plus three public MSigDB Hallmark gene sets, for scoring myeloid cell states
# against glial (OPC/Astrocyte) MAPK exposure in 2_MAPK_adj.qmd. See the plan
# rationale for why each signature was chosen and how it was filtered.
import pandas as pd

SRC_DIR = "/mnt/storage/dept/medonc/beroukhim/youyun/plgg/code/sc_integration/ENVI/2_spatial_imputation/myeloid_signature_sources/"
OUT = "/mnt/storage/dept/medonc/beroukhim/youyun/plgg/code/sc_integration/ENVI/2_spatial_imputation/myeloid_signatures.csv"

rows = []


def add(signature, genes):
    genes = sorted(set(genes))
    rows.extend([(signature, g) for g in genes])
    print(f"{signature}: {len(genes)} genes")


# --- 1. Full MSigDB Hallmark collection (all 50 gene sets), via msigdbr ---
# Broadened from the original 3 immune-focused Hallmark sets to the full
# collection as a specificity check: do myeloid/immune-relevant gene sets show
# larger effects (R^2) than generic/unrelated pathways (cell cycle, metabolism,
# DNA repair, etc.), or does everything look similar given ~2,000-9,000 cells
# per patient (high power inflates significance regardless of effect size --
# see conversation)?
import subprocess
hallmark_csv = SRC_DIR + "../hallmark_all50.csv"
subprocess.run([
    "Rscript", "-e",
    f"""
    library(msigdbr); library(data.table)
    h = as.data.table(msigdbr(species = 'Homo sapiens', collection = 'H'))
    fwrite(unique(h[, .(gs_name, gene_symbol)]), '{hallmark_csv}')
    """
], check=True)
hallmark = pd.read_csv(hallmark_csv)
for gs_name, grp in hallmark.groupby("gs_name"):
    add(gs_name, grp.gene_symbol.tolist())

# --- 2. MG-Act (Tripathi et al. JCI 2024, PMC11444160) -- BRAF-fusion PA specific ---
add("MG_Act_Tripathi", """C1QC, APOE, C1QB, C1QA, CST3, APOC1, HLA-DRA, C3, AIF1, CD74, MARCKS, FTL, TREM2, SPP1,
APOC2, CD68, TYROBP, HLA-DRB5, HLA-DPA1, SPI1, NPC2, CTSB, TMEM176B, SERPINA1, HLA-DRB1, FCER1G,
IFI30, GSN, MS4A6A, CSF1R, GPR34, LY86, CD14, VSIG4, HLA-DPB1, TUBA1B, SCIN""".replace("\n", " ").split(", "))

# --- 3/4. Andrade/Jabado (Nat Immunol 2025, PMID 40954250) Supplementary Table 5:
#     TIM3+ vs TIM3- microglia DEGs, user-uploaded. Filter to genes significantly
#     (BH q<0.05) and substantially (|log2FC|>0.5) different between the groups.
deg = pd.read_excel(SRC_DIR + "41590_2025_2268_MOESM6_ESM.xlsx", sheet_name="DEG")
sig = deg[deg.p_val_adj < 0.05]
add("PA_TIM3pos_Microglia", sig[sig.avg_log2FC > 0.5].gene.tolist())
add("PA_TIM3neg_Microglia", sig[sig.avg_log2FC < -0.5].gene.tolist())

# --- 5. Bernstein/Miller (Nature 2025, PMID 40011771) myeloid cNMF programs,
#     user-uploaded. Top-100-gene lists for the two programs explicitly named
#     "Immunosuppressive".
bern = pd.read_excel(SRC_DIR + "41586_2025_8633_MOESM4_ESM.xlsx", sheet_name="Myeloid Top 100 genes per prog.")
add("Bernstein_Complement_Immunosuppressive", bern["Complement Immunosuppressive"].dropna().tolist())
add("Bernstein_Scavenger_Immunosuppressive", bern["Scavenger Immunosuppressive"].dropna().tolist())

# --- 6. Yu/Tabar (Cancer Cell, in press) Table S2, user-uploaded: FOSL2 regulon
#     (small/high-confidence subset) downstream of their proposed pathogenic
#     glioma-associated-macrophage driver TF. Adult glioma progression context.
yu = pd.read_excel(SRC_DIR + "mmc2.xlsx", header=1)
add("FOSL2_Regulon_Small", yu.iloc[:, 1].dropna().tolist())

out_df = pd.DataFrame(rows, columns=["Signature", "Gene"])
out_df.to_csv(OUT, index=False)
print(f"\nWrote {len(out_df)} rows, {out_df.Signature.nunique()} signatures, to {OUT}")
print(out_df.Signature.value_counts())
