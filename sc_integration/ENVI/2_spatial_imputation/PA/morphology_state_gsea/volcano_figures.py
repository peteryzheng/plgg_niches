#!/usr/bin/env python3
"""Volcano plots for the morphology-state contrast, after the scale/filter fix.
Marks which genes are on the measured 266-gene panel (directly validated) vs
imputation-dependent, and shows the before/after effect of the pipeline fix."""
from pathlib import Path
import numpy as np, pandas as pd, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns
sns.set_theme(style="white", context="talk", font_scale=1.0)
plt.rcParams.update({"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"savefig.dpi":150})

B=Path("/mnt/storage/dept/medonc/beroukhim/youyun/plgg/code/sc_integration/ENVI/2_spatial_imputation")
SPB=Path("/ncache0/3638588.1.all.q/claude-68628/-mnt-storage-dept-medonc-beroukhim-youyun-plgg-code/992160f9-a942-451f-8e36-4a0f6cfb4f2d/scratchpad/prefix_baseline")
DISP={"astro_glial":"Astrocytes","opc":"OPC","myeloid":"Myeloid","neuronal":"Neuronal",
      "radial_glia":"Radial glia","oligo_mature":"Oligodendrocytes"}
CFG={"PA":dict(col="log2FC_loose_vs_compact",old="logFC_loose_vs_compact",cts=["astro_glial","opc","myeloid"],
               lo="compact",hi="myxoid",h5="188570_pilocytic_astrocytoma_spatial_envi.h5ad"),
     "CL":dict(col="log2FC_lipid_vs_nonlipid",old="logFC_lipid_vs_nonlipid",
               cts=["neuronal","astro_glial","opc","radial_glia","oligo_mature","myeloid"],
               lo="non-lipidized",hi="lipidized",h5="42112_cerebellar_liponeurocytoma_spatial_envi.h5ad")}
ENVI=Path("/mnt/storage/dept/medonc/beroukhim/youyun/plgg/data/sc_integration/ENVI_results_091126")

def panel(h5):
    with h5py.File(ENVI/h5,"r") as h:
        vx=h["var"].attrs.get("_index","_index"); vx=vx.decode() if isinstance(vx,bytes) else vx
        return {(g.decode() if isinstance(g,bytes) else str(g)).upper() for g in h["var"][vx][:]}

def volcano(coh, out):
    c=CFG[coh]; de=pd.read_csv(B/coh/"morphology_state_gsea/morphology_state_limma_de.csv")
    meas=panel(c["h5"]); cts=[x for x in c["cts"] if x in set(de.cell_type)]
    # one x-limit for every panel: lineages are only comparable on a shared scale
    xlim=max(0.4, de[de.cell_type.isin(cts)][c["col"]].abs().quantile(.9995)*1.3)
    ncol=3; nrow=int(np.ceil(len(cts)/ncol))
    fig,axes=plt.subplots(nrow,ncol,figsize=(7.2*ncol,6.2*nrow),squeeze=False)
    for ax,ct in zip(axes.flat,cts):
        d=de[de.cell_type==ct].copy(); d["nlp"]=-np.log10(d.pvalue.clip(lower=1e-300))
        sig=d[d.fdr<0.25]; ns=d[d.fdr>=0.25]
        ax.scatter(ns[c["col"]],ns.nlp,s=6,c="#d9d9d9",lw=0,alpha=.7)
        if len(sig):
            up=sig[sig[c["col"]]>0]; dn=sig[sig[c["col"]]<0]
            ax.scatter(dn[c["col"]],dn.nlp,s=44,c="#4C72B0",lw=.5,ec="black",zorder=3,label=f"down in {c['hi']} (FDR<0.25)")
            ax.scatter(up[c["col"]],up.nlp,s=44,c="#C44E52",lw=.5,ec="black",zorder=3,label=f"up in {c['hi']} (FDR<0.25)")
            mm=sig[sig.gene.isin(meas)]
            ax.scatter(mm[c["col"]],mm.nlp,s=150,facecolors="none",edgecolors="#F0A800",lw=2.4,zorder=4,
                       label="on measured 266-gene panel")
            # stagger labels so crowded regions stay readable
            lab=sig.nsmallest(8,"fdr").sort_values("nlp",ascending=False).reset_index(drop=True)
            for i,r in lab.iterrows():
                dx = 8 if r[c["col"]] >= 0 else -8
                dy = 10 if i % 2 == 0 else -12
                ax.annotate(r.gene,(r[c["col"]],r.nlp),fontsize=10,fontweight="bold",
                            xytext=(dx,dy),textcoords="offset points",
                            ha="left" if dx>0 else "right",
                            arrowprops=dict(arrowstyle="-",lw=.6,color="0.55",shrinkA=0,shrinkB=2))
        # FDR<0.25 threshold line, if any gene crosses it
        if len(sig):
            thr=-np.log10(sig.pvalue.max())
            ax.axhline(thr,color="0.55",ls="--",lw=1)
            ax.text(0.985,thr,"FDR 0.25",fontsize=10,va="bottom",ha="right",color="0.45",
                    transform=ax.get_yaxis_transform())
        ax.axvline(0,color="0.35",lw=1)
        ax.set_xlim(-xlim,xlim)
        ax.set_title(f"{DISP.get(ct,ct)}\n{len(sig)} genes at FDR<0.25  ·  {len(d)} tested",fontsize=14)
        ax.set_xlabel(f"log2FC  ({c['lo']} → {c['hi']})"); ax.set_ylabel("-log10 p")
    for ax in axes.flat[len(cts):]: ax.axis("off")
    h,l=axes.flat[0].get_legend_handles_labels()
    if h: fig.legend(h,l,loc="lower center",ncol=3,frameon=False,fontsize=13,bbox_to_anchor=(0.5,-0.04))
    fig.suptitle(f"{coh}: gene-level DE after the pipeline fix  ·  log2 transform + expression filter",
                 fontsize=19,fontweight="bold",y=1.0)
    fig.tight_layout(); fig.savefig(out,bbox_inches="tight"); print("wrote",out)

def before_after(out):
    fig,axes=plt.subplots(2,3,figsize=(20,11))
    for row,(coh,lbl) in enumerate([("PA","PA  (myxoid vs compact)"),("CL","CL  (lipidized vs non)")]):
        c=CFG[coh]
        old=pd.read_csv(SPB/f"{coh}_de_old.csv"); new=pd.read_csv(B/coh/"morphology_state_gsea/morphology_state_limma_de.csv")
        ct=c["cts"][2] if coh=="PA" else "myeloid"
        for j,(t,d,col) in enumerate([("BEFORE (linear input)",old,c["old"]),("AFTER (log2 + filter)",new,c["col"])]):
            ax=axes[row][j]; dd=d[d.cell_type==ct].copy(); dd["nlp"]=-np.log10(dd.pvalue.clip(lower=1e-300))
            sig=dd[dd.fdr<0.25]
            ax.scatter(dd[col],dd.nlp,s=5,c="#d9d9d9",lw=0,alpha=.7)
            if len(sig): ax.scatter(sig[col],sig.nlp,s=34,c="#C44E52",lw=.4,ec="black",zorder=3)
            ax.axvline(0,color="0.35",lw=1)
            ax.set_title(f"{lbl} · {DISP[ct]}\n{t}\n{len(sig)} genes FDR<0.25 · min p {dd.pvalue.min():.1e}",fontsize=12)
            ax.set_xlabel("effect (linear diff)" if j==0 else "log2FC"); ax.set_ylabel("-log10 p")
        # p-value histogram
        ax=axes[row][2]
        for t,d,cl in [("before",old,"#999999"),("after",new,"#C44E52")]:
            ax.hist(d[d.cell_type==ct].pvalue,bins=40,density=True,alpha=.55,color=cl,label=t)
        ax.axhline(1,color="black",ls="--",lw=1.4,label="uniform null")
        ax.set_xlabel("raw p-value"); ax.set_ylabel("density")
        ax.set_title(f"{lbl} · {DISP[ct]}\np-value distribution",fontsize=12); ax.legend(fontsize=11)
    fig.suptitle("Pipeline fix: linear-scale input → log2 + expression filter",fontsize=19,fontweight="bold",y=1.0)
    fig.tight_layout(); fig.savefig(out,bbox_inches="tight"); print("wrote",out)

volcano("PA", B/"PA/morphology_state_gsea/figures/fig_deg_volcano.png")
volcano("CL", B/"CL/morphology_state_gsea/figures/fig_deg_volcano.png")
before_after(B/"PA/morphology_state_gsea/figures/fig_pipeline_fix_before_after.png")
