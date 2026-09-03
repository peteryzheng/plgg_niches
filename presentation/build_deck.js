// Deck generator v2 — Spatial Architecture of Pediatric Brain Tumors.
// Bigger fonts, restructured intro, embedded H&E + cartoons, figure-drop frames for qmd figures.
const pptxgen = require("pptxgenjs");
const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.33 x 7.5
pres.author = "Peter Youyun Zheng";
pres.title = "Spatial Architecture of Pediatric Brain Tumors";

const NAVY="1F2A44", INK="1F2A37", ACCENT="A4243B", STEEL="2E5A87", MUTED="6B7280",
      CARD="F3F5F8", BORDER="D5DBE3", WHITE="FFFFFF", ICE="C7D2E4",
      GREEN="2C7A4B", AMBER="C98A2B", NUC="3A3F4A";
const SERIF="Cambria", SANS="Calibri";
const PAL=[ACCENT,STEEL,GREEN,AMBER];
const GENECOLS=["A4243B","2E5A87","2C7A4B","C98A2B","6D5199","1C7293","B5651D","4C8577","8A2846"];
const path=require("path");
const A=path.join(__dirname,"assets")+path.sep;   // assets/ next to this script
const R=(a,b)=>a+Math.random()*(b-a);

// ---------- shared helpers ----------
const dot=(s,x,y,d,c,tr=0)=>s.addShape("ellipse",{x:x-d/2,y:y-d/2,w:d,h:d,fill:{color:c,transparency:tr},line:{color:c,width:0.4}});
function line2(s,x1,y1,x2,y2,o={}){
  const x=Math.min(x1,x2),y=Math.min(y1,y2),w=Math.abs(x2-x1),h=Math.abs(y2-y1);
  s.addShape("line",{x,y,w,h,flipV:((x2-x1)*(y2-y1))<0,line:Object.assign({color:MUTED,width:1},o)});
}
function motifDots(s,n,x0,x1,y0,y1){ for(let i=0;i<n;i++){const c=[ACCENT,STEEL,"8892A6"][i%3],d=R(0.05,0.14);
  s.addShape("ellipse",{x:R(x0,x1),y:R(y0,y1),w:d,h:d,fill:{color:c,transparency:R(20,62)},line:{color:c,width:0.5,transparency:40}});}}
function head(s,kicker,title){
  let ty=0.62;
  if(kicker){ s.addShape("ellipse",{x:0.6,y:0.66,w:0.14,h:0.14,fill:{color:ACCENT},line:{color:ACCENT,width:0.5}});
    s.addText(kicker.toUpperCase(),{x:0.84,y:0.56,w:11.6,h:0.34,fontFace:SANS,fontSize:14,bold:true,color:ACCENT,charSpacing:2,margin:0,valign:"middle"});
    ty=0.98; }
  s.addText(title,{x:0.57,y:ty,w:12.2,h:0.9,fontFace:SERIF,fontSize:34,bold:true,color:INK,margin:0,valign:"top"});
}
function foot(s){ s.slideNumber={x:12.5,y:7.02,w:0.5,h:0.3,fontFace:SANS,fontSize:10,color:MUTED,align:"right"}; }
function bullets(s,items,x,y,w,h,sz){
  s.addText(items.map(t=>({text:t,options:{bullet:{code:"2022",indent:20},breakLine:true,paraSpaceAfter:14,fontFace:SANS,fontSize:sz||22,color:INK}})),
    {x,y,w,h,margin:0,valign:"top",lineSpacingMultiple:1.04});
}
function figPh(s,x,y,w,h,label,src){
  s.addShape("roundRect",{x,y,w,h,rectRadius:0.04,fill:{color:CARD},line:{color:BORDER,width:1.25,dashType:"dash"}});
  s.addText([{text:"▸ paste figure\n",options:{fontSize:14,bold:true,color:MUTED,fontFace:SANS}},
    {text:label+"\n",options:{fontSize:16,bold:true,color:STEEL,fontFace:SANS}},
    {text:src,options:{fontSize:12,italic:true,color:MUTED,fontFace:SANS}}],
    {x,y,w,h,align:"center",valign:"middle",margin:8,lineSpacingMultiple:1.15});
}
function img(s,file,x,y,w,h){
  s.addShape("roundRect",{x:x-0.03,y:y-0.03,w:w+0.06,h:h+0.06,rectRadius:0.03,fill:{color:WHITE},line:{color:BORDER,width:1}});
  s.addImage({path:A+file,x,y,w,h});
}
// center a figure inside a box, preserving its native aspect ratio
function imgFit(s,file,natW,natH,bx,by,bw,bh){
  const ar=natW/natH; let w=bw, h=bw/ar; if(h>bh){h=bh; w=bh*ar;}
  s.addImage({path:A+file,x:bx+(bw-w)/2,y:by+(bh-h)/2,w,h});
}

// ---------- pipeline cartoons (ported from preview.js) ----------
function cProseg(s,ax,ay,aw,ah){
  const nuclei=[[0.5,0.55],[1.28,0.5],[0.92,1.12]];
  const own=[{c:ACCENT,pts:[[0.30,0.42],[0.62,0.46],[0.42,0.78],[0.70,0.70]]},
    {c:STEEL,pts:[[1.10,0.36],[1.42,0.40],[1.20,0.66],[1.46,0.62]]},
    {c:GREEN,pts:[[0.76,1.00],[1.04,1.02],[0.82,1.30],[1.08,1.24]]}];
  const spur=[[0.92,0.54,AMBER],[1.04,0.90,ACCENT],[0.64,0.92,STEEL]];
  const sc=aw/3.05; // scale local (~3.05 wide design) to area width portion
  const P=(ox,u,v)=>[ox+u*sc, ay+v*sc];
  function panel(ox,mode){
    if(mode==="b") nuclei.forEach(n=>{const c=P(ox,n[0],n[1]);s.addShape("ellipse",{x:c[0]-0.45*sc,y:c[1]-0.45*sc,w:0.9*sc,h:0.9*sc,fill:{color:"E2E5EA",transparency:25},line:{color:"9AA1AC",width:1}});});
    else {const el=[[0.5,0.57,0.64,0.46,-18],[1.28,0.5,0.6,0.44,14],[0.92,1.14,0.64,0.48,8]];
      el.forEach(e=>{const c=P(ox,e[0],e[1]);s.addShape("ellipse",{x:c[0]-e[2]/2*sc,y:c[1]-e[3]/2*sc,w:e[2]*sc,h:e[3]*sc,rotate:e[4],fill:{color:"DCE7F2",transparency:15},line:{color:STEEL,width:1.25}});});}
    nuclei.forEach(n=>{const c=P(ox,n[0],n[1]);dot(s,c[0],c[1],0.1*sc,NUC);});
    own.forEach(g=>g.pts.forEach(p=>{const c=P(ox,p[0],p[1]);dot(s,c[0],c[1],0.08*sc,g.c);}));
    spur.forEach(p=>{const c=P(ox,p[0],p[1]);dot(s,c[0],c[1],0.08*sc,p[2]);});
  }
  panel(ax,"b");
  line2(s,ax+aw*0.45,ay+ah*0.5,ax+aw*0.55,ay+ah*0.5,{color:STEEL,width:2,endArrowType:"triangle"});
  panel(ax+aw*0.55,"a");
}
function cBanksy(s,ax,ay,aw,ah){
  const sx=aw/5.9, sy=ah/1.6, x0=ax+0.34*sx, y0=ay+0.28*sy, cols=8, rows=3, dd=0.15*sx;
  const bg=(i,j)=>PAL[(i*3+j*2+(i%2))%4];
  const motif={"0,0":ACCENT,"1,0":GREEN,"0,1":STEEL,"1,1":AMBER};
  const anchors=[[1,0],[4,1],[6,0]];
  const mc=(i,j)=>{for(const[ai,aj]of anchors){const di=i-ai,dj=j-aj;if(di>=0&&di<=1&&dj>=0&&dj<=1)return motif[`${di},${dj}`];}return null;};
  for(let i=0;i<cols;i++)for(let j=0;j<rows;j++)dot(s,x0+i*0.62*sx,y0+j*0.5*sy,dd,mc(i,j)||bg(i,j));
  anchors.forEach(a=>{const cx=x0+(a[0]+0.5)*0.62*sx,cy=y0+(a[1]+0.5)*0.5*sy;
    s.addShape("ellipse",{x:cx-0.5*sx,y:cy-0.46*sy,w:1.0*sx,h:0.92*sy,fill:{color:WHITE,transparency:100},line:{color:NAVY,width:1.4,dashType:"dash"}});});
}
function cEnvi(s,ax,ay,aw,ah){
  const sc=Math.min(aw/5.0,ah/1.6), rowTr=[12,40,58,28];
  const mat=(ox,cols,cw)=>{for(let i=0;i<cols;i++)for(let j=0;j<4;j++)
    s.addShape("rect",{x:ox+i*(cw+0.05*sc),y:ay+0.15*sc+j*(0.16*sc+0.05*sc),w:cw,h:0.16*sc,fill:{color:GENECOLS[i%GENECOLS.length],transparency:rowTr[j]},line:{color:WHITE,width:0.75}});};
  mat(ax+0.15*sc,3,0.2*sc);
  line2(s,ax+1.35*sc,ay+0.62*sc,ax+1.95*sc,ay+0.62*sc,{color:STEEL,width:2,endArrowType:"triangle"});
  mat(ax+2.15*sc,9,0.2*sc);
}
function cStalign(s,ax,ay,aw,ah){
  const th=40,rad=th*Math.PI/180,c=Math.cos(rad),sn=Math.sin(rad);
  const sc=Math.min(aw/5.0,ah/1.55), hw=0.5*sc, hh=0.5*sc;
  const Hc=[ax+0.95*sc,ay+ah*0.55], Xc=[ax+3.7*sc,ay+ah*0.55];
  const rot=(ct,u,v)=>[ct[0]+u*c-v*sn, ct[1]+u*sn+v*c];
  const ann=[[-0.26*sc,-0.16*sc],[0.28*sc,-0.04*sc],[0.02*sc,0.24*sc]];
  s.addShape("roundRect",{x:Hc[0]-hw,y:Hc[1]-hh,w:2*hw,h:2*hh,rectRadius:0.04,rotate:th,fill:{color:"F7EEF0"},line:{color:ACCENT,width:1.25}});
  s.addShape("roundRect",{x:Xc[0]-hw,y:Xc[1]-hh,w:2*hw,h:2*hh,rectRadius:0.04,rotate:th,fill:{color:WHITE},line:{color:STEEL,width:1.25}});
  const gn=6;
  for(let gi=0;gi<gn;gi++)for(let gj=0;gj<gn;gj++){const u=(-hw+0.09*sc)+(2*(hw-0.09*sc))*gi/(gn-1),v=(-hh+0.09*sc)+(2*(hh-0.09*sc))*gj/(gn-1);
    const p=rot(Xc,u,v);dot(s,p[0],p[1],0.06*sc,PAL[(gi*2+gj)%4],20);}
  ann.forEach(a=>{const ph=rot(Hc,a[0],a[1]),px=rot(Xc,a[0],a[1]);
    [ph,px].forEach(q=>s.addShape("ellipse",{x:q[0]-0.09*sc,y:q[1]-0.09*sc,w:0.18*sc,h:0.18*sc,fill:{color:WHITE,transparency:100},line:{color:ACCENT,width:1.5}}));});
  [[-hw,-hh],[hw,-hh],[-hw,hh],[hw,hh]].forEach(cn=>{const ph=rot(Hc,cn[0],cn[1]),px=rot(Xc,cn[0],cn[1]);
    dot(s,ph[0],ph[1],0.09*sc,NUC);dot(s,px[0],px[1],0.09*sc,NUC);
    line2(s,ph[0],ph[1],px[0],px[1],{color:MUTED,width:1,dashType:"dash",endArrowType:"triangle"});});
}

// ============================================================ 1. TITLE
(()=>{const s=pres.addSlide(); s.background={color:NAVY};
  motifDots(s,52,8.6,13.0,3.4,7.2);
  s.addText("Spatial Architecture of\nPediatric Brain Tumors",
    {x:0.8,y:2.3,w:9.6,h:2.2,fontFace:SERIF,fontSize:48,bold:true,color:WHITE,margin:0,lineSpacingMultiple:1.02});
  s.addText([{text:"Peter Youyun Zheng",options:{bold:true,color:WHITE,fontSize:16,fontFace:SANS,breakLine:true}},
    {text:"Beroukhim Lab · Dana-Farber Cancer Institute / Broad Institute",options:{color:ICE,fontSize:14,fontFace:SANS,breakLine:true}},
    {text:"September 2026",options:{color:MUTED,fontSize:12.5,fontFace:SANS}}],
    {x:0.8,y:5.5,w:9,h:1.3,margin:0,lineSpacingMultiple:1.3});
})();

// ============================================================ 2. MOTIVATION
(()=>{const s=pres.addSlide(); head(s,null,"Pediatric brain tumors are spatially structured");
  bullets(s,[
    "Pilocytic astrocytoma is biphasic — dense compact-fibrillary tissue abuts loose myxoid tissue (H&E, right).",
    "Across these tumors cells are arranged, not mixed — glia, neurons, and immune cells hold distinct neighborhoods.",
    "Bulk and dissociated single-cell assays discard that arrangement; in-situ transcriptomics preserves it.",
  ],0.6,1.9,6.3,4.4,22);
  img(s,"pa_biphasic_slide_annot.jpg",7.75,1.75,4.6,4.6);  // square crop, room for caption below
  s.addText("Pilocytic astrocytoma H&E — green = compact fibrillary, blue = loose myxoid (pathologist-annotated).",
    {x:7.75,y:6.5,w:4.6,h:0.6,fontFace:SANS,fontSize:12,italic:true,color:MUTED,align:"center",margin:0,valign:"top",lineSpacingMultiple:1.05});
  foot(s);
})();

// ============================================================ 3. COHORT
(()=>{const s=pres.addSlide(); head(s,null,"The cohort");
  const stats=[["11","tumors"],["4","tumor histologies"],["266","gene brain panel"]];
  stats.forEach((st,i)=>{const x=0.6+i*4.05;
    s.addText(st[0],{x,y:1.8,w:2.4,h:1.0,fontFace:SERIF,fontSize:56,bold:true,color:ACCENT,margin:0});
    s.addText(st[1],{x:x+0.06,y:2.85,w:3.8,h:0.4,fontFace:SANS,fontSize:15,color:MUTED,margin:0});});
  const rows=[
    [{text:"Histology",options:{bold:true}},{text:"n",options:{bold:true,align:"center"}},{text:"Driver alteration(s)",options:{bold:true}},{text:"Location",options:{bold:true}}],
    ["Pilocytic astrocytoma","4","KIAA1549::BRAF, FAM131B::BRAF, RAF1::TOP2B, FGFR1-ITD","Cerebellum / temporal"],
    ["Ganglioglioma","2","BRAF-V600E, FGFR1-ITD","Cortical"],
    ["Central neurocytoma","2","Unknown","Intraventricular"],
    ["Cerebellar liponeurocytoma","3","Unknown","Cerebellum"],
  ];
  const styled=rows.map((r,ri)=>r.map(cn=>{const b=typeof cn==="string"?{text:cn}:cn;
    return Object.assign({options:{}},b,{options:Object.assign({fontFace:SANS,fontSize:15,color:ri===0?WHITE:INK,
      fill:{color:ri===0?NAVY:(ri%2?WHITE:CARD)},valign:"middle",margin:[4,6,4,6]},(typeof cn==="object"?cn.options:{}))});}));
  s.addTable(styled,{x:0.6,y:3.6,w:12.13,colW:[3.1,0.6,5.63,2.8],border:{type:"solid",color:BORDER,pt:0.5},rowH:0.52});
  s.addText("Panel: 10x Xenium Human Brain panel (266 genes).",{x:0.6,y:6.5,w:12,h:0.4,fontFace:SANS,fontSize:13,italic:true,color:MUTED,margin:0});
  foot(s);
})();

// ============================================================ 4. FOUR HISTOLOGIES (H&E thumbnails)
(()=>{const s=pres.addSlide(); head(s,null,"Four tumor histologies");
  const cols=[["thumb_PA.jpg","Pilocytic astrocytoma","WHO 1 · biphasic · BRAF"],
    ["thumb_GG.jpg","Ganglioglioma","WHO 1 · neurons + glia · BRAF-V600E"],
    ["thumb_CN.jpg","Central neurocytoma","WHO 2 · intraventricular"],
    ["thumb_CLN.jpg","Cerebellar liponeurocytoma","WHO 2 · lipidized"]];
  const tw=2.86, gap=0.23, x0=0.6, iy=2.05, ih=2.86;
  cols.forEach((c,i)=>{const x=x0+i*(tw+gap);
    img(s,c[0],x,iy,tw,ih);
    s.addText(c[1],{x:x,y:iy+ih+0.12,w:tw,h:0.7,fontFace:SANS,fontSize:15,bold:true,color:INK,align:"center",margin:0,valign:"top",lineSpacingMultiple:1.0});
    s.addText(c[2],{x:x,y:iy+ih+0.82,w:tw,h:0.5,fontFace:SANS,fontSize:12.5,color:MUTED,align:"center",margin:0,valign:"top"});});
  foot(s);
})();

// ============================================================ 5. GOALS
(()=>{const s=pres.addSlide(); head(s,null,"Questions");
  const items=[["01","What cell types and spatial niches build each tumor?"],
    ["02","Do those spatial niches match what a pathologist sees?"],
    ["03","In PA, how does MAPK activation shape glial–immune interaction?"]];
  let y=1.85;
  items.forEach(it=>{
    s.addText(it[0],{x:0.7,y,w:1.4,h:1.3,fontFace:SERIF,fontSize:46,bold:true,color:ACCENT,margin:0,valign:"top"});
    s.addText(it[1],{x:2.25,y:y+0.08,w:10.4,h:1.4,fontFace:SANS,fontSize:26,bold:true,color:INK,margin:0,valign:"top",lineSpacingMultiple:1.02});
    y+=1.62;});
  foot(s);
})();

// ============================================================ 6. PIPELINE (flow + cartoons)
(()=>{const s=pres.addSlide(); head(s,null,"One resegmentation, feeding three analyses");
  // ProSeg foundation (left) — cartoon centered in the rectangle
  s.addText("Xenium segmentation",{x:0.5,y:1.72,w:4.5,h:0.38,fontFace:SANS,fontSize:14,color:MUTED,align:"center",margin:0});
  line2(s,2.75,2.12,2.75,2.46,{color:STEEL,width:2,endArrowType:"triangle"});
  s.addShape("roundRect",{x:0.5,y:2.52,w:4.5,h:3.95,rectRadius:0.06,fill:{color:"FBF2F3"},line:{color:ACCENT,width:1.25}});
  s.addText("ProSeg — resegment cells",{x:0.65,y:2.66,w:4.2,h:0.4,fontFace:SANS,fontSize:16,bold:true,color:ACCENT,align:"center",margin:0});
  cProseg(s,0.8,3.25,3.9,2.4);
  s.addText("default segmentation was poor → shape cells to their transcripts",
    {x:0.65,y:5.88,w:4.2,h:0.55,fontFace:SANS,fontSize:13,color:INK,align:"center",margin:0,lineSpacingMultiple:1.05});
  // three downstream modules (right) — bigger cartoons, wrapped text
  const mods=[["BANKSY",cBanksy,"find spatial niches — recurrent local neighborhoods"],
    ["ENVI",cEnvi,"impute the whole transcriptome from paired single-cell data"],
    ["STalign",cStalign,"overlay pathologist annotations onto the data"]];
  const mx=5.5, mw=7.3, mh=1.62, y0=1.72, gap=0.22, hub=[5.0,4.49];
  mods.forEach((m,i)=>{const y=y0+i*(mh+gap);
    line2(s,hub[0],hub[1],mx-0.02,y+mh/2,{color:STEEL,width:1.5,endArrowType:"triangle"});
    s.addShape("roundRect",{x:mx,y,w:mw,h:mh,rectRadius:0.06,fill:{color:CARD},line:{color:BORDER,width:1}});
    m[1](s,mx+0.2,y+0.14,3.9,mh-0.28);
    s.addText([{text:m[0]+"\n",options:{fontSize:16,bold:true,color:STEEL,fontFace:SANS,breakLine:true}},
      {text:m[2],options:{fontSize:13,color:INK,fontFace:SANS}}],
      {x:mx+4.25,y,w:mw-4.4,h:mh,valign:"middle",margin:0,lineSpacingMultiple:1.05});});
  foot(s);
})();

// ============================================================ CELL TYPING
(()=>{const s=pres.addSlide(); head(s,"cell typing","Cell types across the cohort");
  imgFit(s,"fig_celltype_comp.png",2496,1036,0.6,1.85,12.13,4.7);
  foot(s);
})();

// ============================================================ 9. SPATIAL MOTIFS
(()=>{const s=pres.addSlide(); head(s,"spatial motifs","Recurrent spatial niches");
  imgFit(s,"fig_niche_comp.png",2496,1075,0.6,1.85,12.13,4.7);
  foot(s);
})();

// ============================================================ 10. ATLAS — cell types by histology
(()=>{const s=pres.addSlide(); head(s,"spatial atlas","Each histology has its own cell-type makeup");
  imgFit(s,"fig_celltype_forest.png",2496,998,0.6,1.9,12.13,4.6);
  foot(s);
})();

// ============================================================ 11. ATLAS — niches by histology
(()=>{const s=pres.addSlide(); head(s,"spatial atlas","…and its own spatial niches");
  imgFit(s,"fig_niche_forest.png",2496,1248,0.6,1.9,12.13,4.65);
  foot(s);
})();

// ============================================================ 12. ATLAS — niches are ecosystems
(()=>{const s=pres.addSlide(); head(s,"spatial atlas","Signature niches are ecosystems of signature cell types");
  imgFit(s,"fig_ecosystems.png",2112,960,0.9,1.95,11.5,4.55);
  foot(s);
})();

// ============================================================ 13. DIVIDER — PA deep dive
(()=>{const s=pres.addSlide(); s.background={color:NAVY};
  motifDots(s,30,9.4,13.0,0.4,7.1);
  s.addText("Pilocytic astrocytoma:\nthe MAPK driver organizes the microenvironment",
    {x:0.78,y:2.7,w:11.4,h:2.0,fontFace:SERIF,fontSize:34,bold:true,color:WHITE,margin:0,lineSpacingMultiple:1.04});
})();

// ============================================================ 14. ENVI integration
(()=>{const s=pres.addSlide(); head(s,null,"Imputing a whole transcriptome in situ with ENVI");
  bullets(s,[
    "Xenium measures ~266 genes; ENVI co-embeds matched single-nucleus whole-transcriptome PA data.",
    "It imputes genome-wide expression onto each spatial cell, so we can score pathway activity (e.g. MAPK) in place.",
    "Four PA patients, four driver alterations.",
  ],0.6,2.0,6.3,4.2,22);
  cStalign; // no-op guard
  figPh(s,7.2,1.95,5.5,4.4,"ENVI co-embedding schematic","sc_integration/ENVI");
  foot(s);
})();

// ============================================================ 15. MAPK in glia
(()=>{const s=pres.addSlide(); head(s,null,"MAPK activity is concentrated in the tumor glial lineage");
  bullets(s,[
    "Imputed MAPK activity concentrates in OPCs and astrocytes across all four drivers.",
    "Myeloid, lymphoid, and stromal cells are consistently MAPK-low.",
  ],0.6,2.0,6.1,3.0,22);
  figPh(s,7.0,1.95,5.7,4.4,"MAPK activity: patient × cell type (heatmap)","2_MAPK_adj.qmd · mapk_patient_cell_type_significance.tsv");
  foot(s);
})();

// ============================================================ 16. MAPK -> myeloid
(()=>{const s=pres.addSlide(); head(s,null,"MAPK-high glia sit among more myeloid cells");
  figPh(s,7.0,1.95,5.7,4.4,"Neighbor composition vs. MAPK (volcano)","2_MAPK_adj.qmd · mapk_adj_results.tsv");
  bullets(s,[
    "As a glial cell's MAPK score rises, its neighborhood gains more myeloid cells (Astro→Myeloid +0.030, z≈19).",
    "Modest but highly consistent across cells and patients.",
    "Echoes the PA axis where microglia/myeloid recruitment enables bypass of BRAF-fusion senescence.",
  ],0.6,2.0,6.1,4.3,20);
  foot(s);
})();

// ============================================================ 17. Myeloid reprogramming
(()=>{const s=pres.addSlide(); head(s,null,"Proximity to MAPK-high glia reprograms myeloid programs");
  figPh(s,7.0,1.95,5.7,4.4,"Glial-exposure effect on myeloid signatures (heatmap)","2_MAPK_adj.qmd · myeloid_glial_exposure_part[A/B]_results.tsv");
  bullets(s,[
    "Myeloid cells near MAPK-high glia show a shifted immunoactivation ↔ immunosuppression balance.",
    "Signatures span Hallmark IFN / inflammatory, TIM3± microglia, complement / scavenger, FOSL2 regulon.",
    "Consistent within patients; small effect (R²≈1%). Emerging.",
  ],0.6,2.0,6.1,4.3,20);
  foot(s);
})();

// ============================================================ 18. Pathology (STalign)
(()=>{const s=pres.addSlide(); head(s,null,"Spatial niches are molecular, not morphological");
  figPh(s,0.6,1.95,6.0,4.4,"STalign registration + niche × region (forest)","annotations/STalign · linear_models_niches.qmd");
  bullets(s,[
    "Pathologist H&E regions registered to Xenium via STalign.",
    "Only the vascular Blood region robustly marks niches (PA niches 9 & 12).",
    "Diagnostic architecture regions do not — so niches add information beyond the H&E.",
  ],6.95,2.05,5.75,4.2,20);
  foot(s);
})();

// ============================================================ 19. Synthesis
(()=>{const s=pres.addSlide(); head(s,null,"The driver shapes the spatial immune architecture");
  const boxes=[["MAPK driver","activation"],["Glial lineage","OPC / astrocyte state"],["Immune organization","myeloid recruitment + reprogramming"]];
  const bw=3.7,bh=1.5,y=1.95,gap=0.7;
  boxes.forEach((b,i)=>{const x=0.6+i*(bw+gap);
    s.addShape("roundRect",{x,y,w:bw,h:bh,rectRadius:0.06,fill:{color:i===0?ACCENT:CARD},line:{color:i===0?ACCENT:BORDER,width:1}});
    s.addText([{text:b[0]+"\n",options:{fontSize:17,bold:true,color:i===0?WHITE:STEEL,fontFace:SANS,breakLine:true}},
      {text:b[1],options:{fontSize:13,color:i===0?"F3D9DD":INK,fontFace:SANS}}],
      {x:x+0.2,y,w:bw-0.4,h:bh,valign:"middle",align:"center",margin:0,lineSpacingMultiple:1.05});
    if(i<2)line2(s,x+bw+0.08,y+bh/2,x+bw+gap-0.08,y+bh/2,{color:MUTED,width:2,endArrowType:"triangle"});});
  bullets(s,[
    "Every histology's dominant spatial signature matches its expected cell-of-origin.",
    "In PA, the driver is glial-restricted and organizes the myeloid compartment.",
    "Molecular niches carry information beyond histopathology.",
  ],0.7,4.0,12.0,2.7,22);
  foot(s);
})();

// ============================================================ 20. Conclusions & limitations
(()=>{const s=pres.addSlide(); head(s,null,"Where this leaves us");
  s.addShape("roundRect",{x:0.6,y:2.0,w:5.95,h:4.5,rectRadius:0.06,fill:{color:CARD},line:{color:BORDER,width:1}});
  s.addText("What we show",{x:0.9,y:2.2,w:5.4,h:0.4,fontFace:SANS,fontSize:17,bold:true,color:STEEL,margin:0});
  bullets(s,["An integrated in-situ pipeline across four rare tumors.",
    "Histology-specific cellular & niche architecture matching known biology.",
    "In PA, MAPK-driven glial–immune crosstalk, seen spatially."],0.9,2.75,5.4,3.6,17);
  s.addShape("roundRect",{x:6.78,y:2.0,w:5.95,h:4.5,rectRadius:0.06,fill:{color:"FBF2F3"},line:{color:ACCENT,width:1}});
  s.addText("Limitations & next steps",{x:7.08,y:2.2,w:5.4,h:0.4,fontFace:SANS,fontSize:17,bold:true,color:ACCENT,margin:0});
  bullets(s,["n = 11, confounded and underpowered — read direction, not p-values.",
    "Sample-private niches unconfirmed (inferCNV).",
    "Niche ↔ pathology correspondence weak beyond vasculature.",
    "ENVI / MAPK deep dive is PA-only so far.",
    "Next: inferCNV, cohort expansion, extend to GG / CN / CLN."],7.08,2.75,5.4,3.6,16);
  foot(s);
})();

// ============================================================ 21. Acknowledgements
(()=>{const s=pres.addSlide(); s.background={color:NAVY};
  motifDots(s,40,9.0,13.0,0.4,7.1);
  s.addText("Acknowledgements",{x:0.8,y:1.6,w:10,h:0.9,fontFace:SERIF,fontSize:36,bold:true,color:WHITE,margin:0});
  const g=[["Collaborators","[names / labs]"],["Neuropathology","[pathologist(s)]"],["Single-nucleus data","[contributors]"],["Funding","[grants]"]];
  g.forEach((it,i)=>{const y=2.95+i*0.92;
    s.addShape("ellipse",{x:0.85,y:y+0.08,w:0.15,h:0.15,fill:{color:ACCENT},line:{color:ACCENT,width:0.5}});
    s.addText([{text:it[0]+"   ",options:{bold:true,color:WHITE,fontSize:17,fontFace:SANS}},
      {text:it[1],options:{color:ICE,fontSize:15,fontFace:SANS}}],{x:1.2,y,w:10.5,h:0.5,valign:"middle",margin:0});});
  s.addText("Thank you",{x:0.8,y:6.5,w:6,h:0.6,fontFace:SERIF,fontSize:22,italic:true,color:ACCENT,margin:0});
})();

pres.writeFile({fileName:path.join(__dirname,"plgg_xenium_deck.pptx")}).then(f=>console.log("WROTE",f));
