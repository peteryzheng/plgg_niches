// Deck generator v4 — Spatial Architecture of Pediatric Brain Tumors.
// Two-act spine: Act I describes composition (the BANKSY atlas); a pivot slide then poses three
// edge-questions linking BANKSY (composition) / ENVI (state) / STalign (morphology), answered in
// Act II (Q1–Q3) and resolved in the synthesis. 20 slides; broad spatial framing (MAPK enters only
// in Q1). Figures embedded from assets/.
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
// mini "you-are-here" triangle locator for the Q opener slides: BANKSY apex, ENVI bottom-left,
// STalign bottom-right; the active edge is drawn in ACCENT, the others muted. active: "Q1"|"Q2"|"Q3".
function triLocator(s,ox,oy,sc,active){
  const B=[ox+0.85*sc,oy], E=[ox,oy+1.0*sc], Sv=[ox+1.7*sc,oy+1.0*sc];
  const edge=(p,q,on)=>line2(s,p[0],p[1],q[0],q[1],{color:on?INK:"CBD2DB",width:on?2.75:1});
  edge(B,E,active==="Q1"); edge(B,Sv,active==="Q2"); edge(E,Sv,active==="Q3");
  [[B,STEEL],[E,ACCENT],[Sv,AMBER]].forEach(v=>dot(s,v[0][0],v[0][1],0.13*sc,v[1]));
  s.addText("BANKSY",{x:B[0]-0.75*sc,y:B[1]-0.34*sc,w:1.5*sc,h:0.2,fontFace:SANS,fontSize:8,bold:true,color:MUTED,align:"center",margin:0});
  s.addText("ENVI",{x:E[0]-0.92*sc,y:E[1]-0.06,w:0.85*sc,h:0.2,fontFace:SANS,fontSize:8,bold:true,color:MUTED,align:"right",margin:0});
  s.addText("STalign",{x:Sv[0]+0.12*sc,y:Sv[1]-0.06,w:1.05*sc,h:0.2,fontFace:SANS,fontSize:8,bold:true,color:MUTED,align:"left",margin:0});
}

// ---------- pipeline cartoons (ported from preview.js) ----------
// cProseg is kept for a potential backup / methods slide even though the main deck no longer
// draws it inline (the triangle roadmap reuses BANKSY/ENVI/STalign only).
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
  // label the measured panel (small matrix) and the ENVI-imputed genome-wide matrix (big matrix)
  s.addText("266 genes",{x:ax,y:ay+1.02*sc,w:1.0*sc,h:0.3,align:"center",valign:"top",margin:0,fontFace:SANS,fontSize:8.5,italic:true,color:MUTED});
  s.addText("full transcriptome",{x:ax+2.0*sc,y:ay+1.02*sc,w:2.5*sc,h:0.3,align:"center",valign:"top",margin:0,fontFace:SANS,fontSize:8.5,italic:true,color:MUTED});
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
    {x:0.8,y:1.95,w:9.8,h:2.0,fontFace:SERIF,fontSize:46,bold:true,color:WHITE,margin:0,lineSpacingMultiple:1.02});
  s.addText([{text:"Peter Youyun Zheng",options:{bold:true,color:WHITE,fontSize:16,fontFace:SANS,breakLine:true}},
    {text:"Beroukhim Lab · Dana-Farber Cancer Institute / Broad Institute",options:{color:ICE,fontSize:14,fontFace:SANS,breakLine:true}},
    {text:"September 2026",options:{color:MUTED,fontSize:12.5,fontFace:SANS}}],
    {x:0.8,y:5.6,w:9,h:1.3,margin:0,lineSpacingMultiple:1.3});
})();

// ============================================================ 2. HOOK (broad — why study tumors spatially)
(()=>{const s=pres.addSlide(); head(s,null,"Tumors are spatially organized — on H&E, you can see it");
  bullets(s,[
    "A tumor is not a bag of cells: on H&E its structure is spatially organized — pilocytic astrocytoma, for instance, is biphasic, with dense compact-fibrillary tissue abutting loose myxoid tissue (right).",
    "Dissociation-based bulk and single-cell methods discard that spatial context.",
    "Single-cell spatial transcriptomics maps every cell in place — so we can ask what recurrent tissue organization exists across tumors, and what it reflects.",
  ],0.6,1.95,6.3,4.5,21);
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

// ============================================================ 5. BANKSY MODULE (composition lens)
// Introduces the first lens right before its results (the atlas, slides 6–10). The triangle roadmap
// now lives at the post-atlas pivot (slide 11), where ENVI + STalign are introduced alongside it.
(()=>{const s=pres.addSlide(); head(s,"method · composition","BANKSY groups cells into recurrent spatial niches");
  bullets(s,[
    "First we assign every cell a type from its expression profile.",
    "BANKSY then clusters each cell on its own profile plus its neighborhood — so cells in similar surroundings group together.",
    "Recurrent multicellular neighborhoods emerge as spatial niches (dashed) — this is how we read composition: which cell types co-occur, and where.",
  ],0.6,2.0,6.2,4.4,20);
  // BANKSY cartoon in a card on the right
  s.addShape("roundRect",{x:7.0,y:2.15,w:5.6,h:3.7,rectRadius:0.06,fill:{color:CARD},line:{color:BORDER,width:1}});
  cBanksy(s,7.35,2.75,4.9,2.3);
  s.addText("BANKSY: cell profiles + neighborhoods → spatial niches",{x:7.0,y:5.5,w:5.6,h:0.3,fontFace:SANS,fontSize:12,italic:true,color:MUTED,align:"center",margin:0});
  foot(s);
})();

// ============================================================ 6. Cell-type composition (bars)
(()=>{const s=pres.addSlide(); head(s,"cell typing","Cell types across the cohort");
  imgFit(s,"fig_celltype_comp.png",2496,1036,0.6,1.85,12.13,4.7);
  foot(s);
})();

// ============================================================ 7. Cell-type enrichment by histology (forest)
(()=>{const s=pres.addSlide(); head(s,"spatial atlas","Each histology has its own cell-type makeup");
  imgFit(s,"fig_celltype_forest.png",2496,998,0.6,1.9,12.13,4.6);
  foot(s);
})();

// ============================================================ 8. Niche composition (bars)
(()=>{const s=pres.addSlide(); head(s,"spatial motifs","Recurrent spatial niches");
  imgFit(s,"fig_niche_comp.png",2496,1075,0.6,1.85,12.13,4.7);
  foot(s);
})();

// ============================================================ 9. Niche enrichment by histology (forest)
(()=>{const s=pres.addSlide(); head(s,"spatial atlas","Each histology has its own spatial niches");
  imgFit(s,"fig_niche_forest.png",2496,1248,0.6,1.9,12.13,4.65);
  foot(s);
})();

// ============================================================ 10. Niches are ecosystems (heatmap) — Act I payoff
(()=>{const s=pres.addSlide(); head(s,"spatial atlas","Signature niches are ecosystems of signature cell types");
  imgFit(s,"fig_ecosystems.png",2112,960,0.9,1.8,11.5,4.4);
  // transition into Act II: niches are made of histology-specific cells — but are those cells also
  // high in the biologically relevant pathway? (sets up Q1 / the ENVI state lens)
  s.addText("Histology-specific niches are built from histology-specific cells — but are those cells also high in the biologically relevant pathway (MAPK in PA)?",
    {x:0.9,y:6.45,w:11.5,h:0.6,fontFace:SANS,fontSize:14,italic:true,color:MUTED,align:"center",margin:0,valign:"top",lineSpacingMultiple:1.05});
  foot(s);
})();

// ============================================================ 11. PIVOT — the triangle + three questions
// Act I described composition (BANKSY). Introduce the other two lenses (ENVI state, STalign morphology)
// and pose the three edge-questions Q1–Q3 — immediately before Act II answers them. Edges are neutral
// here; the synthesis reveals each verdict.
(()=>{const s=pres.addSlide(); head(s,"the questions","Three views of one tissue");
  s.addText("Act I described composition — histology-specific niches of histology-specific cells. Now ENVI (cell state) and STalign (morphology) describe the same cells; every pairing is a question.",
    {x:0.6,y:1.45,w:12.2,h:0.55,fontFace:SANS,fontSize:14.5,italic:true,color:MUTED,margin:0,lineSpacingMultiple:1.06});
  // vertex card = cartoon + "MODULE · layer" label; layerColor distinguishes the three layers
  function vertex(x,y,w,h,name,layer,layerColor,draw){
    s.addShape("roundRect",{x,y,w,h,rectRadius:0.06,fill:{color:CARD},line:{color:BORDER,width:1}});
    s.addText([{text:name,options:{bold:true,color:layerColor,fontSize:15,fontFace:SANS}},
      {text:"   ·  "+layer,options:{color:MUTED,fontSize:13,fontFace:SANS}}],
      {x:x+0.16,y:y+0.08,w:w-0.32,h:0.3,margin:0,valign:"middle"});
    draw(s,x+0.2,y+0.44,w-0.4,h-0.58);
  }
  // vertex geometry: top center = composition; bottom-left = state; bottom-right = morphology
  const BX=4.82,BY=2.25,BW=3.7,BH=1.75;   // BANKSY (composition) — apex
  const EX=0.7,EY=5.0,EW=3.7,EH=1.75;     // ENVI  (cell state)   — bottom-left
  const SX=8.93,SY=5.0,SW=3.7,SH=1.75;    // STalign (morphology) — bottom-right
  // edges neutral here — the questions are posed, not yet answered (verdicts revealed in the synthesis)
  line2(s,BX+0.55,BY+BH-0.15,EX+EW-0.65,EY+0.15,{color:STEEL,width:2.25});   // Q1 BANKSY–ENVI
  line2(s,BX+BW-0.55,BY+BH-0.15,SX+0.65,SY+0.15,{color:STEEL,width:2.25});   // Q2 BANKSY–STalign
  line2(s,EX+EW,EY+EH/2,SX,SY+SH/2,{color:STEEL,width:2.25});                // Q3 ENVI–STalign
  vertex(BX,BY,BW,BH,"BANKSY","composition",STEEL,cBanksy);
  vertex(EX,EY,EW,EH,"ENVI","cell state",ACCENT,cEnvi);
  vertex(SX,SY,SW,SH,"STalign","morphology",AMBER,cStalign);
  // each edge: a short correlational question + one-line expansion (this slide only — later kickers stay short)
  const q=(x,y,w,qt,d)=>s.addText([{text:qt+"\n",options:{bold:true,fontSize:14,color:STEEL,fontFace:SANS,breakLine:true}},
    {text:d,options:{fontSize:11.5,color:MUTED,fontFace:SANS}}],
    {x,y,w,h:1.4,align:"center",valign:"top",margin:0,lineSpacingMultiple:1.06});
  q(0.3,3.05,3.85,"Q1 · Does niche relate to cell state?","Are a niche's cells also high in the biologically relevant pathway?");
  q(9.0,3.05,4.0,"Q2 · Does niche match morphology?","Do BANKSY niches line up with the pathologist's annotations?");
  q(4.4,5.95,4.1,"Q3 · Does cell state relate to morphology?","Do cells in different morphology regions carry different states?");
  foot(s);
})();

// ============================================================ 12. Q1 — niche → state (ENVI beat + result)
(()=>{const s=pres.addSlide(); head(s,"Q1 · does niche shape cell state?","Niches organize glial state");
  triLocator(s,11.3,0.5,0.68,"Q1");   // you-are-here: edge 1 highlighted
  // niche three-group MAPK figure enlarged on the left; two succinct points on the right
  imgFit(s,"fig_mapk_by_niche.png",1920,1152,0.3,1.8,8.2,5.2);
  bullets(s,[
    "ENVI imputes genome-wide expression per Xenium cell — so MAPK activity is scored in place (4 PA patients).",
    "Within glia (OPC & astrocytes), MAPK runs higher in niches 3 and 5 than elsewhere: the niche a glial cell occupies tracks its state.",
  ],8.75,2.85,4.15,3.8,17);
  foot(s);
})();

// ============================================================ 13. Q1 — state → immune (recruitment, patient-adjusted)
// Honesty: the model adjusts for patient only (neighbor_prop ~ MAPK + sample_id), NOT niche. Keep the
// circularity caveat — MAPK-high glia sit in niches 3/5 (slide 12), so a niche-adjustment is still owed.
(()=>{const s=pres.addSlide(); head(s,"Q1 · does glial state reshape its neighborhood?","MAPK-high glia gain myeloid, lose endothelial neighbors");
  // per-patient neighbor-composition figure enlarged on the left; three succinct points on the right
  imgFit(s,"fig_mapk_neighbor_by_patient.png",2496,1440,0.3,1.8,8.2,5.2);
  bullets(s,[
    "As glial MAPK rises, the neighborhood gains myeloid cells — enriched across patients (the recruitment signal).",
    "Endothelial neighbors move the opposite way — depleted as MAPK rises — a shift toward immune, away from vascular cells.",
    "Model adjusts for patient, not niche; a planned niche-adjustment will confirm this is state-driven, not just niche composition.",
  ],8.75,2.4,4.15,4.7,16);
  foot(s);
})();

// ============================================================ 14. Q1 — state → immune (reprogramming: immunosuppressive-TAM axis)
// Focused 4-program figure from PA_myeloid_mapk_neighborhood.qmd: Complement + activated microglia up,
// IFN-γ + inflammatory response down = immunosuppressive-TAM shift, same direction in all 4 PA samples.
(()=>{const s=pres.addSlide(); head(s,"Q1 · does glial state reprogram myeloid cells?","Myeloid cells skew immunosuppressive");
  s.addText("↑ Complement & activated microglia (MG_Act) · ↓ IFN-γ & inflammatory response — an immunosuppressive-TAM shift near MAPK-high glia.",
    {x:0.6,y:1.5,w:12.2,h:0.4,fontFace:SANS,fontSize:15,italic:true,color:MUTED,margin:0,lineSpacingMultiple:1.05});
  imgFit(s,"fig_myeloid_tam_axis.png",2520,1200,0.4,1.9,12.5,4.95);   // as large as fits below the tagline
  foot(s);
})();

// ============================================================ 15. Q2 — setup: STalign registration + the question
// Ease-in for the null: introduce STalign, show the registration checks out (endothelial cells on
// H&E vessels), then pose the compact-vs-myxoid question. The null result lands on the next slide.
(()=>{const s=pres.addSlide(); head(s,"Q2 · does niche match morphology?","Do niches track the pathologist's annotations?");
  triLocator(s,11.3,0.5,0.68,"Q2");   // you-are-here: edge 2 highlighted
  bullets(s,[
    "STalign registers the H&E + the pathologist's annotations on the H&E onto the spatial data, cell by cell.",
    "The alignment checks out: endothelial cells land squarely on H&E blood vessels (right).",
    "So we can ask directly — do BANKSY niches line up with compact-fibrillary vs. loose-myxoid regions?",
  ],0.6,2.15,6.6,4.4,19);
  imgFit(s,"fig_endothelial_vessel.png",1038,1121,7.6,2.2,4.0,4.3);
  s.addText("Registration control · PA258482 H&E",{x:7.6,y:6.55,w:4.0,h:0.3,fontFace:SANS,fontSize:11,italic:true,color:MUTED,align:"center",margin:0});
  foot(s);
})();

// ============================================================ 16. Q2 — result: the null
(()=>{const s=pres.addSlide(); head(s,"Q2 · does niche match morphology?","Niches don't track the pathologist's annotations");
  imgFit(s,"fig_niche_pa_morphology.png",1920,844,0.6,2.15,12.13,3.9);
  s.addText("Compact-fibrillary vs. loose-myxoid regions show no enriched or depleted niches (all n.s. at q < 0.25).",
    {x:0.6,y:6.35,w:12.13,h:0.5,fontFace:SANS,fontSize:15,italic:true,color:MUTED,align:"center",margin:0,valign:"top"});
  foot(s);
})();

// ============================================================ 17. Q3 — setup: do cell states differ by region?
// Bridge from Q2's null: composition didn't separate the two regions, so ask state directly. Reuse the
// hook H&E crop; a central "Do cell states differ?" box points (green/blue arrows) at the compact and
// loose-myxoid regions, which are labeled on either side.
(()=>{const s=pres.addSlide(); head(s,"Q3 · does cell state relate to morphology?","Composition didn't explain it — does cell state?");
  triLocator(s,11.3,0.5,0.68,"Q3");   // you-are-here: edge 3 highlighted
  img(s,"pa_biphasic_slide_annot.jpg",4.72,3.0,3.9,3.9);   // hook crop, moved down to clear the question box
  // central question box (below the title), connected by two arrows to the compact (left) and loose (right) regions
  s.addShape("roundRect",{x:5.0,y:1.95,w:3.5,h:0.62,rectRadius:0.06,fill:{color:CARD},line:{color:INK,width:1.25}});
  s.addText("Do cell states differ?",{x:5.0,y:1.95,w:3.5,h:0.62,align:"center",valign:"middle",margin:0,fontFace:SANS,fontSize:16,bold:true,color:INK});
  s.addShape("line",{x:5.25,y:2.6,w:0.45,h:0.95,flipV:true,line:{color:GREEN,width:3,beginArrowType:"triangle"}});  // → compact (left)
  s.addShape("line",{x:7.85,y:2.6,w:0.45,h:0.95,line:{color:STEEL,width:3,endArrowType:"triangle"}});                // → loose (right)
  // region labels beside their sides (no boundary subtext)
  s.addText("COMPACT",{x:0.5,y:4.6,w:3.6,h:0.5,align:"center",valign:"middle",margin:0,fontFace:SANS,fontSize:19,bold:true,color:GREEN});
  s.addText("LOOSE MYXOID",{x:9.2,y:4.6,w:3.6,h:0.5,align:"center",valign:"middle",margin:0,fontFace:SANS,fontSize:19,bold:true,color:STEEL});
  foot(s);
})();

// ============================================================ 18. Q3 — result: glycan-ECM state
(()=>{const s=pres.addSlide(); head(s,"Q3 · does state shape morphology?","Myxoid morphology carries a glycan-ECM state");
  // ~1.5:1 heatmap enlarged on the left; one concise point in a narrow right column
  imgFit(s,"fig_glycan_ecm.png",3252,2171,0.3,1.75,7.9,5.3);
  bullets(s,[
    "Across OPC, astrocytes, and myeloid, myxoid regions turn up GAG / proteoglycan synthesis and down its breakdown — a coordinated matrix-accumulation state.",
  ],8.45,3.1,4.35,3.0,17);
  foot(s);
})();

// ============================================================ 19. SYNTHESIS (triangle answered)
// Echoes the pivot triangle, now with each edge's verdict (↔ = associated, ≠ = null). No takeaway
// strip, no MAPK-driver box — the framing stays broad and correlational.
(()=>{const s=pres.addSlide(); head(s,null,"Synthesis — cell state is the common thread");
  function vertex(x,y,w,h,name,layer,layerColor,draw){
    s.addShape("roundRect",{x,y,w,h,rectRadius:0.06,fill:{color:CARD},line:{color:BORDER,width:1}});
    s.addText([{text:name,options:{bold:true,color:layerColor,fontSize:15,fontFace:SANS}},
      {text:"   ·  "+layer,options:{color:MUTED,fontSize:13,fontFace:SANS}}],
      {x:x+0.16,y:y+0.08,w:w-0.32,h:0.3,margin:0,valign:"middle"});
    draw(s,x+0.2,y+0.44,w-0.4,h-0.58);
  }
  // triangle enlarged to fill the slide (no takeaway strip); edges carry each verdict
  const BX=4.82,BY=2.05,BW=3.7,BH=1.75;   // BANKSY apex
  const EX=0.7,EY=4.85,EW=3.7,EH=1.75;    // ENVI bottom-left
  const SX=8.93,SY=4.85,SW=3.7,SH=1.75;   // STalign bottom-right
  line2(s,BX+0.55,BY+BH-0.15,EX+EW-0.65,EY+0.15,{color:GREEN,width:2.75});                 // Q1 holds
  line2(s,BX+BW-0.55,BY+BH-0.15,SX+0.65,SY+0.15,{color:ACCENT,width:2,dashType:"dash"});    // Q2 null
  line2(s,EX+EW,EY+EH/2,SX,SY+SH/2,{color:GREEN,width:2.75});                              // Q3 holds
  vertex(BX,BY,BW,BH,"BANKSY","composition",STEEL,cBanksy);
  vertex(EX,EY,EW,EH,"ENVI","cell state",ACCENT,cEnvi);
  vertex(SX,SY,SW,SH,"STalign","morphology",AMBER,cStalign);
  // verdict labels — ↔ = associated with (correlational, no causal direction); ≠ = no association
  const v=(x,y,w,t,col)=>s.addText(t,{x,y,w,h:0.7,align:"center",valign:"middle",margin:0,fontFace:SANS,fontSize:15,bold:true,color:col,lineSpacingMultiple:1.04});
  v(0.35,3.5,4.0,"Q1 · niche ↔ state ↔ immune infiltration",GREEN);
  v(8.95,3.5,4.05,"Q2 · niche ≠ morphology",ACCENT);
  v(4.4,4.15,4.1,"Q3 · state ↔ morphology",GREEN);
  foot(s);
})();

// ============================================================ 20. Acknowledgements (headshot bubbles)
(()=>{const s=pres.addSlide(); s.background={color:NAVY};
  motifDots(s,40,9.0,13.0,0.4,7.1);
  s.addText("Acknowledgements",{x:0.8,y:0.55,w:11.7,h:0.9,fontFace:SERIF,fontSize:36,bold:true,color:WHITE,margin:0});
  // circular headshot (rounding=ellipse mask) + themed ring behind + name + affiliation
  function bubble(file,cx,cy,d,ring,name,sub){
    const x=cx-d/2,y=cy-d/2;
    s.addShape("ellipse",{x:x-0.05,y:y-0.05,w:d+0.10,h:d+0.10,fill:{color:ring}});
    s.addImage({path:A+file,x,y,w:d,h:d,rounding:true});
    s.addText(name,{x:cx-1.5,y:y+d+0.10,w:3.0,h:0.30,align:"center",fontFace:SANS,fontSize:15,bold:true,color:WHITE,margin:0});
    s.addText(sub ,{x:cx-1.5,y:y+d+0.42,w:3.0,h:0.34,align:"center",fontFace:SANS,fontSize:11,color:ICE,margin:0,valign:"top",lineSpacingMultiple:1.02});
  }
  // 2x3 grid — row 1 = DKFZ / Heidelberg, row 2 = US institutions
  const PPL=[
    ["ack_savran.jpg","Dilan Savran","DKFZ / Univ. Hospital Heidelberg"],
    ["ack_sievers.jpg","Philipp Sievers","DKFZ / Univ. Hospital Heidelberg"],
    ["ack_sahm.jpg","Felix Sahm","DKFZ / Univ. Hospital Heidelberg"],
    ["ack_collins.jpg","Jared Collins","MIT"],
    ["ack_bandopadhayay.jpg","Pratiti Bandopadhayay","DFCI / Broad"],
    ["ack_beroukhim.jpg","Rameen Beroukhim","DFCI / Broad"]];
  const colX=[2.95,6.67,10.39], rowY=[2.5,4.85], D=1.4;
  PPL.forEach((p,i)=>bubble(p[0],colX[i%3],rowY[Math.floor(i/3)],D,PAL[i%4],p[1],p[2]));
  // additional non-bubble thanks
  s.addText("WITH THANKS ALSO TO",{x:0,y:6.55,w:13.33,h:0.28,align:"center",fontFace:SANS,fontSize:11,bold:true,color:ACCENT,charSpacing:2,margin:0});
  s.addText("Natalie Schoebe    ·    Sabrina Sprengart",{x:0,y:6.83,w:13.33,h:0.4,align:"center",fontFace:SANS,fontSize:14,color:WHITE,margin:0});
})();

pres.writeFile({fileName:path.join(__dirname,"plgg_xenium_deck.pptx")}).then(f=>console.log("WROTE",f));
