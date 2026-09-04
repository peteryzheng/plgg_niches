// mockup_spine.js — standalone mockup of the 3 conceptually-new "spine" slides for the
// deck restructure: (A) reworded hook, (B) the triangle roadmap, (C) hub-and-shadows synthesis.
// It is NOT the real deck; it exists so we can iterate on these three slides before rewriting
// build_deck.js. Palette + shared helpers + the four pipeline cartoons are copied verbatim from
// build_deck.js (the parameterized c* versions that scale to arbitrary boxes); only the triangle
// and synthesis *assembly* below is new. Reads assets from ../assets, writes ../mockup_spine.pptx.
//   run:  cd presentation && NODE_PATH="$(npm root -g)" node src/mockup_spine.js
const pptxgen = require("pptxgenjs");
const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.33 x 7.5
pres.title = "Spine mockup — Spatial Architecture of Pediatric Brain Tumors";

const NAVY="1F2A44", INK="1F2A37", ACCENT="A4243B", STEEL="2E5A87", MUTED="6B7280",
      CARD="F3F5F8", BORDER="D5DBE3", WHITE="FFFFFF", ICE="C7D2E4",
      GREEN="2C7A4B", AMBER="C98A2B", NUC="3A3F4A";
const SERIF="Cambria", SANS="Calibri";
const PAL=[ACCENT,STEEL,GREEN,AMBER];
const GENECOLS=["A4243B","2E5A87","2C7A4B","C98A2B","6D5199","1C7293","B5651D","4C8577","8A2846"];
const path=require("path");
const A=path.join(__dirname,"..","assets")+path.sep;   // assets/ is one level up from src/
const R=(a,b)=>a+Math.random()*(b-a);

// ---------- shared helpers (copied from build_deck.js) ----------
const dot=(s,x,y,d,c,tr=0)=>s.addShape("ellipse",{x:x-d/2,y:y-d/2,w:d,h:d,fill:{color:c,transparency:tr},line:{color:c,width:0.4}});
function line2(s,x1,y1,x2,y2,o={}){
  const x=Math.min(x1,x2),y=Math.min(y1,y2),w=Math.abs(x2-x1),h=Math.abs(y2-y1);
  s.addShape("line",{x,y,w,h,flipV:((x2-x1)*(y2-y1))<0,line:Object.assign({color:MUTED,width:1},o)});
}
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
function img(s,file,x,y,w,h){
  s.addShape("roundRect",{x:x-0.03,y:y-0.03,w:w+0.06,h:h+0.06,rectRadius:0.03,fill:{color:WHITE},line:{color:BORDER,width:1}});
  s.addImage({path:A+file,x,y,w,h});
}

// ---------- pipeline cartoons (copied verbatim from build_deck.js) ----------
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

// ============================================================ A. HOOK (reworded)
// Plant the driver + the biphasic morphology, then pose the two spine questions.
(()=>{const s=pres.addSlide(); head(s,null,"One MAPK driver — does it organize the tissue?");
  bullets(s,[
    "Pilocytic astrocytoma is driven by a single MAPK alteration (BRAF fusion / V600E / FGFR1).",
    "It is biphasic — dense compact-fibrillary tissue abuts loose myxoid tissue (H&E, right).",
    "Q1  Does that one driver organize the tissue's spatial architecture?",
    "Q2  And does that architecture match what the pathologist sees?",
  ],0.6,1.9,6.3,4.5,22);
  img(s,"pa_biphasic_slide_annot.jpg",7.75,1.75,4.6,4.6);
  s.addText("Pilocytic astrocytoma H&E — green = compact fibrillary, blue = loose myxoid (pathologist-annotated).",
    {x:7.75,y:6.5,w:4.6,h:0.6,fontFace:SANS,fontSize:12,italic:true,color:MUTED,align:"center",margin:0,valign:"top",lineSpacingMultiple:1.05});
  foot(s);
})();

// ============================================================ B. THE TRIANGLE (roadmap)
// Reuse the four cartoons: ProSeg = foundation slab; BANKSY / STalign / ENVI = the 3 vertices
// (composition / morphology / state); each pairing is a scientific question with a status.
(()=>{const s=pres.addSlide(); head(s,"approach & roadmap","One tissue, three descriptions of it");
  s.addText("ProSeg gives one set of cells; BANKSY, STalign, and ENVI each describe them — and every pairing raises a question (Q1–Q3).",
    {x:0.6,y:1.5,w:12.1,h:0.34,fontFace:SANS,fontSize:14.5,italic:true,color:MUTED,margin:0});

  // vertex card = cartoon + "MODULE · layer" label; layerColor distinguishes the three layers
  function vertex(x,y,w,h,name,layer,layerColor,draw){
    s.addShape("roundRect",{x,y,w,h,rectRadius:0.06,fill:{color:CARD},line:{color:BORDER,width:1}});
    s.addText([{text:name,options:{bold:true,color:layerColor,fontSize:15,fontFace:SANS}},
      {text:"   ·  "+layer,options:{color:MUTED,fontSize:13,fontFace:SANS}}],
      {x:x+0.16,y:y+0.08,w:w-0.32,h:0.3,margin:0,valign:"middle"});
    draw(s,x+0.2,y+0.44,w-0.4,h-0.58);
  }
  // vertex geometry (top center = composition; bottom-left = state; bottom-right = morphology).
  // Enlarged + pushed down to fill the full slide now that the ProSeg slab is gone.
  const BX=4.82,BY=2.25,BW=3.7,BH=1.75;   // BANKSY (composition) — apex
  const EX=0.7,EY=5.0,EW=3.7,EH=1.75;     // ENVI  (cell state)   — bottom-left
  const SX=8.93,SY=5.0,SW=3.7,SH=1.75;    // STalign (morphology) — bottom-right
  // edges drawn BEHIND the cards: solid = answered, dashed = open/absent; color+status per question
  line2(s,BX+0.55,BY+BH-0.15,EX+EW-0.65,EY+0.15,{color:GREEN,width:2.75});                    // Q1 BANKSY–ENVI (answered)
  line2(s,BX+BW-0.55,BY+BH-0.15,SX+0.65,SY+0.15,{color:ACCENT,width:2.5,dashType:"dash"});     // Q2 BANKSY–STalign (no)
  line2(s,EX+EW,EY+EH/2,SX,SY+SH/2,{color:AMBER,width:2.25,dashType:"dash"});                  // Q3 ENVI–STalign (open)
  // vertices on top
  vertex(BX,BY,BW,BH,"BANKSY","composition",STEEL,cBanksy);
  vertex(EX,EY,EW,EH,"ENVI","cell STATE",ACCENT,cEnvi);
  vertex(SX,SY,SW,SH,"STalign","morphology",AMBER,cStalign);
  // edge labels (question + status), placed in the open space beside each edge
  const elabel=(x,y,w,t,st,col)=>s.addText([{text:t+"\n",options:{bold:true,fontSize:14,color:col,fontFace:SANS,breakLine:true}},
    {text:st,options:{italic:true,fontSize:12,color:MUTED,fontFace:SANS}}],
    {x,y,w,h:0.75,align:"center",valign:"middle",margin:0,lineSpacingMultiple:1.03});
  elabel(0.85,3.9,3.3,"Q1 · Does niche reflect cell state?","answered",GREEN);
  elabel(9.15,3.9,3.3,"Q2 · Does niche correlate with morphology?","no",ACCENT);
  elabel(4.82,5.05,3.7,"Q3 · Does cell state correlate with morphology?","open",AMBER);
  foot(s);
})();

// ============================================================ C. SYNTHESIS (hub-and-shadows)
// State is the generative hub; composition + morphology are two SEPARABLE downstream shadows
// reading different axes of state (which is why Q2 is null). Full-circle to biphasic PA.
(()=>{const s=pres.addSlide(); head(s,null,"Synthesis — state is the generative hub");
  // driver -> state hub
  const DX=0.7,DY=2.75,DW=2.15,DH=1.2;
  s.addShape("roundRect",{x:DX,y:DY,w:DW,h:DH,rectRadius:0.06,fill:{color:ACCENT},line:{color:ACCENT,width:1}});
  s.addText([{text:"MAPK\n",options:{bold:true,fontSize:18,color:WHITE,fontFace:SANS,breakLine:true}},
    {text:"driver",options:{fontSize:14,color:"F3D9DD",fontFace:SANS}}],
    {x:DX,y:DY,w:DW,h:DH,align:"center",valign:"middle",margin:0,lineSpacingMultiple:1.02});
  // the hub
  const HX=4.35,HY=2.45,HW=3.45,HH=1.8;
  // driver -> hub arrow spans the full gap so the two boxes visibly connect
  line2(s,DX+DW,DY+DH/2,HX,DY+DH/2,{color:MUTED,width:2.5,endArrowType:"triangle"});
  s.addShape("roundRect",{x:HX,y:HY,w:HW,h:HH,rectRadius:0.08,fill:{color:"F7EEF0"},line:{color:ACCENT,width:1.75}});
  s.addText([{text:"CELL STATE\n",options:{bold:true,fontSize:20,color:ACCENT,fontFace:SERIF,breakLine:true}},
    {text:"ENVI-imputed · set by the driver",options:{fontSize:13,color:INK,fontFace:SANS}}],
    {x:HX,y:HY,w:HW,h:HH,align:"center",valign:"middle",margin:0,lineSpacingMultiple:1.1});
  // two separable downstream shadows
  const shadow=(x,y,w,h,name,sub,col)=>{s.addShape("roundRect",{x,y,w,h,rectRadius:0.06,fill:{color:CARD},line:{color:col,width:1.5}});
    s.addText([{text:name+"\n",options:{bold:true,fontSize:15,color:col,fontFace:SANS,breakLine:true}},
      {text:sub,options:{fontSize:12.5,color:INK,fontFace:SANS}}],{x:x+0.15,y,w:w-0.3,h,align:"center",valign:"middle",margin:0,lineSpacingMultiple:1.05});};
  const CX=9.0,CY=1.75,CW=3.55,CH=1.3;   // composition shadow
  const MX=9.0,MY=3.55,MW=3.55,MH=1.3;   // morphology shadow
  shadow(CX,CY,CW,CH,"COMPOSITION","BANKSY niches · reads the immune axis",STEEL);
  shadow(MX,MY,MW,MH,"MORPHOLOGY","STalign regions · reads the ECM / secretory axis",AMBER);
  line2(s,HX+HW,HY+0.45,CX,CY+CH/2,{color:STEEL,width:2,endArrowType:"triangle"});
  line2(s,HX+HW,HY+HH-0.45,MX,MY+MH/2,{color:AMBER,width:2,endArrowType:"triangle"});
  s.addText("two separable shadows",{x:HX+HW+0.05,y:3.15,w:1.9,h:0.3,fontFace:SANS,fontSize:11,italic:true,color:MUTED,align:"center",margin:0});
  // takeaways: decoupling + full-circle
  s.addShape("roundRect",{x:0.7,y:5.55,w:11.9,h:1.35,rectRadius:0.06,fill:{color:CARD},line:{color:BORDER,width:1}});
  s.addText([{text:"Composition and morphology read ",options:{fontSize:16,color:INK,fontFace:SANS}},
    {text:"different axes of state",options:{fontSize:16,bold:true,color:STEEL,fontFace:SANS}},
    {text:" — so they need not agree (why Q2 is null).\n",options:{fontSize:16,color:INK,fontFace:SANS,breakLine:true}},
    {text:"Full circle:  ",options:{fontSize:16,bold:true,color:ACCENT,fontFace:SANS}},
    {text:"biphasic PA = the state depositing glycan-rich myxoid ECM — what the pathologist sees.",options:{fontSize:16,color:INK,fontFace:SANS}}],
    {x:1.0,y:5.62,w:11.3,h:1.2,valign:"middle",margin:0,lineSpacingMultiple:1.15});
  foot(s);
})();

pres.writeFile({fileName:path.join(__dirname,"..","mockup_spine.pptx")}).then(f=>console.log("WROTE",f));
