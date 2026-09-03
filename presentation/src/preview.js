// PREVIEW v2 — cartoons redone per feedback. Not the real deck.
const pptxgen = require("pptxgenjs");
const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
const NAVY="1F2A44", INK="1F2A37", ACCENT="A4243B", STEEL="2E5A87", MUTED="6B7280",
      CARD="F3F5F8", BORDER="D5DBE3", WHITE="FFFFFF", GREEN="2C7A4B", AMBER="C98A2B", NUC="3A3F4A";
const SERIF="Cambria", SANS="Calibri";
const GENECOLS=["A4243B","2E5A87","2C7A4B","C98A2B","6D5199","1C7293","B5651D","4C8577","8A2846"];
const PAL=[ACCENT,STEEL,GREEN,AMBER];
const R=(a,b)=>a+Math.random()*(b-a);
const dot=(s,x,y,d,c,tr=0)=>s.addShape("ellipse",{x:x-d/2,y:y-d/2,w:d,h:d,fill:{color:c,transparency:tr},line:{color:c,width:0.4}});
function line2(s,x1,y1,x2,y2,o={}){
  const x=Math.min(x1,x2),y=Math.min(y1,y2),w=Math.abs(x2-x1),h=Math.abs(y2-y1);
  s.addShape("line",{x,y,w,h,flipV:((x2-x1)*(y2-y1))<0,line:Object.assign({color:MUTED,width:1},o)});
}

(()=>{
  const s=pres.addSlide(); s.background={color:WHITE};
  s.addText("Pipeline module cartoons — preview v3",{x:0.5,y:0.3,w:12.3,h:0.6,fontFace:SERIF,fontSize:25,bold:true,color:INK,margin:0});
  const qw=6.0,qh=2.8,cols=[0.5,6.83],rows=[1.2,4.28];
  const quads=[
    {c:0,r:0,label:"ProSeg",cap:"DAPI + 5 µm grabs stray transcripts; ProSeg fits cell shape to its transcripts",draw:drawProseg},
    {c:1,r:0,label:"BANKSY",cap:"the same local cell pattern recurs across tissue → a spatial niche",draw:drawBanksy},
    {c:0,r:1,label:"ENVI",cap:"each column = one gene; impute the whole transcriptome from paired single-cell data",draw:drawEnvi},
    {c:1,r:1,label:"STalign",cap:"align H&E and Xenium on one tissue via shared corner landmarks",draw:drawStalign},
  ];
  quads.forEach(q=>{
    const qx=cols[q.c],qy=rows[q.r];
    s.addShape("roundRect",{x:qx,y:qy,w:qw,h:qh,rectRadius:0.06,fill:{color:CARD},line:{color:BORDER,width:1}});
    s.addText(q.label,{x:qx+0.25,y:qy+0.12,w:qw-0.5,h:0.35,fontFace:SANS,fontSize:16,bold:true,color:ACCENT,margin:0});
    q.draw(s,qx+0.35,qy+0.6,qw-0.7,qh-1.15);
    s.addText(q.cap,{x:qx+0.25,y:qy+qh-0.44,w:qw-0.5,h:0.4,fontFace:SANS,fontSize:11,color:MUTED,margin:0});
  });

  // ---- ProSeg: same transcripts, different segmentation ----
  function drawProseg(s,ax,ay,aw,ah){
    const nuclei=[[0.5,0.55],[1.28,0.5],[0.92,1.12]];
    const own=[
      {c:ACCENT,pts:[[0.30,0.42],[0.62,0.46],[0.42,0.78],[0.70,0.70]]},
      {c:STEEL, pts:[[1.10,0.36],[1.42,0.40],[1.20,0.66],[1.46,0.62]]},
      {c:GREEN, pts:[[0.76,1.00],[1.04,1.02],[0.82,1.30],[1.08,1.24]]},
    ];
    const spurious=[[0.92,0.54,AMBER],[1.04,0.90,ACCENT],[0.64,0.92,STEEL]];
    function panel(ox,oy,mode){
      if(mode==="before"){ // uniform round DAPI-expansion cells (grab strays)
        nuclei.forEach(n=>s.addShape("ellipse",{x:ox+n[0]-0.45,y:oy+n[1]-0.45,w:0.9,h:0.9,fill:{color:"E2E5EA",transparency:25},line:{color:"9AA1AC",width:1}}));
      } else { // irregular transcript-aware cells (exclude strays)
        const el=[[0.5,0.57,0.64,0.46,-18],[1.28,0.5,0.6,0.44,14],[0.92,1.14,0.64,0.48,8]];
        el.forEach(e=>s.addShape("ellipse",{x:ox+e[0]-e[2]/2,y:oy+e[1]-e[3]/2,w:e[2],h:e[3],rotate:e[4],fill:{color:"DCE7F2",transparency:15},line:{color:STEEL,width:1.25}}));
      }
      nuclei.forEach(n=>dot(s,ox+n[0],oy+n[1],0.11,NUC));
      own.forEach(g=>g.pts.forEach(p=>dot(s,ox+p[0],oy+p[1],0.085,g.c)));
      spurious.forEach(p=>dot(s,ox+p[0],oy+p[1],0.085,p[2]));
    }
    panel(ax,ay,"before");
    line2(s,ax+aw*0.45,ay+ah/2,ax+aw*0.55,ay+ah/2,{color:STEEL,width:2.25,endArrowType:"triangle"});
    panel(ax+aw*0.55,ay,"after");
    s.addText("default: DAPI + 5 µm",{x:ax,y:ay+ah-0.02,w:aw*0.42,h:0.24,fontFace:SANS,fontSize:9,italic:true,color:MUTED,align:"center",margin:0});
    s.addText("ProSeg",{x:ax+aw*0.55,y:ay+ah-0.02,w:aw*0.42,h:0.24,fontFace:SANS,fontSize:9,italic:true,color:MUTED,align:"center",margin:0});
  }

  // ---- BANKSY: recurrent motif = niche ----
  function drawBanksy(s,ax,ay,aw,ah){
    // one EVEN lattice; the "motif" is a recurring 2x2 COLOR pattern (same spacing everywhere)
    const x0=ax+0.34, y0=ay+0.3, sx=0.53, sy=0.5, cols=9, rows=3;
    const bg=(i,j)=>PAL[(i*3 + j*2 + (i%2)) % 4];       // varied background colors
    const motif={"0,0":ACCENT,"1,0":GREEN,"0,1":STEEL,"1,1":AMBER};   // the recurring pattern
    const anchors=[[1,0],[4,1],[7,0]];                   // top-left cell of each motif block
    const mColor=(i,j)=>{ for(const [ai,aj] of anchors){ const di=i-ai,dj=j-aj;
      if(di>=0&&di<=1&&dj>=0&&dj<=1) return motif[`${di},${dj}`]; } return null; };
    for(let i=0;i<cols;i++)for(let j=0;j<rows;j++)
      dot(s, x0+i*sx, y0+j*sy, 0.15, mColor(i,j) || bg(i,j));
    anchors.forEach(a=>{ const cx=x0+(a[0]+0.5)*sx, cy=y0+(a[1]+0.5)*sy;
      s.addShape("ellipse",{x:cx-0.49,y:cy-0.47,w:0.98,h:0.94,fill:{color:WHITE,transparency:100},line:{color:NAVY,width:1.5,dashType:"dash"}}); });
    const la=anchors[2], lx=x0+(la[0]+0.5)*sx, ly=y0+(la[1]+0.5)*sy;
    s.addText("= niche",{x:lx-0.5,y:ly+0.52,w:1.0,h:0.24,fontFace:SANS,fontSize:9.5,italic:true,color:NAVY,align:"center",margin:0});
  }

  // ---- ENVI: columns = genes ----
  function drawEnvi(s,ax,ay,aw,ah){
    const rowTr=[12,40,58,28];
    function mat(ox,cols,rows,cw,ch,gap){
      for(let i=0;i<cols;i++)for(let j=0;j<rows;j++)
        s.addShape("rect",{x:ox+i*(cw+gap),y:ay+0.2+j*(ch+gap),w:cw,h:ch,fill:{color:GENECOLS[i%GENECOLS.length],transparency:rowTr[j%4]},line:{color:"FFFFFF",width:0.75}});
    }
    mat(ax+0.2,3,4,0.2,0.16,0.05);
    s.addText("266 genes",{x:ax+0.05,y:ay+1.28,w:1.3,h:0.24,fontFace:SANS,fontSize:9.5,color:MUTED,align:"center",margin:0});
    line2(s,ax+1.42,ay+0.7,ax+2.02,ay+0.7,{color:STEEL,width:2.25,endArrowType:"triangle"});
    mat(ax+2.2,9,4,0.2,0.16,0.04);
    s.addText("whole transcriptome",{x:ax+2.2,y:ay+1.28,w:2.3,h:0.24,fontFace:SANS,fontSize:9.5,color:MUTED,align:"center",margin:0});
  }

  // ---- STalign: two ROTATED tissue squares, corner-arrows, annotations overlaid ----
  function drawStalign(s,ax,ay,aw,ah){
    const th=40, rad=th*Math.PI/180, c=Math.cos(rad), sn=Math.sin(rad);
    const hw=0.52, hh=0.52;                                // square
    const Hc=[ax+1.0, ay+1.0], Xc=[ax+3.95, ay+1.0];
    const rot=(ctr,u,v)=>[ctr[0]+u*c - v*sn, ctr[1]+u*sn + v*c];
    const ann=[[-0.26,-0.16],[0.28,-0.04],[0.02,0.24]];    // annotation-circle local positions
    s.addShape("roundRect",{x:Hc[0]-hw,y:Hc[1]-hh,w:2*hw,h:2*hh,rectRadius:0.05,rotate:th,fill:{color:"F7EEF0"},line:{color:ACCENT,width:1.25}});
    s.addShape("roundRect",{x:Xc[0]-hw,y:Xc[1]-hh,w:2*hw,h:2*hh,rectRadius:0.05,rotate:th,fill:{color:WHITE},line:{color:STEEL,width:1.25}});
    // Xenium: roughly even grid of multicolor dots (rotated to fill the square)
    const gn=6;
    for(let gi=0; gi<gn; gi++)for(let gj=0; gj<gn; gj++){
      const u=(-hw+0.09)+(2*(hw-0.09))*gi/(gn-1)+R(-0.02,0.02);
      const v=(-hh+0.09)+(2*(hh-0.09))*gj/(gn-1)+R(-0.02,0.02);
      const p=rot(Xc,u,v); dot(s,p[0],p[1],0.06,PAL[(gi*2+gj)%4],20);
    }
    // annotation circles: in H&E, and overlaid at the SAME relative spot in Xenium
    ann.forEach(a=>{ const ph=rot(Hc,a[0],a[1]), px=rot(Xc,a[0],a[1]);
      [ph,px].forEach(q=>s.addShape("ellipse",{x:q[0]-0.09,y:q[1]-0.09,w:0.18,h:0.18,fill:{color:WHITE,transparency:100},line:{color:ACCENT,width:1.75}})); });
    // corner landmarks + ARROWS from H&E -> Xenium
    [[-hw,-hh],[hw,-hh],[-hw,hh],[hw,hh]].forEach(cn=>{ const ph=rot(Hc,cn[0],cn[1]), px=rot(Xc,cn[0],cn[1]);
      dot(s,ph[0],ph[1],0.1,NUC); dot(s,px[0],px[1],0.1,NUC);
      line2(s,ph[0],ph[1],px[0],px[1],{color:MUTED,width:1,dashType:"dash",endArrowType:"triangle"}); });
    s.addText("H&E",{x:Hc[0]-0.9,y:ay+0.02,w:0.85,h:0.24,fontFace:SANS,fontSize:10,bold:true,color:ACCENT,align:"right",margin:0});
    s.addText("Xenium",{x:Xc[0]+0.05,y:ay+0.02,w:0.95,h:0.24,fontFace:SANS,fontSize:10,bold:true,color:STEEL,align:"left",margin:0});
  }
})();

// goals slide (01/02 wording pending user) — kept for continuity
(()=>{
  const s=pres.addSlide(); s.background={color:WHITE};
  s.addText("GOALS  (01/02 wording pending)",{x:0.7,y:0.55,w:9,h:0.4,fontFace:SANS,fontSize:13,bold:true,color:ACCENT,charSpacing:2,margin:0});
  const items=[["01","What cell types and spatial niches build each tumor?"],
    ["02","Do those spatial niches match what a pathologist sees?"],
    ["03","In PA, how does MAPK activation shape glial–immune interaction?"]];
  let y=1.55;
  items.forEach(it=>{
    s.addText(it[0],{x:0.7,y,w:1.35,h:1.3,fontFace:SERIF,fontSize:44,bold:true,color:ACCENT,margin:0,valign:"top"});
    s.addText(it[1],{x:2.2,y:y+0.06,w:10.4,h:1.4,fontFace:SANS,fontSize:27,bold:true,color:INK,margin:0,valign:"top"});
    y+=1.75;
  });
})();

pres.writeFile({fileName:"/Users/youyun/plgg_deck/preview.pptx"}).then(f=>console.log("WROTE",f));
