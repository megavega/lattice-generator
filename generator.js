/* Lattice Generator — Engine
   Computational geometry for SVG lattice patterns.
   Shared between EN and RU pages. */

const S3=Math.sqrt(3),
  C6=Array.from({length:6},(_,i)=>Math.cos((-90+60*i)*Math.PI/180)),
  S6=Array.from({length:6},(_,i)=>Math.sin((-90+60*i)*Math.PI/180));

/* geometry */
function hexV(cx,cy,R){return C6.map((c,i)=>[cx+R*c,cy+R*S6[i]])}
function cubeRh(cx,cy,R){
  const v=hexV(cx,cy,R),c=[cx,cy];
  return[[c,v[5],v[0],v[1]],[c,v[1],v[2],v[3]],[c,v[3],v[4],v[5]]];
}
function inset(vs,d){
  const n=vs.length;if(n<3||d<=0)return vs.slice();
  const sg=[];
  for(let i=0;i<n;i++){
    const[ax,ay]=vs[i],[bx,by]=vs[(i+1)%n],dx=bx-ax,dy=by-ay,L=Math.hypot(dx,dy);
    if(L<1e-10)continue;
    const nx=-dy/L*d,ny=dx/L*d;
    sg.push([[ax+nx,ay+ny],[bx+nx,by+ny]]);
  }
  const p=[];
  for(let i=0;i<sg.length;i++){const r=xs(sg[i],sg[(i+1)%sg.length]);if(r)p.push(r);}
  return p.length>=3?p:[];
}
function xs([[x1,y1],[x2,y2]],[[x3,y3],[x4,y4]]){
  const d=(x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);if(Math.abs(d)<1e-10)return null;
  const t=((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/d;
  return[x1+t*(x2-x1),y1+t*(y2-y1)];
}
function clipR(poly,x0,y0,x1,y1){
  function cl(pts,ok,ix){
    const o=[];
    for(let i=0;i<pts.length;i++){
      const c=pts[i],p=pts[(i-1+pts.length)%pts.length];
      if(ok(c)){if(!ok(p))o.push(ix(p,c));o.push(c);}
      else if(ok(p))o.push(ix(p,c));
    }return o;
  }
  const ix=(a,b,x)=>{if(Math.abs(b[0]-a[0])<1e-10)return[x,a[1]];const t=(x-a[0])/(b[0]-a[0]);return[x,a[1]+t*(b[1]-a[1])]};
  const iy=(a,b,y)=>{if(Math.abs(b[1]-a[1])<1e-10)return[a[0],y];const t=(y-a[1])/(b[1]-a[1]);return[a[0]+t*(b[0]-a[0]),y]};
  let r=poly.slice();
  r=cl(r,v=>v[0]>=x0,(a,b)=>ix(a,b,x0));
  r=cl(r,v=>v[0]<=x1,(a,b)=>ix(a,b,x1));
  r=cl(r,v=>v[1]>=y0,(a,b)=>iy(a,b,y0));
  r=cl(r,v=>v[1]<=y1,(a,b)=>iy(a,b,y1));
  return r;
}
function clipCirc(poly,cx,cy,r){
  const r2=r*r,n=poly.length;if(n<3)return[];
  const out=[];
  for(let i=0;i<n;i++){
    const c=poly[i],p=poly[(i-1+n)%n];
    const cIn=(c[0]-cx)**2+(c[1]-cy)**2<=r2;
    const pIn=(p[0]-cx)**2+(p[1]-cy)**2<=r2;
    if(cIn){
      if(!pIn){const ix=lcIx(p,c,cx,cy,r);for(const pt of ix)out.push(pt)}
      out.push(c);
    }else if(pIn){const ix=lcIx(p,c,cx,cy,r);for(const pt of ix)out.push(pt)}
  }
  return out.length>=3?out:[];
}
function lcIx(a,b,cx,cy,r){
  const dx=b[0]-a[0],dy=b[1]-a[1],fx=a[0]-cx,fy=a[1]-cy;
  const A=dx*dx+dy*dy,B=2*(fx*dx+fy*dy),C=fx*fx+fy*fy-r*r;
  const disc=B*B-4*A*C;if(disc<0)return[];
  const sq=Math.sqrt(disc),res=[];
  for(const t of[(-B-sq)/(2*A),(-B+sq)/(2*A)])
    if(t>=0&&t<=1)res.push([a[0]+t*dx,a[1]+t*dy]);
  return res;
}
function pA(v){let a=0;for(let i=0;i<v.length;i++){const j=(i+1)%v.length;a+=v[i][0]*v[j][1]-v[j][0]*v[i][1]}return Math.abs(a)/2}
function pD(v){let s=`M${v[0][0].toFixed(3)},${v[0][1].toFixed(3)}`;for(let i=1;i<v.length;i++)s+=`L${v[i][0].toFixed(3)},${v[i][1].toFixed(3)}`;return s+'Z'}
function pL(a,b){return`M${a[0].toFixed(3)},${a[1].toFixed(3)}L${b[0].toFixed(3)},${b[1].toFixed(3)}`}
function clipLine(a,b,x0,y0,x1,y1){
  let ax=a[0],ay=a[1],bx=b[0],by=b[1];
  function code(x,y){return(x<x0?1:0)|(x>x1?2:0)|(y<y0?4:0)|(y>y1?8:0)}
  let ca=code(ax,ay),cb=code(bx,by);
  for(let i=0;i<20;i++){
    if(!(ca|cb))return[[ax,ay],[bx,by]];
    if(ca&cb)return null;
    const c=ca||cb;let x,y;
    if(c&8){x=ax+(bx-ax)*(y1-ay)/(by-ay);y=y1}
    else if(c&4){x=ax+(bx-ax)*(y0-ay)/(by-ay);y=y0}
    else if(c&2){y=ay+(by-ay)*(x1-ax)/(bx-ax);x=x1}
    else{y=ay+(by-ay)*(x0-ax)/(bx-ax);x=x0}
    if(c===ca){ax=x;ay=y;ca=code(ax,ay)}else{bx=x;by=y;cb=code(bx,by)}
  }
  return null;
}

/* gradient */
function grT(hy,hx,H,W,mode){
  let t;
  switch(mode){
    case'uniform':t=1;break;
    case'center-v':t=1-2*Math.abs(hy/H-.5);t=(1-Math.cos(t*Math.PI))/2;break;
    case'edges-v':t=2*Math.abs(hy/H-.5);t=(1-Math.cos(t*Math.PI))/2;break;
    case'top-bottom':t=1-hy/H;break;
    case'bottom-top':t=hy/H;break;
    case'center-r':{const a=2*(hx/W-.5),b=2*(hy/H-.5);t=Math.max(0,1-Math.sqrt(a*a+b*b));t=(1-Math.cos(t*Math.PI))/2;break}
    case'edges-r':{const a=2*(hx/W-.5),b=2*(hy/H-.5);t=Math.min(1,Math.sqrt(a*a+b*b));t=(1-Math.cos(t*Math.PI))/2;break}
    default:t=1;
  }
  return t;
}
function br(hy,hx,H,W,mn,mx,mode){
  return mn+grT(hy,hx,H,W,mode)*(mx-mn);
}

/* generator */
let LP=[],LL=[],LW=0,LH=0,LR=0,LCirc=false;

function gen(W,H,R,mn,mx,gr,brd,rad,kerf,snap,pat,hflip,hhr,hbyR,htyR,hspR,hrsM,hfrM,hszMin,hdynR,hdynF,spts,srl,csn,csr,csm,gsz,csw,circ){
  if(mn>mx){const t=mn;mn=mx;mx=t}
  const b0=brd,b1=W-brd,b2=brd,b3=H-brd;
  const cx0=W/2,cy0=H/2;
  const cR=circ?Math.min(W,H)/2-brd:0;
  const clip=circ?(p=>clipCirc(p,cx0,cy0,cR)):(p=>clipR(p,b0,b2,b1,b3));
  const snapIn=circ?((x,y,r)=>Math.hypot(x-cx0,y-cy0)+r<=cR):((x,y,r)=>y-r>=b2&&y+r<=b3&&x-r>=b0&&x+r<=b1);
  const paths=[],lines=[];

  if(pat==='parquet'){
    const cs=R*2;
    const nSlots=Math.max(2,Math.min(8,Math.round(R/3)));
    const nc=Math.ceil(W/cs)+2, nr=Math.ceil(H/cs)+2;
    const ox=W/2, oy=H/2;
    for(let ri=-nr;ri<=nr;ri++){
      for(let ci=-nc;ci<=nc;ci++){
        const cx=ox+ci*cs, cy=oy+ri*cs;
        if(snap&&!snapIn(cx,cy,cs/2))continue;
        const bw=br(cy,cx,H,W,mn,mx,gr);
        const off=bw/2+kerf/2;
        const horiz=((ci%2+2)%2+((ri%2+2)%2))%2===0;
        for(let s=0;s<nSlots;s++){
          let rect;
          if(horiz){
            const sy=cy-cs/2+s*cs/nSlots;
            const x0=cx-cs/2+off, x1=cx+cs/2-off;
            const y0=sy+off, y1=sy+cs/nSlots-off;
            if(x1<=x0||y1<=y0)continue;
            rect=[[x0,y0],[x1,y0],[x1,y1],[x0,y1]];
          }else{
            const sx=cx-cs/2+s*cs/nSlots;
            const x0=sx+off, x1=sx+cs/nSlots-off;
            const y0=cy-cs/2+off, y1=cy+cs/2-off;
            if(x1<=x0||y1<=y0)continue;
            rect=[[x0,y0],[x1,y0],[x1,y1],[x0,y1]];
          }
          const cl=clip(rect);
          if(cl.length>=3&&pA(cl)>.5)paths.push(pD(cl));
        }
      }
    }
  }else if(pat==='stars'){
    const cs=R*2;
    const nc=Math.ceil(W/cs)+2,nr=Math.ceil(H/cs)+2;
    const N=spts||5;
    const rl=srl||0.5;
    const rOut=R*.92;
    const rIn=rOut*(1-rl);
    const step=Math.PI/N;
    const szMin=hszMin!=null?hszMin:0.3;
    for(let ri=-nr;ri<=nr;ri++){
      for(let ci=-nc;ci<=nc;ci++){
        const cx=cx0+ci*cs+((ri%2+2)%2)*cs/2;
        const cy=cy0+ri*cs;
        if(snap&&!snapIn(cx,cy,R))continue;
        const t=grT(cy,cx,H,W,gr);
        const sc=szMin+(1-szMin)*t;
        if(sc<0.01)continue;
        const off=kerf/2;
        const v=[];
        for(let i=0;i<2*N;i++){
          const a=-Math.PI/2+i*step;
          const r=i%2===0?rOut:rIn;
          v.push([cx+r*Math.cos(a), cy+r*Math.sin(a)]);
        }
        if(sc<0.999) for(let k=0;k<v.length;k++){v[k][0]=cx+(v[k][0]-cx)*sc;v[k][1]=cy+(v[k][1]-cy)*sc;}
        const ins=off>0?inset(v,off):v;
        if(!ins.length)continue;
        const cl=clip(ins);
        if(cl.length>=3&&pA(cl)>.5)paths.push(pD(cl));
      }
    }
  }else if(pat==='hearts'){
    const csx=R*2*(hfrM||1), csy=R*2*(hrsM||1), Na=20;
    const nc=Math.ceil(W/csx)+2, nr=Math.ceil(H/csy)+2;
    const hr=R*(hhr||.48);
    const sp=hr*(hspR||1);
    const byOff=R*(hbyR||-0.15);
    const tyOff=R*(htyR||0.72);
    const tdx=-sp, tdy=tyOff-byOff;
    const td=Math.hypot(tdx,tdy);
    const angT=Math.atan2(tdy,tdx);
    const halfC=Math.acos(hr/td);
    const tpA=angT-halfC;
    const tpL=Math.PI-tpA;
    const dipD=sp<hr?Math.sqrt(hr*hr-sp*sp):0;
    const dipAngR=(Math.atan2(-dipD,-sp)+2*Math.PI)%(2*Math.PI);
    const dipAngL=Math.atan2(-dipD,sp);
    const swR2=(2*Math.PI+tpA)-dipAngR;
    const swL2=(dipAngL+2*Math.PI)-tpL;
    function heartPoly(cx,cy,flip){
      const by=cy+byOff, ty=cy+tyOff;
      const v=[];
      for(let i=0;i<=Na;i++){
        const a=dipAngR+i*swR2/Na;
        v.push([cx+sp+hr*Math.cos(a), by+hr*Math.sin(a)]);
      }
      v.push([cx, ty]);
      for(let i=0;i<=Na;i++){
        const a=tpL+i*swL2/Na;
        v.push([cx-sp+hr*Math.cos(a), by+hr*Math.sin(a)]);
      }
      if(flip){
        for(let k=0;k<v.length;k++) v[k][1]=2*cy-v[k][1];
        v.reverse();
      }
      return v;
    }
    const szMin=hszMin!=null?hszMin:0.3;
    const dynRows=[];
    if(hdynR){
      dynRows.push({y:cy0,idx:0});
      let yU=cy0,iU=0;
      for(;;){const t=grT(yU,cx0,H,W,gr);const sc=szMin+(1-szMin)*t;
        const step=csy*Math.max(sc,0.12);yU-=step;if(yU<-R*2)break;iU--;dynRows.push({y:yU,idx:iU});}
      let yD=cy0,iD=0;
      for(;;){const t=grT(yD,cx0,H,W,gr);const sc=szMin+(1-szMin)*t;
        const step=csy*Math.max(sc,0.12);yD+=step;if(yD>H+R*2)break;iD++;dynRows.push({y:yD,idx:iD});}
    }else{
      for(let ri=-nr;ri<=nr;ri++) dynRows.push({y:cy0+ri*csy,idx:ri});
    }
    for(const row of dynRows){
      const cy=row.y,ri=row.idx;
      let csxD=csx,ncD=nc;
      if(hdynF){const tRow=grT(cy,cx0,H,W,gr);const scRow=szMin+(1-szMin)*tRow;
        csxD=csx*Math.max(scRow,0.12);ncD=Math.ceil(W/csxD)+2;}
      for(let ci=-ncD;ci<=ncD;ci++){
        const cx=cx0+ci*csxD+((ri%2+2)%2)*csxD/2;
        if(snap&&!snapIn(cx,cy,R))continue;
        const t=grT(cy,cx,H,W,gr);
        const sc=szMin+(1-szMin)*t;
        if(sc<0.01)continue;
        const off=kerf/2;
        let flip=false;
        if(hflip==='up') flip=true;
        else if(hflip==='alt') flip=((ri+ci)%2+2)%2===1;
        else if(hflip==='row') flip=((ri%2)+2)%2===1;
        else if(hflip==='half') flip=cy<H/2;
        const v=heartPoly(cx,cy,flip);
        if(sc<0.999) for(let k=0;k<v.length;k++){v[k][0]=cx+(v[k][0]-cx)*sc;v[k][1]=cy+(v[k][1]-cy)*sc;}
        const ins=off>0?inset(v,off):v;
        if(!ins.length)continue;
        const cl=clip(ins);
        if(cl.length>=3&&pA(cl)>.5)paths.push(pD(cl));
      }
    }
  }else if(pat==='gyroid'){
    const period=R*2, kk=2*Math.PI/period;
    const zAng=(gsz!=null?gsz:0.13)*2*Math.PI;
    const sZ=Math.sin(zAng), cZ=Math.cos(zAng);
    const res=Math.max(0.4, period/25);
    const nx=Math.ceil(W/res), ny=Math.ceil(H/res), nxp=nx+1;
    function gf(x,y){return Math.sin(kk*x)*Math.cos(kk*y)+Math.sin(kk*y)*cZ+sZ*Math.cos(kk*x)}
    const gp=new Float32Array(nxp*(ny+1)), gn=new Float32Array(nxp*(ny+1));
    for(let j=0;j<=ny;j++)for(let i=0;i<=nx;i++){
      const x=i*res, y=j*res, f=gf(x,y), t=br(y,x,H,W,mn,mx,gr)/R;
      const idx=j*nxp+i; gp[idx]=f-t; gn[idx]=f+t;
    }
    function gyLoops(gv){
      const segs=[];
      for(let j=0;j<ny;j++){const r0=j*nxp,r1=r0+nxp;
        for(let i=0;i<nx;i++){
          const v0=gv[r0+i],v1=gv[r0+i+1],v2=gv[r1+i+1],v3=gv[r1+i];
          const c=(v0>0?1:0)|(v1>0?2:0)|(v2>0?4:0)|(v3>0?8:0);
          if(!c||c===15)continue;
          const x0=i*res,y0=j*res,x1=x0+res,y1=y0+res;
          const lp=(va,vb,xa,ya,xb,yb)=>{const u=-va/(vb-va);return[xa+u*(xb-xa),ya+u*(yb-ya)]};
          const eT=lp(v0,v1,x0,y0,x1,y0),eR=lp(v1,v2,x1,y0,x1,y1);
          const eB=lp(v3,v2,x0,y1,x1,y1),eL=lp(v0,v3,x0,y0,x0,y1);
          switch(c){
            case 1:case 14:segs.push([eT,eL]);break;
            case 2:case 13:segs.push([eT,eR]);break;
            case 3:case 12:segs.push([eR,eL]);break;
            case 4:case 11:segs.push([eR,eB]);break;
            case 6:case 9:segs.push([eT,eB]);break;
            case 7:case 8:segs.push([eL,eB]);break;
            case 5:{const m=(v0+v1+v2+v3)/4;
              if(m>0){segs.push([eT,eR]);segs.push([eL,eB])}
              else{segs.push([eT,eL]);segs.push([eR,eB])}break}
            case 10:{const m=(v0+v1+v2+v3)/4;
              if(m>0){segs.push([eT,eL]);segs.push([eB,eR])}
              else{segs.push([eT,eR]);segs.push([eB,eL])}break}
          }
        }
      }
      const ep=res*0.001;
      const pk=p=>Math.round(p[0]/ep)+','+Math.round(p[1]/ep);
      const adj=new Map();
      for(let i=0;i<segs.length;i++){
        const ka=pk(segs[i][0]),kb=pk(segs[i][1]);
        if(!adj.has(ka))adj.set(ka,[]);if(!adj.has(kb))adj.set(kb,[]);
        adj.get(ka).push({i,e:0});adj.get(kb).push({i,e:1});
      }
      const used=new Uint8Array(segs.length);
      const loops=[];
      for(let si=0;si<segs.length;si++){
        if(used[si])continue;used[si]=1;
        const pts=[segs[si][0],segs[si][1]];
        const sk=pk(segs[si][0]);let ck=pk(segs[si][1]);
        for(let it=0;it<segs.length;it++){
          if(ck===sk&&pts.length>2)break;
          const nb=adj.get(ck);if(!nb)break;
          let fd=false;
          for(const n of nb){if(used[n.i])continue;used[n.i]=1;
            pts.push(n.e===0?segs[n.i][1]:segs[n.i][0]);
            ck=pk(pts[pts.length-1]);fd=true;break;}
          if(!fd)break;
        }
        if(pts.length>=4&&ck===sk)loops.push(pts);
      }
      return loops;
    }
    const kOff=kerf/2;
    for(const grid of[gp,gn]){
      const isPos=grid===gp;
      const loops=gyLoops(grid);
      for(const lp of loops){
        let sx=0,sy=0;
        for(const p of lp){sx+=p[0];sy+=p[1]}
        sx/=lp.length;sy/=lp.length;
        const fc=gf(sx,sy), t=br(sy,sx,H,W,mn,mx,gr)/R;
        if(isPos?fc<=t:fc>=-t)continue;
        if(snap){let out=false;for(const p of lp){if(circ?Math.hypot(p[0]-cx0,p[1]-cy0)>cR:(p[0]<b0||p[0]>b1||p[1]<b2||p[1]>b3)){out=true;break}}if(out)continue}
        let sa=0;
        for(let i=0;i<lp.length;i++){const j=(i+1)%lp.length;sa+=lp[i][0]*lp[j][1]-lp[j][0]*lp[i][1]}
        if(sa<0)lp.reverse();
        let poly=lp;
        if(kOff>0){poly=inset(poly,kOff);if(!poly.length)continue}
        const cl=clip(poly);
        if(cl.length>=3&&pA(cl)>.3)paths.push(pD(cl));
      }
    }
  }else if(pat==='hilbert'){
    /* Hilbert curve — space-filling fractal for laser cutting.
       R controls cell size; order is derived from panel dimensions.
       mn/mx control bridge width; gsz controls slot width ratio. */
    const iW=W-2*brd, iH=H-2*brd;
    const side=Math.min(iW,iH);
    const order=Math.max(2,Math.min(6,Math.round(Math.log2(side/(R*2)))));
    const n=1<<order;
    const cell=side/n;
    const ox=circ?cx0-side/2:brd+(iW-n*cell)/2;
    const oy=circ?cy0-side/2:brd+(iH-n*cell)/2;

    function d2xy(nn,d){
      let x=0,y=0,s=1;
      while(s<nn){
        const rx=1&(d/2|0),ry=1&(d^rx);
        if(ry===0){if(rx===1){x=s-1-x;y=s-1-y;}const t=x;x=y;y=t;}
        x+=s*rx;y+=s*ry;
        d=Math.floor(d/4);s*=2;
      }
      return[x,y];
    }
    const total=n*n,hpts=[];
    for(let d=0;d<total;d++){
      const[gx,gy]=d2xy(n,d);
      hpts.push([ox+(gx+0.5)*cell, oy+(gy+0.5)*cell]);
    }

    /* Chaikin corner-cutting: smooth 90° turns into curves */
    function chaikin(pts,iter){
      let p=pts;
      for(let it=0;it<iter;it++){
        const q=[p[0]];
        for(let i=0;i<p.length-1;i++){
          const[ax,ay]=p[i],[bx,by]=p[i+1];
          q.push([ax*0.75+bx*0.25, ay*0.75+by*0.25]);
          q.push([ax*0.25+bx*0.75, ay*0.25+by*0.75]);
        }
        q.push(p[p.length-1]);
        p=q;
      }
      return p;
    }
    const smooth=chaikin(hpts,3);

    /* slot width ratio from gsz (0→thin 0.25, 1→wide 0.45) */
    const slotR=gsz!=null?0.25+gsz*0.2:0.35;
    const perpW=cell*slotR-kerf/2;
    if(perpW>=0.15){
      /* bridge at structural connectors: every 4th segment midpoint
         (aligns with fractal sub-curve boundaries).
         Bridge = rectangular band perpendicular to segment direction. */
      const stride=4;
      const bridges=[];
      const bridgeHalf=cell*0.12;
      for(let i=0;i<hpts.length-1;i++){
        if((i+1)%stride===0){
          const mx=(hpts[i][0]+hpts[i+1][0])/2, my=(hpts[i][1]+hpts[i+1][1])/2;
          const dx=hpts[i+1][0]-hpts[i][0], dy=hpts[i+1][1]-hpts[i][1];
          const l=Math.hypot(dx,dy)||1;
          bridges.push({mx,my,ux:dx/l,uy:dy/l});
        }
      }

      /* compute normals along smooth curve */
      const norms=[];
      for(let i=0;i<smooth.length;i++){
        let dx,dy;
        if(i===0){dx=smooth[1][0]-smooth[0][0];dy=smooth[1][1]-smooth[0][1];}
        else if(i===smooth.length-1){dx=smooth[i][0]-smooth[i-1][0];dy=smooth[i][1]-smooth[i-1][1];}
        else{dx=smooth[i+1][0]-smooth[i-1][0];dy=smooth[i+1][1]-smooth[i-1][1];}
        const l=Math.hypot(dx,dy)||1;
        norms.push([-dy/l,dx/l]);
      }

      /* build slot polygons, bridging at original segment midpoints */
      let left=[],right=[];
      let wasInBridge=true;
      function addCap(idx){
        /* add straight perpendicular endcap at point idx */
        const sx=smooth[idx][0],sy=smooth[idx][1];
        const nx=norms[idx][0],ny=norms[idx][1];
        left.push([sx+nx*perpW, sy+ny*perpW]);
        right.push([sx-nx*perpW, sy-ny*perpW]);
      }
      function flushSlot(){
        if(left.length<2){left=[];right=[];return;}
        const poly=left.concat(right.reverse());
        left=[];right=[];
        const cl=clip(poly);
        if(cl.length>=3&&pA(cl)>0.3)paths.push(pD(cl));
      }

      for(let i=0;i<smooth.length;i++){
        const sx=smooth[i][0],sy=smooth[i][1];
        let inBridge=i===0||i===smooth.length-1;
        if(!inBridge){
          for(let m=0;m<bridges.length;m++){
            const b=bridges[m];
            const along=(sx-b.mx)*b.ux+(sy-b.my)*b.uy;
            const perp=(sx-b.mx)*(-b.uy)+(sy-b.my)*b.ux;
            if(Math.abs(along)<bridgeHalf&&Math.abs(perp)<cell*0.6){inBridge=true;break;}
          }
        }
        if(inBridge){
          if(!wasInBridge){
            /* entering bridge — cap the slot end */
            addCap(i>0?i-1:i);
            flushSlot();
          }
          wasInBridge=true;
        }else{
          if(wasInBridge){
            /* exiting bridge — cap the slot start */
            addCap(i);
          }
          left.push([sx+norms[i][0]*perpW, sy+norms[i][1]*perpW]);
          right.push([sx-norms[i][0]*perpW, sy-norms[i][1]*perpW]);
          wasInBridge=false;
        }
      }
      flushSlot();
    }
  }else{
    const cs=R*S3,rs=R*1.5;
    const nc=Math.ceil(W/cs)+2,nr=Math.ceil(H/rs)+2;
    const hmx=R*S3/2;
    for(let ri=-nr;ri<=nr;ri++){
      for(let ci=-nc;ci<=nc;ci++){
        const hx=cx0+ci*cs+((ri%2+2)%2)*cs/2;
        const hy=cy0+ri*rs;
        if(snap&&!snapIn(hx,hy,R))continue;
        const bw=br(hy,hx,H,W,mn,mx,gr);
        if(pat==='cube'){
          const off=bw/2+kerf/2;
          if(off>=hmx*.95)continue;
          for(const rh of cubeRh(hx,hy,R)){
            const ins=inset(rh,off);
            if(!ins.length)continue;
            const cl=clip(ins);
            if(cl.length>=3&&pA(cl)>.5)paths.push(pD(cl));
          }
        }else if(pat==='honeycomb'){
          const off=bw/2+kerf/2;
          if(off>=R*.9)continue;
          const hex=hexV(hx,hy,R);
          const ins=inset(hex,off);
          if(!ins.length)continue;
          const cl=clip(ins);
          if(cl.length>=3&&pA(cl)>.5)paths.push(pD(cl));
        }else if(pat==='cubeslats'){
          const nS=csn||3;
          const margin=csm!=null?csm:0.1;
          const ratio=csr||0.8;
          const bwVal=br(hy,hx,H,W,mn,mx,gr);
          const bridgeT=bwVal/R;
          const avail=1-2*margin;
          const maxSlotT=(avail-bridgeT*(nS-1))/nS;
          const slotT=maxSlotT*ratio;
          const pad=(maxSlotT-slotT)/2;
          if(slotT>=0.02){
            const faces=cubeRh(hx,hy,R);
            const kOff=kerf/2;
            for(const face of faces){
              const[C0,A,B,D]=face;
              for(let s=0;s<nS;s++){
                const t0=margin+pad+s*(maxSlotT+bridgeT);
                const t1=t0+slotT;
                const qa=[C0[0]+(A[0]-C0[0])*t0, C0[1]+(A[1]-C0[1])*t0];
                const qb=[C0[0]+(A[0]-C0[0])*t1, C0[1]+(A[1]-C0[1])*t1];
                const qc=[D[0]+(B[0]-D[0])*t1, D[1]+(B[1]-D[1])*t1];
                const qd=[D[0]+(B[0]-D[0])*t0, D[1]+(B[1]-D[1])*t0];
                let slot=[qa,qb,qc,qd];
                if(kOff>0){slot=inset(slot,kOff);if(!slot.length)continue;}
                const cl=clip(slot);
                if(cl.length>=3&&pA(cl)>.3)paths.push(pD(cl));
              }
            }
          }
        }
      }
    }
  }
  LP=paths;LL=lines;LW=W;LH=H;LR=rad;
  return{paths,lines,csw:csw||0.5};
}

function svg(paths,lines,W,H,rad,mode,lw,circ){
  lines=lines||[];lw=lw||0.5;
  let s='';
  if(mode!=='inline')s+='<?xml version="1.0" encoding="UTF-8"?>\n';
  const wa=mode==='inline'?'':`width="${W}mm" height="${H}mm" `;
  s+=`<svg xmlns="http://www.w3.org/2000/svg" ${wa}viewBox="0 0 ${W} ${H}">\n`;
  if(mode==='cut'){
    s+=`<g id="cuts">\n`;
    for(const d of paths)s+=`<path d="${d}" fill="none" stroke="red" stroke-width="0.1"/>\n`;
    for(const d of lines)s+=`<path d="${d}" fill="none" stroke="red" stroke-width="0.1"/>\n`;
    if(circ) s+=`</g>\n<g id="outline"><circle cx="${W/2}" cy="${H/2}" r="${Math.min(W,H)/2}" fill="none" stroke="red" stroke-width="0.1"/></g>\n`;
    else s+=`</g>\n<g id="outline"><rect x="0" y="0" width="${W}" height="${H}" rx="${rad}" ry="${rad}" fill="none" stroke="red" stroke-width="0.1"/></g>\n`;
  }else{
    if(circ) s+=`<circle cx="${W/2}" cy="${H/2}" r="${Math.min(W,H)/2}" fill="#333"/>\n`;
    else s+=`<rect width="${W}" height="${H}" rx="${rad}" fill="#333"/>\n`;
    if(paths.length)s+=`<path d="${paths.join('')}" fill="#eee"/>\n`;
    if(lines.length)s+=`<path d="${lines.join('')}" fill="none" stroke="#eee" stroke-width="${lw}"/>\n`;
  }
  s+=`</svg>`;return s;
}
