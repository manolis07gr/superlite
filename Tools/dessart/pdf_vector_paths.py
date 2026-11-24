import re, sys, json
def load(svg):
    s=open(svg).read(); out=[]
    for p in re.findall(r'<path ([^>]*)/>', s):
        m=re.search(r'd="([^"]*)"', p)
        if not m: continue
        toks=re.findall(r'[MLCZmlcz]|-?\d+\.?\d*(?:e-?\d+)?', m.group(1))
        P=[]; i=0; cmd=None
        while i<len(toks):
            t=toks[i]
            if t in 'MLCZmlcz': cmd=t; i+=1; continue
            if cmd=='C': P.append((float(toks[i+4]),float(toks[i+5]))); i+=6
            else: P.append((float(t),float(toks[i+1]))); i+=2
        tr=re.search(r'transform="matrix\(([^)]*)\)"', p)
        if tr:
            a,b,c,d,e,f=[float(x) for x in tr.group(1).split(',')]
            P=[(a*x+c*y+e, b*x+d*y+f) for x,y in P]
        st=re.search(r'style="([^"]*)"', p); st=st.group(1) if st else ''
        col=re.search(r'stroke:([^;]*)', st); w=re.search(r'stroke-width:([^;]*)', st)
        out.append(dict(P=P, stroke=col.group(1) if col else None, w=float(w.group(1)) if w else None, fill='fill:none' not in st))
    return out
if __name__=='__main__':
    L=load(sys.argv[1])
    for k,o in enumerate(L):
        P=o['P']; xs=[a for a,b in P]; ys=[b for a,b in P]
        w=max(xs)-min(xs); h=max(ys)-min(ys)
        if len(P) in (4,5) and w>100 and h>40: print('FRAME', k, '%.2f %.2f %.2f %.2f' % (min(xs),max(xs),min(ys),max(ys)), o['stroke'])
