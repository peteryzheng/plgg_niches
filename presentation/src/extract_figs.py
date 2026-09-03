import re, base64, os, sys
html_path = sys.argv[1]; outdir = sys.argv[2]
os.makedirs(outdir, exist_ok=True)
data = open(html_path, encoding="utf-8", errors="ignore").read()
# embedded figures in quarto HTML: <img src="data:image/png;base64,....">
imgs = re.findall(r'data:image/(png|jpe?g);base64,([A-Za-z0-9+/=]+)', data)
kept=0
for i,(ext,b64) in enumerate(imgs):
    raw = base64.b64decode(b64)
    if len(raw) < 6000:   # skip tiny icons/logos
        continue
    e = "jpg" if ext.startswith("jp") else "png"
    path=f"{outdir}/img_{i:02d}.{e}"
    open(path,"wb").write(raw)
    print(f"img_{i:02d}.{e}  {len(raw)//1024} KB")
    kept+=1
print("kept", kept, "of", len(imgs), "->", outdir)
