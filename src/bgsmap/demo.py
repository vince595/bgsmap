"""Small, labeled SYNTHETIC examples. Not FinnGen/AZ observations."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd


def generate_demo(out,seed=19):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(seed)
    nx,ny=12,8; nspot=nx*ny; ngene=160
    x,y=np.meshgrid(np.arange(nx),np.arange(ny));x=x.ravel();y=y.ravel()
    domains=(x//4+1).astype(str)
    coords=pd.DataFrame({"spot":[f"spot_{i:03d}" for i in range(nspot)],
                         "x":x,"y":y,"domain":domains})
    coords.to_csv(out/"spatial_coords.tsv",sep="\t",index=False)
    genes=[f"GENE{i:04d}" for i in range(ngene)]
    values=rng.lognormal(0,0.7,size=(ngene,nspot)).astype(np.float32)
    planted=np.arange(25)
    values[np.ix_(planted,np.flatnonzero(domains=="3"))]*=9
    gss=pd.DataFrame(values,columns=coords.spot)
    gss.insert(0,"HUMAN_GENE_SYM",genes)
    gss.to_csv(out/"GSS.tsv.gz",sep="\t",index=False,compression="gzip")
    p=rng.uniform(0.03,1,ngene)
    p[planted]=10**(-rng.uniform(2.1,7.5,len(planted)))
    burden=pd.DataFrame({"gene":genes,"p":p,"mask":"Mask1.0.01",
                         "gene_length":rng.integers(1000,250_000,ngene)})
    burden.to_csv(out/"burden.tsv.gz",sep="\t",index=False,compression="gzip")
    (out/"README.txt").write_text("SYNTHETIC demonstration only. Domain 3 has planted trait-linked GSS.\n",encoding="utf-8")
    return out
