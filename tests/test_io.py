import numpy as np
import pandas as pd
import pytest

from bgsmap.io import load_burden, load_gss, load_coords, read_gene_map


def test_burden_masks_duplicates_and_bonferroni(tmp_path):
    p=tmp_path/"burden.tsv"
    pd.DataFrame({"Gene":["A","A","B"],"P":[0.001,0.01,0.5],"mask":["1","2","1"]}).to_csv(p,sep="\t",index=False)
    with pytest.raises(ValueError, match="Duplicate burden genes"):
        load_burden(p)
    a,_=load_burden(p,mask="1")
    assert set(a.gene)=={"A","B"}
    b,_=load_burden(p,duplicate_policy="bonferroni")
    assert np.isclose(float(b.set_index("gene").loc["A","p"]),0.002)


def test_gss_gene_ids_and_ortholog_map(tmp_path):
    p=tmp_path/"gss.csv"
    pd.DataFrame({"HUMAN_GENE_SYM":["ENSG000001.12","abc"],
                  **{f"s{i}":[float(i),float(i+1)] for i in range(5)}}).to_csv(p,index=False)
    genes,spots,gss=load_gss(p)
    assert list(genes)==["ENSG000001","ABC"]
    assert gss.shape==(2,5)
    m=tmp_path/"map.tsv"
    pd.DataFrame({"st_gene":["a","b"],"burden_gene":["A","A"]}).to_csv(m,sep="\t",index=False)
    with pytest.raises(ValueError,match="one-to-one"):
        read_gene_map(m)


def test_coords_table(tmp_path):
    p=tmp_path/"coords.tsv"
    pd.DataFrame({"spot":["a","b"],"x":[0,1],"y":[2,3],"region":["r1","r2"]}).to_csv(p,sep="\t",index=False)
    c=load_coords(coords=p,annotation="region")
    assert list(c.domain)==["r1","r2"]


def test_anndata_h5ad_lightweight_reader(tmp_path):
    h5py=pytest.importorskip("h5py")
    p=tmp_path/"st.h5ad"
    with h5py.File(p,"w") as h:
        obs=h.create_group("obs")
        obs.attrs["_index"]="_index"
        obs.create_dataset("_index",data=np.array([b"a",b"b",b"c"]))
        group=obs.create_group("domain")
        group.create_dataset("categories",data=np.array([b"D1",b"D2"]))
        group.create_dataset("codes",data=np.array([0,1,-1],dtype=np.int8))
        h.create_group("obsm").create_dataset("spatial",data=np.array([[1,2],[3,4],[5,6]],dtype=float))
    c=load_coords(h5ad=p,annotation="domain")
    assert list(c.domain)==["D1","D2","unassigned"]
