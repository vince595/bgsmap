import numpy as np
import pandas as pd
import pytest

from bgsmap.statistics import run_enrichment, bh, make_strata


def _toy(seed=11):
    rng=np.random.default_rng(seed)
    n_gene,n_spot=120,60
    domains=np.repeat(["1","2","3"],20)
    gss=rng.lognormal(0,0.4,(n_gene,n_spot)).astype(np.float32)
    gss[:20,40:60]*=8
    p=rng.uniform(0.05,1,n_gene);p[:20]=1e-7
    burden=pd.DataFrame({"gene":[f"g{i}" for i in range(n_gene)],"p":p,
                         "length":rng.uniform(1,1e5,n_gene)})
    coords=pd.DataFrame({"spot":[f"s{i}" for i in range(n_spot)],
                          "x":np.arange(n_spot),"y":np.zeros(n_spot),
                          "domain":domains}).set_index("spot",drop=False)
    return gss,burden,coords


def test_bh():
    assert np.allclose(bh([0.01,0.5,0.03]),[0.03,0.5,0.045])


def test_enrichment_recovers_planted_domain_and_reproduces():
    gss,b,coords=_toy()
    kw=dict(gss=gss,genes=b.gene.to_list(),spots=coords.spot.to_list(),burden=b,coord=coords,
            permutations=199,seed=23,chunk_size=37)
    r=run_enrichment(**kw)
    s=r.region.set_index("domain")
    assert s.loc["3","p_competitive_analytic"]<=0.01
    assert s.loc["3","logp_competitive_analytic"]>2
    assert np.all(r.spot.p_holm_analytic>=r.spot.p_analytic-1e-12)
    assert np.all(r.region.p_holm_competitive_analytic>=r.region.p_competitive_analytic-1e-12)
    assert "p_empirical" not in r.spot.columns
    assert "p_empirical" not in r.region.columns
    assert r.region_permutation.set_index('domain').loc['3','p_competitive_empirical'] <= 0.01
    v=run_enrichment(**kw)
    assert np.array_equal(r.spot_permutation.p_empirical,v.spot_permutation.p_empirical)
    assert np.array_equal(r.spot.p_analytic,v.spot.p_analytic)
    assert r.meta["permutations"]==199


def test_matching_and_no_variation():
    gss,b,coords=_toy()
    labels,info=make_strata(b,["length"],bins=3)
    assert len(np.unique(labels))==3
    assert info["matched"]
    with pytest.raises(ValueError,match="Matching bins too small"):
        make_strata(b,["length"],bins=20)
    with pytest.raises(ValueError,match="no variation"):
        b["p"]=0.5
        run_enrichment(gss,b.gene.to_list(),coords.spot.to_list(),b,coords,permutations=19)


def test_exact_finite_population_variance_and_underflow_free_logs():
    import itertools
    from bgsmap.statistics import _permutation_variance, _p_columns, _log10_adjust
    # Enumerate every within-stratum permutation and compare with the exact formula.
    a=np.array([[2,3,7,0,4,1],[-5,3,6,1,2,4]],dtype=np.float32)
    w=np.array([1.,-2.,1.,3.,-3.,0.],dtype=float)
    groups=[np.array([0,1,2]),np.array([3,4,5])]
    var=_permutation_variance(a,w,groups)
    samples=[]
    for p in itertools.permutations(groups[0]):
        for q in itertools.permutations(groups[1]):
            ww=np.r_[w[list(p)],w[list(q)]]
            samples.append(a @ ww / 6.)
    assert np.allclose(var,np.var(samples,axis=0),rtol=1e-12,atol=1e-12)
    t=pd.DataFrame({'id':[0,1,2]})
    _p_columns(t,'analytic',np.array([50.,5.,-1.]),np.ones(3))
    assert t.loc[0,'p_analytic']==0.0  # IEEE underflow is expected
    assert t.loc[0,'logp_analytic']>500  # no artificial -log10 P cap
    assert t.loc[0,'logq_bh_analytic']>500
    assert np.allclose(_log10_adjust([3,2,1],'bh'),[3-np.log10(3),2-np.log10(1.5),1])


def test_analytic_p_invariant_to_number_of_diagnostic_permutations():
    gss,b,coords=_toy()
    kw=dict(gss=gss,genes=b.gene.to_list(),spots=coords.spot.to_list(),burden=b,coord=coords,seed=21)
    zero=run_enrichment(**kw,permutations=0)
    many=run_enrichment(**kw,permutations=39)
    assert zero.spot_permutation is None and zero.region_permutation is None
    assert np.array_equal(zero.spot.logp_analytic,many.spot.logp_analytic)
    assert np.array_equal(zero.region.logp_competitive_analytic,many.region.logp_competitive_analytic)
