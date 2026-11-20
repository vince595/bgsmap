"""Offline interactive HTML and publication-quality, editable-vector PDF plots."""
from __future__ import annotations

from html import escape
from pathlib import Path
import json

import numpy as np
import pandas as pd


def _minmax(z):
    vals = pd.to_numeric(z, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(vals).any():
        raise ValueError("Cannot plot a non-finite signal")
    lo,hi = float(np.nanmin(vals)),float(np.nanmax(vals))
    if hi <= lo:
        hi = lo + 1e-9
    return lo,hi


def write_plots(spot, region, gene, out, trait="Trait", invert_y=True, max_gene_bars=25):
    out = Path(out)
    out.mkdir(parents=True,exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    title = str(trait)
    coords = spot[["x","y"]].to_numpy(dtype=float)
    width = np.ptp(coords[:,0]); height = np.ptp(coords[:,1])
    from scipy.spatial import cKDTree
    near=cKDTree(coords).query(coords,k=2)[0][:,1]
    near=near[np.isfinite(near) & (near>0)]
    typical_spacing=float(np.median(near)) if len(near) else max(width,height)/50
    scale_ratio=typical_spacing/max(width,height,1e-9)
    # Marker diameters adapt to Visium spot density; remain vector in PDF.
    point_size_pdf=float(np.clip((0.77*325*scale_ratio)**2,9,650))
    point_size_html=float(np.clip(0.72*510*scale_ratio,4,35))
    plt.rcParams.update({"font.family":"DejaVu Sans","pdf.fonttype":42,"ps.fonttype":42,
                          "axes.linewidth":0.8,"font.size":9,"savefig.bbox":"tight"})
    # Plotly offline, hover details, stable numeric limits.
    fig = make_subplots(rows=1, cols=2, subplot_titles=["Analytic enrichment (-log10 p)","Burden–GSS score"], horizontal_spacing=0.08)
    customdata = np.column_stack([spot.spot.astype(str),spot.domain.astype(str),
                                   spot.p_analytic.to_numpy(),spot.q_bh_analytic.to_numpy(),
                                   spot.z_analytic.to_numpy()])
    for col,signal,scale,color_title in ((1,"logp_analytic","YlOrRd","-log10(analytic p)"),(2,"score","RdBu_r","Score")):
        fig.add_trace(go.Scattergl(
            x=spot.x,y=spot.y,mode="markers",name=color_title,
            customdata=customdata,
            marker=dict(size=point_size_html, color=spot[signal],colorscale=scale,
                        showscale=True, colorbar=dict(title=color_title,len=0.72,
                                                      x=0.445 if col==1 else 1.02)),
            hovertemplate="Spot %{customdata[0]}<br>Domain %{customdata[1]}<br>Analytic p %{customdata[2]:.3g}<br>Analytic BH q %{customdata[3]:.3g}<br>Analytic Z %{customdata[4]:.3g}<extra></extra>",
        ),row=1,col=col)
        fig.update_xaxes(showgrid=False,zeroline=False,visible=False,scaleanchor=f"y{'' if col==1 else '2'}",row=1,col=col)
        fig.update_yaxes(showgrid=False,zeroline=False,visible=False,autorange="reversed" if invert_y else True,row=1,col=col)
    fig.update_layout(title=f"bgsmap | {title}",template="plotly_white",height=660,
                      width=1250,showlegend=False,margin=dict(t=95,l=20,b=40,r=90))
    fig.write_html(str(out/"spatial_enrichment.html"),include_plotlyjs=True,full_html=True)
    # PDF artists remain vector circles and vector text, not raster screenshots.
    f,axs = plt.subplots(1,2,figsize=(12.4,5.4),layout="constrained")
    for ax,signal,cmap,label in ((axs[0],"logp_analytic","YlOrRd",r"$-\log_{10}(p_{analytic})$"),
                                  (axs[1],"score","RdBu_r","Burden–GSS score")):
        low,high = _minmax(spot[signal]);
        if signal=="score":
            low,high = -max(abs(low),abs(high)),max(abs(low),abs(high))
        s=ax.scatter(spot.x,spot.y,c=spot[signal],s=point_size_pdf,cmap=cmap,
                     vmin=low,vmax=high,edgecolors="none",rasterized=False)
        if invert_y: ax.invert_yaxis()
        ax.set_aspect("equal",adjustable="box"); ax.axis("off")
        ax.set_title(label,fontsize=11,pad=10)
        f.colorbar(s,ax=ax,shrink=0.77,fraction=0.05,pad=0.025)
    f.suptitle(f"{title} · Spatial burden enrichment",fontsize=15)
    f.savefig(out/"spatial_enrichment.pdf",metadata={"Title":f"{title} · bgsmap spatial map"})
    plt.close(f)
    # Domain graph: competitive contrast is primary for spatial specificity.
    table=region.sort_values("p_competitive_analytic").copy()
    fig2=go.Figure(go.Bar(y=table.domain.astype(str), x=table.logp_competitive_analytic,
                          orientation="h",marker_color="#3976A4",
                          customdata=np.stack([table.competitive_score,table.p_competitive_analytic,table.q_bh_competitive_analytic,table.z_competitive_analytic],axis=-1),
                          hovertemplate="%{y}<br>Competitive score: %{customdata[0]:.3g}<br>Analytic p: %{customdata[1]:.3g}<br>Analytic BH q: %{customdata[2]:.3g}<br>Analytic Z: %{customdata[3]:.3g}<extra></extra>"))
    fig2.update_layout(title=f"{title} | Competitive domain enrichment",template="plotly_white",
                       xaxis_title="-log10(competitive analytic p)",height=max(400,60*len(table)))
    fig2.update_yaxes(autorange="reversed")
    fig2.write_html(str(out/"region_enrichment.html"),include_plotlyjs=True)
    f,ax=plt.subplots(figsize=(8,max(3.2,0.45*len(table)+1.6)),layout="constrained")
    ax.barh(table.domain.astype(str),table.logp_competitive_analytic,color="#3976A4")
    ax.invert_yaxis();ax.set_xlabel("-log10(competitive analytic p)")
    ax.set_title(f"{title} · Spatial domain specificity");ax.spines[["right","top"]].set_visible(False)
    f.savefig(out/"region_enrichment.pdf");plt.close(f)
    good=gene.head(max_gene_bars).iloc[::-1]
    f,ax=plt.subplots(figsize=(9,max(4.5,0.26*len(good)+1.4)),layout="constrained")
    ax.barh(good.gene.astype(str),good.contribution.astype(float),color="#447EA6")
    ax.set_xlabel("Contribution to top-domain contrast (descriptive)")
    ax.set_title(f"{title} · Top gene contributions")
    ax.spines[["right","top"]].set_visible(False)
    f.savefig(out/"gene_contributions.pdf");plt.close(f)
    # A lightweight report references the companion interactive plot; a directory
    # can be copied intact and viewed offline.
    region_show = table[["domain","n_spots","mean_score","p_analytic","logp_analytic",
                         "competitive_score","z_competitive_analytic","p_competitive_analytic",
                         "logp_competitive_analytic","q_bh_competitive_analytic","p_holm_competitive_analytic"]]
    html_table = region_show.to_html(index=False,float_format=lambda x:f"{x:.4g}",border=0,escape=True,classes="data")
    gene_table = gene[["gene","p","contribution"]].head(30).to_html(index=False,float_format=lambda x:f"{x:.4g}",border=0,escape=True,classes="data")
    doc = """<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>bgsmap report</title><style>body{font:15px/1.6 system-ui,Arial,sans-serif;max-width:1200px;margin:20px auto;color:#243047;padding:0 16px}h1,h2{color:#153a5b}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:8px;border-bottom:1px solid #e1e6e9;text-align:left}thead{background:#eef3f7}.note{background:#f4f6f9;border-left:4px solid #708198;padding:12px}iframe{width:100%;height:690px;border:0}a{color:#17649c}</style></head><body>"""
    doc += f"<h1>bgsmap · {escape(title)}</h1><p>Gene-level rare-variant burden and spatial transcriptomics enrichment.</p>"
    doc += "<div class='note'><strong>Scientific caution:</strong> Primary P values use a normal approximation with the exact within-stratum label-permutation variance, not gsMap spatial S-LDSC. They are approximate and require independent calibration; optional finite permutations are saved only as diagnostics. Calibration assumes exchangeability within matching strata. Correlated genes, variable testing power, and unmatched gene covariates can invalidate nominal error rates. Gene contributions are descriptive, not causal.</div>"
    doc += '<h2>Spatial associations</h2><iframe title="spatial map" src="spatial_enrichment.html"></iframe>'
    doc += '<h2>Domain enrichment</h2><p>Competitive tests compare each domain to all other spots. Plots and this report use analytic (Gaussian-tail) P values; BH and Holm are within-trait adjustments. Optional permutation diagnostics are saved separately, not used for plotting.</p>' + html_table
    doc += '<h2>Top gene contributions</h2>' + gene_table
    doc += '<h2>Files</h2><p><a href="region_enrichment.html">Interactive region plot</a> · <a href="../spot_results.tsv.gz">Spot table</a> · <a href="../region_results.tsv">Region table</a> · <a href="../gene_contributions.tsv.gz">Gene contributions</a></p>'
    doc += '</body></html>'
    (out/"report.html").write_text(doc,encoding="utf-8")
