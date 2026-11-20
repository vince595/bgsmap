"""Regression tests for MOSTA-style GSS spot-ID alignment."""
import json

import numpy as np
import pandas as pd
import pytest

from bgsmap.io import load_coords, _decode_h5
from bgsmap.pipeline import _align_spot_ids


def _coords(ids):
    return pd.DataFrame({"spot": ids, "x": np.arange(len(ids)),
                         "y": np.arange(len(ids)), "domain": "brain"}).set_index("spot", drop=False)


def test_mosta_style_underscore_spots_match_exactly(tmp_path):
    h5py = pytest.importorskip("h5py")
    ids = ["190_194", "190_195", "190_196", "700_239", "700_240"]
    h5ad = tmp_path / "MOSTA.h5ad"
    with h5py.File(h5ad, "w") as h:
        obs = h.create_group("obs")
        obs.attrs["_index"] = "_index"
        obs.create_dataset("_index", data=np.array([x.encode() for x in ids]))
        obs.create_dataset("annotation", data=np.array([b"forebrain"] * 5))
        h.create_group("obsm").create_dataset("spatial", data=np.arange(10).reshape(5, 2))
    coord = load_coords(h5ad=h5ad, annotation="annotation")
    keep, aligned = _align_spot_ids(pd.Index(ids), coord)
    assert keep.tolist() == [0, 1, 2, 3, 4]
    assert list(aligned.index) == ids
    assert _decode_h5(np.array([np.bytes_("190_194")])).tolist() == ["190_194"]


def test_anndata_recovery_uses_exact_ids_not_row_order(monkeypatch):
    ids = pd.Index(["190_194", "190_195", "190_196", "700_239", "700_240"])
    # Mimic a broken low-level H5 index reader with malformed identifiers.
    fast = _coords(["b'" + x + "'" for x in ids])
    correct = _coords(list(reversed(ids)))  # H5AD spatial order differs from GSS order.
    monkeypatch.setattr("bgsmap.pipeline.load_coords_anndata",
                        lambda *a, **kw: correct)
    keep, aligned = _align_spot_ids(ids, fast, h5ad="MOSTA.h5ad")
    assert keep.tolist() == list(range(5))
    # The downstream lookup is by IDs; never assume same row ordering.
    assert list(aligned.loc[list(ids), "spot"]) == list(ids)


def test_unmatched_spots_show_two_id_samples():
    spots = pd.Index(["190_194", "190_195", "190_196", "190_197", "190_198"])
    with pytest.raises(ValueError, match="First GSS spot IDs") as err:
        _align_spot_ids(spots, _coords(["X1", "X2", "X3", "X4", "X5"]))
    assert "first spatial IDs" in str(err.value)
    assert "190_194" in str(err.value)


def test_inspect_spot_alignment_without_loading_expression(tmp_path, capsys):
    h5py = pytest.importorskip("h5py")
    from bgsmap.cli import main
    h5ad = tmp_path / "st.h5ad"
    with h5py.File(h5ad, "w") as h:
        obs = h.create_group("obs")
        obs.attrs["_index"] = "_index"
        obs.create_dataset("_index", data=np.array([b"190_194", b"190_195", b"190_196", b"190_197", b"190_198"]))
        obs.create_dataset("annotation", data=np.array([b"A"] * 5))
        h.create_group("obsm").create_dataset("spatial", data=np.arange(10).reshape(5, 2))
    gss = tmp_path / "gss.tsv"
    pd.DataFrame({"HUMAN_GENE_SYM": ["A", "B"],
                  **{x: [1,2] for x in ["190_194", "190_195", "190_196", "190_197", "190_198"]}}).to_csv(gss,sep="\t",index=False)
    main(["inspect", "--gss", str(gss), "--h5ad", str(h5ad),
          "--check-spots", "--annotation", "annotation"])
    doc = json.loads(capsys.readouterr().out)
    assert doc["spot_alignment"]["matched"] == 5
    assert doc["spot_alignment"]["reader"] == "h5py"
