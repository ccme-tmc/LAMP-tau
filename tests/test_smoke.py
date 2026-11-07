"""Smoke test for the shipped DyAC example.

The test keeps the full 201-mode CFP file but uses only two temperatures and
disables figure generation, so it runs in well under a minute on a laptop.
"""
from pathlib import Path

from mag_anal.Anal_Tau import CAL_RELAX


def test_example_runs(tmp_path, monkeypatch):
    root = Path(__file__).resolve().parents[1]
    monkeypatch.chdir(root)

    cfg = (root / "examples/relaxation/CFPs_v1.yml").read_text(encoding="utf-8")
    cfg = cfg.replace("examples/relaxation/TPs_file", str(tmp_path / "out"))
    cfg = cfg.replace("examples/relaxation/TPs", str(tmp_path / "TPs"))
    cfg = cfg.replace("Trange1: 10,50,5", "Trange1: 10,20,10")
    cfg = cfg.replace("If_get_picture_Tt: True", "If_get_picture_Tt: False")

    yml = tmp_path / "CFPs_test.yml"
    yml.write_text(cfg, encoding="utf-8")

    tau = CAL_RELAX(str(yml))
    tau.kernel()

    assert (tmp_path / "out" / "OB_v1.dat").is_file()
