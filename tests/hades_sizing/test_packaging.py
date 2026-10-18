"""Packaging smoke test: data files must ship in the wheel, legacy files must not.

Builds the wheel with ``uv build`` from a *pristine temporary copy* of the
source tree (no ``build/``, ``*.egg-info``, ``.venv``, ``__pycache__``, so no
stale setuptools state can leak in), checks the wheel contents, installs it
(no deps) into a fresh venv outside the checkout, exposes the already
installed third-party dependencies through a plain path entry (which does
NOT activate the repo's editable install), and runs a sizing-model call that
reads the bundled xlsx.
"""

from __future__ import annotations

import json
import os
import shutil
import site
import subprocess
import sys
import venv
import zipfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]

pytestmark = pytest.mark.skipif(shutil.which("uv") is None, reason="uv not available")

REQUIRED = [
    "HADES/sizing_model/HADES_size_cellsizedata.xlsx",
    "HADES/feed_sizing/flow_requirements.csv",
    "HADES/feed_sizing/flow_requirements_manifold_input.csv",
    "HADES/feed_sizing/orifice_blowdown_input.csv",
    "HADES/data_analysis/sample_thermocouple_data.csv",
    "MODELING/moc/LICENSE",
]
LEGACY_DIRS = ("HADES/SIZING MODEL/", "HADES/FEED SIZING/", "HADES/DATA ANALYSIS/", "MODELING/MODELING_size_")
LEGACY_SUFFIXES = (".m", ".exe", ".asv")

SCRIPT = r"""
import contextlib, io, json
import HADES.sizing_model.HADES_size_annulusgap_fillheight as m
from HADES.sizing_model.HADES_size_main import HADES_size_main
with contextlib.redirect_stdout(io.StringIO()):
    r = HADES_size_main(plot=False)
print("RESULT" + json.dumps({"file": m.__file__, "rows": len(m.read_cellsize_data()),
                             "thrust": r.thrust, "pc": r.avg_chamber_p}))
"""


@pytest.fixture(scope="module")
def wheel(tmp_path_factory) -> Path:
    tmp = tmp_path_factory.mktemp("wheel")
    src = tmp / "src"
    shutil.copytree(
        REPO, src,
        ignore=shutil.ignore_patterns("build", "dist", "*.egg-info", ".venv", "__pycache__", ".git",
                                      ".pytest_cache", "*.pyc", ".DS_Store"),
    )
    dist = tmp / "dist"
    subprocess.run(["uv", "build", "--wheel", "-q", "-o", str(dist), str(src)], check=True,
                   capture_output=True, timeout=300)
    return next(dist.glob("rde-*.whl"))


def test_wheel_contains_required_assets_and_no_legacy_files(wheel):
    names = zipfile.ZipFile(wheel).namelist()
    for req in REQUIRED:
        assert req in names, req
    assert any(n.startswith("HADES/data_analysis/bos_print_assets/") and n.endswith(".svg") for n in names)
    assert any(n.startswith("HADES/data_analysis/bos_print_assets/") and n.endswith(".png") for n in names)
    legacy = [n for n in names if n.startswith(LEGACY_DIRS) or n.endswith(LEGACY_SUFFIXES)]
    assert legacy == []


def test_wheel_runs_sizing_model_outside_checkout(wheel, tmp_path):
    env_dir = tmp_path / "venv"
    venv.create(env_dir, with_pip=False, symlinks=(sys.platform != "win32"))
    py = env_dir / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
    subprocess.run(["uv", "pip", "install", "-q", "--no-deps", "--offline", "--python", str(py), str(wheel)],
                   check=True, capture_output=True, timeout=300)
    # the site-packages entry (on Windows getsitepackages() also lists the venv root first)
    env_site = subprocess.run(
        [str(py), "-c", "import site; print([p for p in site.getsitepackages() if p.endswith('site-packages')][0])"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    (Path(env_site) / "_dev_deps.pth").write_text("\n".join(site.getsitepackages()) + "\n")
    run_dir = tmp_path / "elsewhere"
    run_dir.mkdir()
    # Minimal environment; Windows needs SYSTEMROOT & co. for the interpreter to start, and an
    # explicit MPLCONFIGDIR avoids depending on a home directory.
    env = {"MPLBACKEND": "Agg", "MPLCONFIGDIR": str(tmp_path / "mplconfig"), "PATH": os.defpath}
    if sys.platform == "win32":
        env.update({k: os.environ[k] for k in ("SYSTEMROOT", "SYSTEMDRIVE", "WINDIR", "TEMP", "TMP", "USERPROFILE")
                    if k in os.environ})
    p = subprocess.run([str(py), "-I", "-c", SCRIPT], cwd=run_dir, capture_output=True, text=True,
                       timeout=300, env=env)
    assert p.returncode == 0, p.stderr[-3000:]
    res = json.loads([ln for ln in p.stdout.splitlines() if ln.startswith("RESULT")][-1][6:])
    assert Path(res["file"]).resolve().is_relative_to(env_dir.resolve())  # code came from the wheel
    assert res["rows"] == 12
    assert res["thrust"] == pytest.approx(863.7703526250729, rel=1e-12)
    assert res["pc"] == pytest.approx(142.74514850616087, rel=1e-12)
