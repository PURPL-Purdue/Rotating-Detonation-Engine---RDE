"""Live cross-check against the original MATLAB scripts run in GNU Octave.

The originals were removed from the working tree after the port; they are
looked up under ``$RDE_MATLAB_SOURCE_ROOT/DEIMOS`` (see ``tests/matlab_sources.py``
for how to extract them from git history). Skipped automatically, with the
recovery instructions, when Octave or the original .m files are unavailable.
The originals are never modified: Octave's missing ``optimoptions`` is provided
by a shim written to a temp directory.
"""

import os
import re
import subprocess
import warnings
from pathlib import Path

import numpy as np
import pytest

from matlab_sources import OCTAVE, original, skip_reason

from DEIMOS.DEIMOS_sys_CdCalc_v1 import DEIMOS_sys_CdCalc_v1
from DEIMOS.DEIMOS_sys_FannoFlowCalc_v2 import DEIMOS_sys_FannoFlowCalc_v2

DEIMOS_DIR = original("DEIMOS")


def _run_octave(tmp_path: Path, script: Path, varnames, shim: bool) -> dict:
    reason = skip_reason([script])
    if reason:
        pytest.skip(reason)
    lines = []
    if shim:
        (tmp_path / "optimoptions.m").write_text(
            "function o = optimoptions(solver, varargin)\n"
            "  o = optimset(varargin{:}, 'TolX', 1e-15, 'TolFun', 1e-15, 'MaxIter', 1000);\nend\n"
        )
        lines.append(f"addpath('{tmp_path}');")
    lines += [
        "set(0, 'defaultfigurevisible', 'off');",
        f"source('{script}');",
        "printf('DUMP_BEGIN\\n');",
    ]
    for v in varnames:
        lines.append(f"printf('{v} ='); printf(' %.17g', {v}); printf('\\n');")
    driver = tmp_path / "driver.m"
    driver.write_text("\n".join(lines) + "\n")
    proc = subprocess.run(
        [OCTAVE, "--no-gui", "--norc", "--quiet", str(driver)],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=120,
        env={**os.environ, "QT_QPA_PLATFORM": "offscreen"},
    )
    assert proc.returncode == 0, proc.stderr
    dump = proc.stdout.split("DUMP_BEGIN", 1)[1]
    out = {}
    for line in dump.strip().splitlines():
        name, vals = line.split("=", 1)
        out[name.strip()] = np.array([float(t) for t in vals.split()])
    printed = re.sub(r"\x1b\[[0-9;]*[A-Za-z]", "", proc.stdout.split("DUMP_BEGIN", 1)[0])
    return out, printed


def test_cdcalc_live_octave(tmp_path):
    names = ["Cd_inc", "Cd_comp", "mdot_comp", "Re_D1", "P0", "f", "Kn", "iter"]
    ref, printed = _run_octave(tmp_path, DEIMOS_DIR / "DEIMOS_sys_CdCalc_v1.m", names, shim=False)
    res = DEIMOS_sys_CdCalc_v1()
    for n in names:
        np.testing.assert_allclose(getattr(res, n), ref[n][0], rtol=1e-12, err_msg=n)
    assert "Cd (incompressible) = 0.613529" in printed


def test_fanno_v2_live_octave(tmp_path):
    names = ["m_dotT", "m_dot", "A_gap", "Rinner", "Router", "iter", "mu", "Re", "f", "M_TubeIN", "P01", "P2", "T_stag"]
    ref, printed = _run_octave(tmp_path, DEIMOS_DIR / "DEIMOS_sys_FannoFlowCalc_v2.m", names, shim=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = DEIMOS_sys_FannoFlowCalc_v2()
    for n in names:
        np.testing.assert_allclose(np.atleast_1d(getattr(res, n)), ref[n], rtol=1e-10, err_msg=n)
    assert "P01 Methane Upstream Pressure   : 32.25 bar | 467.77 psi" in printed
