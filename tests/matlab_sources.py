"""Locate GNU Octave and the original MATLAB sources for the live cross-checks.

The MATLAB originals were deleted from the working tree once the Python ports
were validated. They are still in git history: commit ``2a30c95`` (``main``
before the rewrite) has all of them. To run the live Octave tests::

    mkdir -p /tmp/rde_src
    git archive 2a30c95 | tar -x -C /tmp/rde_src
    RDE_MATLAB_SOURCE_ROOT=/tmp/rde_src uv run pytest

``RDE_MATLAB_SOURCE_ROOT`` is the root of the extracted tree (the folder that
holds ``DEIMOS/``, ``HADES/`` and ``MODELING/``). Without it the repository root
is searched, where the originals no longer exist, so the live tests skip.
"""

from __future__ import annotations

import functools
import os
import shutil
import subprocess
from pathlib import Path

ENV_VAR = "RDE_MATLAB_SOURCE_ROOT"
RECOVERY_COMMIT = "2a30c95"
REPO = Path(__file__).resolve().parents[1]
RECOVERY_HINT = (
    f"recover the originals with `mkdir -p /tmp/rde_src && git archive {RECOVERY_COMMIT} | "
    f"tar -x -C /tmp/rde_src` and set {ENV_VAR}=/tmp/rde_src"
)

_HOMEBREW_OCTAVE = Path("/opt/homebrew/bin/octave")
OCTAVE: str | None = shutil.which("octave") or (str(_HOMEBREW_OCTAVE) if _HOMEBREW_OCTAVE.exists() else None)


def source_root() -> Path:
    """Root of the tree holding the original MATLAB sources."""
    env = os.environ.get(ENV_VAR)
    return Path(env).expanduser().resolve() if env else REPO


def original(*parts: str) -> Path:
    """Path of an original file or folder, e.g. ``original("HADES", "DATA ANALYSIS")``."""
    return source_root().joinpath(*parts)


def missing(paths) -> list[Path]:
    return [p for p in paths if not Path(p).exists()]


@functools.cache
def octave_has_package(name: str) -> bool:
    """True if ``pkg load <name>`` succeeds in the Octave found above."""
    if OCTAVE is None:
        return False
    try:
        proc = subprocess.run([OCTAVE, "--no-gui", "--norc", "--quiet", "--eval", f"pkg load {name}"],
                              stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return False
    return proc.returncode == 0


def skip_reason(paths, *, octave: bool = True, octave_packages=()) -> str | None:
    """``None`` if the live check can run, else an actionable skip reason.

    ``octave_packages`` are Octave Forge packages the check loads (e.g. ``signal``);
    they are only probed when Octave and the sources are both present.
    """
    problems, fixes = [], []
    if octave and OCTAVE is None:
        problems.append("GNU Octave not found on PATH")
        fixes.append("install GNU Octave")
    lost = missing(paths)
    if lost:
        shown = str(lost[0]) + (f" and {len(lost) - 1} more" if len(lost) > 1 else "")
        problems.append(f"original MATLAB sources not found ({shown})")
        fixes.append(RECOVERY_HINT)
    if not problems and octave:
        for pkg in octave_packages:
            if not octave_has_package(pkg):
                problems.append(f"Octave package '{pkg}' not installed")
                fixes.append(f"install it with `octave --eval \"pkg install -forge {pkg}\"` "
                             f"(or the OS package, e.g. `apt install octave-{pkg}`)")
    if not problems:
        return None
    return f"live check skipped: {'; '.join(problems)}. To run it, {' and '.join(fixes)}"
