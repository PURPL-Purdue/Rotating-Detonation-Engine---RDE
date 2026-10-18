"""Run the ORIGINAL HADES/DATA ANALYSIS .m files in GNU Octave for cross-validation.

Nothing in the original folder is modified.  Temporary shims are written to a
scratch directory and put first on the Octave path:

* ``rand`` / ``randn`` / ``randi`` replay pre-generated number streams, so the
  Python port (driven by :class:`ReplayRNG` with the same streams) and the
  MATLAB code see identical random numbers.
* graphics calls (``figure``, ``plot``, ``title``, ...) are no-ops; ``plot``
  records its X/Y data so function outputs that MATLAB only plotted can be
  compared.
* ``readmatrix`` / ``writematrix`` / ``subtitle`` / ``sgtitle`` fill Octave gaps
  (``readmatrix`` also redirects the hard-coded ``/Users/Abby/...`` path).

The originals were removed from the working tree after the port. They are
looked up under ``$RDE_MATLAB_SOURCE_ROOT/HADES/DATA ANALYSIS`` (see
``tests/matlab_sources.py`` for how to extract them from git history). These
tests are skipped when Octave or the original .m files are unavailable; the
hard-coded reference numbers in the test files keep the regression check.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pytest
from matlab_sources import OCTAVE, original, skip_reason

ORIGINAL_DIR = original("HADES", "DATA ANALYSIS")


def live_skip_reason(*m_files: str, octave: bool = True, packages=()) -> str | None:
    """``None`` if the live check can run, else the actionable skip reason."""
    return skip_reason([ORIGINAL_DIR / f for f in m_files] or [ORIGINAL_DIR], octave=octave,
                       octave_packages=packages)


def requires_originals(*m_files: str, octave: bool = True, packages=()):
    """Skip marker for tests that need Octave (and Octave ``packages``) and/or the original files."""
    reason = live_skip_reason(*m_files, octave=octave, packages=packages)
    return pytest.mark.skipif(reason is not None, reason=reason or "")


class ReplayRNG:
    """numpy-Generator-like object replaying fixed uniform/normal/integer streams."""

    def __init__(self, uniform, normal, integer_uniform):
        self.u = np.asarray(uniform, dtype=float)
        self.n = np.asarray(normal, dtype=float)
        self.iu = np.asarray(integer_uniform, dtype=float)
        self.pu = self.pn = self.pi = 0

    @staticmethod
    def _take(seq, pos, size):
        k = 1 if size is None else int(np.prod(size))
        out = seq[pos : pos + k]
        if len(out) < k:
            raise RuntimeError("replay stream exhausted")
        return (out[0] if size is None else out.reshape(size)), pos + k

    def random(self, size=None):
        out, self.pu = self._take(self.u, self.pu, size)
        return out

    def standard_normal(self, size=None):
        out, self.pn = self._take(self.n, self.pn, size)
        return out

    def integers(self, low, high, size=None, endpoint=False):
        hi = high if endpoint else high - 1
        u, self.pi = self._take(self.iu, self.pi, size)
        return low + np.floor(u * (hi - low + 1)).astype(int)


def make_streams(seed: int, n_uniform: int = 100, n_normal: int = 200_000, n_int: int = 1000):
    rs = np.random.RandomState(seed)
    return rs.random_sample(n_uniform), rs.standard_normal(n_normal), rs.random_sample(n_int)


_SHIMS = {
    "rand.m": """function r = rand(varargin)
  r = replay_take__('U', varargin{:});
end
""",
    "randn.m": """function r = randn(varargin)
  r = replay_take__('N', varargin{:});
end
""",
    "randi.m": """function r = randi(lim, varargin)
  if isscalar(lim), lo = 1; hi = lim; else, lo = lim(1); hi = lim(2); end
  u = replay_take__('I', varargin{:});
  r = lo + floor(u * (hi - lo + 1));
end
""",
    "replay_take__.m": """function r = replay_take__(which, varargin)
  global REPLAY_U REPLAY_N REPLAY_I REPLAY_POS
  args = varargin(cellfun(@isnumeric, varargin));
  if isempty(args), dims = [1 1];
  elseif numel(args) == 1 && isscalar(args{1}), dims = [args{1} args{1}];
  elseif numel(args) == 1, dims = args{1};
  else, dims = [args{:}]; end
  k = prod(dims);
  switch which
    case 'U', seq = REPLAY_U; p = REPLAY_POS(1);
    case 'N', seq = REPLAY_N; p = REPLAY_POS(2);
    case 'I', seq = REPLAY_I; p = REPLAY_POS(3);
  end
  r = reshape(seq(p+1:p+k), dims);
  switch which
    case 'U', REPLAY_POS(1) = p + k;
    case 'N', REPLAY_POS(2) = p + k;
    case 'I', REPLAY_POS(3) = p + k;
  end
end
""",
    "plot.m": """function varargout = plot(varargin)
  global PLOT_CAPTURE
  nums = varargin(cellfun(@isnumeric, varargin));
  if numel(nums) >= 2, PLOT_CAPTURE{end+1} = {nums{1}(:), nums{2}(:)};
  else, PLOT_CAPTURE{end+1} = {(1:numel(nums{1}))', nums{1}(:)}; end
  varargout = cell(1, nargout);
end
""",
    "readmatrix.m": """function data = readmatrix(varargin)
  global READMATRIX_FILE
  data = dlmread(READMATRIX_FILE, ',', 1, 0);
end
""",
    "writematrix.m": """function writematrix(data, filename)
  dlmwrite(filename, data, 'precision', '%.15g');
end
""",
}
for _name in ("figure", "hold", "xlabel", "ylabel", "title", "legend", "grid", "subplot",
              "subtitle", "sgtitle", "imshow", "set", "clc"):  # fmt: skip
    _SHIMS[f"{_name}.m"] = f"function varargout = {_name}(varargin)\n  varargout = cell(1, nargout);\nend\n"
_SHIMS["gca.m"] = "function h = gca()\n  h = 0;\nend\n"


def _write_vec(path: Path, x) -> None:
    np.savetxt(path, np.ravel(np.asarray(x, dtype=float)), fmt="%.17g")


def run_octave(workdir: Path, statements: str, outputs: list[str], streams=None, readmatrix_file=None,
               packages=(), timeout=600):  # fmt: skip
    """Run ``statements`` in Octave and return ``{name: array}`` for each output.

    ``outputs`` are workspace variable names, or ``plot<k>x`` / ``plot<k>y`` for
    the X/Y data of the k-th (1-based) ``plot`` call.
    """
    workdir = Path(workdir)
    shim_dir = workdir / "shims"
    out_dir = workdir / "out"
    shim_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, text in _SHIMS.items():
        (shim_dir / name).write_text(text)

    setup = ["global REPLAY_U REPLAY_N REPLAY_I REPLAY_POS PLOT_CAPTURE READMATRIX_FILE",
             "PLOT_CAPTURE = {}; REPLAY_POS = [0 0 0];"]  # fmt: skip
    if streams is not None:
        for tag, seq in zip("UNI", streams):
            _write_vec(workdir / f"replay_{tag}.txt", seq)
            setup.append(f"REPLAY_{tag} = load('-ascii', '{workdir / f'replay_{tag}.txt'}');")
    if readmatrix_file is not None:
        setup.append(f"READMATRIX_FILE = '{readmatrix_file}';")
    for pkg in packages:
        setup.append(f"pkg load {pkg};")

    save = ["global PLOT_CAPTURE"]
    for name in outputs:
        fname = out_dir / f"{name}.txt"
        if name.startswith("plot") and name[-1] in "xy":
            k, axis = int(name[4:-1]), {"x": 1, "y": 2}[name[-1]]
            expr = f"PLOT_CAPTURE{{{k}}}{{{axis}}}"
        else:
            expr = name
        save.append(f"dlmwrite('{fname}', double({expr}(:)), 'precision', '%.17g');")

    # shims first (shadow built-ins), then the untouched original folder
    script = "\n".join(
        [f"addpath('{ORIGINAL_DIR}');", f"addpath('{shim_dir}');", "warning('off', 'all');",
         *setup, statements, *save]
    )  # fmt: skip
    (workdir / "driver.m").write_text(script)
    proc = subprocess.run(
        [OCTAVE, "--no-gui", "--no-window-system", "--quiet", "--norc", str(workdir / "driver.m")],
        cwd=workdir, capture_output=True, text=True, timeout=timeout,
    )  # fmt: skip
    if proc.returncode != 0:
        raise RuntimeError(f"Octave failed:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    result = {name: np.loadtxt(out_dir / f"{name}.txt", delimiter=",", ndmin=1) for name in outputs}
    result["_stdout"] = proc.stdout
    return result
