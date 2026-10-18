"""BOS_Background_Generation.m -> HADES.data_analysis.bos_background_generation.

The oracles are the files MATLAB wrote with ``rng(1)``. They are stored
read-only under ``reference/bos/`` and pinned by SHA-256; the 15 MB 60 mm-FOV
SVG is pinned by hash only. The port must reproduce them byte-for-byte and
pixel-for-pixel. No test writes to the fixtures or to the package.
"""

import hashlib
import math
import warnings
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

from HADES.data_analysis import bos_background_generation as bos

from ._octave_harness import ORIGINAL_DIR

FIXTURES = Path(__file__).parent / "reference" / "bos"
Image.MAX_IMAGE_PIXELS = None  # 11811 x 8268 poster

# SHA-256 of the files MATLAB wrote (HADES/DATA ANALYSIS, commit daa9d8b)
SHA_SVG = "b4abeb0a29915733ae139be506a913d2995584b1b2c68138be267fe83585a058"
SHA_PNG_FILE = "09888f327f2ad67d12671c17e5a60b9a4bc3e2fceaaed5111e941c92d432e226"
SHA_PNG_PIXELS = "cc84477a8045463878704126bdc208e00d591c69fef64cd06d5448cb9ee09bbf"  # uint8 (8268, 11811), C order
SHA_SVG_FOV60 = "b7c3dc82ce31ba584a42309726ce35fe70cb2b4d135cd1c19bb5fde9f70ad8f9"


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _pixels_sha(img) -> str:
    return _sha(np.ascontiguousarray(img, dtype=np.uint8).tobytes())


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    out = tmp_path_factory.mktemp("bos")
    return bos.BOS_Background_Generation(output_folder=out, plot=False, verbose=False)


def test_fixtures_are_the_matlab_files():
    assert _sha((FIXTURES / "BOS_random_dot_background.svg").read_bytes()) == SHA_SVG
    assert _sha((FIXTURES / "BOS_random_dot_background.png").read_bytes()) == SHA_PNG_FILE
    assert _pixels_sha(np.array(Image.open(FIXTURES / "BOS_random_dot_background.png"))) == SHA_PNG_PIXELS


def test_svg_byte_identical_to_matlab_output(generated):
    data = generated.svg_file.read_bytes()
    assert _sha(data) == SHA_SVG
    assert data == (FIXTURES / "BOS_random_dot_background.svg").read_bytes()


def test_png_pixel_identical_to_matlab_output(generated):
    new = Image.open(generated.png_file)
    assert new.mode == "L" and new.size == (11811, 8268)
    pixels = np.array(new)
    assert _pixels_sha(generated.img) == SHA_PNG_PIXELS
    if _pixels_sha(pixels) != SHA_PNG_PIXELS:  # the hash is the oracle; count differences to help debugging
        ref = np.array(Image.open(FIXTURES / "BOS_random_dot_background.png"))
        pytest.fail(f"{np.count_nonzero(pixels != ref)} pixels differ from the MATLAB PNG")


def test_fov60mm_svg_matches_older_matlab_run(tmp_path):
    # The stray top-level SVG in the MATLAB folder came from a run with a 60 mm FOV.
    r = bos.BOS_Background_Generation(FOV_width_mm=60, FOV_height_mm=60, output_folder=tmp_path,
                                      make_png=False, plot=False, verbose=False)  # fmt: skip
    assert r.number_of_dots == 231483
    assert _sha(r.svg_file.read_bytes()) == SHA_SVG_FOV60


def test_parameters_hand_calculation(generated):
    g = generated
    assert g.mm_per_camera_pixel == pytest.approx(300 / 256)
    assert g.nominal_dot_diameter_mm == pytest.approx(2.5 * 300 / 256)  # 2.9297 mm
    n = -math.log(0.7) * 500 * 350 / (math.pi * (2.5 * 300 / 256 / 2) ** 2)
    assert g.number_of_dots == round(n) == 9259
    assert (g.poster_width_px, g.poster_height_px) == (11811, 8268)
    assert g.dot_diameter_mm.min() >= 0.8 * g.nominal_dot_diameter_mm
    assert g.dot_diameter_mm.max() <= 1.2 * g.nominal_dot_diameter_mm


def test_black_coverage_statistics(generated):
    # Poisson coverage with jittered diameters: E[A] = A_nom * (1 + j^2/3)
    expected = 1 - math.exp(-generated.number_of_dots * generated.nominal_dot_area_mm2 * (1 + 0.2**2 / 3) / (500 * 350))
    black = np.mean(generated.img == 0)
    assert black == pytest.approx(expected, abs=0.003)
    assert black == pytest.approx(0.30, abs=0.01)


def test_printed_summary(tmp_path, capsys):
    bos.BOS_Background_Generation(output_folder=tmp_path, make_png=False, plot=False)
    out = capsys.readouterr().out
    expected_lines = [
        "                 BOS BACKGROUND PARAMETERS",
        "Camera frame                  : 256 x 256 pixels",
        "Background FOV                : 300.000 x 300.000 mm",
        "Horizontal physical scale     : 1.17188 mm/pixel",
        "Vertical physical scale       : 1.17188 mm/pixel",
        "Target dot diameter           : 2.50 camera pixels",
        "Nominal physical dot diameter : 2.9297 mm",
        "Target black coverage         : 30.0 %",
        "Number of dots                : 9259",
        "Poster size                   : 500.0 x 350.0 mm",
        "PNG resolution                : 600 dpi",
        "PNG raster dimensions:",
        "11811 x 8268 pixels",
        "1. Print at 100% / ACTUAL SIZE.",
        '2. Do NOT use "Fit to Page".',
        "4. PNG should be interpreted at 600 dpi.",
    ]
    for line in expected_lines:
        assert line in out.splitlines()
    assert str(tmp_path / "BOS_random_dot_background.svg") in out


def test_warnings():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        bos.BOS_Background_Generation(FOV_width_mm=600, FOV_height_mm=300, output_folder=None,
                                      make_png=False, plot=False, verbose=False)  # fmt: skip
    msgs = [str(x.message) for x in w]
    assert any("scales differ by 66.67%" in m for m in msgs)
    assert "Poster width is smaller than camera FOV width." in msgs


def test_preview_plot_smoke():
    bos.BOS_Background_Generation(poster_width_mm=50, poster_height_mm=35, FOV_width_mm=30, FOV_height_mm=30,
                                  print_resolution_dpi=100, output_folder=None, plot=True, verbose=False)  # fmt: skip


def test_default_output_goes_to_cwd_not_package(tmp_path, monkeypatch, no_show):
    monkeypatch.chdir(tmp_path)
    before = {p.name: _sha(p.read_bytes()) for p in bos.PRINT_ASSETS.iterdir()}
    r = bos.BOS_Background_Generation(make_png=False, plot=False, verbose=False)
    assert r.svg_file == tmp_path / "BOS Backgrounds" / "BOS_random_dot_background.svg"
    assert _sha(r.svg_file.read_bytes()) == SHA_SVG
    bos.main(["--output-dir", str(tmp_path / "cli"), "--no-plot"])
    assert _sha((tmp_path / "cli" / "BOS_random_dot_background.svg").read_bytes()) == SHA_SVG
    assert _pixels_sha(np.array(Image.open(tmp_path / "cli" / "BOS_random_dot_background.png"))) == SHA_PNG_PIXELS
    assert {p.name: _sha(p.read_bytes()) for p in bos.PRINT_ASSETS.iterdir()} == before


def test_print_assets_are_the_matlab_files():
    assert _sha((bos.PRINT_ASSETS / "BOS_random_dot_background.svg").read_bytes()) == SHA_SVG
    assert _sha((bos.PRINT_ASSETS / "BOS_random_dot_background.png").read_bytes()) == SHA_PNG_FILE
    assert _sha((bos.PRINT_ASSETS / "BOS_random_dot_background_FOV60mm.svg").read_bytes()) == SHA_SVG_FOV60


@pytest.mark.skipif(not (ORIGINAL_DIR / "BOS Backgrounds").exists(), reason="original MATLAB folder removed")
def test_hashes_match_original_matlab_folder():
    assert _sha((ORIGINAL_DIR / "BOS Backgrounds" / "BOS_random_dot_background.svg").read_bytes()) == SHA_SVG
    assert _sha((ORIGINAL_DIR / "BOS Backgrounds" / "BOS_random_dot_background.png").read_bytes()) == SHA_PNG_FILE
    assert _sha((ORIGINAL_DIR / "BOS_random_dot_background.svg").read_bytes()) == SHA_SVG_FOV60
