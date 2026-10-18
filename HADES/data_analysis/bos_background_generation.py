"""BOS random-dot background generator (port of ``BOS_Background_Generation.m``).

Generates one random black-dot pattern for background-oriented schlieren (BOS)
and saves the SAME pattern as

    1. SVG (vector, preferred for professional printing)
    2. PNG (8-bit grayscale raster at ``print_resolution_dpi``)

into a ``BOS Backgrounds`` folder in the CURRENT WORKING DIRECTORY (as MATLAB
did with ``fullfile(pwd, 'BOS Backgrounds')``), or into ``output_folder`` /
``--output-dir``.  Dot size is chosen so
that each dot spans ``target_dot_pixels`` pixels of the 256 x 256 camera frame,
and the dot count so that the expected black coverage of the overlapping random
dots is ``target_black_fraction`` (``f = 1 - exp(-N*A_dot/A_total)``).

Random stream: with the default ``random_seed = 1`` the dot positions and sizes
are bit-identical to MATLAB's ``rng(1)`` output (see ``_matlab_compat``), so the
MATLAB-generated SVG/PNG are reproduced exactly.  Copies of those files for
printing are in ``bos_print_assets/`` next to this module (never written by
this code); the test oracles are separate, read-only fixtures under
``tests/hades_data_analysis/reference/bos/``.

Run as a script::

    uv run python -m HADES.data_analysis.bos_background_generation [--output-dir DIR]
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ._matlab_compat import matlab_rng, matlab_round

# ======================================================================
#  USER INPUTS (same defaults as the MATLAB script)
# ======================================================================

# Camera frame used during BOS testing
CAMERA_WIDTH_PX = 256
CAMERA_HEIGHT_PX = 256

# Physical field of view at the BACKGROUND plane: the physical width/height of
# the background visible in the 256 x 256 camera frame (e.g. 300 mm x 300 mm).
FOV_WIDTH_MM = 300.0
FOV_HEIGHT_MM = 300.0

# Desired apparent dot diameter on camera (typical BOS target: ~2-3 pixels)
TARGET_DOT_PIXELS = 2.5

# Physical printed background size (only needs to be larger than the camera FOV)
POSTER_WIDTH_MM = 500.0
POSTER_HEIGHT_MM = 350.0

# Approximate fraction of the background covered by black dots (0.30 = ~30 %)
TARGET_BLACK_FRACTION = 0.30

# Random dot size variation (0.20 = +/- 20 % variation in dot diameter)
DIAMETER_JITTER = 0.20

# PNG print resolution (300 dpi = good, 600 dpi = very good)
PRINT_RESOLUTION_DPI = 600

# Random seed (same seed = same pattern every time)
RANDOM_SEED = 1

# Output filename (without extension) and folder (relative paths resolve against
# the current working directory at call time, like MATLAB's pwd)
BASE_FILENAME = "BOS_random_dot_background"
OUTPUT_FOLDER = "BOS Backgrounds"

# MATLAB-generated backgrounds kept for printing (read-only; not an output location)
PRINT_ASSETS = Path(__file__).resolve().parent / "bos_print_assets"

_RULE = "============================================================"


@dataclass
class BOSBackground:
    """Everything the MATLAB script computed (workspace variables)."""

    mm_per_camera_pixel_x: float
    mm_per_camera_pixel_y: float
    scale_difference_percent: float
    mm_per_camera_pixel: float
    nominal_dot_diameter_mm: float
    nominal_dot_area_mm2: float
    poster_area_mm2: float
    number_of_dots: int
    x_mm: np.ndarray
    y_mm: np.ndarray
    dot_diameter_mm: np.ndarray
    dot_radius_mm: np.ndarray
    poster_width_px: int
    poster_height_px: int
    px_per_mm: float
    img: np.ndarray | None
    svg_file: Path | None
    png_file: Path | None


def svg_text(poster_width_mm, poster_height_mm, x_mm, y_mm, dot_radius_mm) -> str:
    """SVG document text, formatted exactly like the MATLAB ``fprintf`` calls."""
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" '
        'width="%.4fmm" height="%.4fmm" '
        'viewBox="0 0 %.4f %.4f">\n'
        % (poster_width_mm, poster_height_mm, poster_width_mm, poster_height_mm),
        '<rect x="0" y="0" width="%.4f" height="%.4f" fill="white"/>\n'
        % (poster_width_mm, poster_height_mm),
    ]
    lines.extend(
        '<circle cx="%.5f" cy="%.5f" r="%.5f" fill="black"/>\n' % (x, y, r)
        for x, y, r in zip(x_mm.tolist(), y_mm.tolist(), dot_radius_mm.tolist())
    )
    lines.append("</svg>\n")
    return "".join(lines)


def rasterize_dots(x_mm, y_mm, dot_radius_mm, poster_width_px, poster_height_px, px_per_mm):
    """Draw the dots into a white uint8 image (255 = white, 0 = black).

    Faithful to the MATLAB loop: pixel *centres* sit at MATLAB's 1-based integer
    coordinates, and a pixel is black when its centre is inside the circle
    (``(X-xc)^2 + (Y-yc)^2 <= r^2``).  Python index ``k`` is MATLAB index ``k+1``.
    """
    img = np.full((poster_height_px, poster_width_px), 255, dtype=np.uint8)
    xc_all = x_mm * px_per_mm
    yc_all = y_mm * px_per_mm
    r_all = dot_radius_mm * px_per_mm
    for x_center_px, y_center_px, radius_px in zip(xc_all.tolist(), yc_all.tolist(), r_all.tolist()):
        # Bounding box in MATLAB 1-based pixel coordinates
        x_min = max(1, math.floor(x_center_px - radius_px))
        x_max = min(poster_width_px, math.ceil(x_center_px + radius_px))
        y_min = max(1, math.floor(y_center_px - radius_px))
        y_max = min(poster_height_px, math.ceil(y_center_px + radius_px))
        if x_min > x_max or y_min > y_max:
            continue
        X = np.arange(x_min, x_max + 1, dtype=float)
        Y = np.arange(y_min, y_max + 1, dtype=float)
        mask = (X[None, :] - x_center_px) ** 2 + (Y[:, None] - y_center_px) ** 2 <= radius_px**2
        img[y_min - 1 : y_max, x_min - 1 : x_max][mask] = 0
    return img


def BOS_Background_Generation(
    camera_width_px: int = CAMERA_WIDTH_PX,
    camera_height_px: int = CAMERA_HEIGHT_PX,
    FOV_width_mm: float = FOV_WIDTH_MM,
    FOV_height_mm: float = FOV_HEIGHT_MM,
    target_dot_pixels: float = TARGET_DOT_PIXELS,
    poster_width_mm: float = POSTER_WIDTH_MM,
    poster_height_mm: float = POSTER_HEIGHT_MM,
    target_black_fraction: float = TARGET_BLACK_FRACTION,
    diameter_jitter: float = DIAMETER_JITTER,
    print_resolution_dpi: int = PRINT_RESOLUTION_DPI,
    random_seed: int = RANDOM_SEED,
    base_filename: str = BASE_FILENAME,
    output_folder: str | Path | None = OUTPUT_FOLDER,
    *,
    make_png: bool = True,
    plot: bool = True,
    verbose: bool = True,
) -> BOSBackground:
    """Generate the BOS background, write SVG + PNG, print the summary.

    Arguments mirror the MATLAB "USER INPUTS" block.  ``output_folder=None``
    skips writing files (the pattern and image are still returned);
    ``make_png=False`` skips the (large) raster entirely.  A relative
    ``output_folder`` (default ``"BOS Backgrounds"``) is resolved against the
    current working directory.
    """
    # rng(random_seed): bit-identical to MATLAB's twister stream
    rs = matlab_rng(random_seed)

    svg_file = png_file = None
    if output_folder is not None:
        output_folder = Path(output_folder)
        if not output_folder.is_absolute():
            output_folder = Path.cwd() / output_folder
        output_folder.mkdir(parents=True, exist_ok=True)
        svg_file = output_folder / f"{base_filename}.svg"
        png_file = output_folder / f"{base_filename}.png"

    # ---------------- CHECK CAMERA / FOV SCALE ----------------
    mm_per_camera_pixel_x = FOV_width_mm / camera_width_px
    mm_per_camera_pixel_y = FOV_height_mm / camera_height_px
    # For square camera pixels and no distortion these should be ~equal.
    mean_scale = (mm_per_camera_pixel_x + mm_per_camera_pixel_y) / 2
    scale_difference_percent = 100 * abs(mm_per_camera_pixel_x - mm_per_camera_pixel_y) / mean_scale
    if scale_difference_percent > 1:
        warnings.warn(
            "Horizontal and vertical physical scales differ by %.2f%%. Check FOV dimensions."
            % scale_difference_percent,
            stacklevel=2,
        )
    # Use average physical scale to define dot diameter.
    mm_per_camera_pixel = mean_scale

    # ---------------- CALCULATE PHYSICAL DOT SIZE ----------------
    nominal_dot_diameter_mm = target_dot_pixels * mm_per_camera_pixel
    nominal_dot_area_mm2 = math.pi * (nominal_dot_diameter_mm / 2) ** 2
    poster_area_mm2 = poster_width_mm * poster_height_mm

    # ---------------- CHECK POSTER SIZE ----------------
    if poster_width_mm < FOV_width_mm:
        warnings.warn("Poster width is smaller than camera FOV width.", stacklevel=2)
    if poster_height_mm < FOV_height_mm:
        warnings.warn("Poster height is smaller than camera FOV height.", stacklevel=2)

    # ---------------- NUMBER OF RANDOM DOTS ----------------
    # Overlapping random dots: f = 1 - exp(-N*A_dot/A_total)  =>  N = -ln(1-f)*A_total/A_dot
    number_of_dots = int(
        matlab_round(-math.log(1 - target_black_fraction) * poster_area_mm2 / nominal_dot_area_mm2)
    )

    # ---------------- RANDOM DOT LOCATIONS AND DIAMETERS ----------------
    x_mm = poster_width_mm * rs.random_sample(number_of_dots)
    y_mm = poster_height_mm * rs.random_sample(number_of_dots)
    diameter_factor = 1 + diameter_jitter * (2 * rs.random_sample(number_of_dots) - 1)
    dot_diameter_mm = nominal_dot_diameter_mm * diameter_factor
    dot_radius_mm = dot_diameter_mm / 2

    # ---------------- DISPLAY PARAMETERS ----------------
    if verbose:
        print()
        print(_RULE)
        print("                 BOS BACKGROUND PARAMETERS")
        print(_RULE)
        print("Camera frame                  : %d x %d pixels" % (camera_width_px, camera_height_px))
        print("Background FOV                : %.3f x %.3f mm" % (FOV_width_mm, FOV_height_mm))
        print("Horizontal physical scale     : %.5f mm/pixel" % mm_per_camera_pixel_x)
        print("Vertical physical scale       : %.5f mm/pixel" % mm_per_camera_pixel_y)
        print("Target dot diameter           : %.2f camera pixels" % target_dot_pixels)
        print("Nominal physical dot diameter : %.4f mm" % nominal_dot_diameter_mm)
        print("Target black coverage         : %.1f %%" % (target_black_fraction * 100))
        print("Number of dots                : %d" % number_of_dots)
        print("Poster size                   : %.1f x %.1f mm" % (poster_width_mm, poster_height_mm))
        print("PNG resolution                : %d dpi" % print_resolution_dpi)
        print(_RULE)

    # ---------------- SAVE SVG ----------------
    if svg_file is not None:
        svg_file.write_text(
            svg_text(poster_width_mm, poster_height_mm, x_mm, y_mm, dot_radius_mm), encoding="ascii"
        )

    # ---------------- CREATE PNG ----------------
    poster_width_px = int(matlab_round((poster_width_mm / 25.4) * print_resolution_dpi))
    poster_height_px = int(matlab_round((poster_height_mm / 25.4) * print_resolution_dpi))
    px_per_mm = print_resolution_dpi / 25.4

    img = None
    if make_png:
        img = rasterize_dots(x_mm, y_mm, dot_radius_mm, poster_width_px, poster_height_px, px_per_mm)
        if png_file is not None:
            from PIL import Image  # Pillow ships with matplotlib

            Image.fromarray(img, mode="L").save(png_file, format="PNG")

        if plot:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(facecolor="w")
            ax.imshow(img, cmap="gray", vmin=0, vmax=255, interpolation="antialiased")
            ax.set_axis_off()
            ax.set_title("BOS Random Dot Background | %.2f Camera Pixels per Dot" % target_dot_pixels)

    # ---------------- FINISHED ----------------
    if verbose and svg_file is not None:
        print()
        print("Files successfully created:\n")
        print("SVG:\n%s\n" % svg_file)
        if png_file is not None and make_png:
            print("PNG:\n%s\n" % png_file)
        print("PNG raster dimensions:")
        print("%d x %d pixels\n" % (poster_width_px, poster_height_px))
        print("IMPORTANT FOR PRINTING:")
        print("1. Print at 100% / ACTUAL SIZE.")
        print('2. Do NOT use "Fit to Page".')
        print("3. SVG is preferred for professional printing.")
        print("4. PNG should be interpreted at %d dpi." % print_resolution_dpi)
        print()

    return BOSBackground(
        mm_per_camera_pixel_x=mm_per_camera_pixel_x,
        mm_per_camera_pixel_y=mm_per_camera_pixel_y,
        scale_difference_percent=scale_difference_percent,
        mm_per_camera_pixel=mm_per_camera_pixel,
        nominal_dot_diameter_mm=nominal_dot_diameter_mm,
        nominal_dot_area_mm2=nominal_dot_area_mm2,
        poster_area_mm2=poster_area_mm2,
        number_of_dots=number_of_dots,
        x_mm=x_mm,
        y_mm=y_mm,
        dot_diameter_mm=dot_diameter_mm,
        dot_radius_mm=dot_radius_mm,
        poster_width_px=poster_width_px,
        poster_height_px=poster_height_px,
        px_per_mm=px_per_mm,
        img=img,
        svg_file=svg_file,
        png_file=png_file if make_png else None,
    )


def main(argv: list[str] | None = None) -> None:
    import argparse

    import matplotlib.pyplot as plt

    p = argparse.ArgumentParser(prog="python -m HADES.data_analysis.bos_background_generation")
    p.add_argument("--output-dir", type=Path, default=Path(OUTPUT_FOLDER),
                   help="output folder (default: ./BOS Backgrounds in the current directory)")  # fmt: skip
    p.add_argument("--no-plot", action="store_true", help="skip the preview figure")
    args = p.parse_args(argv)
    BOS_Background_Generation(output_folder=args.output_dir, plot=not args.no_plot)
    plt.show()


if __name__ == "__main__":
    main()
