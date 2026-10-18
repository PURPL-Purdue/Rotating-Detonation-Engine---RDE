"""Empty plotting domain for the RDE graph.

Python translation of ``RDE_Toolbox_2/post_processing/Setup_Graph_Domain.m``.
"""

import matplotlib.pyplot as plt


def Setup_Graph_Domain(xMax, yMax):
    """Create the 'RDE Flow Field' figure and return its axes (x in [0, xMax], y in [0, yMax])."""
    fig = plt.figure("RDE Flow Field")
    ax = fig.add_subplot(111)

    ax.set_xlim(0, xMax)
    ax.set_ylim(0, yMax)

    # Lock the plot box shape to the data ratio (xMax:yMax) so one unit in x and
    # one unit in y have the same length (MATLAB pbaspect, without changing limits).
    # Characteristic/wave angles only look correct under 1:1 scaling.
    ax.set_box_aspect(yMax / xMax)

    ax.set_xlabel("x (mm)")  # TODO: set real units once decided (m, mm, deg...)
    ax.set_ylabel("y (mm)")
    ax.set_title("RDE Domain")

    ax.grid(True)  # later functions add to this axes (MATLAB 'hold on')
    return ax
