"""Entry point for the RDE post-processing / plotting steps.

Python translation of ``RDE_Toolbox_2/post_processing/Post_Process_Main.m``.
"""

from .Setup_Graph_Domain import Setup_Graph_Domain


def Post_Process_Main(PostProc):
    """Run the post-processing steps in order.

    ``PostProc`` fields used so far: ``xMax``, ``yMax`` (domain extents).
    """
    ax = Setup_Graph_Domain(PostProc.xMax, PostProc.yMax)  # noqa: F841

    # Future stuff goes here
