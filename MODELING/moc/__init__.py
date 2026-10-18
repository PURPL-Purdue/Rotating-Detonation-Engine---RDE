"""Method-of-characteristics (MoC) nozzle solver and RDE MoC model.

Python translation of the GPL-3.0 project *Nozzle_MOC_supersonic* by Xavier
Dechamps (https://github.com/xavierdechamps/Nozzle_MOC_supersonic), 2-D steady
irrotational MoC after Zucrow & Hoffman, *Gas Dynamics* Vol. II (1977), as
modified by the PURPL RDE team, plus the team's RDE MoC model
(``MODELING.moc.rde``).  Licensed GPL-3.0, see ``LICENSE`` in this folder.

One module per original MATLAB function, with the same name, argument order
and outputs, e.g.::

    from MODELING.moc.MOC_2D_steady_irrotational_main import default_inputs, run
    geom, params, plots = default_inputs()
    sol = run(geom, params, plots=None)        # no plotting
    sol.X, sol.Y, sol.U, sol.V                 # characteristic-net nodes

Index-valued outputs (``i_tri``, ``j_tri``, ``i_quad``, ``j_quad``) keep the
MATLAB 1-based values; subtract 1 to index the NumPy node arrays.
"""
