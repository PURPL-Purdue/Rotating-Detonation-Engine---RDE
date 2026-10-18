"""HADES feed-system sizing: sonic nozzle/orifice sizing and blowdown orifice Cd calibration.

Modules
-------
sonic_nozzle_sizing
    Port of ``SonicNozzleSizingV1.m`` ... ``V5.m`` (V5 is the canonical workflow).
blowdown_orifice_cd
    Orifice discharge-coefficient calibration from a supply-tank blowdown test
    (consolidates ``blowdown_orifice_cd.py`` and ``blowdown_orifice_cd_chat.py``).
choked_flow
    Ideal-gas choked-flow relations shared by both.
"""
