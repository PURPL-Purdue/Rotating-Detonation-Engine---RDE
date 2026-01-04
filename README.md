Hello, This is the official RDE GITHUB Repo. ALL USERS MUST FOLLOW THE NAMING CONVENTION SPECIFIED BELOW. NOTE: THIS CONVENTION IS SUBJECT TO CHANGE


FOR HADES - The following Subteams are given as 
HADES_<subteam>_<category>_<descriptor>_v<number>.

| Subteam                   | Code    |
| ------------------------- | ------- |
| RDE Sizing Model          | `size`  |
| CAD                       | `cad`   |
| CFD                       | `cfd`   |
| FEA                       | `fea`   |
| Integration               | `int`   |
| Testing / Instrumentation | `test`  |
| Pre-Detonator             | `pd`    |
| Logistics / Systems       | `sys`   |
| Manufacturing             | `mfg`   |
| Media                     | `media` |

EXAMPLES
1. HADES_size_thermal_heat_flux_v1.m
2. HADES_size_thermal_heat_flux_v2.m
3. HADES_cfd_solver_axial_v3.m
4. HADES_fea_struct_wall_stress_v1.m
5. HADES_pd_flow_ignition_delay_v1.m
6. HADES_test_val_hotfire_tc_v2.m

FOR DEVIL - The following Subteams are given as 
DEVIL_<subteam>_<category>_<descriptor>_v<number>.m

| Subteam             | Code    |
| ------------------- | ------- |
| Sizing              | `size`  |
| CAD                 | `cad`   |
| CFD                 | `cfd`   |
| FEA                 | `fea`   |
| Integration         | `int`   |
| Testing             | `test`  |
| Pre-Detonator       | `pd`    |
| Systems / Logistics | `sys`   |
| Manufacturing       | `mfg`   |
| Media               | `media` |

EXAMPLES
1. DEVIL_size_thermal_heat_flux_v1.m
2. DEVIL_size_thermal_heat_flux_v2.m
3. DEVIL_cfd_solver_axial_v3.m
4. DEVIL_fea_struct_wall_stress_v1.m
5. DEVIL_pd_flow_ignition_delay_v1.m
6. DEVIL_test_val_hotfire_tc_v2.m
7. DEVIL_int_perf_mass_budget_v1.m
8. DEVIL_mfg_util_tolerance_stackup_v1.m

---------------------------------------------------------------------------------------------------------------------------------------------------------
FOR GENERAL BOOKEEPING DO NOT NAME VERSIONS as FINAL until the code is used in CDR. 
for the <category> tag, please make sure to be consistent and appropriate. HERE are some general catagory tags that you can use: 
| Category     | Code      |
| ------------ | --------- |
| Geometry     | `geom`    |
| Thermal      | `thermal` |
| Structural   | `struct`  |
| Flow         | `flow`    |
| Performance  | `perf`    |
| Optimization | `opt`     |
| Solver       | `solver`  |
| Validation   | `val`     |
| Plotting     | `plot`    |
| Utilities    | `util`    |
| Run Script   | `run`     |

