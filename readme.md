================================================================================
                        MULTI-STAGE ROCKET SIMULATION & OPTIMIZATION TOOL
================================================================================

DESCRIPTION
-----------
This software provides a modular Python-based tool for simulating, sizing, and 
optimizing multi-stage orbital launch vehicles. It supports customized staging, 
thrust-to-weight modeling, delta-v breakdown (including boostback maneuvers), 
and tank sizing based on propellant densities and engine performance.

The tool reads input parameters from a structured JSON file and generates 
detailed output in a human-readable format. It is particularly suited for 
preliminary mission design and performance analysis of staged rockets.

--------------------------------------------------------------------------------

FILES & DIRECTORY STRUCTURE
---------------------------
The project is organized as follows:

- evaluate.py              : Main driver script to run simulations
- delta_v.py               : Module for computing required delta-v (with/without boostback)
- stage_optimizer.py       : Module for staging optimization and mass allocation
- rocket.py                : Rocket model including tank sizing and simulation
- __init__.py              : Common imports and shared functions
- input files/             : Contains JSON files with rocket configurations
- output files/            : Stores the output text files from simulation runs

--------------------------------------------------------------------------------

INPUT FILE FORMAT
-----------------
The input JSON file (e.g., falcon 9.json) should include the following fields:

{
  "target orbit"         : Target orbital altitude in kilometers (e.g., 300),
  "average isp"          : Average specific impulse (used for estimating total delta-v),
  "thrust by weight"     : List of thrust-to-weight ratios per stage (e.g., [1.5, 0.9]),
  "number of stages"     : Integer count of rocket stages,
  "isps"                 : List of ISPs (specific impulse) per stage (e.g., [311, 348]),
  "structural factors"   : List of structural mass fractions per stage,
  "payload"              : Payload mass in kilograms,
  "Boostback"            : "True" or "False" (string) to enable or disable boostback,
  "step diameters"       : List of rocket stage diameters in meters,
  "fuel density"         : Density of fuel (kg/m³),
  "oxidiser density"     : Density of oxidizer (kg/m³),
  "tank diameters"       : Outer tank diameter in meters,
  "o/f"                  : Oxidizer-to-fuel ratio,
  "engine height"        : Height of each engine in meters,
  "engine thrust"        : List of thrust values for each engine (N)
}

--------------------------------------------------------------------------------

USAGE INSTRUCTIONS
------------------

1. Place your JSON configuration file in the `input files/` directory.
   Example: input files/falcon 9.json

2. Open `evaluate.py` and set the rocket name at the top of the file:
   rocket_name = "falcon 9"

3. Run the script using the following command:
   > python evaluate.py

4. Output will be saved automatically in the `output files/` directory
   with the filename matching the rocket name (e.g., falcon 9.txt).

--------------------------------------------------------------------------------

OUTPUT
------
The output text file contains:

- Delta-v breakdown (including boostback if applicable)
- Mass optimization results per stage
- Thrust-to-weight values per stage
- Tabulated stage-by-stage breakdown
- Rocket tank sizing (based on densities and O/F ratio)



NOTES TO SELF
- TRC HEAVY WITH TANK DIA SWEEP FILES/7.txt OR 4.txt will propably work


**********************************************************************************TRC Heavy tank dia
sweep**********************************************************************************
---Input Parameters---
rocket name : TRC Heavy sweep
target orbit : 500
average isp : 360
thrust by weight : [1.9, 1.4]
number of stages : 2
isps : [350, 380]
structural factors : [0.0505 0.1746]
modify factors : True
payload : 70000.0
Boostback : [True, True]
step diameters : [8.0, 8.0]
fuel density : 820
oxidiser density : 1141
tank diameters : 7.5
o/f : 2.6
engine height : 3.1
engine thrust : [2800000.0, 3285000.0]

###################################   DELTA-V BREAKDOWN   ##############################

|                                     |   Delta-v (m/s) |
|:------------------------------------|----------------:|
| Required Delta-v (m/s)              |        12908.1  |
| Initial Velocity (m/s)              |          405.22 |
| Orbital Velocity (m/s)              |         7901.32 |
| Gravitational Loss (m/s)            |         1862    |
| Drag Loss (m/s)                     |           50    |
| Boostback Delta-v for stage 1 (m/s) |         1500    |
| Boostback Delta-v for stage 2 (m/s) |         2000    |

################################# ROCKET MASS & Delta-v BREAKDOWN ###################################

=====================================================================================
|                              |                1 |      2 |
|:-----------------------------|-----------------:|-------:|
| Step Total Mass (kg)         |      3.61887e+06 | 164561 |
| Step Structural Mass (kg)    | 182752           |  28732 |
| Step Propellant Mass (kg)    |      3.43611e+06 | 135829 |
| Fuel Mass for Ascent (kg)    |      3.33586e+06 | 115401 |
| Fuel Mass for Boostback (kg) | 100248           |  20427 |
| Delta-v Ascent (m/s)         |   6886           |   2522 |
| Delta-v Boostback (m/s)      |   1499           |   1842 |
| Stage Total Mass (kg)        |      3.78343e+06 | 164561 |

===================================================================================
Delta-v Achieved (excluding boostback): 9408.1 m/s
Boostback Delta-v Achieved for stage 1: 1500.0 m/s
Boostback Delta-v Achieved for stage 2: 2000.0 m/s
Total Delta-v Achieved: 12908.1 m/s
===================================================================================

###################################   ROCKET DIMENSIONS   ##############################

|    | Stage   |   Fuel Mass (kg) |   Oxidiser Mass (kg) |   Fuel Tank Height (m) |   Oxidiser Tank Height (m) |   Intertank Height (m) |   Engine Count |   Total Stage Height (m) |
|---:|:--------|-----------------:|---------------------:|-----------------------:|---------------------------:|-----------------------:|---------------:|-------------------------:|
|  0 | Stage 1 |         954476   |          2.48164e+06 |                  28.85 |                      51.73 |                      0 |             26 |                    83.68 |
|  1 | Stage 2 |          37730.3 |      98098.8         |                   3.54 |                       4.45 |                      0 |              1 |                    11.09 |

Total Rocket Height: 96.77 m
Note that the total rocket height does not include the height of the payload fairing
Updated structural factors are [0.0505 0.1746]
********************************************************************************************************************************

$$$$$$$$$ THIS CAN ALSO WORK