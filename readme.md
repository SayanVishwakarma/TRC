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

--------------------------------------------------------------------------------

REQUIREMENTS
------------
- Python 3.8 or later
- Numpy

Install missing dependencies with:
> pip install numpy

--------------------------------------------------------------------------------

NOTES
-----
- The code supports both manual and automatic staging.
- Boostback modeling assumes a delta-v penalty and adjusts the optimization accordingly.
- Thrust and tank parameters must be physically consistent for valid results.

--------------------------------------------------------------------------------

CONTACT
-------
Developed by [Your Name / Team Name]  
[Institution or Organization]  
Email: [your.email@example.com]

--------------------------------------------------------------------------------

LICENSE
-------
This code is provided under the MIT License. You may modify and distribute it
with appropriate attribution.
