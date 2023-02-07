# HALOS - Heliostat Aimpoint and Layout Optimization Software 

## Introduction
HALOS is an open-source software decision tool implemented in Python that uses mixed-integer programming models to determine the best aimpoint strategy for the solar collection field of a concentrating solar power (CSP) central receiver plant. Design and operations decisions addressed by this tool include: (i) the location of each heliostat in the solar field, and (ii) the intended aimpoint of each heliostat to the receiver for each hour, across a representative collection of days. Given weather and location data, heliostat specifications, a solar field layout, and a receiver’s size, location, and geometry as input, HALOS outputs an aiming strategy that maximizes thermal power delivery to the receiver while ensuring that the thermal flux profile of the receiver falls within design specifications.  A key feature of the tool includes a method to subdivide the solar field into subfields, for which the aiming strategy can be optimized separately and in parallel.  This allows for aiming strategies to be obtainable using integer programming methods for commercial-scale plants within a matter of minutes.

The tool includes a module that directly interfaces with SolarPILOT, an NREL-developed solar field performance characterization tool, to obtain high-fidelity flux maps and solar field layouts.  HALOS also accepts flux images in flat-file format.

The software development is carried out at [National Renewable Energy Laboratory (NREL)](https://www.nrel.gov/) and was previously funded by the U.S. Department of Energy under the award number DE-EE00035930. 

<!---## Access
Access to the repository is currently limited to project development team and to gain access contact [Alex Zolan](mailto://alexander.zolan@nrel.gov). --->

## SolarPILOT Integration/Support
SolarPILOT is a design and optimization tool for concentrating solar power (CSP) central receiver plant developed by NREL and is available open source. HALOS provides access to SolarPILOT through python API and SolarPILOT can be used to populate flux/field into HALOS optimization model.
* The [SolarPILOT Python API](https://github.com/NREL/SolarPILOT/tree/copilot/deploy/api) is currently available within the [copilot](https://github.com/NREL/SolarPILOT/tree/copilot) branch of [SolarPILOT](https://github.com/NREL/SolarPILOT). 
* To integrate SolarPILOT into HALOS, [solarpilot.dll](https://github.com/NREL/SolarPILOT/blob/copilot/deploy/api/solarpilot.dll) and [copylot.py](https://github.com/NREL/SolarPILOT/blob/copilot/deploy/api/copylot.py) have to be in the HALOS directory at  ".\HALOS\code".  Within HALOS, [sp_module.py](https://github.com/NREL/HALOS/blob/master/code/sp_module.py) can be used to interact with SolarPilot.
* Microsoft Visual Studio 2019 or above is required for SolarPILOT API to work. 

## Solver Binary Setup
HALOS is implmemented in Python 3.7 and uses Pyomo, an algebraic modeling language implemented in Python for its optimization model.  While a large collection of solvers are accessible to Pyomo, the default solver used in HALOS is Cbc.  The following steps will allow Cbc to be accessible within HALOS on Windows systems:
* Download the Cbc binary from its [location at COIN-OR](https://www.coin-or.org/download/binary/Cbc/).
* Add the folder location of the binary to the PATH environment variable
After adding this folder, placing binaries of other solvers that work with Pyomo in the same folder should allow Pyomo models to access them.  HALOS includes some settings for Cbc, GLPK, and CPLEX in the source code.  

## Running Case Studies 

Once SolarPILOT is setup to work, follow the following steps to run different cases.
1. Create CSV files for inputs 
    * Create a main CSV file that has the links to specific CSV input files, here's an [Example](https://github.com/NREL/HALOS/blob/master/case_inputs/radial_50_ca_case.csv). 
    * Now create CSV files for receiver, settings and provide the path to weather file and field file (if using existing field). 
    * Examples: [Receiver input CSV](https://github.com/NREL/HALOS/blob/master/case_inputs/radial_dagget_50/extcyl_50MW.csv), [Settings](https://github.com/NREL/HALOS/blob/master/case_inputs/radial_dagget_50/case_settings.csv)
2. If running for single hour with DNI provided in receiver file, run [main.py](https://github.com/NREL/HALOS/blob/master/code/main.py). For annual simulation of aimpoint strategy use [ao_annual_case_study.py](https://github.com/NREL/HALOS/blob/master/code/ao_annual_case_study.py). 
3. Follow step 2 for running HALOS only, SolarPILOT only or both. 
4. [sp_module.py](https://github.com/NREL/HALOS/blob/master/code/sp_module.py) can be used to obtain flux maps and field from SolarPILOT using the Python API. 

## Tool Organization

The software's module files are described as follows: 

### [main.py](https://github.com/NREL/HALOS/blob/master/code/main.py)
This is the main aimpoint optimization module and is used to run a single case simulation. User can decide whether to decompose the problem or not. 

### [annual_layout_optimize.py](https://github.com/NREL/HALOS/blob/master/code/annual_layout_optimize.py)
Runs the layout optimization for a number of hours. 

### [ao_annual_case_study.py](https://github.com/NREL/HALOS/blob/master/code/ao_annual_case_study.py)
Runs the aimpoint optimization for a number of selected Hours. Creates a summary CSV file.

### [case_study_flux_tests.py](https://github.com/NREL/HALOS/blob/master/code/case_study_flux_tests.py)
This module is used to investigate performance and accuracy of image-shift approximation annually. 

### [field.py](https://github.com/NREL/HALOS/blob/master/code/field.py)
Reads field from file or uses SolarPilot to generate field. Divides field in sections using angle or distance.

### [flux_method.py](https://github.com/NREL/HALOS/blob/master/code/flux_method.py)
HALOS' flux calculator, uses Hermite expansion for flux characterization. 

### [flux_model.py](https://github.com/NREL/HALOS/blob/master/code/flux_model.py)
This module is used for creating models for flux modeling and outputting images to either a data file or an optimization model. It can also use SolarPilot to generate Flux images in parallel.

### [geometry.py](https://github.com/NREL/HALOS/blob/master/code/geometry.py)
Builds different type of receiver geometries and generates coordinates of aimpoints and measurement points. 

### [heliostat.py](https://github.com/NREL/HALOS/blob/master/code/heliostat.py)
Contains methods for modelling heliostats. 

### [hermite_calculation.py](https://github.com/NREL/HALOS/blob/master/code/hermite_calculation.py)
Calculates Hermite expansion for flux characterization.

### [inputs.py](https://github.com/NREL/HALOS/blob/master/code/inputs.py)
Inputs reading module. Creates receiver and flux model from input csv data. 

### [mirror_model.py](https://github.com/NREL/HALOS/blob/master/code/mirror_model.py)
Mirror and slope error characterization model. 

### [opt_heuristic.py](https://github.com/NREL/HALOS/blob/master/code/opt_heuristic.py)
This module includes heuristics that can generate an initial feasible solution to the aimpoint strategy optimization model.

### [optimize_aimpoint.py](https://github.com/NREL/HALOS/blob/master/code/optimize_aimpoint.py)
This module contains the aimpoint optimization model, which is implemented in the Pyomo modeling language.

### [plotting.py](https://github.com/NREL/HALOS/blob/master/code/plotting.py)
Contains code for creating different types of plots. 

### [process_aimpoint_outputs.py](https://github.com/NREL/HALOS/blob/master/code/process_aimpoint_outputs.py)
Serves as aimpoint optimization model's output processing module.

### [run_hourly_case.py](https://github.com/NREL/HALOS/blob/master/code/run_hourly_case.py)
Runs a single hour simulation for layout optimization. 

### [single_hour_for_annual.py](https://github.com/NREL/HALOS/blob/master/code/single_hour_for_annual.py)
Runs a single hour at a time for annual aimpoint optimization simulation. 

### [sol_pos.py](https://github.com/NREL/HALOS/blob/master/code/sol_pos.py)
Solar Positioning module. Calculates solar azimuth and zenith. 

### [solve_aim_model.py](https://github.com/NREL/HALOS/blob/master/code/solve_aim_model.py)
Solves aimpoint optimization model, either directly or using decomposition approach. 

### [sp_aimpoint_heuristic.py](https://github.com/NREL/HALOS/blob/master/code/sp_aimpoint_heuristic.py)
This module contains aimpoint heuristics using SolarPILOT flux calculations. Removes heliostats to overcome flux violation at the receiver for SolarPilot. 

### [sp_module.py](https://github.com/NREL/HALOS/blob/master/code/sp_module.py)
HALOS and SolarPilot (python API) integration module. Interacts with Solarpilot API and generates field and flux images for different cases.

### [sun_shape.py](https://github.com/NREL/HALOS/blob/master/code/sun_shape.py)
This module contains source code for building sun shape objects and their related error models.

### [Toolbox.py](https://github.com/NREL/HALOS/blob/master/code/Toolbox.py) 
Classes for points, vectors, and their operations. 

### [WELL512.py](https://github.com/NREL/HALOS/blob/master/code/WELL512.py)
Random number generator. An implementation of the WELL512a RNG created by L'Ecuyer and Matsumoto, originally programmed in C.
 

## [Inputs](https://github.com/NREL/HALOS/tree/master/case_inputs)

This section defines the inputs for the optimization model and SolarPILOT Python API.

## 1. Main Case File
The main input CSV file is used to assign paths to CSV files and the parameters are defined as:

| Parameter | Data Type | Description | 
| --- | --- | --- | 
| field_filename | string | Path to Solar Field if using field from file |
| receiver_filename | string | Path to Receiver settings CSV |
| mirror_filename | string | Path to mirror settings CSV |
| settings | string | Path to case settings CSV |
| weather_filename | string | Path to Weather File such as (TMY3) |
| heliostat_file_dir | string | (OPTIONAL) Path to directory containing heliostat-specific flux images in flat-file format |
| flux_limit_filename | string | (OPTIONAL) Path to receiver flux limits table CSV |
| rec_obj_filename | string | (OPTIONAL) Path to objective value table CSV |

Note that the tables in "flux_limit_filename" and "rec_obj_filename" do not include any headers, and must match the dimensions of the receiver given by "pts_per_ht_dim" and "pts_per_len_dim" for rows and columns, respectively.  The table in flux_limit_filename should provide inputs in kW, while the table in "rec_obj_filename" should consist of unitless multipliers, generally 1 for measurement points on the receiver surface and 0 for measurement points on the heat shields, if used.  If "rec_obj_filename" is not provided, a default of 1 wil be applied to every measurement point.  

## 2. Case Settings
The case_setting file contains the following parameters. 

| Parameter | Data Type | Description  | Default Value (if any) | 
| --- | --- | --- | --- |
| mirror_model | String | Select mirror model  | SinglePointGaussian |
| method | String | Flux Calculation method for HALOS | SimpleNormalFluxCalc |
| num_sections | Integer | Number of field subsections | 8 |
| section_method | String | Method for field division into subsections (angle/distance) | angle |
| use_sp_flux | Integer | Boolean (1/0): Use SolarPILOT to calculate Flux | 1 |
| use_sp_field | Integer | Boolean (1/0): Use SolarPILOT to generate field | 1 |
| hour_idx | Integer | Simulation hour index from weather file | 4311 |
| heliostat_group_size | Integer | Heliostat group size to focus on same aimpoint | 4 |

## 3. Receiver Settings
The receiver input csv has the following parameters  

| Parameter | Data Type | Description  | Default Value (if any) | 
| --- | --- | --- | --- |
| tow_height | Floating-point number | Tower Height  | 150 |
| length | Floating-point number | Receiver horizontal Length | 21 |
| height | Floating-point number | Receiver height | 17 |
| diameter | Floating-point number | Receiver diameter | 10.38 |
| pts_per_ht_dim | Integer | Number of rows in receiver measurement point grid | 20 |
| pts_per_len_dim | Integer | Number of columns in receiver measurement point grid | 20 |
| zenith_deg | Floating-point number | Receiver zenith angle - degrees | 90 |
| azimuth_deg | Floating-point number | Receiver azimuth angle - degrees | 180 |
| rec_cent_offset_x | Integer | Receiver offset from center - x_axis | 0 |
| rec_cent_offset_y | Integer | Receiver offset from center - y_axis | 0 |
| rec_cent_offset_z | Integer | Receiver offset from center - z_axis | 0 |
| aim_rows | Integer | Number of Aimpoint rows or number of aimpoints per column | 7 |
| aim_cols | Integer | Number of Aimpoint columns or number of aimpoints per row | 7 |
| power_rating | Floating-point number | Solar field design power - Receiver power rating (MWt) | 5.00E+07 |
| receiver_type | String | Receiver geometry (Flat plate/External cylindrical) | External cylindrical |
| field_dni | Floating-point number | Design point DNI for field generation (W/m^2) | 950 |
| flux_dni | Floating-point number | Design point DNI for flux calculation (W/m^2) | 950 |
| mirror_len | Floating-point number | Length of Heliostat (m) | 12.2 |
| mirror_ht | Floating-point number | height of Heliostat (m) | 12.2 |
| mirror_area | Floating-point number | Area of Heliostat (m^2) | 148.84|
| flux_ub | Floating-point number | Maximum allowable flux (kW/m^2) | 1000 |
| flux_lb | Floating-point number | Flux lower bound (kW/m^2) | 0 |
| n_circulation | Integer | Number of flow circulations (Assumes flow enters from top ) | 5 |
| min_col_fraction | Float | Ratio of flux from lowest to highest column on receiver | 0.0 |

Note that "flux_ub," "flux_lb," and "n_circulation" are overridden if a filepath is provided to a user-specific flux limit file via "flux_limit_filename" in the Case Settings file.

## External libraries

Open source packages used in HALOS:

| Project | Version | Usage |
|---------|---------|-------|
| [Pysolar](https://github.com/pingswept/pysolar)  	| v0.8 	| Pysolar is an open source module that provides code for calculating solar irradiance and solar position based on time and location (Solar angles).|
| [Pyomo](http://www.pyomo.org/)  	| v5.7	| An open-source optimization modeling language for linear, non-linear, quadratic, and mixed integer programming. |
| [GLPK](https://www.gnu.org/software/glpk/)  	| v4.65	| A solver for large-scale linear,mixed integer, and related programming problems. |


