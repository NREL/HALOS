# -*- coding: utf-8 -*-
"""
Layout optimization methods
"""
import inputs
import solve_aim_model
import flux_model


def hourlyHalosRun(case_name, case_filename, h_id):
    #TODO: include an optional layout filename or full year run.
    solve_aim_model.runHourlyCase(case_name,case_name+str(h_id), case_filename, 
        hour_id = h_id, decomp = True, parallel = True,
        plots = False, print_outputs = True, append = True
        )
        # append = True


def hourlySPRun(case_name, case_filename, h_id):
    filenames = inputs.readCaseFile(case_filename)
    weather_data = flux_model.ReadWeatherFile(filenames["weather_filename"])
    solve_aim_model.runSPHourlyCase(case_name,case_name+str(h_id), case_filename, hour_id = h_id, 
                                    append = True, weather_data=weather_data)
    
def singleHourlyRun(case_name, case_filename,hour_id, sp_case_name=None):

    ##Perform a single-hour run.  Need to adjust to use layout file if possible
    #TODO: look into memory issue when running full-year analysis
    hourlyHalosRun(case_name, case_filename, hour_id)
    #if comparing to solarpilot, include this.
    if sp_case_name != None:
        hourlySPRun(sp_case_name, case_filename, hour_id)
    #add in method that analyzes utilization by section for all hours
    #add in method of highest-potential helios from each potential section
    #redo full year run
    
    
    
    
if __name__ == "__main__":
    import sys
    
    args = sys.argv
    hour_id = int(args[1])
    
    case_name = args[2]
    case_filename = args[3]
    if len(args) > 4: 
        sp_case_name = args[4]
    else:
        sp_case_name = None
    
    singleHourlyRun(case_name, case_filename,hour_id, sp_case_name=sp_case_name)
    