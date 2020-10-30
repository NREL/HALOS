# -*- coding: utf-8 -*-
"""
Layout optimization methods
"""
import inputs
import solve_aim_model
import flux_model

def getHourIDs(case_filename):
    filenames = inputs.readCaseFile(case_filename)
    weather_data = flux_model.ReadWeatherFile(filenames["weather_filename"])
    hour_ids = []
    for i in range(len(weather_data["dni"])):
        if (weather_data["month"][i], weather_data["day"][i]) in [(2,3), (5,5), (8,4), (11,4)] and weather_data["dni"][i] >= 500:
            hour_ids.append(i)
    return hour_ids


def fullYearRun(case_name, case_filename):
    #TODO: include an optional layout filename or full year run.
    hour_ids = getHourIDs(case_filename)
    append = False
    for h_id in hour_ids:
        solve_aim_model.runHourlyCase(case_name+str(h_id), case_filename, 
            hour_id = h_id, decomp = True, parallel = True,
            plots = False, print_outputs = True, append = append
            )
        append = True


def fullYearSPRun(case_name, case_filename):
    filenames = inputs.readCaseFile(case_filename)
    weather_data = flux_model.ReadWeatherFile(filenames["weather_filename"])
    hour_ids = getHourIDs(case_filename)
    append = False
    for h_id in hour_ids:
        solve_aim_model.runSPHourlyCase(case_name+str(h_id), case_filename, hour_id = h_id, append = append, weather_data=weather_data)
        append = True
        
        
def optimizeLayout(case_name, case_filename, sp_case_name=None):
    import sp_module
    filenames = inputs.readCaseFile(case_filename)
    weather_data = flux_model.ReadWeatherFile(filenames["weather_filename"])
    sp_field = sp_module.SP_Field(filenames)
    ##Starting layout from SP
    sp_field.get_sp_layouts( weather_data, case_name+"-layout.csv")
    ##Perform a full-year run.  Need to adjust to use layout file if possible
    #TODO: look into memory issue when running full-year analysis
    fullYearRun(case_name, case_filename)
    #if comparing to solarpilot, include this.
    if sp_case_name != None:
        fullYearSPRun(sp_case_name, case_filename)
    #add in method that analyzes utilization by section for all hours
    #add in method of highest-potential helios from each potential section
    #redo full year run
    
    
    
    
if __name__ == "__main__":
    case_name = "radial-250-daggett"
    sp_case_name = "radial-250-daggett-sp"
    case_filename = "./../case_inputs/radial_250_ca_case.csv"
    hour_ids = getHourIDs(case_filename)
    print(hour_ids)
    optimizeLayout(case_name, case_filename, sp_case_name=None)
    