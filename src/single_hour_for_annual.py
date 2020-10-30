# -*- coding: utf-8 -*-
"""
Generate case study for all periods.
"""

import solve_aim_model

def runCaseAnnual(case_name, case_filename, hour_id = None, decomp = False, parallel = True,
            plots = True, print_outputs = True):
    if decomp:
        results = solve_aim_model.solveDecomposedModel(case_name, case_filename, hour_id, parallel)
    else: 
        results = solve_aim_model.solveModelDirectly(case_name, case_filename)
    if plots:
        results.plotOutputs(case_name)
    if print_outputs:
        results.printOutput(case_name)
    return results.obj_value,results.num_defocused_heliostats, results.solve_time



if __name__ == "__main__":
    
    import flux_model
    import inputs
    import sp_flux
    import field
    import sys
    
    args = sys.argv
    hour_id = int(args[1])
    
    main_case_name = "flat_daggett_250MW_SP_Field"
    case_filename = "./../case_inputs/flat_250_ca_case.csv"
    filenames = inputs.readCaseFile(case_filename)
    settings = inputs.readSettingsFile(filenames["settings"])
    solar_field = field.Field(filenames,params=settings, use_sp_field = settings["use_sp_field"])
    
    weather_data = flux_model.ReadWeatherFile(filenames["weather_filename"])
    
    #for i in range(2): 
    case_name = str(hour_id) + main_case_name
    print(case_name)
    obj_value,num_defocused_heliostats, solve_time = runCaseAnnual(case_name, case_filename, hour_id = hour_id, decomp = True, parallel = True)
    outputs = sp_flux.full_field_flux_CaseStudy(filenames,solar_field,weather_data,case_name, use_sp_field = True,hour_id = hour_id)
    ofile = open(main_case_name +"_Halos_annual_summary.csv", 'a')
    ofile.write(case_name+","+str(weather_data["month"][hour_id])+","+str(weather_data["day"][hour_id])+","+str(weather_data["hour"][hour_id])+","+
                str(weather_data["dni"][hour_id])+","+str(obj_value)+","+str(num_defocused_heliostats)+","+str(solve_time)+","+str(outputs["obj_value"])+","+str(outputs["num_defocus"])+","+str(outputs["max_flux"])+"\n")
    ofile.close()