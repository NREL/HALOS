# -*- coding: utf-8 -*-
"""
Generate case study for all periods.
"""

import solve_aim_model

def runCaseAnnual(case_name, case_filename, hour_id = None, decomp = False, parallel = True,
            plots = True, print_outputs = True):
    """
    Run Aimpoint optimization model for any period

    Returns
    -------
    obj_value
        Power delivered to receiver 
    num_defocused_heliostats
        Number of Heliostats defocused 
    solve_time
        Solution time 

    """
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
    import sp_module
    import field
    import sys
    
    args = sys.argv
    hour_id = int(args[1])

    
    main_case_name = "radial_daggett_50MW"
    case_filename = "./../case_inputs/radial_50_ca_case.csv"
    filenames = inputs.readCaseFile(case_filename)
    settings = inputs.readSettingsFile(filenames["settings"])
    solar_field = field.Field(filenames,params=settings, use_sp_field = settings["use_sp_field"])
    
    weather_data = flux_model.ReadWeatherFile(filenames["weather_filename"])
    sp_flux = sp_module.SP_Flux(filenames, field = solar_field)
    #for i in range(2): 
    case_name = str(hour_id) + main_case_name
    print(case_name)
    obj_value,num_defocused_heliostats, solve_time = runCaseAnnual(case_name, case_filename, hour_id = hour_id, decomp = True, parallel = True)
    outputs = sp_flux.run_sp_case(weather_data, case_name, hour_id = hour_id, sp_aimpoint_heur=True)
    ofile = open(main_case_name +"_Halos_annual_summary.csv", 'a')
    ofile.write(case_name+","+str(weather_data["month"][hour_id])+","+str(weather_data["day"][hour_id])+","+str(weather_data["hour"][hour_id])+","+
                str(weather_data["dni"][hour_id])+","+str(obj_value)+","+str(num_defocused_heliostats)+","+str(solve_time)+","+str(outputs["obj_value"])
                +","+str(outputs["num_defocus"])+","+str(outputs["max_flux"])+","+str(outputs["before_flux"])+","+str(outputs["post_heur_obj"])+","+str(outputs["post_num_defocus"])+"\n")
    ofile.close()