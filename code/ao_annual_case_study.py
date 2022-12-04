"""
Aimpoint Optimization Annual Simulation Module

Runs the aimpoint optimization for a number of selected Hours. 
Creates a summary CSV file. 
"""


if __name__ == "__main__":
    
    import os
    
    import single_hour_for_annual
    import flux_model
    import inputs
    import field
    
    
    main_case_name = "radial_daggett_50MW_SP_Field"
    case_filename = "./../case_inputs/flat_50_ca_case.csv"
    filenames = inputs.readCaseFile(case_filename)
    settings = inputs.readSettingsFile(filenames["settings"])
    solar_field = field.Field(filenames,params=settings, use_sp_field = settings["use_sp_field"])
    
    weather_data = flux_model.ReadWeatherFile(filenames["weather_filename"])
    hour_id = []
    for i in range(len(weather_data["dni"])):
        if weather_data["dni"][i] >= 500:
            if weather_data["month"][i] == 2 and weather_data["day"][i] == 3:
                hour_id.append(i)
            elif weather_data["month"][i] == 5 and weather_data["day"][i] == 5:
                hour_id.append(i)
            elif weather_data["month"][i] == 8 and weather_data["day"][i] == 4:
                hour_id.append(i)
            elif weather_data["month"][i] == 11 and weather_data["day"][i] == 4:
                hour_id.append(i)
                    
    ofile = open(main_case_name +"_Halos_annual_summary.csv", 'w')
    ofile.write("case_name,month,day,hour,dni,obj_value,num_defocused_heliostats,solve_time,SP_obj_value,SP_num_defocus,SP_max_flux\n")
    ofile.close()
    
    for i in range(len(hour_id)):
        command = "python single_hour_for_annual.py "+str(hour_id[i])
        os.system(command)
    