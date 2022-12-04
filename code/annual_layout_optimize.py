# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 17:41:21 2020

@author: kliaqat
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

def layout_from_SP(case_name, case_filename):
    import sp_module
    filenames = inputs.readCaseFile(case_filename)
    weather_data = flux_model.ReadWeatherFile(filenames["weather_filename"])
    sp_field = sp_module.SP_Field(filenames)
    ##Starting layout from SP
    sp_field.get_sp_layouts( weather_data, case_name+"-layout.csv")
    
    
def getAnnualSummary(case_name,case_filename):
    hour_ids = getHourIDs(case_filename)
    ofile_1 = open("halos_"+case_name+"_summary.csv",'w')
    ofile_2 = open("utilization_"+case_name+"_summary.csv",'w')
    suffix = "_summary.csv"
    first = True
    for hour in hour_ids:
        f1 = open(case_name+str(hour)+suffix,'r')
        f2 = open(case_name+str(hour)+suffix,'r')
        l1 = f1.readlines()
        l2 = f2.readlines()
        if first:
            ofile_1.write(l1[0])
            ofile_2.write(l2[0])
            ofile_1.write(l1[1])
            ofile_2.write(l2[1])
            first = False
        else:
            if l1[0].startswith("case"):
                ofile_1.write(l1[1])
            else: 
                ofile_1.write(l1[0])
            ofile_2.write(l2[0])
        f1.close()
        f2.close()
    ofile_1.close()
    ofile_2.close()
    
def getAnnualSPSummary(sp_case_name,case_filename):
    hour_ids = getHourIDs(case_filename)
    ofile = open("./../outputs/sp_"+sp_case_name+"_summary.csv",'w')
    first = True
    for hour in hour_ids:
        f = open("./../outputs/"+sp_case_name+str(hour)+"_summary.csv",'r')
        l = f.readlines()
        if first:
            ofile.write(l[0])
            first = False
        else:
            ofile.write(l[0])
        f.close()
    ofile.close()
    
def runAnnualStudy(case_name, case_filename, sp_case_name):
    """
    Runs a simulated annual estimate of aimpoint strategies, using the 
    collection of time periods used as the default in SolarPILOT, and reporting
    results.  Results are written once for each hour.
    
    arguments: 
        String case_name -- identifier of case
        String case_filename -- path to case settings file
        String sp_case_name -- path to SolarPILOT case settings file; set to 
            None if not running SolarPILOT
    returns: 
        None
    """
    hour_ids = getHourIDs(case_filename)
    ofile = open(case_name+"_summary.csv",'w')
    ofile.write("case_name,obj_value,num_defocused_heliostats,solve_time\n")
    ofile.close()
    if sp_case_name != None:
        spfile = open("SolarPilot_" + sp_case_name +"_summary.csv",'w')
        spfile.write("case_name,obj_value,num_defocus,max_flux,obj_pre_heur, obj_post_heur, post_num_defocused\n")
        spfile.close()
    for i in range(len(hour_ids)):
        command = "python run_hourly_case.py "+str(hour_ids[i])+" "+case_name+" "+case_filename
        if sp_case_name != None: 
            command = command + " " + sp_case_name
            
        os.system(command)
        print(hour_ids[i], " Hour is Done!")

if __name__ == "__main__":
    import os
    
    case_name = "radial-250-daggett"
    sp_case_name = "radial-250-daggett-sp"
    case_filename = "./../case_inputs/radial_250_ca_case.csv"
    
    layout_from_SP(case_name, case_filename)
    runAnnualStudy(case_name, case_filename, sp_case_name)
    
        

        
    
    

