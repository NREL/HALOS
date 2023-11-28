#!/usr/bin/env python
# coding: utf-8

# In[1]:


#min_col_fraction = 0
# -*- coding: utf-8 -*-
"""
HALOS main module

This module creates a flux_model object according to user inputs, then creates
an instance of the aimpoint strategy optimization model.

"""
#Edit solve_aim_model.runHourlyCase to return optimization results in order to store aimpoint info
import solve_aim_model
import inputs
import flux_model
import sp_module
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    case_filename = "./../case_inputs/flat_1mw_ca_case.csv"
    case_name = "flat_daggett_1mw"
    filenames = inputs.readCaseFile(case_filename)
    settings = inputs.readSettingsFile(filenames["settings"])
    results = solve_aim_model.runHourlyCase(case_name, case_name, case_filename, hour_id=settings["hour_idx"], decomp = True, parallel = True) #Save results of optimization model
    # results.outputAimpointsToFile(case_name)
    #solve_aim_model.runSPHourlyCase(case_name, case_name, case_filename, hour_id=settings["hour_idx"],weather_data = filenames["weather_filename"],sp_aimpoint_heur = True)
    sp_flux = sp_module.SP_Flux(filenames,use_sp_field = settings["use_sp_field"]) #Create SolarPILOT flux model

    #Input HALOS aimpoints into SolarPILOT, plots SP fluxmap and saves values to flat file
    array = sp_flux.sp_flux_with_opt_aim(case_filename, results, case_name=case_name, hour_id=settings["hour_idx"]) #Get SolarPILOT results with HALOS aimpoints


    path = "./../outputs/" #Insert path to HALOS outputs folder

    halos_flux = np.loadtxt(path + "flat_daggett_1mw_fluxmap_values.csv", delimiter = ',') #Load halos fluxmap values
    sp_flux = np.loadtxt(path + 'flat_daggett_1mw_opt_sp_flux_values.csv', delimiter = ',') #Load SolarPILOT fluxmap values
    diff = sp_flux-halos_flux #Difference between HALOS and SolarPILOT fluxmaps

    #Plot fluxmap of difference between SP and HALOS
    plt.imshow(diff, cmap = 'hot', aspect = 'auto', extent = (-21.6/2, 21.6/2, 0, 12))
    plt.xlabel('Receiver horizontal postion [m]')
    plt.ylabel('Receiver vertical position [m]')
    plt.colorbar()

    #Find HALOS & SP & difference in sum & avg fluxes 
    halos_sum_of_flux = halos_flux.sum()
    sp_sum_of_flux = sp_flux.sum()
    halos_max_flux = halos_flux.max()
    sp_max_flux = sp_flux.max()
    halos_avg_flux = np.average(halos_flux)
    sp_avg_flux = np.average(sp_flux)
    sum_of_flux_diff = diff.sum()
    max_flux_diff = diff.max()
    avg_flux_diff = np.average(diff)
    print(f'Difference in sum of flux: {sum_of_flux_diff}\nDifference in max flux: {max_flux_diff}\nDifference in avg flux: {avg_flux_diff}')




