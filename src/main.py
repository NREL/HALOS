# -*- coding: utf-8 -*-
"""
OHLASTS main module

This module creates a flux_model object according to user inputs, then creates
an instance of the aimpoint strategy optimization model.

"""
import solve_aim_model


if __name__ == "__main__":
    case_filename = "./../case_inputs/radial_50_ca_case.csv"
    #case_filename = "./../case_inputs/radial_tonopah.csv"
    case_name = "radial_daggett_50"
    #case_name = "p_flat_daggett_50"
    
    
    solve_aim_model.runHourlyCase(case_name, case_filename, decomp = True, parallel = True)
    solve_aim_model.runSPHourlyCase(case_name, case_filename)


        
    