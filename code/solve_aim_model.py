# -*- coding: utf-8 -*-
"""
Module that solves the aimpoint optimization model, whether by decomposition 
or directly.
"""

import time

import inputs
import optimize_aimpoint
import process_aimpoint_outputs

def solveModelDirectly(case_name, case_filename):
    """
    Solve complete problem without decomp

    """
    t_start= time.time()
    fm = inputs.getFullFluxModelFromFiles(case_filename)
    ao_params = {"section_id":1,"num_sections":1,"hour_idx":fm.settings["hour_idx"],"aimpoint_cons_only":False}
    ao = optimize_aimpoint.AimpointOptimizer(fm,ao_params)
    elapsed_1 = time.time()-t_start
    print("time to set up objects: ",elapsed_1)
    t_start= time.time()
    ao.createFullProblem()
    elapsed_2 = time.time()-t_start
    print("time to build model: ",elapsed_2)
    t_start= time.time()
    ao.optimize()
    elapsed_3 = time.time()-t_start
    print("time to solve:",elapsed_3)
    t_start = time.time()
    outputs = process_aimpoint_outputs.AimpointOptOutputs(ao,fm)
    elapsed_4 = time.time() - t_start
    print("Output processing time: ",elapsed_4)
    build_and_solve_time = elapsed_1+elapsed_2
    print("Total computing time: ",(elapsed_1+elapsed_2+elapsed_3))
    outputs.setSolveTime(build_and_solve_time)
    return outputs
    
def subproblemSolve(kwargs):
    """
    Solve each section as subproblem

    Parameters
    ----------
    kwargs : dict
        Input parameters

    Returns
    -------
    d : dict
        aimpoint optimization results

    """
    fm = kwargs["flux_model"]
    ao_params = {"section_id":kwargs["section_id"], 
                 "num_sections":kwargs["num_sections"],
                 "ordered_defocus":kwargs["ordered_defocus"],
                 "order_method":kwargs["order_method"],
                 "aimpoint_cons_only":False}
    t_start_1 = time.time()
    ao = optimize_aimpoint.AimpointOptimizer(fm,ao_params)
    print("Start")
    ao.createFullProblem()
    print("Problem Created")
    elapsed_build = time.time() - t_start_1 
    print("Time to build: ", elapsed_build)
    ao.optimize()
    print("Processing results")
    #d = ao.processResults()
    d = ao.processOutputs()
    return d

def solveDecomposedModel(case_name, case_filename,hour_id = None, 
                         parallel=False, ordered_defocus=False,
                         order_method="eff"):
    """
    Decompose complete problem into subproblems 

    Parameters
    ----------
    hour_id : 
        Hour index or hour number of the year 
    ordered_defocus : 
        orderly defocus Heliostats
    order_method : 
        order defocus method, "eff" for efficiency based defocus and "dist" for
        distance based ordered defocus.

    Returns
    -------
    outputs : Combined outputs of subproblems

    """
    t_start= time.time()
    fm = inputs.getFullFluxModelFromFiles(case_filename,hour_id)
    elapsed_1 = time.time()-t_start
    print("time to set up objects: ",elapsed_1)
    t_start= time.time()
    sub_inputs = []
    kwargs = {"flux_model":fm, "num_sections":fm.field.num_sections,
              "ordered_defocus":ordered_defocus, "order_method":order_method}
    for idx in range(fm.field.num_sections):
        kwargs["section_id"] = idx
        sub_inputs.append(kwargs.copy())
    if parallel:
        import multiprocessing
        p = multiprocessing.Pool(multiprocessing.cpu_count())
        results = p.map(subproblemSolve, sub_inputs)
        del(p)
    else: 
        results = []
        for idx in range(len(sub_inputs)):
            results.append(subproblemSolve(sub_inputs[idx]))
    
    elapsed_2 = time.time()-t_start
    print("time to build and solve:",elapsed_2)
    t_start = time.time()
    outputs = process_aimpoint_outputs.AimpointOptOutputs(results, fm)
    elapsed_3 = time.time() - t_start
    print("Output processing time: ",elapsed_3)
    build_and_solve_time = elapsed_1+elapsed_2
    print("Total computing time: ",(elapsed_1+elapsed_2+elapsed_3))
    outputs.setSolveTime(build_and_solve_time)
    return outputs
            
        
def runHourlyCase(main_case_name, case_name, case_filename, hour_id = None, decomp = False, parallel = True,
            plots = True, print_outputs = True, append=False, bigfile = None):
    """
    Run aimpoint optimization for any specific hour 

    Parameters
    ----------
    main_case_name : 
        Name for the main result CSV file
    case_name : 
        Name for individual hour case CSV and plots
    case_filename : 
        CSV with input filenames path
    hour_id : 
        The default is None, optional hour_id to be simulated
    decomp : 
        Option to use solve using decomposition 
    plots : 
        plottin of results, optional      
    append : 
        If True: Appends hourly results into the main csv summary file
    bigfile:
        If passed in, to create summary csv with a row for each case run.
        Useful when running a batch file or many cases in a row.

    Returns
    -------
    None.Writes results into CSV and creates plots
    """
    if decomp:
        results = solveDecomposedModel(case_name, case_filename,hour_id, parallel)
    else: 
        results = solveModelDirectly(case_name, case_filename)
    if plots:
        results.plotOutputs(case_name)
    if print_outputs:
        results.printOutput(case_name)
    rw = 'a' if append else 'w'
    ofile = open("./../outputs/"+main_case_name+"_summary.csv",rw)
    if not append:
        #ofile.write("case_name,obj_value,num_defocused_heliostats,solve_time,init_feas_obj,init_feas_time,post_obj,post_num_defoc\n") # if want to record all initial and post-refocusing values (obj val, time, number defocused)
        ofile.write("case_name,obj_value,num_defocused_heliostats,solve_time,post_obj\n")
    #ofile.write(case_name+","+str(results.obj_value)+","+str(results.num_defocused_heliostats)+","+str(results.solve_time)+","+str(results.obj_val_feas_add)+","+str(results.time_add)+","+str(results.new_obj)+","+str(results.post_num_defocused)+"\n") # if recording all info
    ofile.write(case_name+","+str(results.obj_value)+","+str(results.num_defocused_heliostats)+","+str(results.solve_time)+","+str(results.new_obj)+"\n")
    ofile.close()
    from csv import writer
    if bigfile != None:
        with open(bigfile,'a+',newline='') as write_obj:
            csv_writer = writer(write_obj)
            row1 = [case_name,str(results.obj_value),str(results.num_defocused_heliostats), str(results.solve_time),str(results.new_obj)] # change according to which values you want to return
            csv_writer.writerow(row1)
    
def runSPHourlyCase(main_case_name,case_name, case_filename, hour_id = None, append = False, weather_data = None, sp_aimpoint_heur = True):
    """
    Run full SolarPilot Case for any specific hour
    
    Parameters
    ----------
    main_case_name : 
        Name for the main result CSV file
    case_name : 
        Name for individual hour case CSV and plots
    case_filename : 
        CSV with input filenames path
    hour_id : 
        The default is None, optional hour_id to be simulated
    append : 
        If True: Appends hourly results into the main csv summary file
    weather_data : Dataframe, optional
        The default is None.

    Returns
    -------
    None. Writes a CSV file with results 

    """
    import sp_module
    import inputs
   
    filenames = inputs.readCaseFile(case_filename)
    sp_flux = sp_module.SP_Flux(filenames,use_sp_field = True)
    outputs = sp_flux.run_sp_case(case_name = case_name, sp_aimpoint_heur = sp_aimpoint_heur, 
                                  saveCSV = True, hour_id = hour_id, weather_data = weather_data)
    rw = 'a' if append else 'w'
    ofile = open("./../outputs/SolarPilot_" + main_case_name +"_summary.csv",rw)
    if sp_aimpoint_heur:
        if not append:
            ofile.write("case_name,obj_value,num_defocus,max_flux,obj_pre_heur,obj_post_heur,post_num_defocused\n")
        ofile.write(case_name+","+str(outputs["obj_value"])+","+str(outputs["num_defocus"])
                    +","+str(outputs["max_flux"])+","+str(outputs["before_flux"])+","+
                    str(outputs["post_heur_obj"])+","+str(outputs["post_num_defocus"])+"\n")
        ofile.close()
    if sp_aimpoint_heur != True:
        if not append:
            ofile.write("case_name,obj_value,num_defocus,max_flux\n")
        ofile.write(case_name+","+str(outputs["obj_value"])+","+str(outputs["num_defocus"])
                    +","+str(outputs["max_flux"])+"\n")
        ofile.close()