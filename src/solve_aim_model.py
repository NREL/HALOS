# -*- coding: utf-8 -*-
"""
Module that solves the aimpoint optimization model, whether by decomposition 
or directly.
"""
import inputs
import optimize_aimpoint
import process_aimpoint_outputs
import time

def solveModelDirectly(case_name, case_filename):
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
            
        
def runHourlyCase(case_name, case_filename, hour_id = None, decomp = False, parallel = True,
            plots = True, print_outputs = True, append=False):
    if decomp:
        results = solveDecomposedModel(case_name, case_filename,hour_id, parallel)
    else: 
        results = solveModelDirectly(case_name, case_filename)
    if plots:
        results.plotOutputs(case_name)
    if print_outputs:
        results.printOutput(case_name)
    rw = 'a' if append else 'w'
    ofile = open(case_name+"_summary.csv",rw)
    if not append:
        ofile.write("case_name,obj_value,num_defocused_heliostats,solve_time\n")
    ofile.write(case_name+","+str(results.obj_value)+","+str(results.num_defocused_heliostats)+","+str(results.solve_time)+"\n")
    ofile.close()
    
def runSPHourlyCase(case_name, case_filename, hour_id = None, append = False, weather_data = None):
    import sp_module
    import inputs
   
    filenames = inputs.readCaseFile(case_filename)
    sp_flux = sp_module.SP_Flux(filenames,use_sp_field = True)
    outputs = sp_flux.run_sp_case(case_name = case_name, sp_aimpoint_heur = True, saveCSV = True, hour_id = hour_id, weather_data = weather_data)
    rw = 'a' if append else 'w'
    ofile = open("SolarPilot_" + case_name +"_summary.csv",rw)
    if not append:
        ofile.write("case_name,obj_value,num_defocus,max_flux,obj_pre_heur, obj_post_heur, post_num_defocused\n")
    ofile.write(case_name+","+str(outputs["obj_value"])+","+str(outputs["num_defocus"])+","+str(outputs["max_flux"])+","+str(outputs["before_flux"])+","+str(outputs["post_heur_obj"])+","+str(outputs["post_num_defocus"])+"\n")
    ofile.close()