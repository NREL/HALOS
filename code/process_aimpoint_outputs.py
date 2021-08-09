# -*- coding: utf-8 -*-
"""
HALOS Aimpoint Optimization Model's Output Processing Module 

"""
import numpy

class AimpointOptOutputs(object):
    def __init__(self, results, flux_model,post_refocus = True):
        '''
        Parameters
        ------------
        post_refocus :
            If post_refocus = True: after running optimization model, loop through each 
            defocused heliostat to see whether it can point at an aimpoint without violating a flux limit.
            Set post_refocus to False in __init__ of AimpointOptOutputs in process_aimpoint_outputs
            to exclude this step.
            Note: Not included as an input to solve_aim_model because of associated pyomo errors when tried that.
        '''
        self.flux_model = flux_model
        # to print zeros for post_num_defocused and new_obj if post_refocus = False
        self.post_num_defocused = 0
        self.new_obj = 0
        if type(results) is list:
            self.processSubproblemOutputs(results)
            self.removeFluxViolation()
            if post_refocus:
                self.refocusHeliostats()
                self.newObj()
        else:
            self.processFullProblemOutputs(results)
    
    def setSolveTime(self,solve_time):
        self.solve_time = solve_time
        
    def processFullProblemOutputs(self, d):
        """
        Processes output of full problem without decomposition 

        Parameters
        ----------
        d : Dict
        Aimpoint Optimization results of Full problem 

        """
        self.flux_map = d["flux_map"]
        self.aimpoint_select_map = d["aimpoint_select_map"]
        self.obj_value = d["obj_value"]
        self.num_defocused_heliostats = d["num_defocused_heliostats"]
        self.utilization_by_section = float(d["num_defocused_heliostats"]) / self.flux_model.field.num_heliostats
        self.contribution_by_heliostat = d["contribution_by_heliostat"]
        self.measurement_points = d["measurement_pts"]
        self.flux_ham = d["flux_ham"]
        self.surface_area = d["surface_area"]
        self.aimpoints = d["aimpoints"]
        self.obj_val_feas_add =  d["obj_val_feas_add"]
        self.time_add = d["time_add"]
        self.init_defocused = d["init_defocused"]
                    
    def processSubproblemOutputs(self, ds):
        """
        Processes output of decomposed problem 

        Parameters
        ----------
        ds : List
            Contains Dict of each section outputs from Aimpoint Optimization

        """
        self.flux_map = numpy.zeros(self.flux_model.receiver.x.size,dtype=float)
        self.aimpoint_select_map = numpy.zeros(self.flux_model.field.x.size)
        self.obj_value = 0.0
        self.num_defocused_heliostats = 0
        self.utilization_by_section = numpy.zeros(self.flux_model.field.num_sections)
        self.contribution_by_heliostat = numpy.zeros([self.flux_model.field.x.size,self.flux_model.receiver.x.size])
        self.aimpoints = []
        self.measurement_points = []
        self.flux_ham = {}
        self.surface_area = {}
        self.heliostats = []
        self.obj_val_feas_add = 0
        self.time_add = 0
        self.init_defocused = 0
        self.ub = []


        for d in ds:
            self.flux_map += d["flux_map"]
            self.aimpoint_select_map += d["aimpoint_select_map"]
            self.obj_value += d["obj_value"]
            self.num_defocused_heliostats += d["num_defocused_heliostats"]
            self.contribution_by_heliostat += d["contribution_by_heliostat"]
            self.utilization_by_section[d["section_id"]] = d["utilization"]
            if d["section_id"] == 1:
                self.aimpoints += d["aimpoints"]
                self.measurement_points = d["measurement_pts"]
            self.flux_ham.update(d["flux_ham"])
            self.surface_area.update(d["surface_area"])
            self.obj_val_feas_add += d["obj_val_feas_add"]
            self.time_add += d["time_add"]
            self.init_defocused += d["init_defocused"]

        # make file with obj val and ubs per sxn, as well as overall obj val - for 8 sxns
        sxn_file = open("section_file_highermipagain.csv","a+",newline="")
        i = 1
        for d in ds:
            sxn_file.write('S'+str(i)+'_ub: '+str(d["upper_bound"])+"\n")
            sxn_file.write('S'+str(i)+'_obj: '+str(d["obj_value"])+"\n")
            i = i + 1
    
    def plotOutputs(self, case_name=""):
        """
        Plot Aimpoint Optimization Outputs

        """
        import plotting
        plotting.plot_optimal_flux_heatmap(self,case_name+"_fluxmap")
        plotting.plot_optimal_aimpoint_allocation(self,case_name+"_aimpoint_map")
        plotting.plot_optimal_aimpoint_guide(self,case_name+"_aimpoint_select_map")
        plotting.plot_measurement_and_aimpoints(self,case_name+"_measurement_and_aimpoints")
        plotting.plot_defocused(self,case_name+"_defocused_Heliostats")
        if self.flux_model.settings["num_sections"] > 1:
            plotting.plot_field_sections(self,case_name+"_field_sections")
        else:
            plotting.plot_field(self,case_name+"_field")
        plotting.plot_flux_violation(self, case_name+"_violation")
        plotting.plot_flux_slack(self, case_name+"_slack")
        
    def getFluxViolation(self):
        """
        calculates flux violation at each measurement point on the receiver
        after aimpoints have been optimized by section.
        """
        flux_ub = self.flux_model.receiver.flux_upper_limits
        flux = self.flux_map
        flux_violation = numpy.zeros_like(flux)
        for m in range(len(flux)):
            flux_violation[m] = max(0.0, flux[m]-flux_ub[m])
        self.flux_violation = numpy.array(flux_violation)
    
    def getMaxViolationContributor(self, violation, m):
        """
        returns the heliostat contibuting either (i) the greatest flux 
        contribution to a measurement point, or (ii) the lowest flux required
        to remove the violation.
        
        Parameters
        ==========
        violation -- amount of flux violation at measurement point
        m -- measurement point identifier
        
        Returns
        =======
        best_h -- Heliostat ID for selected heliostat
        """
        over_v = False  #indicator that flux> violation found
        best_f = 0.0
        best_h = -1
        for h in range(len(self.aimpoint_select_map)):
            if self.contribution_by_heliostat[h,m] > best_f: #
                if not over_v:
                    best_h = h
                    best_f = self.contribution_by_heliostat[h,m]
                    if best_f > violation:
                        over_v = True
            elif self.contribution_by_heliostat[h,m] >= violation and over_v:  #closer to violation than best_f
                best_f = self.contribution_by_heliostat[h,m]
                best_h = h
        return best_h
        
    
    def removeFluxViolation(self, threshold=1.0e-3):
        """
        removes heliostats from aimpoint strategy optimization model until
        flux is feasible.
        """
        self.getFluxViolation()
        while self.flux_violation.max() > threshold:
            print(self.flux_violation.max())
            #obtain best candidate heliostat to remove
            m = self.flux_violation.argmax()
            violation = self.flux_violation[m]
            h = self.getMaxViolationContributor(violation, m)
            self.flux_map -= self.contribution_by_heliostat[h,:]
            self.flux_violation -= self.contribution_by_heliostat[h,:]
            self.flux_violation = numpy.maximum(0,self.flux_violation)
            self.contribution_by_heliostat[h,:] *= 0.0
        m_rows = self.flux_model.receiver.params["pts_per_ht_dim"]
        m_cols = self.flux_model.receiver.params["pts_per_len_dim"]
        self.flux_violation.reshape(m_rows,m_cols)


    def refocusHeliostats(self, threshold=1.0e-3):
        '''
        Adds previously defocused heliostats from the model if doing so does not violate flux limits.
        For each previously defocused heliostat, cycles through the possible aimpoints to see 
        if for one of them the additional flux added to each of the measurement points does not 
        violate any flux limit. Each heliostat refocuses on the first aimpoint for which it 
        can focus without violating flux limits.
        '''
        
        # obtain list of heliostat ids that are defocused
        defoc_dict = {}
        defoc_list = []
        for h in range(len(self.aimpoint_select_map)):
            defoc_dict[h] = self.aimpoint_select_map[h]
            if defoc_dict[h] == 0:
                defoc_list.append(h)
        self.post_num_defocused = self.num_defocused_heliostats
        self.h_refocused = []
        self.a_refocused = []
        for h in defoc_list:
            for a in self.aimpoints:
                for m in self.measurement_points:
                    self.flux_map += self.flux_ham[h,m,a]
                # calculate flux violation with this new flux map
                self.getFluxViolation()
                # if violates flux, delete the flux contribution of h-a combo on each m
                if self.flux_violation.max() > threshold:
                    for m in self.measurement_points:
                        self.flux_map -= self.flux_ham[h,m,a]
                # otherwise no flux limit violated and the heliostat is refocused to that aimpoint
                else: 
                    self.post_num_defocused -= 1
                    self.h_refocused.append(h)
                    self.a_refocused.append(a)
                    # so heliostat refocuses only at the first aimpoint for which it can  
                    break

    def newObj(self):
        '''
        Adds contributions of any refocused heliostats to the previously calculated obj val
        '''
        post_add_obj = 0
        for i in range(len(self.h_refocused)):
            h = self.h_refocused[i]
            a = self.a_refocused[i]
            print('heliostat ',h,'refocuses at aimpoint ',a)
            for m in self.measurement_points:
                post_add_obj += self.surface_area[m]* self.flux_ham[h,m,a]

        self.new_obj = post_add_obj + self.obj_value

    def printOutput(self, case_name, console = False):
        """
        Write/print outputs as CSV

        Parameters
        ----------
        case_name : TYPE
            DESCRIPTION.
        console :  optional
            Prints output on Console if True

        Returns
        -------
        None. Writes CSV & prints outputs on Console 

        """
        ofile = open(case_name+"_solution.csv",'w')
        ofile.write("obj_value: "+str(self.obj_value)+"\n\n")
        if self.new_obj != 0:
            ofile.write("new_obj: "+str(self.new_obj)+"\n\n")
        ofile.write("defocused heliostats: "+str(self.num_defocused_heliostats)+"\n\naimpoint summary:\n")
        if self.post_num_defocused != 0:
            ofile.write("post defocused heliostats: "+str(self.post_num_defocused)+"\n\naimpoint summary:\n")
        ofile.write("Initial feasible objective: "+str(self.obj_val_feas_add))
        ofile.write("Initial Heuristic Time: "+str(self.time_add))
        ofile.write("heliostat_id,aimpoint_id\n")
        for h in range(len(self.aimpoint_select_map)):
            ofile.write(str(h)+","+str(self.aimpoint_select_map[h])+"\n")
        ofile.close()
        if console:
            for h in range(len(self.aimpoint_select_map)):
                print("aim heliostat ",h," at aimpoint ",self.aimpoint_select_map[h])
            print("obj value: ",self.obj_value)
        ofile2 = open(case_name+"_utilization.csv","w")
        for i in range(len(self.utilization_by_section)):
            ofile2.write(str(self.utilization_by_section[i])+",")
        ofile2.write("\n")
        ofile2.close()
