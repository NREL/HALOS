# -*- coding: utf-8 -*-
"""
HALOS Aimpoint Optimization Model's Output Processing Module 

"""
import numpy
import pyomo.environ as pe


class AimpointOptOutputs(object):
    def __init__(self, results, flux_model):
        self.flux_model = flux_model
        if type(results) is list:
            self.processSubproblemOutputs(results)
            self.removeFluxViolation()
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
        for d in ds:
            self.flux_map += d["flux_map"]
            self.aimpoint_select_map += d["aimpoint_select_map"]
            self.obj_value += d["obj_value"]
            self.num_defocused_heliostats += d["num_defocused_heliostats"]
            self.contribution_by_heliostat += d["contribution_by_heliostat"]
            self.utilization_by_section[d["section_id"]] = d["utilization"]
    
    def plotOutputs(self, case_name="",fm = None):
        """
        Plot Aimpoint Optimization Outputs

        """
        import plotting
        plotting.plot_optimal_flux_heatmap(self,case_name+"_fluxmap")
        plotting.plot_optimal_aimpoint_allocation(self,case_name+"_aimpoint_map")
        plotting.plot_optimal_aimpoint_guide(self,case_name+"_aimpoint_select_map", fm)
        plotting.plot_defocused(self,case_name+"_defocused_Heliostats")
        if self.flux_model.settings["num_sections"] > 1:
            plotting.plot_field_sections(self,case_name+"_field_sections")
        else:
            plotting.plot_field(self,case_name+"_field")
        plotting.plot_flux_violation(self, case_name+"_violation")
        
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
        ofile.write("defocused heliostats: "+str(self.num_defocused_heliostats)+"\n\naimpoint summary:\n")
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
