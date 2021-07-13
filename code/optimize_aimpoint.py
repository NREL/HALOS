# -*- coding: utf-8 -*-
"""
HALOS optimize module

This module contains the optimization model, which is implemented in the 
Pyomo modeling language.

"""

import pyomo
import pyomo.environ as pe
import sol_pos
import numpy

### Methods or objective and constraint rules
EPSILON = 1e-8
def objectiveRule(model):
        return sum(  #m 
                    sum( #h
                        sum(#a
                            model.obj_by_point[m] * model.surface_area[m]* model.flux[h,m,a] * model.select_aimpoint[h,a]
                            for a in model.aimpoints
                        ) for h in model.heliostats
                    ) for m in model.measurement_points
                )


def aimSelectRule(model,h):
    return sum(model.select_aimpoint[h,a] for a in model.aimpoints) + model.defocus[h] == 1


def fluxLBRule(model,m):
    if sum(sum(model.flux[h,m,a] for a in model.aimpoints) for h in model.heliostats) < EPSILON:
        return pe.Constraint.Feasible
    return sum(sum(model.flux[h,m,a] * model.select_aimpoint[h,a] for a in model.aimpoints) for h in model.heliostats) >= model.flux_lbs[m]


def fluxUBRule(model,m):
    # If the flux limit for the section is lower than the min provided as input,
    # then do not create a constraint.  This is intended to remove columns on
    # cylindrical receivers with very low flux for a given section that might
    # yield infeasibilities for the subsection.
    # TODO check for infeasibilities caused by this relaxation in post-processing
    if model.flux_ubs[m] < model.flux_constraint_limit:
        return pe.Constraint.Feasible
    if sum(sum(model.flux[h,m,a] for a in model.aimpoints) for h in model.heliostats) < EPSILON:
        return pe.Constraint.Feasible
    return sum(sum(model.flux[h,m,a] * model.select_aimpoint[h,a] for a in model.aimpoints) for h in model.heliostats) <= model.flux_ubs[m]


def fluxDiffRule(model,m,mp):
    if m not in model.neighboring_points[mp]:
        return pe.Constraint.Feasible
    return (
            sum(sum(model.flux[h,m,a] * model.select_aimpoint[h,a] for a in model.aimpoints) for h in model.heliostats) -
            sum(sum(model.flux[h,mp,a] * model.select_aimpoint[h,a] for a in model.aimpoints) for h in model.heliostats) 
            <= model.flux_diff
            )


def orderedDefocusingRule(model, h, hp):
    if model.order[h]+1 == model.order[hp]:
        return model.defocus[h] >= model.defocus[hp]
    return pe.Constraint.Feasible

def groupingRule(model, h, hp, a):
    if h < hp and model.group_assignment[h] == model.group_assignment[hp]:
        return model.select_aimpoint[h,a] == model.select_aimpoint[hp,a]
    return pe.Constraint.Feasible

def neighborRule(model, m):
    """
    determines neighboring measurement points.  in this case, we assume that the neighbors are all vertical.
    """
    if m+model.num_cols <= model.num_measurement_points:
        yield m+model.num_cols
    if m-model.num_cols <= 0:
        yield m-model.num_cols

class AimpointOptimizer(object):
    def __init__(self, flux_model, params={"section_id":1,"num_sections":1,
                                           "aimpoint_cons_only":False,
                                           "ordered_defocus":True,
                                           "order_method":"eff"}):
        """
        Creates the basic pyomo object for the aimpoint optimization model.
        
        ===== Parameters =====
        flux_model -- flux_model object containing all model inputs
        params -- dictionary with additional parameters for the model, 
            including deomposition settings.
        
        ===== Returns =====
        None
        """
        self.params = params
        self.flux_model = flux_model
        self.model = pe.ConcreteModel()
        self.params["flux_constraint_limit"] = 500 /params["num_sections"]
        self.solver = params.get("solver")
        if self.solver is None:
            self.solver = "cbc"

    def generateMeasurementSubset(self):
        measurement_points = []
        m_rows = int(self.flux_model.receiver.params["pts_per_ht_dim"])
        m_cols = int(self.flux_model.receiver.params["pts_per_len_dim"])
        if self.flux_model.receiver.params["receiver_type"] == "Flat plate":
            aim_cols = int(self.flux_model.receiver.params["aim_cols"])

        elif self.flux_model.receiver.params["receiver_type"] == "External cylindrical":
            aim_cols = 1
  
        aim_rows = int(self.flux_model.receiver.params["aim_rows"])
        row_spacing = m_cols // aim_cols
        col_spacing = m_rows // aim_rows
        first_row = row_spacing // 2
        first_col = col_spacing // 2
        for r in range(aim_rows):
            for c in range(aim_cols):
                pt = (first_row + (r*row_spacing)) * m_rows + (first_col + (c*col_spacing)) + 1
                measurement_points.append(pt)
        self.model.check_measurement_points = pe.Set(initialize = measurement_points)


    def generateSets(self):
        """
        Generate sets for the pyomomodel object.
        """
        self.model.num_measurement_points = self.flux_model.receiver.x.size
        self.model.num_rows = self.flux_model.receiver.params["pts_per_ht_dim"]
        self.model.num_cols = self.flux_model.receiver.params["pts_per_len_dim"]
        self.model.measurement_points = pe.Set(
                initialize = range(1,self.model.num_measurement_points+1)
                )  #M
        
        self.generateMeasurementSubset()

        if self.flux_model.receiver.params["use_flux_gradient"] == 1:
            self.model.neighboring_points = pe.Set(self.model.measurement_points, initialize=neighborRule)
        
        self.num_aimpoints = self.flux_model.receiver.num_aimpoints
        self.model.aimpoints = pe.Set(initialize = range(1,self.num_aimpoints+1)) #A
        
        if self.params["num_sections"] == 1:    
            self.num_heliostats = self.flux_model.field.x.size
            self.model.heliostats = pe.Set(initialize = range(self.num_heliostats)) #H
        else:
            heliostats = self.flux_model.field.helios_by_section[self.params["section_id"]]
            self.num_heliostats = len(heliostats)
            self.model.heliostats = pe.Set(initialize = heliostats)
        self.helio_list = list(self.model.heliostats)
        if self.flux_model.settings["heliostat_group_size"] > 1:
            self.num_heliostat_groups = ((self.num_heliostats // 
                    self.flux_model.settings["heliostat_group_size"])  + (
                    1 if (self.num_heliostats % 
                    self.flux_model.settings["heliostat_group_size"] > 0) else 0
                    ))
            self.heliostat_groups = pe.Set(initialize = range(self.num_heliostat_groups))
            group_d = {}
            for i in range(self.num_heliostat_groups):
                for hidx in range(i*self.flux_model.settings["heliostat_group_size"],min(self.num_heliostats,(i+1)*self.flux_model.settings["heliostat_group_size"])):
                    group_d[self.helio_list[hidx]] = i
            self.model.group_assignment = pe.Param(self.model.heliostats, initialize = group_d)
       
        
    def getFluxParameters(self, approx=True):
        """
        calculates flux for each heliostat-aimpoint pairing, at each masurement
        point.
        [az] note: this is done in for loops for now, but I believe we should
        be able to do some of these by passing arrays into GetFlux() 
        in flux_method; this will require rework of GetFlux().
        
        Parameters
        ==========
        t -- time index (from which weather data will be taken)
        Returns 
        =======
        """
        dni = self.flux_model.dni
        zenith, azimuth = sol_pos.get_single_solar_position(self.flux_model.weather_data, self.flux_model.settings["hour_idx"])
        solar_vector = self.flux_model.getSolarVector(azimuth, zenith)
        ## Get flux, gradient limits from receiver object
        flux_lbs = {}
        flux_ubs = {}
        surface_area = {}
        obj_by_point = {}
        if self.params["num_sections"] == 1:
            fraction = 1.0
        else: 
            if self.flux_model.receiver.params["receiver_type"] == "External cylindrical":
                fraction = self.flux_model.fraction_maps[self.params["section_id"]].flatten()
            else:   
                fraction = self.flux_model.field.section_flux[self.params["section_id"]]
        for i in range(self.model.num_measurement_points):
            if self.flux_model.receiver.params["receiver_type"] == "External cylindrical":
                flux_lbs[i+1] = self.flux_model.receiver.flux_lower_limits[i] * fraction[i]
                flux_ubs[i+1] = self.flux_model.receiver.flux_upper_limits[i] * fraction[i]
            else:
                flux_lbs[i+1] = self.flux_model.receiver.flux_lower_limits[i] * fraction
                flux_ubs[i+1] = self.flux_model.receiver.flux_upper_limits[i] * fraction      
            surface_area[i+1] = self.flux_model.receiver.surface_area[i]
            obj_by_point[i+1] = self.flux_model.receiver.obj_by_point[i]
        #if specifying flux maps from a file, do so. #TODO remove this as the .csv's method will replace
        try: 
            if self.params["flux_from_file"]:
                flux = self.getFluxFromFile()
                return flux, flux_lbs, flux_ubs, surface_area, obj_by_point
        except KeyError:
            pass
        flux = {}
        for h in self.model.heliostats:
            center_idx = self.flux_model.receiver.num_aimpoints//2
            h_map = self.flux_model.ShiftImage_GenerateSingleHeliostatFluxMaps(h,solar_vector,approx)
            for m in self.model.measurement_points:
                if self.flux_model.receiver.params["receiver_type"] == "Flat plate":
                    for a in self.model.aimpoints:
                        flux[h,m,a] = h_map[a-1][m-1]
                elif self.flux_model.receiver.params["receiver_type"] == "External cylindrical":
                    for a in self.model.aimpoints:
                        flux[h,m,a] = h_map[a-1][m-1]
        return flux, flux_lbs, flux_ubs, surface_area, obj_by_point
    
    
    def getFluxFromFile(self):
        """
        Uses CSV file to populate flux

        Returns
        -------
        flux : Dict

        """
        import csv
        f = csv.reader(open(self.params["flux_filename"],'r'))
        next(f)
        flux = {}
        for line in f:
            if len(line) > 3:
                flux[int(line[0]),int(line[1]),int(line[2])] = float(line[3])
        return flux

        
    def generateParameters(self):
        """
        Generate Parameters for the optimization model

        Returns
        -------
        None. Assigns inputs to the pyomo model

        """
        flux, flux_lbs, flux_ubs, surface_area, obj_by_point = self.getFluxParameters()
        self.model.flux = pe.Param(self.model.heliostats, self.model.measurement_points, self.model.aimpoints, initialize = flux)
        self.model.flux_lbs = pe.Param(self.model.measurement_points, initialize = flux_lbs)
        self.model.flux_ubs = pe.Param(self.model.measurement_points, initialize = flux_ubs)
        if self.flux_model.receiver.params["use_flux_gradient"] == 1:
            self.model.flux_diff = pe.Param(mutable=True, initialize = self.flux_model.receiver.params["gradient_limit"]/self.params["num_sections"])
        self.model.surface_area = pe.Param(self.model.measurement_points, initialize = surface_area)
        self.model.obj_by_point = pe.Param(self.model.measurement_points, initialize = obj_by_point)
        if self.params["ordered_defocus"]:
            self.getHeliostatOrdering(self.params["order_method"])
        self.model.flux_constraint_limit = self.params["flux_constraint_limit"]


    def printFluxToFile(self,filename):
        """
        Write Flux to a CSV File

        Parameters
        ----------
        filename :         
            Output filename

        """
        ofile = open(filename,"w")
        ofile.write("h,m,a,flux\n")
        for h in self.model.heliostats:
            for m in self.model.measurement_points:
                for a in self.model.aimpoints:
                    ofile.write(str(h)+","+str(m)+","+str(a)+","+str(self.model.flux[h,m,a])+"\n")
        ofile.close()                    


    def generateVariables(self):
        self.model.select_aimpoint = pe.Var(self.model.heliostats * self.model.aimpoints, domain=pe.Binary)
        self.model.defocus = pe.Var(self.model.heliostats, domain=pe.NonNegativeReals, bounds=(0,1))

    def getHeliostatOrdering(self, method="eff"):
        """
        Orderly defocuses Heliostats based of efficiency or distance from 
        receiver. 

        Parameters
        ----------
        method : 
            Efficiency ("eff") or distance ("dist") method

        Returns
        -------
        Assigns Heliostat ordering to optimization model

        """
        if method == "eff":
            vals = numpy.array([self.flux_model.field.eff[h] for h in self.helio_list])
            ordering = numpy.argsort(vals)
        elif method == "dist":
            vals = numpy.array([self.flux_model.field.x[h]**2 * self.flux_model.field.y[h]**2 for h in self.helio_list])
            ordering = numpy.size(vals)-numpy.argsort(vals)-1
        d_order = {}
        for idx in range(self.num_heliostats):
            d_order[self.helio_list[ordering[idx]]] = idx
        self.model.order = pe.Param(self.model.heliostats, initialize=d_order)

             
    def setObjective(self):
        self.model.OBJ = pe.Objective(rule=objectiveRule, sense = pe.maximize)
      
    def genConstraintsBinOnly(self): 
        self.model.select_con = pe.Constraint(self.model.heliostats, rule=aimSelectRule)
        if self.params["aimpoint_cons_only"]:
            self.model.flux_ub_con = pe.Constraint(self.model.check_measurement_points, rule=fluxUBRule)    
        else:
            self.model.flux_ub_con = pe.Constraint(self.model.measurement_points, rule=fluxUBRule)
        if self.flux_model.receiver.params["use_flux_gradient"] == 1:
            self.model.flux_diff_con = pe.Constraint(self.model.measurement_points * self.model.measurement_points, rule=fluxDiffRule)
        if self.params["ordered_defocus"]:
            self.model.ordered_defocusing_con = pe.Constraint(self.model.heliostats * self.model.heliostats, rule=orderedDefocusingRule)
        if self.flux_model.settings["heliostat_group_size"] > 1: 
            self.model.ordered_defocusing_con = pe.Constraint(self.model.heliostats * self.model.heliostats * self.model.aimpoints, rule=groupingRule)
            print("group cons made")

    def createFullProblem(self):
        """
        Creates full optimization problem for solver

        Returns
        -------
        None.

        """
        self.generateSets()
        self.generateParameters()
        self.generateVariables()
        self.setObjective()
        self.genConstraintsBinOnly()


    def optimize(self, mipgap=0.001, timelimit=300, tee=False, keepfiles=False, warmstart=False):
        """
        Solves the optimization model 

        Parameters
        ----------
        mipgap : 
            optimality gap. The default is 0.001.
        timelimit : 
            Time limit for the solve in seconds. The default is 300.
        solver : 
            Open Source: 'cbc', 'glpk'
            Commercial: 'cplex'
        warmstart : 
            Uses opt_heuristic to generate an initial feasible solution.  

        Raises
        ------
        Solver 

        Returns
        -------
        Assigns results to self.opt_results

        """
        if warmstart: 
            import opt_heuristic
            heuristic = opt_heuristic.AcceptRejectHeuristic(3)
            heuristic.getIFS(self.model)
        opt = pe.SolverFactory(self.solver)
        if self.solver == "cbc":
            opt.options["ratioGap"] = mipgap
            opt.options["seconds"] = timelimit
        elif self.solver == "glpk":
            opt.options["mipgap"] = mipgap
            opt.options["tmlim"] = timelimit
            opt.options["cuts"] = 1
        elif self.solver == "cplex":
            opt.options["mipgap"] = mipgap
            opt.options["timelimit"] = timelimit
        else:
            raise Exception("invalid solver.")
        self.opt_results = opt.solve(self.model, tee=tee, keepfiles=keepfiles, warmstart=warmstart, load_solutions=False)
        self.gap = self.opt_results.solution[0].gap
        self.model.solutions.load_from(self.opt_results)
        
            
    def processOutputs(self):
        """
        Processess outputs of aimpoint optimization model

        Returns
        -------
        d : Dict 
            Contains results from optimization model

        """
        self.flux_map = numpy.zeros(self.model.num_measurement_points,dtype=float)
        for m in range(1,self.model.num_measurement_points+1):
            try: 
                self.flux_map[m-1] = sum( #h
                        sum(#a
                            self.model.flux[h,m,a] * self.model.select_aimpoint[h,a].value
                            for a in self.model.aimpoints
                        ) for h in self.model.heliostats
                    )
            except TypeError:
                assert(False)
        self.aimpoint_select_map = numpy.zeros(self.flux_model.field.x.size)
        self.contribution_by_heliostat = numpy.zeros([self.flux_model.field.x.size, self.model.num_measurement_points])
        self.num_defocused_heliostats = 0
        for h in self.model.heliostats:
            if self.model.defocus[h].value > 0.5:
                self.num_defocused_heliostats += 1
            else:
                for a in self.model.aimpoints:
                    if self.model.select_aimpoint[h,a].value > 0.5:
                        self.aimpoint_select_map[h-1] = a
                        for m in range(self.model.num_measurement_points):
                            self.contribution_by_heliostat[h,m] = self.model.flux[h,m+1,a]
                        break
        if self.gap == None:
            ub = 1.0e10  
            # TODO parse the results text to obtain the gap when pyomo times out
        else:
            ub = pe.value(self.model.OBJ) * (self.gap + 1)
        d = {
            "flux_map":self.flux_map, 
             "aimpoint_select_map":self.aimpoint_select_map,
             "obj_value":pe.value(self.model.OBJ), 
             "num_defocused_heliostats":self.num_defocused_heliostats,
             "upper_bound":ub,
             "contribution_by_heliostat":self.contribution_by_heliostat,
             "section_id":self.params["section_id"],
             "utilization": 1.0 - (float(self.num_defocused_heliostats) / self.num_heliostats),
             "flux_ham":self.model.flux,
             "measurement_pts":self.model.measurement_points,
             "aimpoints":self.model.aimpoints,
             "surface_area":self.model.surface_area,
             }
        return d


    def plotOutputs(self):
        """
        Generates plots using aimpoint optimization results

        """
        import plotting
        plotting.plot_optimal_flux_heatmap(self,"fluxmap")
        plotting.plot_optimal_aimpoint_allocation(self,"aimpoint_map")
        plotting.plot_optimal_aimpoint_guide(self,"aimpoint_select_map")


    def printOutput(self):
        """
        Prints Aimpoint optimization results on Console

        """
        for h in self.model.heliostats:
            for a in self.model.aimpoints: 
                if self.model.select_aimpoint[h,a].value > 0.5:
                    print("aim heliostat ",h," at aimpoint ",a)
                    

if __name__ == "__main__":
    import flux_model
    import flux_method
    import geometry
    import field
    
    import numpy
    
    weather_filename = "./../weather_files/USA NV Tonopah Airport (TMY3).csv"
    aimpoints = numpy.array([[-3,0,3],[0,0,0],[150,150,150]])
    receiver_params = {"pts_per_dim":5,"length":20,"width":0,"height":20, "receiver_type": "Flat plate"}
    t_index = int(8760/2)
    flux_diff_max = 1e6
    r = geometry.FlatPlateReceiver(150,receiver_params)
    r.build_measurement_points()
#    r.generateAimpointsGrid()
#    r.getCoords()
##    print(r.coords)
#    r.getNormals()
    r.generateAimpointsGrid(
            3,3,0,0,1e-6
            )
    r.generateFixedFluxLimits(0,1000)
    r.getSurfaceArea()
    f = field.Field()
    f.GenerateCoords(numpy.array([-200,-100,0,100,200]),numpy.array([-300,-300,-300,-300,-300]),numpy.array([0,0,0,0,0]))
    flux_model = flux_model.FluxModel(None,None,r,flux_method.SimpleNormalFluxCalc(5),weather_filename,f)
    ao = AimpointOptimizer(flux_model)
    
    ao.generateSets()
    ao.model.pprint()
    ao.generateParameters(t_index,flux_diff_max)
    ao.generateVariables()
    ao.setObjective()
    ao.genConstraintsBinOnly()
    
    ao.model.pprint()
    
    
#    print(ao.getFluxParameters(0))
