# -*- coding: utf-8 -*-
"""
OHLASTS optimize module

This module contains the optimization model, which is implemented in the 
Pyomo modeling language.

"""
#import flux_model
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
                            model.surface_area[m]* model.flux[h,m,a] * model.select_aimpoint[h,a]
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
    if (sum(sum(model.flux[h,m,a] for a in model.aimpoints) for h in model.heliostats) < EPSILON and
            sum(sum(model.flux[h,mp,a] for a in model.aimpoints) for h in model.heliostats) < EPSILON):
        return pe.Constraint.Feasible
    return (
            sum(sum(model.flux[h,m,a] * model.select_aimpoint[h,a] for a in model.aimpoints) for h in model.heliostats) -
            sum(sum(model.flux[h,mp,a] * model.select_aimpoint[h,a] for a in model.aimpoints) for h in model.heliostats) 
            <= model.flux_diff[m,mp]
            )


def orderedDefocusingRule(model, h, hp):
    if model.order[h]+1 == model.order[hp]:
        return model.defocus[h] >= model.defocus[hp]
    return pe.Constraint.Feasible

def groupingRule(model, h, hp, a):
    if h < hp and model.group_assignment[h] == model.group_assignment[hp]:
        return model.select_aimpoint[h,a] == model.select_aimpoint[hp,a]
    return pe.Constraint.Feasible

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
#        self.aimpoints = aimpoints

    def generateMeasurementSubset(self):
        measurement_points = []
        pts_per_dim = int(self.flux_model.receiver.params["pts_per_dim"])
        if self.flux_model.receiver.params["receiver_type"] == "Flat plate":
            aim_cols = int(self.flux_model.receiver.params["aim_cols"])

        elif self.flux_model.receiver.params["receiver_type"] == "External cylindrical":
            aim_cols = 1
  
        aim_rows = int(self.flux_model.receiver.params["aim_rows"])
        row_spacing = pts_per_dim // aim_cols
        col_spacing = pts_per_dim // aim_rows
        first_row = row_spacing // 2
        first_col = col_spacing // 2
        for r in range(aim_rows):
            for c in range(aim_cols):
                pt = (first_row + (r*row_spacing)) * pts_per_dim + (first_col + (c*col_spacing)) + 1
                measurement_points.append(pt)
        self.model.check_measurement_points = pe.Set(initialize = measurement_points)


    def generateSets(self):
        """
        Generate sets for the pyomomodel object.
        """
        self.num_measurement_points = self.flux_model.receiver.x.size
        self.model.measurement_points = pe.Set(
                initialize = range(1,self.num_measurement_points+1)
                )  #M
        
        self.generateMeasurementSubset()
        
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
       
    # def getSectionFluxFraction(self):    THIS CODE GENERATES COLUMN SUM OF EACH SECTION IN PARALLEL 
        ##THE FUNCTION IN FLUX MODEL NEEDS TO BE ADJUSTED 
    #     self.flux_sum_cols = {}
    #     print(self.params["section_id"])
    #     self.flux_sum_cols[self.params["section_id"]] = self.flux_model.SectionFluxMap(self.params["section_id"])
        
    # def col_sum(self):
    #     self.getSectionFluxFraction()
    #     col_sums = []  # Stores sum as: 1st list is sum of 1st column of each section
    #     print(self.flux_sum_cols)
    #     for ncol in range(self.flux_model.receiver.params["pts_per_dim"]):
    #         col_sum_each_column = []   #Stores sum of one specific column from each section
    #         print("goes into the loop")
    #         for s in range(self.params["num_sections"]):
    #             col_sum_each_column.append((self.flux_sum_cols[s])[ncol])
    #         col_sums.append(col_sum_each_column)
    #     print(col_sums)
        
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
        if self.params["num_sections"] == 1:
            fraction = 1.0
        else: 
            if self.flux_model.receiver.params["receiver_type"] == "External cylindrical":
                fraction = self.flux_model.fraction_maps[self.params["section_id"]].flatten()
            else:   
                fraction = self.flux_model.field.section_flux[self.params["section_id"]]
        for i in range(self.num_measurement_points):
            if self.flux_model.receiver.params["receiver_type"] == "External cylindrical":
                flux_lbs[i+1] = self.flux_model.receiver.flux_lower_limits[i] * fraction[i]
                flux_ubs[i+1] = self.flux_model.receiver.flux_upper_limits[i] * fraction[i]
            else:
                flux_lbs[i+1] = self.flux_model.receiver.flux_lower_limits[i] * fraction
                flux_ubs[i+1] = self.flux_model.receiver.flux_upper_limits[i] * fraction
            # flux_ub_no_frac[i+1] = self.flux_model.receiver.flux_upper_limits[i]
            # flux_lbs[i+1] = 0.0 * fraction 
            # flux_ubs[i+1] = 200 * fraction      
            surface_area[i+1] = self.flux_model.receiver.surface_area[i]
        
        #if specifying flux maps from a file, do so.
        try: 
            if self.params["flux_from_file"]:
                flux = self.getFluxFromFile()
                return flux, flux_lbs, flux_ubs, surface_area
        except KeyError:
            pass
        flux = {}
        # ofile = open(str(self.params["section_id"])+"_flagged_helios.csv",'w')
        # ofile.write("Section_id,Helio_ID,m_point,flux_at_m,flux_ub_with_fraction,flux_ub_no_fraction\n")
        
        for h in self.model.heliostats:
            center_idx = self.flux_model.receiver.num_aimpoints//2
            h_map = self.flux_model.ShiftImage_GenerateSingleHeliostatFluxMaps(h,solar_vector,approx)
            for m in self.model.measurement_points:
                # if (h_map[center_idx][m-1]) > flux_ubs[m]:
                #     print("Heliostat: ", h," exceeds flux limit")
                #     ofile.write(str(self.params["section_id"])+","+str(h-1)+","+str(m-1)+","+str(h_map[center_idx][m-1])+","+str(flux_ubs[m])+","+str(flux_ub_no_frac[m])+"\n")
                if self.flux_model.receiver.params["receiver_type"] == "Flat plate":
                    for a in self.model.aimpoints:
                        #print("Flat a: ", a)
                        flux[h,m,a] = h_map[a-1][m-1]
                elif self.flux_model.receiver.params["receiver_type"] == "External cylindrical":
                    for a in self.model.aimpoints:
                        flux[h,m,a] = h_map[a-1][m-1]
        return flux, flux_lbs, flux_ubs, surface_area
    
    
    def getFluxFromFile(self):
        import csv
        f = csv.reader(open(self.params["flux_filename"],'r'))
        next(f)
        flux = {}
        for line in f:
            if len(line) > 3:
                flux[int(line[0]),int(line[1]),int(line[2])] = float(line[3])
        return flux

        
    def generateParameters(self):
        flux, flux_lbs, flux_ubs, surface_area = self.getFluxParameters()
        self.model.flux = pe.Param(self.model.heliostats, self.model.measurement_points, self.model.aimpoints, initialize = flux)
        self.model.flux_lbs = pe.Param(self.model.measurement_points, initialize = flux_lbs)
        self.model.flux_ubs = pe.Param(self.model.measurement_points, initialize = flux_ubs) 
        self.model.flux_diff = pe.Param(self.model.measurement_points, self.model.measurement_points, mutable=True, initialize = 1e9)
        self.model.surface_area = pe.Param(self.model.measurement_points, initialize = surface_area)
        if self.params["ordered_defocus"]:
            self.getHeliostatOrdering(self.params["order_method"])
        self.model.flux_constraint_limit = self.params["flux_constraint_limit"]


    def printFluxToFile(self,filename):
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
        if self.params["num_sections"] == 1:  #relax differential for subproblems
            self.model.flux_diff_con = pe.Constraint(self.model.measurement_points * self.model.measurement_points, rule=fluxDiffRule)
        if self.params["ordered_defocus"]:
            self.model.ordered_defocusing_con = pe.Constraint(self.model.heliostats * self.model.heliostats, rule=orderedDefocusingRule)
        if self.flux_model.settings["heliostat_group_size"] > 1: 
            self.model.ordered_defocusing_con = pe.Constraint(self.model.heliostats * self.model.heliostats * self.model.aimpoints, rule=groupingRule)
            print("group cons made")

    def createFullProblem(self):
        self.generateSets()
        self.generateParameters()
        self.generateVariables()
        self.setObjective()
        self.genConstraintsBinOnly()
        #print(self.model.select_aimpoint.value)


    def optimize(self, mipgap=0.01, timelimit=300, solver='cbc', tee=False, keepfiles=False, warmstart=False): 
        if warmstart: 
            import opt_heuristic
            heuristic = opt_heuristic.AcceptRejectHeuristic(3)
            heuristic.getIFS(self.model)
        opt = pe.SolverFactory(solver)
        if solver == "cbc":
            opt.options["ratioGap"] = mipgap
            opt.options["seconds"] = timelimit
        elif solver == "glpk":
            opt.options["mipgap"] = mipgap
            opt.options["tmlim"] = timelimit
            opt.options["cuts"] = 1
        elif solver == "cplex":
            opt.options["mipgap"] = mipgap
            opt.options["timelimit"] = timelimit
        else:
            raise Exception("invalid solver.")
        self.opt_results = opt.solve(self.model, tee=tee, keepfiles=keepfiles, warmstart=warmstart, load_solutions=False)
        self.gap = self.opt_results.solution[0].gap
        self.model.solutions.load_from(self.opt_results)
        
            
    def processOutputs(self):
        self.flux_map = numpy.zeros(self.num_measurement_points,dtype=float)
        for m in range(1,self.num_measurement_points+1):
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
        self.contribution_by_heliostat = numpy.zeros([self.flux_model.field.x.size, self.num_measurement_points])
        self.num_defocused_heliostats = 0
        
        for h in self.model.heliostats:
            if self.model.defocus[h].value > 0.5:
                self.num_defocused_heliostats += 1
            else: 
                for a in self.model.aimpoints:
                    if self.model.select_aimpoint[h,a].value > 0.5:
                        self.aimpoint_select_map[h-1] = a
                        for m in range(self.num_measurement_points):
                            self.contribution_by_heliostat[h,m] = self.model.flux[h,m+1,a]
                        break
        if self.gap == None:
            ub = 1.0e10  
            # TODO parse the results text to obtain the gap when pyomo times out
        else:
            ub = pe.value(self.model.OBJ) * (self.gap + 1)
        d = {"flux_map":self.flux_map, 
             "aimpoint_select_map":self.aimpoint_select_map,
             "obj_value":pe.value(self.model.OBJ), 
             "num_defocused_heliostats":self.num_defocused_heliostats,
             "upper_bound":ub,
             "contribution_by_heliostat":self.contribution_by_heliostat,
             "section_id":self.params["section_id"],
             "utilization": 1.0 - (float(self.num_defocused_heliostats) / self.num_heliostats)}
        return d


    def plotOutputs(self):
        import plotting
        plotting.plot_optimal_flux_heatmap(self,"fluxmap")
        plotting.plot_optimal_aimpoint_allocation(self,"aimpoint_map")
        plotting.plot_optimal_aimpoint_guide(self,"aimpoint_select_map")


    def printOutput(self):
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
