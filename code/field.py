# -*- coding: utf-8 -*-
"""
HALOS Field Module

Reads field from file or uses SolarPilot to generate field 
"""

import numpy
import pandas

class Field(object):
    """
    This object consists of coordinates for all the mirrors in the field.
    For now, it's just a collection of 3 arrays containing x-, y-, and 
    z-coordinates, but may include mirror model information as well. 
    """
    def __init__(self,filenames=None, params={"mirror_area":100}, use_sp_field = False):
        if use_sp_field:
            self.GetFieldFromSP(filenames, params)
        elif filenames == None:
            self.x = []
            self.y = []
            self.z = []
            self.coords = numpy.array([self.x,self.y,self.z]).transpose()
            self.num_heliostats = 0
            self.mirror_area = params.get("mirror_area")
            self.eff = []
            self.rej_x = []
            self.rej_y = []
            self.rej_z = []
            self.rej_coords = numpy.array([self.rej_x, self.rej_y, self.rej_z]).transpose()
            self.utilization_by_section = []
            self.min_utilization_by_section = []
            self.distance = []
            self.polar_angles = []
            self.rej_distance = []
            self.rej_polar_angles = []
            self.rej_annual_power = []
        else: 
            filename = filenames["field_filename"]
            self.GetFieldFromFile(filename,params)
            self.SetMirrorArea(params.get("mirror_area"))
        
    def GetX(self):
        return self.x
    def GetY(self):
        return self.y
    def GetZ(self):
        return self.z
    def GetCoords(self):
        return self.coords
    def GetEfficiencies(self):
        return self.eff
    def GetNumHeliostats(self):
        return self.num_heliostats
    def GetMirrorArea(self):
        return self.mirror_area
    
    def SetMirrorArea(self, area):
        self.mirror_area = area
    
    def GenerateCoords(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
        self.coords = numpy.array([self.x,self.y,self.z]).transpose()
        self.num_heliostats = len(self.x)
        
    def GetFieldFromSP(self, filenames, params):
        """
        Use SolarPILOT to get field

        Parameters
        ----------
        filenames : Paths to input files
        params : Dict
            Contains parameters for field subsection

        Returns
        -------
        None. populates self.coords

        """
        import sp_module
        sp_field = sp_module.SP_Field(filenames)
        df = sp_field.generate_field()
        self.x = df["x_location"].values
        self.y = df["y_location"].values
        self.z = df["z_location"].values
        self.coords = numpy.array([self.x,self.y,self.z]).transpose()
        self.num_heliostats = len(self.x)
        self.eff = df["efficiency"].values
        try:
            self.getPolarAngles()
            if params["num_sections"] > 1:
                if params["section_method"] == "angle":
                    self.getSectionsByPolarAngle(params["num_sections"])
                elif params["section_method"] == "distance":
                    self.getSectionsByDistance(params["num_sections"])
            else:
                self.getSectionsByPolarAngle(params["num_sections"])
        except KeyError:
            print("Warning: no sections mentioned in field parameters.")
            
    def GetFieldFromFile(self, filename, params):
        """
        Reads field CSV file for coordinates

        Parameters
        ----------
        filename : Field Filename
            path to file containing solar field heliostat positions.
        params : Dict
            Contains parameters for field subsection

        Returns
        -------
        None. populates self.coords

        """
        df = pandas.read_csv(filename)
        if params.get("hold_sp_rejects"):
            self.rejected_df = df[df["in_field"] == 0]
            df = df[df["in_field"] == 1]
        self.x = df["Pos-x"].values
        self.y = df["Pos-y"].values
        self.z = df["Pos-z"].values
        self.coords = numpy.array([self.x,self.y,self.z]).transpose()
        self.num_heliostats = len(self.x)
        try: 
            self.eff = df["Total eff"].values
        except KeyError:
            self.eff = numpy.ones_like(self.x)  
            # TODO load the annual efficiency from SolarPILOT when using a
            # solar field from file and SolarPILOT for flux calculation
        try:
            self.getPolarAngles()
            if params["num_sections"] > 1:
                if params["section_method"] == "angle":
                    self.getSectionsByPolarAngle(params["num_sections"])
                elif params["section_method"] == "distance":
                    self.getSectionsByDistance(params["num_sections"])
            else:
                self.getSectionsByPolarAngle(1)
        except KeyError:
            print("Warning: no sections mentioned in field parameters.")
                    
        
    def getPolarAngles(self, rejected = False):
        """
        obtains distance (m) and polar angle coordinate (rad) of each heliostat
        
        Parameters
        ----------
        rejected : Bool
            calculate the polar angle of rejected heliostats if True
    
        Returns
        -------
        None
        """
        #calcualate theta via the arc-cosine, using the x and y coordinates
        self.distance = numpy.sqrt(self.y**2+self.x**2)
        self.polar_angles = numpy.arccos(self.x / self.distance)
        #if y < 0, result is 2pi - theta
        self.polar_angles = (
            (self.y >= 0) *  self.polar_angles + 
            (self.y < 0 ) * (2 * numpy.pi - self.polar_angles) 
        )
        if rejected:
            self.rej_distance = numpy.sqrt(self.y ** 2 + self.x ** 2)
            self.rej_polar_angles = numpy.arccos(self.x / self.distance)
            # if y < 0, result is 2pi - theta
            self.rej_polar_angles = (
                    (self.y >= 0) * self.polar_angles +
                    (self.y < 0) * (2 * numpy.pi - self.polar_angles)
            )

    def getSectionsByPolarAngle(self, num_sections):
        """
        Subdivides the field into sections by sorting the field by angle w.r.t.
        due east.
        
        Parameters
        ----------
        num_sections : number of sections in which the solar field is 
                        subdivided (int)
    
        Returns
        -------
        None; populates self.section_id
        """
        self.num_sections = num_sections
        heliostats_per_section = self.num_heliostats // num_sections
        if self.num_heliostats % num_sections > 0 and num_sections < heliostats_per_section: 
            heliostats_per_section += 1
        sorting = self.polar_angles.argsort()
        num_helios_by_section = [heliostats_per_section]*self.num_sections
        num_helios_by_section[-1] -= (heliostats_per_section * self.num_sections - self.num_heliostats) 
        self.helios_by_section = []
        for section_idx in range(num_sections):
            self.helios_by_section.append( list(sorted([sorting[section_idx*heliostats_per_section+idx] 
            for idx in range(num_helios_by_section[section_idx])])))
        self.section_flux = [
            sum([self.eff[idx] for idx in self.helios_by_section[s]])
            / sum(self.eff) for s in range(self.num_sections)]
        self.min_angles = [min([self.polar_angles[idx] for idx in self.helios_by_section[s]])
                           for s in range(self.num_sections)]


    def getSectionsForRejectedHeliostats(self):
        """
        Assigns a section to each rejected heliostat.
        """
        self.rej_x = self.rejected_df["Pos-x"].values
        self.rej_y = self.rejected_df["Pos-y"].values
        self.rej_z = self.rejected_df["Pos-z"].values
        self.rej_coords = numpy.array([self.rej_x, self.rej_y, self.rej_z]).transpose()
        self.rej_annual_power = self.rejected_df["annual_power"].values
        self.rej_section_ids = numpy.zeros(self.rej_coords.size, dtype=int)
        self.getPolarAngles(rejected=True)
        for ridx in range(self.rej_x.size):
            self.rej_section_ids[ridx] = self.getRejSectionID(ridx)


    def getRejSectionID(self, ridx):
        """
        assigns a section ID for a rejected heliostat, using its polar angle as input.

        Parameters
        ----------
        ridx : Int
            rejected heliostat index
        Returns
        -------
        section_idx : Int
            section identifier
        """
        for section_idx in range(1,self.num_sections-1):
            if self.rej_polar_angles[ridx] <= self.min_angles[section_idx+1]:
                return section_idx
        return self.num_sections-1


    def getUtilizationStats(self, case_name, periods):
        """
        gets utilization stats by obtaining period-specific utilization by section from flat files, then averaging
        them over time.

        Parameters
        ----------
        case_name -- case identifier

        Returns
        -------
        None: assigns self.utilization_by_section
        """
        self.utilization_by_section = numpy.zeros(self.num_sections, dtype=float)
        self.min_utilization_by_section = numpy.ones(self.num_sections, dtype=float)
        for period in periods:
            filename = case_name+str(period)+"_utilization.csv"
            fin = open(filename,'r')
            lines = fin.readlines()
            sline = lines[0].split(",")
            for i in range(self.num_sections):
                self.utilization_by_section[i] += float(sline[i])
                self.min_utilization_by_section[i] = min(self.min_utilization_by_section[i], float(sline[i]) )
        self.utilization_by_section /= len(periods)


    def updateFieldForUtilization(self, threshold = 0.995):
        fully_utilized_sections = self.utilization_by_section > 0.9999
        sections_to_replace = self.utilization_by_section < threshold
        n_replaced = 0
        for sidx in range(self.num_sections):
            if sections_to_replace[sidx]:
                # replace max mirrors defocused
                num_helios = int(len(self.helios_by_section[sidx])*(1-self.min_utilization_by_section[sidx]))
                for idx in range(num_helios):
                    replaced = self.replaceHeliostat(fully_utilized_sections,sidx)
                    n_replaced += 1
                    if not replaced:
                        return
        return


    def replaceHeliostat(self, fully_utilized_sections, sidx):
        #Find existing heliostat to replace, rejected heliostat to introduce
        #existing_idx = self.getMinEffHelioInSection(sidx)
        existing_idx = self.getMaxDistHelioInSection(sidx)
        replace_idx = self.findMaxPowerReplaceableHelio()
        if replace_idx == -1:
            return False
        #update coordinates
        self.x[existing_idx] = self.rej_x[replace_idx]
        self.y[existing_idx] = self.rej_y[replace_idx]
        self.z[existing_idx] = self.rej_z[replace_idx]
        self.coords[existing_idx] = self.coords[existing_idx]
        self.helios_by_section[sidx].remove(existing_idx)
        self.helios_by_section[self.rej_section_ids[replace_idx]].append(existing_idx)
        #update efficiency, dist, annual power so heliostats aren't replaced again
        self.eff[existing_idx] = 1.0
        self.distance[existing_idx] = 0.0
        self.rej_annual_power[replace_idx] = 0.0
        return True

    def getMaxDistHelioInSection(self, sidx):
        max_dist = 0.0
        max_idx = -1
        for hidx in self.helios_by_section[sidx]:
            if self.distance[hidx] > max_dist:
                max_idx = hidx
                max_dist = self.distance[hidx]
        return max_idx

    def getMinEffHelioInSection(self, sidx):
        min_eff = 1.0
        min_idx = -1
        for hidx in self.helios_by_section[sidx]:
            if self.eff[hidx] < min_eff:
                min_idx = hidx
                min_eff = self.eff[hidx]
        return min_idx


    def findMaxPowerReplaceableHelio(self):
        replaced = False
        max_power = 1.0 #don't count any rejected heliostats whose annual power has been set to zero
        replace_idx = -1
        for ridx in range(len(self.rej_x)):
            if self.rej_annual_power[ridx] > max_power and self.utilization_by_section[self.rej_section_ids[ridx]] > 0.995:
                replace_idx = ridx
                max_power = self.rej_annual_power[ridx]
        return replace_idx


    def outputFieldToFile(self, filename):
        outfile = open(filename, 'w')
        outfile.write('id,Pos-x,Pos-y,Pos-z\n')
        for hidx in range(len(self.x)):
            outfile.write(str(hidx)+","+str(self.x[hidx])+","+str(self.y[hidx])+","+str(self.z[hidx])+"\n")
        outfile.close()


    def getSectionsByDistance(self, num_sections):
        """
        Subdivides the field into sections by distance from the receiver 
        center.  The field subdivisions are such that there is an equal number
        of heliostats in each section (subject to a difference of one if the
        number of heliostats is not divisible by the number of sections).
        
        Parameters
        ----------
        num_sections : number of sections in which the solar field is 
                        subdivided (int)
    
        Returns
        -------
        None; populates self.section_id
        """
        self.num_sections = num_sections
        heliostats_per_section = self.num_heliostats // self.num_sections
        if self.num_heliostats % self.num_sections > 0: 
            heliostats_per_section += 1
        sorting = self.distance.argsort()
        num_helios_by_section = [heliostats_per_section]*self.num_sections
        num_helios_by_section[-1] -= (heliostats_per_section * self.num_sections - self.num_heliostats) 
        self.helios_by_section = []
        for section_idx in range(num_sections):
            self.helios_by_section.append( list(sorted([sorting[section_idx*heliostats_per_section+idx] 
            for idx in range(num_helios_by_section[section_idx])]))) 
        self.section_flux = [
            sum([self.eff[idx] for idx in self.helios_by_section[s]])
            / sum(self.eff) for s in range(self.num_sections)]


    def updateSPfield(self, sp_field_filename, case_filename, case_name, new_field_filename, num_sections):
        import inputs
        filenames = {"field_filename": sp_field_filename}
        params = {}
        params["num_sections"] = num_sections
        params["section_method"] = "angle"
        params["mirror_area"] = 100
        params["hold_sp_rejects"] = True
        field = Field(filenames, params, use_sp_field=False)
        field.getSectionsForRejectedHeliostats()
        import annual_layout_optimize
        periods = annual_layout_optimize.getHourIDs(case_filename)
        periods = [5170, 5171]
        field.getUtilizationStats(case_name, periods)
        field.updateFieldForUtilization(0.995)
        field.outputFieldToFile(new_field_filename)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import numpy as np
    # fname = "./../solarpilot_cases/radial-daggett-50.csv"
    import inputs

    # case_filename = "./../case_inputs/flat_50_ca_case.csv"
    filenames = {"field_filename": "radial-250-daggett-layout.csv"}
    case_filename = "./../case_inputs/radial_250_ca_case.csv"
    case_name = "radial-250-daggett"
    params = {}
    params["num_sections"] = 16
    params["section_method"] = "angle"
    params["mirror_area"] = 100
    params["hold_sp_rejects"] = True
    field = Field(filenames, params, use_sp_field = False)
    old_field_filename = case_name + "_old_field.csv"
    field.outputFieldToFile(old_field_filename)
    field.getSectionsForRejectedHeliostats()
    # print(field.rej_section_ids)
    import annual_layout_optimize
    case_filename = "./../case_inputs/radial_250_ca_case.csv"
    case_name = "radial-250-daggett"
    new_field_filename = case_name + "_rev_fielddfsaf.csv"
    periods = annual_layout_optimize.getHourIDs(case_filename)
    field.getUtilizationStats(case_name, periods)
    #field.updateFieldForUtilization(0.995)
    # field.outputFieldToFile(new_field_filename)
    # print(field.eff)
    # print(field.utilization_by_section)
    # print(field.min_utilization_by_section)
    # print(field.min_utilization_by_section.argsort())
    if True:
        replaced_idxs = [2715, 2591,2592,2593,2594,2595,2511,2514,2512,2513,2407,2408,2409,2410,2412,2413,2414,2415,2416,
            2411,2294,2295,2296,2297,2405,2404,2406,2283,2281,2282,2284,2285,2286,2288,2289,2290,2291,2287,2156,2153,2154,
            2155,2157,2159,2160,2161,2162,2163,2713,2714,2588,2589,2590,2504,2500,2501,2502,2505,2506,2507]
        x = [0.,0.]
        y = [0.,0.]
        col = [0.,1.0]
        for idx in range(len(field.helios_by_section)):
            #r = np.random.random_sample()
            for h in field.helios_by_section[idx]:
                x.append(field.x.flatten()[h])
                y.append(field.y.flatten()[h])
                if h in replaced_idxs:
                    col.append(0.0)
                else:
                    col.append(1.0)
        plt.scatter(x,y,s=6,c = col,cmap="Set1")
        plt.ylim([-750.,1250.])
        plt.xlim([-1000.,1000.])
        plt.savefig("old.png")
        plt.cla()
        plt.clf()

