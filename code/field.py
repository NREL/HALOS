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
        if filenames == None:
            self.x = []
            self.y = []
            self.z = []
            self.coords = numpy.array([self.x,self.y,self.z]).transpose()
            self.num_heliostats = 0
            self.mirror_area = params["mirror_area"]
            self.eff = []
            self.rej_x = []
            self.rej_y = []
            self.rej_z = []
            self.rej_coords = numpy.array([self.rej_x, self.rej_y, self.rej_z]).transpose()
        elif use_sp_field:
            self.GetFieldFromSP(filenames,params)
        else: 
            filename = filenames["field_filename"]
            self.GetFieldFromFile(filename,params)
            self.SetMirrorArea(params["mirror_area"])
        
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
        
    def GetFieldFromSP(self,filenames, params):
        """
        Use SolarPolit to get field

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
            if params["num_sections"] > 1:
                self.getPolarAngles()
                if params["section_method"] == "angle":
                    self.getSectionsByPolarAngle(params["num_sections"])
                elif params["section_method"] == "distance":
                    self.getSectionsByDistance(params["num_sections"])
        except KeyError:
            print("Warning: no sections mentioned in field parameters.")
            
    def GetFieldFromFile(self,filename, params):
        """
        Reads field CSV file

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
        if params.get("hold_sp_rejects")    :
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
            if params["num_sections"] > 1:
                self.getPolarAngles()
                if params["section_method"] == "angle":
                    self.getSectionsByPolarAngle(params["num_sections"])
                elif params["section_method"] == "distance":
                    self.getSectionsByDistance(params["num_sections"])
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
        self.rej_section_ids = numpy.zeros(self.rej_coords.size)
        self.getPolarAngles(rejected=True)
        for ridx in range(self.rej_x.size):
            self.rej_section_ids[ridx] = self.getRejSectionID(ridx)


    def getRejSectionID(self, ridx):
        """
        assigns a section ID for a rejected heliostat, using its polar angle as input.

        Parameters
            ridx : Int
                rejected heliostat index
        Returns
            section_idx : Int
                section identifier
        """
        for section_idx in range(1,self.num_sections-1):
            if self.rej_polar_angles[ridx] <= self.min_angles[section_idx+1]:
                return section_idx
        return self.num_sections-1

    def getUtilizationStats(self, util_filename):
        pass
    
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
        
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import numpy as np
    # fname = "./../solarpilot_cases/radial-daggett-50.csv"
    import inputs

    case_filename = "./../case_inputs/flat_50_ca_case.csv"
    filenames = {"field_filename": "radial-250-daggett-layout.csv"}
    params = {}
    params["num_sections"] = 8
    params["section_method"] = "angle"
    params["mirror_area"] = 100
    params["hold_sp_rejects"] = True
    field = Field(filenames, params, use_sp_field = False)
    field.getSectionsForRejectedHeliostats()
    print(field.rej_section_ids)
    if False:
        x = []
        y = []
        col = []
        for idx in range(len(field.helios_by_section)):
            #r = np.random.random_sample()
            for h in field.helios_by_section[idx]:
                x.append(field.x.flatten()[h])
                y.append(field.y.flatten()[h])
                col.append(float(idx))
        plt.scatter(x,y,s=6,c = col,cmap="jet")
        plt.show()
        plt.cla()
        plt.clf()

