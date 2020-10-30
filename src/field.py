# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 11:14:15 2020

@author: azolan
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
        df = pandas.read_csv(filename)
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
                    
        
    def getPolarAngles(self):
        """
        obtains distance (m) and polar angle coordinate (rad) of each heliostat
        
        Parameters
        ----------
        None
    
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
        #TODO: Make a new method, GetSectionFractions(), which provides measurement-point-specific fractions.  For a cylinder these will be identical by vertical colun, and for a flat plate these will be identical for all measurement points.

        # a = (20,20)
        # fraction_maps = []
        # sum_cols = []           #Stores lists - one list contains sum of 1 column of all sections
        # for ncol in range(3):
        #     sums = []           #Stores sum of one specific column of each section
        #     for s in range(self.num_sections):
        #         fraction_map = numpy.zeros(a)+self.section_flux[s]
        #         fraction_map = pandas.DataFrame(fraction_map)
        #         sum_col = sum(fraction_map[ncol])
        #         sums.append(sum_col)
        #     sum_cols.append(sums)       
        
        # for s in range(self.num_sections):
        #     map_fraction = pandas.DataFrame(numpy.zeros(a))
        #     for ncol in range(3): 
        #         map_fraction[ncol] = map_fraction + (sum_cols[ncol][s]/sum(sum_cols[ncol]))
        #         #print(map_fraction)
        #     map_fraction = numpy.array(map_fraction)
        #     fraction_maps.append(map_fraction)
        
        
    # def getSectionFluxColumn(self, receiver):
    #     self.getSectionsByPolarAngle(params["num_sections"])
    #     for s in range(self.num_sections):
            
    #         print(fraction_map)
    #     # for ncol in receiver.params["pts_per_dim"]:
    #     #     for mrow in receiver.param["pts_per_dim"]:
    #     #         fraction_map

    
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
    filenames = inputs.readCaseFile(case_filename)
    params = {}
    params["num_sections"] = 8
    params["section_method"] = "angle"
    params["mirror_area"] = 100
    field = Field(filenames, params, use_sp_field = False)

    if True:
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

        