# -*- coding: utf-8 -*-
"""
This module includes heuristics that can generate an initial feasible solution
to the aimpoint strategy optimization model.
"""

import WELL512

import numpy

class Heuristic(object):
    """generates an initial feasible solution to an existing aimpoint strategy
    optimization model."""
    def __init__(self,name):
        self.name = name
        
    def getName():
        return self.name
        
    def getIFS(self, model):
        raise ValueError("Inherited classes must overwrite 'getIFS'")
        
    def updateVals(self, model, aimpoint_select, flux):
        helio_list = list(model.heliostats)
        for h in range(len(helio_list)):
            if aimpoint_select[h] > 0:
                model.select_aimpoint[helio_list[h], aimpoint_select[h]] = 1
            else: 
                #defocus heliostat
                model.defocus[helio_list[h]] = 1
    
    def attemptToAim(self, model, flux, a, h):
        added_flux = numpy.array([model.flux[h,m,a] for m in model.measurement_points])
        if any(flux+added_flux > model.flux_ubs):
            return flux, -1, False
        return flux+added_flux, a, True
        

class AcceptRejectHeuristic(Heuristic):
    def __init__(self, num_tries):
        Heuristic.__init__(self,"accept-reject")
        self.num_tries = num_tries
        self.gen = WELL512.WELL512("rngstates.csv")
    
    def getIFS(self, model):
        helio_list = list(model.heliostats)
        aim_list = list(model.aimpoints)
        aimpoint_select = numpy.zeros(len(helio_list), dtype=float)
        flux = numpy.zeros(len(model.measurement_points), dtype=float)
        for h in range(len(helio_list)):
            success = False
            for i in range(self.num_tries):
                if success: break
                a = aim_list[int(self.gen.getVariate()*len(model.aimpoints))]
                flux, aim_idx, success = self.attemptToAim(model, flux, a, helio_list[h])
            if success:
                aimpoint_select[h] = aim_idx
        self.updateVals(model, aimpoint_select, flux)
    
    
class SpreadHeuristic(Heuristic):
    def __init__(self):
        Heuristic.__init__(self,"spread")
    
    def getIFS(self, model):
        helio_list = list(model.heliostats)
        aim_list = list(model.aimpoints)
        aimpoint_select = numpy.zeros(len(helio_list), dtype=int)
        flux = numpy.zeros(len(model.measurement_points), dtype=float)
        aim_idx = max(0, (len(aim_list) - len(helio_list)) // 2)
        for h in range(len(helio_list)):
            print("attempting heliostat ",h)
            flux, a, success = self.attemptToAim(model, flux, aim_list[aim_idx], helio_list[h])
            for i in range(len(aim_list)):
                if success: break
                aim_idx += 1
                aim_idx = aim_idx % len(aim_list)
                print("attempting aimpoint ",aim_idx)    
                flux, a, success = self.attemptToAim(model, flux, aim_list[aim_idx], helio_list[h])
            if success:
                aimpoint_select[h] = a
                aim_idx += 1
                aim_idx = aim_idx % len(aim_list)
                print("heliostat #",h," aimed at ",a)
            else: 
                ## 
                print("heuristic terminated at heliostat #",h)
                break
        self.updateVals(model, aimpoint_select, flux)

