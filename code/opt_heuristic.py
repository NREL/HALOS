# -*- coding: utf-8 -*-
"""
This module includes heuristics that can generate an initial feasible solution
to the aimpoint strategy optimization model.
"""

import WELL512

import numpy

import math

import time

import random

class Heuristic(object):
    """generates an initial feasible solution to an existing aimpoint strategy
    optimization model."""
    def __init__(self,name):
        self.name = name
        self.obj_val = 0
        
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
            return flux,self.obj_val, -1, False
        else:
            for m in model.measurement_points:
                self.obj_val += model.flux[h,m,a]*model.surface_area[m]
            return flux+added_flux, self.obj_val, a, True


class AcceptRejectHeuristic(Heuristic):
    def __init__(self, num_tries):
        Heuristic.__init__(self,"accept-reject")
        self.num_tries = num_tries
        self.gen = WELL512.WELL512("rngstates.csv")
    
    def getIFS(self, model):
        t_start = time.time()
        helio_list = list(model.heliostats)
        aim_list = list(model.aimpoints)
        aimpoint_select = numpy.zeros(len(helio_list), dtype=float)
        flux = numpy.zeros(len(model.measurement_points), dtype=float)
        defocused_helio = []
        for h in range(len(helio_list)):
            success = False
            for i in range(self.num_tries):
                if success: break
                a = aim_list[int(self.gen.getVariate()*len(model.aimpoints))]
                flux, obj_val, aim_idx, success = self.attemptToAim(model, flux, a, helio_list[h])
            if success:
                aimpoint_select[h] = aim_idx
            else:
                defocused_helio.append(h)
        self.updateVals(model, aimpoint_select, flux)
        t_end = time.time()
        print('time taken for initital heuristic section: ',t_end-t_start)
        #print('AcceptReject heuristic ran')
        print('obj val: ',obj_val)
        return(obj_val,t_end-t_start,len(defocused_helio))
    
    
class SpreadHeuristic(Heuristic):
    def __init__(self):
        Heuristic.__init__(self,"spread")
    
    def getIFS(self, model):
        t_start = time.time()
        helio_list = list(model.heliostats)
        aim_list = list(model.aimpoints)
        aimpoint_select = numpy.zeros(len(helio_list), dtype=int)
        flux = numpy.zeros(len(model.measurement_points), dtype=float)
        aim_idx = max(0, (len(aim_list) - len(helio_list)) // 2)
        defocused_helio = []

        for h in range(len(helio_list)):
            print("attempting heliostat ",h)
            flux, obj_val, aim_idx, success = self.attemptToAim(model, flux, aim_list[aim_idx], helio_list[h])
            for i in range(len(aim_list)):
                if success: break
                aim_idx += 1
                aim_idx = aim_idx % len(aim_list)
                print("attempting aimpoint ",aim_idx)    
                flux, obj_val, aim_idx, success = self.attemptToAim(model, flux, aim_list[aim_idx], helio_list[h])
            if success:
                aimpoint_select[h] = aim_idx
                aim_idx += 1
                aim_idx = aim_idx % len(aim_list)
                print("heliostat #",h," aimed at ",aim_idx)
            else: 
                ## 
                print("heuristic terminated at heliostat #",h)
                defocused_helio.append(h)
                break
        self.updateVals(model, aimpoint_select, flux)
        print('obj val: ',obj_val)
        t_end = time.time()
        return(obj_val,t_end-t_start,len(defocused_helio))


class FluxStoreSize(Heuristic):
    def __init__(self,num_tries):
        Heuristic.__init__(self,"center_out")
        self.num_tries = num_tries
        self.gen = WELL512.WELL512("rngstates.csv")

    def resetValues(self,model):
        # to reset values before a new method is used to generate new feasible h-a pairings
        self.obj_val = 0
        flux = numpy.zeros(len(model.measurement_points), dtype=float)
        defocused_helio = []
        helio_list = list(model.heliostats)
        aimpoint_select = numpy.zeros(len(helio_list), dtype=int)
        return flux,defocused_helio,helio_list,aimpoint_select

    def calculateFluxSizes(self,model,a_cent):
        # makes list num_zeross for which each element is the number of zeros in the flux 
        # image of the corresponding heliostat when it is aimed at the center aimpoint
        # calculates std and avg of list
        helio_list = list(model.heliostats)
        num_zeross = []
        for h in range(len(helio_list)):
            fluxes = numpy.array([model.flux[helio_list[h],m,a_cent] for m in model.measurement_points])
            num_zeros = numpy.count_nonzero(fluxes == 0)
            num_zeross.append(num_zeros)
        standard_deviation = numpy.std(num_zeross, ddof=0)
        avg_nz = sum(num_zeross)/len(num_zeross)
        return num_zeross,standard_deviation,avg_nz
    
    def getUpdatedAimlist(self,model,a_cent):
        # changes aim_list to a list of length 8 by copying elements, starting from center mostly
        # or by deleting elements if list is longer than num_sxns 8
        aim_list = list(model.aimpoints)
        num_sxns = 8
        aimlist_new = []
        if len(aim_list) < num_sxns:
            aimlist_new = aim_list.copy()
            # if 4 or less num aim, cycle through aim list adding elements until reach length of 8
            if len(aim_list) <= 4:
                i = 0
                while len(aimlist_new) < num_sxns:
                    aimlist_new.append(aim_list[i])
                    i += 1
                    if i > len(aim_list)-1:
                        i = 0
            else: 
                # if 5-7 a, start from center and add aimpoints increasingly far away from center until reach length of 8
                num_from_center = 1
                i = 1
                while len(aimlist_new) < num_sxns:
                    if i == 1:
                        aimlist_new.append(a_cent)
                    else:
                        if i % 2 != 0:
                            aimlist_new.append(a_cent-num_from_center)
                            num_from_center += 1
                        else:
                            aimlist_new.append(a_cent+num_from_center)
                    i += 1
                aimlist_new.sort()
            aim_list = aimlist_new
        elif len(aim_list) > num_sxns:
            y = 2
            count = 1
            # not thoroughly tested for high numbers of aimpoints
            # would want to change to make more spread out and ensure no index errors if using cases with many aimpoints
            while len(aim_list) > num_sxns:
                if count % 2 != 0:
                    aim_list.remove(a_cent-y)
                else:
                    aim_list.remove(a_cent+y)
                    y += 2
                count += 1
                if a_cent + y > len(aim_list):
                    y = 1
        return aim_list
    
    def fluxSizeMethod(self,model,order,flux,defocused_helio,helio_list,aimpoint_select,aim_list,avg_nz,standard_deviation,num_zeross):
        # creates 8 similarly sized/distributed groups of heliostats based on flux size
        # its group determines the aimpoint the heliostat attempts on first try (larger flux size --> closer to center)
        for h in range(len(helio_list)):
            for n in range(self.num_tries):
                list_to_change = aim_list.copy()
                if n == 0:
                    # flux image size increases as go down list, so aimpoints closer to center as go down the list
                    # use weighted list rather than directly choosing a on first try because get higher output that way on preliminary testing
                    # about 5 in 6 chance choose weighted aimpoint - very high chance
                    if num_zeross[h] >= avg_nz+standard_deviation*1.5:
                        for x in range(int(len(aim_list)*5)):
                            list_to_change.append(aim_list[order[0]])
                    if num_zeross[h] >= avg_nz+standard_deviation and num_zeross[h] < avg_nz+standard_deviation*1.5:
                        for x in range(int(len(aim_list)*5)):
                            list_to_change.append(aim_list[order[1]])
                    if num_zeross[h] >= avg_nz+standard_deviation*0.4 and num_zeross[h] < avg_nz+standard_deviation:
                        for x in range(int(len(aim_list)*5)):
                            list_to_change.append(aim_list[order[2]])
                    if num_zeross[h] >= avg_nz and num_zeross[h] < avg_nz+standard_deviation*0.4:
                        for x in range(int(len(aim_list)*5)):
                            list_to_change.append(aim_list[order[3]])
                    if num_zeross[h] >= avg_nz-standard_deviation*0.4 and num_zeross[h] < avg_nz:
                        for x in range(int(len(aim_list)*5)):
                            list_to_change.append(aim_list[order[4]])
                    if num_zeross[h] >= avg_nz-standard_deviation and num_zeross[h] < avg_nz-standard_deviation*0.4:
                        for x in range(int(len(aim_list)*5)):
                            list_to_change.append(aim_list[order[5]])
                    if num_zeross[h] >= avg_nz-standard_deviation*1.5 and num_zeross[h] < avg_nz-standard_deviation:
                        for x in range(int(len(aim_list)*5)):
                            list_to_change.append(aim_list[order[6]])
                    if num_zeross[h] <= avg_nz-standard_deviation*1.5:
                        for x in range(int(len(aim_list)*5)):
                            list_to_change.append(aim_list[order[7]])

                a = list_to_change[int(self.gen.getVariate()*len(list_to_change))]
                flux, obj_val, aim_idx, success = self.attemptToAim(model, flux, a, helio_list[h])
                if success: break
            if success:
                aimpoint_select[h] = aim_idx
            else:
                defocused_helio.append(h)

        return obj_val,defocused_helio,flux,aimpoint_select

    def getIFS(self,model):
        t_start = time.time()
        ub = 0
        aim_list = list(model.aimpoints)
        a_cent = math.ceil(len(aim_list) / 2) # center pt if odd number of aimpoints, center/closer to bottom if even number of aimpoints
        num_zeross,standard_deviation,avg_nz = self.calculateFluxSizes(model,a_cent)
        aim_list = self.getUpdatedAimlist(model,a_cent)

        # aiming method 1
        flux,defocused_helio,helio_list,aimpoint_select = self.resetValues(model)
        order = [7,0,6,1,5,2,4,3]
        obj_val,defocused_helio,flux,aimpoint_select = self.fluxSizeMethod(model,order,flux,defocused_helio,helio_list,aimpoint_select,aim_list,avg_nz,standard_deviation,num_zeross)
        print('method 1 obj val: ',obj_val,' method 1 defocused heliostats: ',len(defocused_helio))
        if obj_val > ub:
            ub = obj_val
            self.updateVals(model, aimpoint_select, flux)
            chosen_defoc_helio = defocused_helio.copy()

        # method 2 - same as method 1 but start alternating pattern from bottom instead of top
        flux,defocused_helio,helio_list,aimpoint_select = self.resetValues(model)
        order = [0,7,1,6,2,5,3,4]
        obj_val,defocused_helio,flux,aimpoint_select = self.fluxSizeMethod(model,order,flux,defocused_helio,helio_list,aimpoint_select,aim_list,avg_nz,standard_deviation,num_zeross)
        print('method 2 obj val: ',obj_val,' method 2 defocused heliostats: ',len(defocused_helio))
        if obj_val > ub:
            ub = obj_val
            print('2nd solution higher than 1st')
            self.updateVals(model, aimpoint_select, flux)
            chosen_defoc_helio = defocused_helio.copy()
        else:
            print('1st solution higher than 2nd')

        # method 3 - now randomly assigning -- most helpful if had some defocused heliostats using the previous method
        # some emphasis on center for first try - 50% chance or more it's chosen
        flux,defocused_helio,helio_list,aimpoint_select = self.resetValues(model)
        for h in range(len(helio_list)):
            for n in range(self.num_tries):
                list_to_change = aim_list.copy()
                if n == 0:
                    for i in range(int(len(aim_list))):
                        list_to_change.append(a_cent)
                a = list_to_change[int(self.gen.getVariate()*len(list_to_change))]
                flux, obj_val, aim_idx, success = self.attemptToAim(model, flux, a, helio_list[h])
                if success: break
            if success:
                aimpoint_select[h] = aim_idx
            else:
                defocused_helio.append(h)

        print('method 3 obj val: ',obj_val,' method 3 defocused heliostats: ',len(defocused_helio))
        if obj_val > ub:
            ub = obj_val
            print('3rd solution highest')
            chosen_defoc_helio = defocused_helio.copy()
            self.updateVals(model, aimpoint_select, flux)

        t_end = time.time()
        return(ub,t_end-t_start,len(chosen_defoc_helio))