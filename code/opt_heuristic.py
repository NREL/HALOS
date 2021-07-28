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
                #a = random.choice(aim_list)
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
                flux, obj_val, a, success = self.attemptToAim(model, flux, aim_list[aim_idx], helio_list[h])
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
        print('obj val: ',obj_val)


class CenterOut(Heuristic):
    def __init__(self,num_tries):
        Heuristic.__init__(self,"center_out")
        self.num_tries = num_tries

    def getIFS(self,model):
        t_start = time.time()
        helio_list = list(model.heliostats)
        aim_list = list(model.aimpoints)
        aimpoint_select = numpy.zeros(len(helio_list), dtype=int)
        flux = numpy.zeros(len(model.measurement_points), dtype=float)
        a_cent = math.ceil(len(aim_list) / 2) # first try the center aimpoint. if odd num aimpts, starts at center. If even aimpts, starts at below center
        #self.num_tries = len(aim_list) # if want to try h at each aimpoint before defocusing

        for h in range(len(helio_list)):
            success = False
            for i in range(self.num_tries):
                if success: break
                if i % 2 != 0:
                    a = a_cent - math.floor(i/2) # 1st round -> a_1, 3rd round -> 1 below a_1, 5th round -> 2 below a_1 and so on
                else:
                    a = a_cent + i//2  # 2nd round -> 1 above a_1, 4th round -> 2 above a_1, and so on
                flux, obj_val, aim_idx, success = self.attemptToAim(model, flux, a, helio_list[h])
            if success:
                aimpoint_select[h] = aim_idx
        self.updateVals(model, aimpoint_select, flux)
        print('CenterOut heuristic ran')
        t_end = time.time()
        print('time taken for initital heuristic: ',t_end-t_start)
        print('obj val: ',obj_val)
        return(obj_val,t_end-t_start)


class WeightedCenter(Heuristic):
    def __init__(self,num_tries):
        Heuristic.__init__(self,"fartherin_centerout")
        self.gen = WELL512.WELL512("rngstates.csv")
        self.num_tries = num_tries

    def getIFS(self,model):
        t_start = time.time()
        helio_list = list(model.heliostats)
        aim_list = list(model.aimpoints)
        flux = numpy.zeros(len(model.measurement_points), dtype=float)
        aimpoint_select = numpy.zeros(len(helio_list), dtype=float)
        
        center_idx = math.ceil(len(aim_list) / 2) 
        orig_aimlist = aim_list.copy()
        # add lots more of the center_idx to aim list, then the ones above and below center if aim_list long enough 
        if len(orig_aimlist) >= 3:
            for i in range(len(aim_list)*10):
                aim_list.append(center_idx)
        if len(orig_aimlist) >= 4:
            for x in range(len(orig_aimlist)*5):
                aim_list.append(center_idx+1)
        if len(orig_aimlist) >= 5:
            for y in range(len(orig_aimlist)*5):
                aim_list.append(center_idx-1)
        print(aim_list)
        for h in range(len(helio_list)):
            success = False
            for n in range(self.num_tries):
                #a = random.choice(aim_list) # if use random.choice instead off WELL512 - seems to have similar times and all
                a = aim_list[int(self.gen.getVariate()*len(model.aimpoints))]
                flux, obj_val, aim_idx, success = self.attemptToAim(model, flux, a, helio_list[h])
                if success: break # put this line at end instead of beginning
            if success:
                aimpoint_select[h] = aim_idx
        self.updateVals(model, aimpoint_select, flux)
        t_end = time.time()
        print('time taken for initital heuristic: ',t_end-t_start)
        print('Weighted Center heuristic ran')
        print('obj val: ',obj_val)
        return(obj_val,t_end-t_start)

class SpecificthenRandom(Heuristic):
    def __init__(self,num_tries):
        Heuristic.__init__(self,"fartherin_centerout")
        self.gen = WELL512.WELL512("rngstates.csv")
        self.num_tries = num_tries

    def getIFS(self,model):
        t_start = time.time()
        helio_list = list(model.heliostats)
        aim_list = list(model.aimpoints)
        flux = numpy.zeros(len(model.measurement_points), dtype=float)
        aimpoint_select = numpy.zeros(len(helio_list), dtype=float)
        defocused_helio = []
        
        center_idx = math.ceil(len(aim_list) / 2)
        for h in range(len(helio_list)):
            success = False
            for n in range(self.num_tries):
                if n == 0:
                    # could also make this weighted random generation
                    a = aim_list[center_idx]
                else:
                    #a = random.choice(aim_list) # if use random.choice instead off WELL512 - seems to have similar times and all
                    a = aim_list[int(self.gen.getVariate()*len(model.aimpoints))]
                flux, obj_val, aim_idx, success = self.attemptToAim(model, flux, a, helio_list[h])
                if success: break # put this line at end instead of beginning
            if success:
                aimpoint_select[h] = aim_idx
                if a == aim_list[center_idx]:
                    print('heliostat aims at ',a)
            if not success:
                defocused_helio.append(defocused_helio)
        self.updateVals(model, aimpoint_select, flux)
        t_end = time.time()
        print('time taken for initital heuristic: ',t_end-t_start)
        print('Middle then Random heuristic ran')
        print('obj val: ',obj_val)
        return(obj_val,t_end-t_start)

class FluxSize(Heuristic):
    def __init__(self,num_tries):
        Heuristic.__init__(self,"center_out")
        self.num_tries = num_tries
        self.gen = WELL512.WELL512("rngstates.csv")

    def getIFS(self,model):
        t_start = time.time()
        helio_list = list(model.heliostats)
        aim_list = list(model.aimpoints)
        aimpoint_select = numpy.zeros(len(helio_list), dtype=int)
        flux = numpy.zeros(len(model.measurement_points), dtype=float)
        a_cent = math.ceil(len(aim_list) / 2) # first try the center aimpoint. if odd num aimpts, starts at center. If even aimpts, starts at below center
        num_zeross = []
        defocused_helio = []

        for h in range(len(helio_list)):
            success = False
            flux = numpy.array([model.flux[helio_list[h],m,a_cent] for m in model.measurement_points])
            num_zeros = numpy.count_nonzero(flux == 0)
            #size = len(model.measurement_points)*len(model.measurement_points) - num_zeros
            #flux_size.append(size) # unnecessary to use flux_sizes -- can just go by num zeros for now
            num_zeross.append(num_zeros)

        for h in range(len(helio_list)):
            for n in range(self.num_tries):
                list_to_change = aim_list.copy()
                if n == 0:
                    # if in smallest 1/8 size of the section of heliostats, random gen weighted at top aimpt
                    if num_zeross[h] >= sum(num_zeross)*7/8:
                        for i in range(int(len(aim_list)*0.5)):
                            list_to_change.append(aim_list[-1])
                    elif num_zeross[h] < sum(num_zeross)*7/8 and num_zeross[h] >= sum(num_zeross)*6/8:
                        for x in range(int(len(aim_list)*0.5)):
                            list_to_change.append(aim_list[0])
                    elif num_zeross[h] < sum(num_zeross)*6/8 and num_zeross[h] >= sum(num_zeross)*5/8:
                        for x in range(int(len(aim_list)*0.5)):
                            list_to_change.append(aim_list[-2])
                    elif num_zeross[h] < sum(num_zeross)*5/8 and num_zeross[h] >= sum(num_zeross)*4/8:
                        for x in range(int(len(aim_list)*0.5)):
                            list_to_change.append(aim_list[1])
                    elif num_zeross[h] < sum(num_zeross)*4/8 and num_zeross[h] >= sum(num_zeross)*3/8:
                        list_to_change = list_to_change
                    elif num_zeross[h] < sum(num_zeross)*4/8 and num_zeross[h] >= sum(num_zeross)*3/8:
                        for x in range(int(len(aim_list)*0.5)):
                            list_to_change.append(a_cent+1)
                    elif num_zeross[h] < sum(num_zeross)*3/8 and num_zeross[h] >= sum(num_zeross)*2/8:
                        for x in range(int(len(aim_list)*0.5)):
                            list_to_change.append(a_cent-1)
                    else:
                        for x in range(int(len(aim_list)*0.5)):
                            list_to_change.append(a_cent-1)
                '''
                if n == 1:
                                        # if in smallest 1/8 size of the section of heliostats, random gen weighted at top aimpt
                    if num_zeross[h] >= sum(num_zeross)*7/8:
                        for i in range(int(len(aim_list)*0.5)):
                            list_to_change.append(aim_list[-2])
                    # if in largest 3/4 size of the section of heliostats, random gen weighted at center
                    elif num_zeross[h] < sum(num_zeross)*6/8:
                        for i in range(int(len(aim_list)*0.5)):
                            list_to_change.append(a_cent+1)
                    else:
                        for x in range(int(len(aim_list)(0.5))):
                            list_to_change.append(aim_list[1])
                        # heliostats in 2nd smallest 1/8 of field point at bottom aimpoint
                '''
                #a = random.choice(list_to_change)
                #print('aim idx: ',int(self.gen.getVariate()*len(list_to_change))) # all updated with using len(list_to_change)
                a = list_to_change[int(self.gen.getVariate()*len(list_to_change))]
                #print('chosen a is ',str(a),' for num zeros ',num_zeross[h])
                flux, obj_val, aim_idx, success = self.attemptToAim(model, flux, a, helio_list[h])
                if success: break
            if success:
                aimpoint_select[h] = aim_idx
            else:
                defocused_helio.append(h)
        self.updateVals(model, aimpoint_select, flux)
        t_end = time.time()
        print('time taken for initital heuristic: ',t_end-t_start)
        print('FluxSize heuristic ran')
        print('obj val: ',obj_val)
        #print('num zeros',num_zeross)
        print('defocused heliostats: ',len(defocused_helio))
        return(obj_val,t_end-t_start,len(defocused_helio))