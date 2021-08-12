# -*- coding: utf-8 -*-
"""
This module includes heuristics that can generate an initial feasible solution
to the aimpoint strategy optimization model.
"""
import math
import time
import random

import numpy
import pandas as pd

import WELL512

class Heuristic(object):
    """generates an initial feasible solution to an existing aimpoint strategy
    optimization model."""
    def __init__(self,name):
        self.name = name
        
    def getName():
        return self.name
        
    def getIFS(self, model):
        raise ValueError("Inherited classes must overwrite 'getIFS'")
        
    def updateVals(self, model, aimpoint_select):
        '''
        Updates select_aimpoint with the initial heliostat-aimpoint pairings. 
        Updates defocus with the heliostats that did not successfully focus.

        Parameters
        ------------
        model -- pyomo model that is defined in optimize_aimpoint.AimpointOptimizer
        (attributes include measurement_points, aimpoints, heliostats, and flux)
        aimpoint_select -- list (with length of the number of heliostats) for which
        each element is the corresponding heliostat's selected aimpoint
        '''
        helio_list = list(model.heliostats)
        for h in range(len(helio_list)):
            if aimpoint_select[h] > 0:
                model.select_aimpoint[helio_list[h], aimpoint_select[h]] = 1
            else: 
                #defocus heliostat
                model.defocus[helio_list[h]] = 1
        
    def attemptToAim(self, model, flux, a, h, obj_val):
        '''
        Tries to aim heliostat h at aimpoint h. Adds flux and updates objective value
        if successfully aims without flux violations. 

        Parameters
        ------------
        flux -- numpy array with flux value at each measurement point for initial solution
        a -- aimpoint at which heliostat h is attempting to aim
        h -- heliostat that is attempting to aim
        obj_val -- existing power delivered by the solution (before h-a pairing attempted)

        Returns
        ------------
        If True (successfully aims): New flux and objective value and aimpoint a
        If False (violates flux): Flux and objective value prior to attempting to aim h at a.
        '''
        added_flux = numpy.array([model.flux[h,m,a] for m in model.measurement_points])
        if any(flux+added_flux > model.flux_ubs):
            return flux,obj_val, -1, False
        else:
            for m in model.measurement_points:
                obj_val += model.flux[h,m,a]*model.surface_area[m]
            return flux+added_flux, obj_val, a, True

    def removeFlux(self,model,flux,a,h,obj_val):
        '''
        Removes flux contributed by heliostat h aiming at aimpoint a.
        Updates objective value accordingly.
        Used for post-processing -- switching pairings around after they are initially selected.

        Returns
        ------------
        Flux and objective value after h-a pairing contributions have been removed.
        '''

        removed_flux = numpy.array([model.flux[h,m,a] for m in model.measurement_points])
        flux -= removed_flux
        for m in model.measurement_points:
            obj_val -= model.flux[h,m,a]*model.surface_area[m]
        return flux,obj_val

    def randomSwitchPairs(self,model,flux,helio_list,aimpoint_select,obj_val):
        ''''
        After initial solution is created, randomly selects heliostats and switches
        their aimpoints if doing so results in a higher objective value and does 
        not violate flux limits.

        Parameters
        ----------
        helio_list -- list of heliostats, extracted from model
        aim_list -- list of aimpoints, extracted from model
        aimpoint_select -- list (with length of the number of heliostats) for which
        each element is the corresponding heliostat's selected aimpoint

        Returns
        --------
        New flux and objective value after switch if successful, or existing flux and
        objective value if switch is unsuccessful.
        '''
        for i in range(50):
            # reset flux and obj val to these values if new aiming strategy doesn't help
            saved_flux = flux
            saved_obj = obj_val
            switched = False
            for n in range(10):
                h1 = helio_list[int(self.gen.getVariate()*len(helio_list))]
                h2 = helio_list[int(self.gen.getVariate()*len(helio_list))]
                h1_idx = helio_list.index(h1)
                h2_idx = helio_list.index(h2)
                a2_old = aimpoint_select[h2_idx]
                a1_old = aimpoint_select[h1_idx]
                if a1_old != a2_old:
                    break
            # remove old fluxes, then attempt to add new switched aimpoint fluxes
            flux, obj_val = self.removeFlux(model,flux,a1_old,h1,obj_val)
            flux,obj_val = self.removeFlux(model,flux,a2_old,h2,obj_val)
            flux, obj_val, aim_idx_new_1, success = self.attemptToAim(model, flux, a2_old,h1,obj_val)
            if success:
                flux, obj_val, aim_idx_new_2, success = self.attemptToAim(model, flux, a1_old,h2,obj_val)
                if success:
                    if obj_val > saved_obj:
                        # change aimpoint_slect and flux if switching aimpoints yielded higher obj val
                        aimpoint_select[h1_idx] = aim_idx_new_1
                        aimpoint_select[h2_idx] = aim_idx_new_2
                        self.updateVals(model, aimpoint_select)
                        #print('switched')
                        switched = True
            if not switched:
                flux = saved_flux
                obj_val = saved_obj
                #print('not better to switch')
        print('obj val after switching: ',obj_val)
        return flux,obj_val

    def createDataFrame(self,model,helio_list,aim_list,aimpoint_select):
        '''
        Creates dataframe with columns of heliostats, aimpoints, and number 
        of zeros in flux image. Then creates dataframes according to the aimpoints
        at which the heliostats aim and sorts these dataframes according to the 
        number of zeros in the heliostats' flux images. For use in post-processing
        switching aimpoints.

        Returns
        ----------
        mid_list -- list of heliostats originally pointed at middle aimpoints in order
        of increasing flux size
        edge_list -- list of heliostats originally pointed at edge aimpoints in order
        of decreasing flux size
        '''
        num_zeross = []
        a_cent = math.ceil(len(aim_list) / 2)
        for h in range(len(helio_list)):
            fluxes = numpy.array([model.flux[helio_list[h],m,a_cent] for m in model.measurement_points])
            num_zeros = numpy.count_nonzero(fluxes == 0)
            num_zeross.append(num_zeros)
        df = pd.DataFrame({'heliostats':helio_list,'aimpoints':aimpoint_select,'zeros':num_zeross})
        top_edge = df.loc[df['aimpoints']==float(len(aim_list))]
        bottom_edge = df.loc[df['aimpoints']==1.0]
        middle = df.loc[df['aimpoints']==float(math.ceil(len(aim_list) / 2))]
        edges = pd.concat([top_edge,bottom_edge])
        if len(aim_list) >= 4:
            mid_up = df.loc[df['aimpoints']==float(1+math.ceil(len(aim_list) / 2))]
            middles = pd.concat([middle,mid_up])

        edges.sort_values(by = 'zeros',ascending=True,inplace=True)
        middles.sort_values(by = 'zeros',ascending=False,inplace=True)
        mid_list = middles['heliostats'].to_list()
        edge_list = edges['heliostats'].to_list()

        return mid_list,edge_list

    def switchEdgeMiddle(self,model,flux,helio_list,aimpoint_select,obj_val,edge_list,middle_list):
        '''
        Switches heliostat-aimpoint pairings around if doing so adheres to flux limits 
        and increases the objective value. Prioritizes switching center-aiming heliostats
        with small flux sizes with edge-aiming heliostats with large flux sizes.

        Parameters
        ----------
        middle_list -- list of heliostats originally pointed at middle aimpoints in order
        of increasing flux size (output of createDataFrame)
        edge_list -- list of heliostats originally pointed at edge aimpoints in order
        of decreasing flux size (output of createDataFrame)

        Returns
        --------
        New flux and objective value after switch if successful, or existing flux and
        objective value if switch is unsuccessful.
        '''

        for h in edge_list:
            saved_flux = flux
            saved_obj = obj_val
            switched = False
            h1 = h
            h1_idx = helio_list.index(h1)
            a1_old = aimpoint_select[h1_idx]
            for i in range(20):
                h2 = middle_list[i] # because middle list ordered now
                h2_idx = helio_list.index(h2)
                a2_old = aimpoint_select[h2_idx]
                # remove old fluxes, then attempt to add new switched aimpoint fluxes
                flux, obj_val = self.removeFlux(model,flux,a1_old,h1,obj_val)
                flux,obj_val = self.removeFlux(model,flux,a2_old,h2,obj_val)
                flux, obj_val, aim_idx_new_1, success = self.attemptToAim(model, flux, a2_old,h1,obj_val)
                if success:
                    flux, obj_val, aim_idx_new_2, success = self.attemptToAim(model, flux, a1_old,h2,obj_val)
                    if success:
                        if obj_val > saved_obj:
                            # change aimpoint_slect and flux if switching aimpoints yielded higher obj val
                            aimpoint_select[h1_idx] = aim_idx_new_1
                            aimpoint_select[h2_idx] = aim_idx_new_2
                            self.updateVals(model, aimpoint_select)
                            middle_list.remove(h2) # since that h now pointed at edge
                            #print('switched')
                            switched = True
                            # if successfully switched, move on to next aimpoint in edge list
                            break
                if not switched:
                    flux = saved_flux
                    obj_val = saved_obj
                    #print('not better to switch')
        print('obj val after switching: ',obj_val)
        return flux,obj_val
        

class AcceptRejectHeuristic(Heuristic):
    def __init__(self, num_tries):
        Heuristic.__init__(self,"accept-reject")
        self.num_tries = num_tries
        self.gen = WELL512.WELL512("rngstates.csv")
    
    def getIFS(self, model):
        t_start = time.time()
        obj_val = 0
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
                flux, obj_val, aim_idx, success = self.attemptToAim(model, flux, a, helio_list[h],obj_val)
            if success:
                aimpoint_select[h] = aim_idx
            else:
                defocused_helio.append(h)
        self.updateVals(model, aimpoint_select)
        t_end = time.time()
        print('time taken for initital heuristic section: ',t_end-t_start)
        #print('AcceptReject heuristic ran')
        print('obj val before switching: ',obj_val)
        #flux,obj_val = self.randomSwitchPairs(model,flux,helio_list,aimpoint_select,obj_val)
        '''
        # doesnt sort lists by flux size --> slightly lower power output than if do sort like I do below
        # not shorter time either so delete very soon
        helio_aim = {'heliostats':helio_list,'aimpoints':aimpoint_select}
        df = pd.DataFrame(helio_aim,columns=['heliostats','aimpoints'])
        top_edge = df.loc[df['aimpoints']==float(len(aim_list))]
        bottom_edge = df.loc[df['aimpoints']==1.0]
        middle = df.loc[df['aimpoints']==float(math.ceil(len(aim_list) / 2))]
        top_list = top_edge['heliostats'].to_list()
        bottom_list = bottom_edge['heliostats'].to_list()
        middle_list = middle['heliostats'].to_list()
        edge_list = top_list + bottom_list
        if len(aim_list) >= 4:
            mid_up = df.loc[df['aimpoints']==float(1+math.ceil(len(aim_list) / 2))]
            mid_up_list = mid_up['heliostats'].to_list()
            middle_list += mid_up_list
        
        #flux,obj_val = self.switchEdgeMiddle(model,flux,helio_list,aimpoint_select,obj_val,edge_list,middle_list)
        '''
        mid_list,edge_list = self.createDataFrame(model,helio_list,aim_list,aimpoint_select)
        flux,obj_val = self.switchEdgeMiddle(model,flux,helio_list,aimpoint_select,obj_val,edge_list,mid_list)
        
        # note: don't record or return flux as of now (or update flux model) but definitely could if wanted to
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
            flux, obj_val, aim_idx, success = self.attemptToAim(model, flux, aim_list[aim_idx], helio_list[h],obj_val)
            for i in range(len(aim_list)):
                if success: break
                aim_idx += 1
                aim_idx = aim_idx % len(aim_list)
                print("attempting aimpoint ",aim_idx)    
                flux, obj_val, aim_idx, success = self.attemptToAim(model, flux, aim_list[aim_idx], helio_list[h],obj_val)
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
        self.updateVals(model, aimpoint_select)
        print('obj val: ',obj_val)
        t_end = time.time()
        return(obj_val,t_end-t_start,len(defocused_helio))


class FluxStoreSize(Heuristic):
    def __init__(self,num_tries):
        Heuristic.__init__(self,"center_out")
        self.num_tries = num_tries
        self.gen = WELL512.WELL512("rngstates.csv")

    def resetValues(self,model):
        """ Reset values before a new method is used to generate new feasible h-a pairings, 
        returns the reset values which are then used in the method."""
        flux = numpy.zeros(len(model.measurement_points), dtype=float)
        defocused_helio = []
        helio_list = list(model.heliostats)
        aimpoint_select = numpy.zeros(len(helio_list), dtype=int)
        obj_val = 0
        return flux,defocused_helio,helio_list,aimpoint_select,obj_val

    def calculateFluxSizes(self,model,a_cent):
        '''
        Makes list containing the number of zeros in the flux image of each heliostat on the 
        center aimpoint a_cent. 

        Returns
        -------
        num_zeross -- list (length of the number of heliostats) for which each element is the
        number of zeros in the corresponding heliostat's flux image
        standard_deviation -- standard deviation of num_zeross
        avg_nz -- average of num_zeross (average number of nonzeros in a heliostat's flux image)
        '''
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
        '''
        Updates aim list to have a length of 8 by adding elements or cutting out elements,
        prioritizing center aimpoints if needed. This is to make the aim list suitable
        for the groups used in FluxSizeMethod.

        Returns
        -------
        Updated aim list.
        '''
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
    
    def fluxSizeMethod(self,model,order,flux,defocused_helio,helio_list,aimpoint_select,aim_list,avg_nz,standard_deviation,num_zeross,obj_val):
        '''
        Creates 8 similarly sized/distributed groups of heliostats based on flux size.
        A heliostat's group heavily determines the aimpoint the heliostat attempts on its first try.
        Heliostats with larger flux sizes aim closer to the center.

        Returns
        --------
        Objective value, list of defocused heliostats, flux, and aimpoint_select that results from method.
        '''
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
                flux, obj_val, aim_idx, success = self.attemptToAim(model, flux, a, helio_list[h],obj_val)
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
        flux,defocused_helio,helio_list,aimpoint_select,obj_val = self.resetValues(model)
        order = [7,0,6,1,5,2,4,3]
        obj_val,defocused_helio,flux,aimpoint_select = self.fluxSizeMethod(model,order,flux,defocused_helio,helio_list,aimpoint_select,aim_list,avg_nz,standard_deviation,num_zeross,obj_val)
        print('method 1 obj val: ',obj_val,' method 1 defocused heliostats: ',len(defocused_helio))
        if obj_val > ub:
            ub = obj_val
            self.updateVals(model, aimpoint_select)
            chosen_defoc_helio = defocused_helio.copy()

        # method 2 - same as method 1 but start alternating pattern from bottom instead of top
        flux,defocused_helio,helio_list,aimpoint_select,obj_val = self.resetValues(model)
        order = [0,7,1,6,2,5,3,4]
        obj_val,defocused_helio,flux,aimpoint_select = self.fluxSizeMethod(model,order,flux,defocused_helio,helio_list,aimpoint_select,aim_list,avg_nz,standard_deviation,num_zeross,obj_val)
        print('method 2 obj val: ',obj_val,' method 2 defocused heliostats: ',len(defocused_helio))
        if obj_val > ub:
            ub = obj_val
            self.updateVals(model, aimpoint_select)
            chosen_defoc_helio = defocused_helio.copy()

        # method 3 - now randomly assigning -- most helpful if had some defocused heliostats using the previous method
        # some emphasis on center for first try - 50% chance or more it's chosen
        flux,defocused_helio,helio_list,aimpoint_select,obj_val = self.resetValues(model)
        for h in range(len(helio_list)):
            for n in range(self.num_tries):
                list_to_change = aim_list.copy()
                if n == 0:
                    for i in range(int(len(aim_list))):
                        list_to_change.append(a_cent)
                a = list_to_change[int(self.gen.getVariate()*len(list_to_change))]
                flux, obj_val, aim_idx, success = self.attemptToAim(model, flux, a, helio_list[h],obj_val)
                if success: break
            if success:
                aimpoint_select[h] = aim_idx
            else:
                defocused_helio.append(h)

        print('method 3 obj val: ',obj_val,' method 3 defocused heliostats: ',len(defocused_helio))
        if obj_val > ub:
            ub = obj_val
            chosen_defoc_helio = defocused_helio.copy()
            self.updateVals(model, aimpoint_select)

        obj_val = ub
        print('obj val before switching: ',obj_val)
        mid_list,edge_list = self.createDataFrame(model,helio_list,aim_list,aimpoint_select)
        flux,obj_val = self.switchEdgeMiddle(model,flux,helio_list,aimpoint_select,obj_val,edge_list,mid_list)
        ub = obj_val
        t_end = time.time()
        return(ub,t_end-t_start,len(chosen_defoc_helio))