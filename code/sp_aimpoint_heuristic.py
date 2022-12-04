# -*- coding: utf-8 -*-
"""
HALOS SolarPILOT API aimpoint heuristic

This module contains aimpoint heuristics using SolarPILOT flux calculations.
"""

import numpy as np

import sp_module

def maximum_violation_heur(helios_flux, flux_limit, print_messages = False):
    """
    Defocuses helistats that create the maximum violation to flux limit
    """
    helio_id = list(helios_flux.keys())
    # Converting over to np.array
    for key in helio_id:
        helios_flux[key] = np.array(helios_flux[key])

    comb_flux = np.zeros((len(helios_flux[helio_id[0]]), len(helios_flux[helio_id[0]][0])))
    for key in helio_id:
        comb_flux += helios_flux[key]

    flux_limit = np.array(flux_limit)
    diff = np.subtract(comb_flux, flux_limit)
    defocused_helios = []
    while np.max(diff) > 0.0:
        # Find maximum violation
        maxV_ind = np.unravel_index(np.argmax(diff, axis=None), diff.shape)
        if print_messages:
            print("Maximum Flux Violation: {:.2f} kW/m2".format(diff[maxV_ind]))
        maxV_dict = {}
        for key in helio_id:
            maxV_dict[key] = helios_flux[key][maxV_ind]
        defocused_helios.append(max(maxV_dict, key=maxV_dict.get))
        if print_messages:
            print("Defocusing Heliostat ID: "+ str(int(defocused_helios[-1])))
        helio_id.remove(defocused_helios[-1])
        comb_flux -= helios_flux[defocused_helios[-1]]
        diff = np.subtract(comb_flux, flux_limit)
    
    return comb_flux, defocused_helios

def flux_limit_heuristic(helios_flux, flux_limit, heur_method = 'maximum violation'):
    """
    Determines the flux profile subject to a flux limit using a heuristic method.

    Parameters
    ----------
    helios_flux : Dictionary 
        Contains flux images for each heliostat (n x m)
    flux_limit : list (n x m) or float
        Either matrix of flux limits with the size of a heliostat image 
            or a single float that will be imposed across the whole receiver

    Returns
    -------
    flux : list 
        Returns a feasiable flux matrix using the choosen heuristic
    defocused_helios : list
        A list of heliostat id that were defocus to meet flux limit
    """
    helio_id = list(helios_flux.keys())
    if type(flux_limit) == list:
        if len(flux_limit) != len(helios_flux[helio_id[0]]) or len(flux_limit[0]) != len(helios_flux[helio_id[0]][0]):
            print("ERROR: if flux_limit is a list, then it must be the same size as the heliostat images")
            return
    elif type(flux_limit) == float:
        flux_list = [flux_limit]*len(helios_flux[helio_id[0]][0])
        flux_limit = [flux_list for _ in range(len(helios_flux[helio_id[0]]))]
    else:
        print("ERROR: flux_limit must be either a list or a float.")
        return

    if heur_method.lower() == 'maximum violation':
        flux, defocused_helios = maximum_violation_heur(helios_flux, flux_limit)
    else:
        print("ERROR: the heuristic method choosen is not supported")
        return

    return flux, defocused_helios

if __name__ == "__main__":
    case_filename = "./../case_inputs/flat_250_ca_case.csv"
    import inputs
    import pickle
    import matplotlib.pyplot as plt

    field_flux_dict_file = "field_flux_image_priority.pkl"
    load_data = False
    if not load_data:
        kwargs = inputs.readCaseFile(case_filename)
        flux_dict = sp_flux.single_helio_flux_dict(kwargs, aim_method='Image size priority')
        with open(field_flux_dict_file, 'wb') as f:
            pickle.dump(flux_dict, f)
    else:
        print("Loading heliostat data in from pickle...")
        with open(field_flux_dict_file, 'rb') as f:
            flux_dict = pickle.load(f)

    # Checking for zero contributions heliostats
    for key in flux_dict.keys():
        if sum([sum(i) for i in zip(*flux_dict[key])]) < 0.01:
            print("Heliostat " + str(key) + " contributes no flux to receiver!")
    
    # Combined flux map
    
    comb_flux = np.zeros((len(flux_dict[0]), len(flux_dict[0][0])))
    for key in flux_dict.keys():
        comb_flux += np.array(flux_dict[key])

    # Before defocusing
    if True:
        im = plt.imshow(comb_flux)
        plt.colorbar(im)
        plt.tight_layout()
        plt.show()
        
        
    flux_before_defocus = np.array(comb_flux)
    if True:   ##Get Flux Limits from CSV
        import csv
        with open('flux_limits.csv', 'r') as data:
            reader = csv.reader(data)
            next(data)
            flux_limit = []
            for row in reader:
                row_value = []
                for cell in row:
                    row_value.append(float(cell))
                flux_limit.append(row_value)
    else:
        flux_limit = [100 + i*15. for i in range(10)]
        flux_limit.extend(flux_limit[::-1])
        flux_limit = [flux_limit for _ in range(20)]

    comb_flux, defocused_helios = flux_limit_heuristic(flux_dict, flux_limit)

    print("Number of Heliostats Defocused: " + str(len(defocused_helios)))

    # After defocusing
    if True:
        im = plt.imshow(comb_flux)
        plt.colorbar(im)
        plt.tight_layout()
        plt.show()
        
    ##The following code is hardcoded for a single example case
    ##This has been generalized in sp_module.py which can be used instead. 
    receiver = "flat"
    if receiver == "cyl":
        area_m_point = (((10.38 * (np.pi))/(20)) * (17 / 20))
    area_m_point = (21/20) * (17/20)

    flux = np.array(comb_flux)
    flux_before_defocus_with_area = flux_before_defocus * area_m_point
    flux_with_area = flux * area_m_point
    obj_value_before = np.sum(flux_before_defocus_with_area)
    obj_value_after = np.sum(flux_with_area)
    print("The Objective value before defocusing is: ", obj_value_before)
    print("The Objective value after defocusing is: ", obj_value_after)