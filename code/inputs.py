# -*- coding: utf-8 -*-
"""
Module that reads in input from files, to be used for input to the 
flux and optimization models.
"""

import flux_model
import flux_method
import sun_shape
import geometry
import mirror_model
import field
import csv
import numpy

def readCaseFile(master_filename):
    """
    reads a csv file with two columns, and returns a dictionary in which the 
    keys and values are strings.  
    
    Parameters:
        master_filename -- path to file
    Returns:
        d -- dictionary mapping file types to paths
    """
    d = {}
    reader = csv.reader(open(master_filename,'r'))
    for line in reader:    
        if len(line)>1:
            d[line[0]] = line[1]
    return d

def readParamFile(param_filename):
    """
    reads a csv file with two columns, and returns a dictionary in which the 
    keys are strings and the values are specific numeric types. 
    
    Parameters:
        param_filename -- path to file
    Returns:
        d -- dictionary mapping parameters to values
    """
    d = {}
    reader = csv.reader(open(param_filename,'r'))
    for line in reader:    
        if len(line)>1:
            if line[0] == "receiver_type":
                d[line[0]] = line[1]
            else:
                d[line[0]] = float(line[1])
    return d

def readSettingsFile(settings_filename):
    """
    reads a csv file with two columns, and returns a dictionary in which the 
    keys and values are strings. 
    
    Parameters:
        param_filename -- path to file
    Returns:
        d -- dictionary mapping parameters to values
    """
    d = {}
    reader = csv.reader(open(settings_filename,'r'))
    for line in reader:    
        if len(line)>1:
            if line[0] in ["num_sections","use_sp_flux","hour_idx",
                    "use_sp_field","heliostat_group_size"]:
                d[line[0]] = int(line[1])
            elif line[0] in ["mirror_area"]:
                d[line[0]] = float(line[1])
            else:
                d[line[0]] = line[1]
    if d["use_sp_flux"] == 1:
        d["use_sp_flux"] = True
    else:
        d["use_sp_flux"] = False
    if d["use_sp_field"] == 1:
        d["use_sp_field"] = True
    else:
        d["use_sp_field"] = False
    if "heliostat_group_size" not in d.keys():
        d["heliostat_group_size"] = 1
    return d

def getReceiverFromFile(filename,solar_field):
    """
    Creates Receiver provided input file and a solar_field 

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    solar_field : Dataframe
        Location of Heliostats

    Returns
    -------
    r : Receiver object

    """
    d = readParamFile(filename)
    if "pts_per_dim" in d.keys():
        d["pts_per_dim"] = int(d["pts_per_dim"])
        num_points = int(d["pts_per_dim"]) * int(d["pts_per_dim"])
    else:
        d["pts_per_len_dim"] = int(d["pts_per_len_dim"])
        d["pts_per_ht_dim"] = int(d["pts_per_ht_dim"])
        num_points = int(d["pts_per_len_dim"]) * int(d["pts_per_ht_dim"])
    d["rec_cent_offset"] = [d["rec_cent_offset_x"],d["rec_cent_offset_y"],d["rec_cent_offset_z"]]
    if d["receiver_type"] == "Flat plate":
        r = geometry.FlatPlateReceiver(d["tow_height"],d)
        print("Flat Plate Receiver")
    elif d["receiver_type"] == "External cylindrical":
        r = geometry.CylindricalPlateReceiver(d["tow_height"],d,solar_field)
    #r.generateFixedFluxLimits(0.0,d["power_rating"] / num_points)
    r.generateDynamicFluxLimits(d["flux_lb"], d["flux_ub"], d["n_circulation"])
    return r
    

def getMirrorModelFromFile(filename,solar_field,settings):
    d = readParamFile(filename)
    if settings["mirror_model"] == "SinglePointGaussian":
        return mirror_model.SinglePointGaussianMirror(solar_field.GetCoords(),d["error"])
    return None

def getMethodFromFiles(settings,mirror):
    """
    obtains a flux calculation method using the a settings dictionary and 
    mirror object provided as input.
    
    Parameters:
        settings (Dict) -- settings, including reference to method type
        mirror (Mirror) -- mirror object that includes error level
    Returns: 
        method -- flux calcuation method object from flux_method module
    """
    if settings["method"] == "SimpleNormalFluxCalc":
        return flux_method.SimpleNormalFluxCalc(mirror.error)
    return None

def getFullFluxModelFromFiles(case_filename, hour_id = None):
    """
    Creates and returns flux model object after creating field and receiver objects

    Parameters
    ----------
    case_filename : CSV with path to filenames
    hour_id : hour index, optional

    Returns
    -------
    fm: Flux model object

    """
    filenames = readCaseFile(case_filename)
    settings = readSettingsFile(filenames["settings"])
    sun = sun_shape.SinglePointSun(0)  #for these case studies, rm sun shape
    solar_field = field.Field(filenames,params=settings, use_sp_field = settings["use_sp_field"])
    receiver = getReceiverFromFile(filenames["receiver_filename"],solar_field)
    mirror = getMirrorModelFromFile(filenames["mirror_filename"],solar_field,settings)
    method = getMethodFromFiles(settings,mirror)
    weather_file = filenames["weather_filename"]
    # change use_sp_flux to True for using SolarPilot API Flux
    fm = flux_model.FluxModel(sun, mirror, receiver, method, weather_file, solar_field,filenames,hour_id, use_sp_flux = settings["use_sp_flux"])
    fm.addSettings(settings)
    return fm

def readFluxMapFromCSV(filename,num_rows,num_cols):
    """
        reads a csv file with a flux map as a matrix.

        Parameters:
            filename -- path to file
            num_rows -- number of rows in flux map matrix
            num_cols -- number of rows in flux map matrix
        Returns:
            arr -- 2-dimensional array with flux map
        """
    arr = numpy.zeros([num_rows,num_cols],dtype=float)
    reader = csv.reader(open(filename, 'r'))
    row_idx = 0
    for line in reader:
        print(len(arr))
        print(len(line))
        arr[row_idx] = [float(line[idx]) for idx in range(num_cols)]
        row_idx += 1
        if row_idx == num_rows:
            break
    return arr



if __name__ == "__main__":
    dirpath = "./../case_inputs/NREL-Tietronix aimpoint case study/flux_maps_21June_noon/"
    helio_idx = 1
    filename = dirpath + "heliostat"+str(helio_idx)+".csv"
    num_rows = 60
    num_cols = 56
    arr = readFluxMapFromCSV(filename,num_rows,num_cols)
    print(arr)