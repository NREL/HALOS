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
                    "use_sp_field","heliostat_group_size","read_flux_from_file"]:
                d[line[0]] = int(line[1])
            elif line[0] in ["mirror_area"]:
                d[line[0]] = float(line[1])
            else:
                d[line[0]] = line[1]
    for key in ["use_sp_flux","use_sp_field","read_flux_from_file"]:
        d[key] =  bool(d.get(key))
    if d.get("heliostat_group_size") is None:
        d["heliostat_group_size"] = 1
    return d

def getReceiverFromFile(filenames,solar_field):
    """
    Creates Receiver provided input file and a solar_field 

    Parameters
    ----------
    filenames : dictionary
        keys include identifiers of filepaths for input; vals are filepaths
    solar_field : Dataframe
        Location of Heliostats

    Returns
    -------
    r : Receiver object

    """
    d = readParamFile(filenames["receiver_filename"])
    if "pts_per_dim" in d.keys():
        d["pts_per_dim"] = int(d["pts_per_dim"])
        num_points = int(d["pts_per_dim"]) * int(d["pts_per_dim"])
        d["pts_per_len_dim"] = d["pts_per_dim"]
        d["pts_per_ht_dim"] = d["pts_per_dim"]
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
    #r.generateFixedFluxLimits(0.0,d["power_rating"] / num_
    if filenames.get("flux_limit_filename") is not None:
        import pandas
        r.flux_upper_limits = readFluxMapFromCSV(filenames["flux_limit_filename"],d["pts_per_ht_dim"],d["pts_per_len_dim"]).flatten()
        r.flux_lower_limits = numpy.ones_like(r.flux_upper_limits) * d["flux_lb"]
        df = pandas.DataFrame(r.flux_upper_limits)
        df.to_csv("flux_limits.csv", index=False)
    else:
        r.generateDynamicFluxLimits(d["flux_lb"], d["flux_ub"], d["n_circulation"])
    if filenames.get("rec_obj_filename") is not None:
        r.obj_by_point = readFluxMapFromCSV(filenames["rec_obj_filename"],d["pts_per_ht_dim"],d["pts_per_len_dim"]).flatten()
    else:
        r.obj_by_point = numpy.ones_like(r.flux_upper_limits)
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
    receiver = getReceiverFromFile(filenames,solar_field)
    mirror = getMirrorModelFromFile(filenames["mirror_filename"],solar_field,settings)
    method = getMethodFromFiles(settings,mirror)
    weather_file = filenames["weather_filename"]
    # change use_sp_flux to True for using SolarPilot API Flux
    fm = flux_model.FluxModel(sun, mirror, receiver, method, weather_file, solar_field,filenames,hour_id,
                              use_sp_flux = settings["use_sp_flux"],
                              read_flux_from_file= settings["read_flux_from_file"])
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
        arr[row_idx] = [float(line[idx]) for idx in range(num_cols)]
        row_idx += 1
        if row_idx == num_rows:
            break
    return arr

def condenseFluxMap(arr, factor=2, use_max=True):
    """
    condenses an array into a smaller one by taking the average of condensed cells.
    """
    num_rows = int(len(arr)/factor)
    num_cols = int(len(arr[0])/factor)
    new_arr = numpy.zeros([num_rows,num_cols],dtype=float)
    if use_max:
        for r in range(num_rows):
            for c in range(num_cols):
                new_arr[r,c] = numpy.max(arr[r*factor:(r+1)*factor,c*factor:(c+1)*factor])
    else:
        for r in range(num_rows):
            for c in range(num_cols):
                new_arr[r,c] = numpy.average(arr[r*factor:(r+1)*factor,c*factor:(c+1)*factor])
    return new_arr

def expandFluxMap(arr, factor=2):
    """
    Expands an array into a larger one by making a copy of cells in arr.
    """
    num_rows = int(len(arr)*factor)
    num_cols = int(len(arr[0])*factor)
    new_arr = numpy.zeros([num_rows,num_cols],dtype=float)
    for r in range(num_rows):
        for c in range(num_cols):
            new_arr[r,c] = arr[r//factor,c//factor]
    return new_arr

if __name__ == "__main__":
    filename = "./../case_inputs/filenames_study/active_surface.csv"
    outfilename = "./../case_inputs/filenames_study/active_surface2.csv"
    num_rows = 30
    num_cols = 28
    arr = readFluxMapFromCSV(filename, num_rows, num_cols)
    print(arr)