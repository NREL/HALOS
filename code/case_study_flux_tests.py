# -*- coding: utf-8 -*-
"""
HALOS Case Studies

This script will be used to run the milestone 2.1 tests.


Current Outputs:
    For a single instance of time: 
    1. Single Flux Map  - (RunSingleMap).
    2. Heatmap to compare the error between shifting image and moving aim (MAPE) of multiple cases 
    with full solar field - (RunTestCase).
    
    Over the time:
    3. Yearly simulation for the hours when dni is more than 500 - (YearlyTestCase).
"""

## HALOS Modules 
import flux_model
import flux_method
import sun_shape
import geometry
import mirror_model
import field
import plotting
import WELL512

# Python Modules 
import numpy
import time
import matplotlib.pyplot as plt

## For every instance of time when the DNI is more than 500 throughout the year
def YealyTestCase(kwargs,num_case):
    """
    For a given number of cases:
    Calculates the yearly error between image shifting (flux approximation) and moving aim flux.
    Condition: DNI > 500, Azimuth and Zenith are based on the time when the condition is satisfied 

    Parameters
    ----------
    kwargs :    Dictionary of inputs, Containes - Solar field, Solar azimuth & zenith, dni, aimpoint rows & cols, pts_per_dim
                     Azimuth and Zenith are numpy arrays 
                     DNI has a single value as the objective of function is to calculate relative error 
    num_case :  Case Number - Defined before calling the function 

    Returns
    -------
    map_shift_image:    FluxMap of image shift approximation. 
    diff_map:           absolute difference between shift image and moving aim flux map.
    MAPES:              Array of MAPE for the hours when DNI > 500 for a single case.
    errs:               Relative error between shift image and moving aim flux map.

    """
    gen = WELL512.WELL512("rngstates.csv")
    sf = field.Field(filenames={"field_filename": kwargs["field_filename"]})
    sun = sun_shape.SinglePointSun(0)
    #    params = {"length": 21, "height": 17, "pts_per_dim": 10}
    #north-facing receiver, no tilt
    params = {"length": 21, "height": 17, "pts_per_dim": kwargs["pts_per_dim"], "zenith_deg": 90, "azimuth_deg": 180,
              "rec_cent_offset": [-20, -10, 0]}
    h = 150      # Receiver Height
    receiver = geometry.FlatPlateReceiver(h, params)
    receiver.generateAimpointsGrid(kwargs["aimpoint_rows"],kwargs["aimpoint_cols"])
    mirror = mirror_model.SinglePointGaussianMirror(numpy.array([300, 300, 0]), 5.0)
    weather_file = "./../weather_files/USA NV Tonopah Airport (TMY3).csv"
    method = flux_method.SimpleNormalFluxCalc(3)
    fm = flux_model.FluxModel(sun, mirror, receiver, method, weather_file, sf)
#    print (fm.field.GetCoords())
    #distribute aimpoints evenly
    num_aimpoints = kwargs["aimpoint_rows"] * kwargs["aimpoint_cols"]
    tstart = time.time()
    print("Case: ",num_case," ",kwargs["field_filename"])
    aimpoints = []
    for i in range(fm.field.num_heliostats):
        u = gen.getVariate()
        aimpoints.append(int(u * num_aimpoints))
    MAPES = []
    errs = []
    """ We have an array of solar angles in this case: 
        Quick Testing - Use Slice - We have more than 2000 hours when DNI is more than 500 -
        TODO: The slice has to be same for ploting where the fucntion has been used """
    #for k in range(len(kwargs["solar_zenith"])):  #Uncomment to iterate for every hour when DNI > 500
    #for k in range(len(kwargs["solar_zenith"][:3])):   #Uncomment and use the slicing if testing for fewer iterations 
    for k in range(0,len(kwargs["solar_zenith"]),15):   #Uncomment to do step wise analysis 
        solar_vec = fm.getSolarVector(kwargs["solar_azimuth"][k],kwargs["solar_zenith"][k])
        map_moving_aim = fm.GenerateFullFieldFluxMap(solar_vec,aimpoints,kwargs["dni"],False)
        map_shift_image = fm.GenerateFullFieldFluxMap(solar_vec,aimpoints,kwargs["dni"],True)
        diff_map = (numpy.absolute(map_moving_aim - map_shift_image) / map_moving_aim)
        MAPE = diff_map.mean() 
        err = (numpy.sum(map_moving_aim) - numpy.sum(map_shift_image)) / numpy.sum(map_moving_aim)
        #print(MAPE, err)
        MAPES.append(MAPE)
        errs.append(err)
        print(k)     #Keeping Count for number of iterations done 
    elapsed = time.time() - tstart
    print ("Case time: ",elapsed)
    return map_shift_image, diff_map, MAPES, errs

## Single instance of time
def RunTestCase(kwargs,num_cases=5):
    """
    For a given number of cases:
    Calculates the error between image shifting (flux approximation) and moving aim flux at single instance of 
    time.
    ----------
    Inputs:
    kwargs : Dictionary of inputs, Containes - Solar field, Solar azimuth & zenith, dni, aimpoint rows & cols, pts_per_dim
                Solar Azimuth, Zenith and DNI have single values.
    num_cases : Number of Case Studies. The default is 5.


    Returns
    -------
    map_shift_image:    FluxMap of image shift approximation. 
    diff_map:           absolute difference between shift image and moving aim flux map.
    MAPES:              Array of MAPE for all the cases.
    errs:               Relative error between shift image and moving aim flux map - All Cases.

    """
    #WELL512 generates random numbers
    gen = WELL512.WELL512("rngstates.csv")
    sf = field.Field(filename=kwargs["field_filename"])
    sun = sun_shape.SinglePointSun(0)
#    params = {"length": 21, "height": 17, "pts_per_dim": 10}
    #north-facing receiver, no tilt
    params = {"length": 21, "height": 17, "pts_per_dim": kwargs["pts_per_dim"], "zenith_deg": 90, "azimuth_deg": 180,
              "rec_cent_offset": [-20, -10, 0]}
    h = 150    #tower height
    receiver = geometry.FlatPlateReceiver(h, params)
    receiver.generateAimpointsGrid(kwargs["aimpoint_rows"],kwargs["aimpoint_cols"])
    mirror = mirror_model.SinglePointGaussianMirror(numpy.array([300, 300, 0]), 5.0)
    weather_file = "./../weather_files/USA NV Tonopah Airport (TMY3).csv"
    method = flux_method.SimpleNormalFluxCalc(3)
    fm = flux_model.FluxModel(sun, mirror, receiver, method, weather_file, sf)
#    print (fm.field.GetCoords())
    solar_vec = fm.getSolarVector(kwargs["solar_azimuth"],kwargs["solar_zenith"])
    #distribute aimpoints evenly
    num_aimpoints = kwargs["aimpoint_rows"] * kwargs["aimpoint_cols"]
    MAPES = []
    errs = []
    for case in range(num_cases):      
        tstart = time.time()
        print("Case: ",case," ",kwargs["field_filename"])
        aimpoints = []
        for i in range(fm.field.num_heliostats):
            u = gen.getVariate()
            aimpoints.append(int(u * num_aimpoints))
        map_moving_aim = fm.GenerateFullFieldFluxMap(solar_vec,aimpoints,kwargs["dni"],False)
        map_shift_image = fm.GenerateFullFieldFluxMap(solar_vec,aimpoints,kwargs["dni"],True)
        diff_map = (numpy.absolute(map_moving_aim - map_shift_image) / map_moving_aim)
        MAPE = diff_map.mean() 
        err = (numpy.sum(map_moving_aim) - numpy.sum(map_shift_image)) / numpy.sum(map_moving_aim)
        MAPES.append(MAPE) 
        errs.append(err)
        elapsed = time.time() - tstart
        print ("Case time: ",elapsed)
    return map_shift_image, diff_map, MAPES, errs


## For single heliostat
def RunSingleMap(kwargs):
    """
    Calculates the error between image shifting (flux approximation) and moving aim flux at single instance of 
    time.  
    Generates Single flux map
    
    Parameters
    ----------
    kwargs : Dictionary of Inputs - same as RunTestCase fucntion

    Returns
    -------
    map_shift_image: Flux map of shifting image - same image is thrown at every aimpoint
    map_moving_aim: Different image for every aim point is calculated
    abs_diff: Absolute difference between moving aim and shift image

    """
    gen = WELL512.WELL512("rngstates.csv")
    sf = field.Field(filename=kwargs["field_filename"])
    sun = sun_shape.SinglePointSun(0)
#    params = {"length": 21, "height": 17, "pts_per_dim": 10}
    params = {"length": 21, "height": 17, "pts_per_dim": kwargs["pts_per_dim"], "zenith_deg": 30, "azimuth_deg": 180,
              "rec_cent_offset": [-20, -10, 0]}
    h = 150
    receiver = geometry.FlatPlateReceiver(h, params)
    receiver.generateAimpointsGrid(kwargs["aimpoint_rows"],kwargs["aimpoint_cols"])
    mirror = mirror_model.SinglePointGaussianMirror(numpy.array([300, 300, 0]), 5.0)
    weather_file = "./../weather_files/USA NV Tonopah Airport (TMY3).csv"
    method = flux_method.SimpleNormalFluxCalc(3)
    fm = flux_model.FluxModel(sun, mirror, receiver, method, weather_file, sf)
#    print (fm.field.GetCoords())
    solar_vec = fm.getSolarVector(kwargs["solar_azimuth"],kwargs["solar_zenith"])
#    print(fm.receiver.aim_x[(3,3)])
    #distribute aimpoints evenly 
    aimpoints = [24]*250
#        for i in range(fm.field.num_heliostats):
#            u = gen.getVariate()
#            aimpoints.append(int(u * num_aimpoints))
#    print(aimpoints)
#    assert(False)
    map_moving_aim = fm.GenerateFullFieldFluxMap(solar_vec,aimpoints,kwargs["dni"],False)
#    print(map_moving_aim)
    map_shift_image = fm.GenerateFullFieldFluxMap(solar_vec,aimpoints,kwargs["dni"],True)
    abs_diff = numpy.absolute(map_moving_aim - map_shift_image)
    diff_map = 100*(numpy.absolute(map_moving_aim - map_shift_image) / map_moving_aim)
    return map_shift_image, map_moving_aim, abs_diff



if __name__ == "__main__":
    cases = []   ##Will store input dictionaries of all cases
    ## Solar Fields obtained from SolarPilot
    names = ["flat-daggett-50","flat-daggett-250","flat-morocco-50","flat-morocco-250"]
    filenames = [
            "./../solarpilot_cases/flat-daggett-50.csv",
            "./../solarpilot_cases/flat-daggett-250.csv",
            "./../solarpilot_cases/flat-morocco-50.csv",
            "./../solarpilot_cases/flat-morocco-250.csv" #,
#            "./../solarpilot_cases/radial-daggett-50.csv",
#            "./../solarpilot_cases/radial-daggett-250.csv",
#            "./../solarpilot_cases/radial-morocco-50.csv",
#            "./../solarpilot_cases/radial-morocco-250.csv"
            ]
    """ 
    OVER THE TIME (DNI>500):
    
    Uncomment the following set of code to run the simulation for a year when DNI > 500. 
    Returns a plot of MAPE for a year whenever the DNI is more than 500 (W/m^2) in a CSV file & a plot. 
    Also plots the ECDF of MAPE
    **********
    #TODO:    SHOULD WE MAKE A FUNCTION OF IT? 
        Given weather file and number of cases gets data, goes through YearlyTestCase and generates plots?
    """ 
    num_cases = 2   #Number of cases to be simulated 
    weather_file = "./../weather_files/USA NV Tonopah Airport (TMY3).csv"
    weather_data = flux_model.ReadWeatherFile(weather_file, get_angles=True)
    # print(weather_data)
    ## Storing DATA in input Dictionary 
    kwargs = {}
    print("Started Collecting weather data")
    for i in range(num_cases):             ## Number of Cases - Using 1 for flat-daggett-50 - Test
        kwargs["solar_zenith"] = []
        kwargs["solar_azimuth"] = []
        kwargs["dni"] = []
        kwargs["field_filename"] = filenames[i]
        kwargs["index"] = []
    # Get zenith and azimuth when dni > 500
        for j in range(len(weather_data["dni"])):
            ##SELECT: One specific day of every month or full year using if statements 
            #if weather_data["dni"][j] >= 500 and weather_data["day"][j] == 15:  #Uncomment to simulate for 15th day of every month
            if weather_data["dni"][j] >= 500:     #Uncomment to simulate for every hour when DNI > 500
                zenith_value = weather_data["solar_zenith"][j]
                zenith_value = zenith_value * numpy.pi / 180.
                kwargs["solar_zenith"].append(zenith_value)
                azimuth_value = weather_data["solar_azimuth"][j]
                azimuth_value = azimuth_value * numpy.pi / 180.
                kwargs["solar_azimuth"].append(azimuth_value)
                kwargs["index"].append(j)
                #dni_value = weather_data["dni"][j]      #get DNI array as well
                #kwargs["dni"].append(dni_value)
                #kwargs["solar_zenith"] = 11.68   #Testing single value
                #kwargs ["solar_azimuth"] = 192.66   #Testing single value
        kwargs["aimpoint_rows"] = 5
        kwargs["aimpoint_cols"] = 5
        kwargs["pts_per_dim"] = 20
        kwargs["dni"] = 950
        kwargs["solar_zenith"] = numpy.array(kwargs["solar_zenith"])
        kwargs["solar_azimuth"] = numpy.array(kwargs["solar_azimuth"])
        #kwargs["dni"] = numpy.array(kwargs["dni"])   #If DNI is an Array
        cases.append(kwargs.copy())
    #Testing Done
    print(range(0, len(cases[0]["solar_zenith"]),15))
    print("total number of iterations: ", len(range(0, len(cases[0]["solar_zenith"]),15)))
    print("Weather Data Stored")
    #print(cases[0]["index"])
    #print(kwargs["field_filename"][0])     
    #print(kwargs["solar_zenith"][0])
    #print(len(cases[0]["dni"]))
    #print("Solar Zenith Range", range(0,len(cases[0]["solar_zenith"],15))
    #print(len(cases[0]["solar_azimuth"]))
    #print(cases[0])
          
    """
    Plotting and Storing the results (MAPE) in a csv file
    """
    for num_case in range(num_cases):
    #num_case = 1 #For using a specific case
        map_shift_image, diff_map, MAPES, errs = YealyTestCase(cases[num_case],num_case)
        #Array of Mapes
        arr = numpy.asarray(MAPES)
        # Saving MAPES of a year for  every case separately in SolarPilot Cases folder:
        f_name = cases[num_case]["field_filename"]
        numpy.savetxt(f_name+"Step_15_Yearly_MAPE.csv", arr)
        ##Plot MAPE of a year for every CASE
        plt.figure(1)
        # """ We have an array of solar angles in this case: 
        # Quick Testing - Use Slice - We have more than 2000 hours when DNI is more than 500 -
        # TODO: The slice has to be same as used in YearlyTestCase function """
        #plt.plot(range(len(cases[num_case]["solar_zenith"][:3])), MAPES)   #Uncomment this to use slicing 
        plt.plot(range(0, len(cases[0]["solar_zenith"]),15), MAPES)         #Uncomment this to use step wise analysis
        #plt.plot(range(len(cases[num_case]["solar_zenith"])), MAPES)        #Uncomment to use every hour when DNI > 500
        plt.title(f_name)
        plt.xlabel('Index for Hour Vector (DNI > 500)')
        plt.ylabel('MAPE')
        #plt.show()
        plt.savefig(f_name+'Step_15.png')
        #plt.savefig(f_name+'Yearly_MAPES.pdf')
        
        ## Plotting Emperical CDF:
        x = numpy.sort(MAPES)
        n = x.size
        y = numpy.arange(1, n+1 ) / n
        plt.figure(2)
        plt.plot(x,y)
        plt.title(f_name)
        plt.xlabel('X')
        plt.ylabel('P(MAPE<=X)')
        plt.savefig(f_name+'Step_15_ECDF_MAPE.png')
        
        
    """ 
    SINGLE INSTANCE OF TIME:
    
    Uncomment the following set of code to run the simulation for a single instance of time. 
    
    Output: Generates a heatmap of MAPE for a single instance of time for every case
            Function RunSingleMap() - The code can be used to get a single flux map -  
            Fuction RunTestCase() - simulates given number of cases with full field - Plots Heatmap of MAPE
    """        
        
   #  zenith_angles = [11.68,11.68,12.17,12.17]
   #  azimuth_angles = [192.66,192.66,142.67,142.67]
   #  kwargs = {}
   #  for i in range(4):  ##Define Number of Cases to run
   #      kwargs["field_filename"] = filenames[i]
   #      kwargs["solar_zenith"] = zenith_angles[i] * numpy.pi / 180.
   #      kwargs["solar_azimuth"] = azimuth_angles[i]* numpy.pi / 180.
   #      kwargs["aimpoint_rows"] = 5
   #      kwargs["aimpoint_cols"] = 5
   #      kwargs["pts_per_dim"] = 20
   #      kwargs["dni"] = 950
   #      cases.append(kwargs.copy())
    
   # # To Run Multiple Cases   
   #  for i in range(1):    ##Define Number of Cases to run
   #      map_shift_image, diff_map, MAPES, errs = RunTestCase(cases[i],num_cases=1)
   #      print(map_shift_image)
   #      print(names[i],numpy.mean(MAPES),numpy.mean(errs),numpy.std(MAPES),numpy.std(errs))
   #      plotting.plot_obj_heatmap(map_shift_image,names[i]+"-objmap.png")
   #      plotting.plot_obj_heatmap(diff_map,names[i]+"-diffmap.png")
        
   #  ##Generate single flux map
   #  # map_shift_image, map_moving_aim, abs_diff = RunSingleMap(cases[0])
   #  # plotting.plot_obj_heatmap(map_shift_image,names[i]+"-shiftmap.png")
   #  # plotting.plot_obj_heatmap(map_moving_aim,names[i]+"-origmap.png")
   #  # plotting.plot_obj_heatmap(abs_diff,names[i]+"-diffmap.png")
    