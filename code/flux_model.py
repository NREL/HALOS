# -*- coding: utf-8 -*-
"""
HALOS
flux_model module
Framework for building out models for flux modeling and outputting images
to either a data file or an optimization model.
"""
import csv
import numpy
import sol_pos
import matplotlib.pyplot as plt
import sp_module

def ReadWeatherFile(weather_file, get_angles=False):
    """
    Reads in a Typical Meteroological Year (TMY) weather file compatible 
    with SolarPILOT and records time-series DNI, plus lat/lon coordinates.
    This is stored as a dictionary under the member name 'weather_data'.
    
    arguments: 
        weather_file -- path to TMY weather file
    
    returns: 
        weather_data, a dictionary with categorical and time-series data
    """
    weather_data = {"dni": numpy.zeros(8760, dtype=float), 
                    "year": numpy.zeros(8760, dtype=int), 
                    "month": numpy.zeros(8760, dtype=int), 
                    "day": numpy.zeros(8760, dtype=int), 
                    "hour": numpy.zeros(8760, dtype=int),
                    "solar_zenith": numpy.zeros(8760, dtype=float),
                    "solar_azimuth": numpy.zeros(8760, dtype=float)}
    reader = csv.reader(open(weather_file, 'r'))
    # First line is column headers for categorical data
    header_keys = next(reader)
    # Second line is data for first line of headers; record lat and lon
    fline = next(reader)
    d = {}
    for i in range(len(fline)):
        d[header_keys[i]] = fline[i]
    weather_data["lat"] = float(d["Latitude"])
    weather_data["lon"] = float(d["Longitude"])
    weather_data["time_zone"] = float(d["Time Zone"])
    # Third line is column headers for time-series data, just track dni
    keys = next(reader)
    year_idx = keys.index("Year")
    month_idx = keys.index("Month")
    day_idx = keys.index("Day")
    hour_idx = keys.index("Hour")
    if d["Source"] == "IWEC" or d["Source"]=="INTL":
        dni_idx = keys.index("Beam")
    else:
        dni_idx = keys.index("DNI")
    #Record time-series data
    t = 0
    for line in reader:
        weather_data["dni"][t] = float(line[dni_idx])
        weather_data["year"][t] = int(line[year_idx])
        weather_data["month"][t] = int(line[month_idx])
        weather_data["day"][t] = int(line[day_idx])
        weather_data["hour"][t] = int(line[hour_idx]) 
        t += 1
        if t == 8760:
            break
    if get_angles:
        import sol_pos
        zenith, azimuth = sol_pos.solar_position(weather_data)
        weather_data["solar_zenith"] = numpy.array(zenith)
        weather_data["solar_azimuth"] = numpy.array(azimuth)
    return weather_data


class FluxModel(object):
    def __init__(self,sun_shape, mirror_shape, receiver, flux_method, 
                 weather_file, field, filenames = None, hour_id = None,
                 use_sp_flux = False, dni = None, read_flux_from_file = False):
        self.sun_shape = sun_shape
        self.mirror_shape = mirror_shape
        self.receiver = receiver
        self.flux_method = flux_method
        self.weather_data = ReadWeatherFile(weather_file)
        self.field = field
        self.filenames = filenames
        self.hour_id = hour_id
        self.use_sp_flux = use_sp_flux
        self.read_flux_from_file = read_flux_from_file
        if self.receiver.dni != None:
            self.dni = self.receiver.dni
        elif self.hour_id != None: 
            self.dni = self.weather_data['dni'][hour_id]
        elif dni != None:
            self.dni = dni
        else: 
            msg = "No inputs given for DNI."
            raise Exception(msg)
        if use_sp_flux:
            self.get_sp_flux_parallel()
        elif self.read_flux_from_file:
            self.get_flux_maps_from_files()
        else:
            self.get_halos_native_flux()
        if self.receiver.params["receiver_type"] == "External cylindrical":
            print("Generating Fraction Maps")
            self.getFractionMap()
        else:
            print("No fraction maps generation - Flat Plate")
            pass


    def addSettings(self,settings):
        self.settings = settings
        
        
    def HermiteConvolution(self, aim_loc, measure_loc):
        pass  # do the hermite evaluation


    def getSolarVector(self,azimuth, zenith):
        """
        returns normalized 3-d vector when given azimuth and zenith as input.
        """
        x = numpy.sin(azimuth)*numpy.cos(zenith)
        y = numpy.sin(zenith)
        z = numpy.cos(azimuth)*numpy.cos(zenith)
        return numpy.array([x,y,z])


    def GenerateSingleFluxMap(self,helio_idx,aimpoint_idx,solar_vector,approx=True):
        """
        for a single heliostat and aimpoint pairing, generates a flux image on 
        the receiver.
                                     "or"
        Given use_sp_flux true and filenames to class FLuxModel generates a flux image on 
        the receiver for a single Heliostat using SolarPilot API
        """
        
        if self.use_sp_flux:
            import sp_flux
            flux_map = sp_flux.single_helio_flux(self.filenames,self.field, helio_idx,self.weather_data,self.hour_id)
            flux_map = numpy.array(flux_map) 
            return flux_map
        else:
            flux_map = numpy.zeros_like(self.receiver.x)
            helio = self.field.coords[helio_idx]
            aim = self.receiver.aimpoints[aimpoint_idx]
            center_aim = self.receiver.aimpoints[len(self.receiver.aimpoints) // 2]
            for mrow in range(len(self.receiver.x)):
                for mcol in range(len(self.receiver.x[mrow])):
                    norm_idx = mrow*len(self.receiver.x[mrow]) + mcol
                    measurement = numpy.array([self.receiver.x[mrow,mcol],self.receiver.y[mrow,mcol],self.receiver.z[mrow,mcol]])
                    flux_map[mrow,mcol] += self.flux_method.GetFlux(helio,aim,measurement,self.receiver.normals[norm_idx],solar_vector,self.dni,approx,center_aim)
            return flux_map


    def GenerateFullFieldFluxMap(self,solar_vector,aimpoints,approx = True):
        """
        given a solar vector and aimpoint indices, generates a flux map of the 
        receiver.
                                     "or"
        Given use_sp_flux true and filenames to class FLuxModel generates a flux image on 
        the receiver using SolarPilot API
        """
        if self.use_sp_flux:
            import sp_flux
            full_map = sp_flux.full_field_flux(self.filenames)
            full_map = numpy.array(full_map) 
            
        else:
            full_map = numpy.zeros_like(self.receiver.x)
            for i in range(self.field.num_heliostats): 
                full_map += self.GenerateSingleFluxMap(i,aimpoints[i],solar_vector,approx)
        return full_map


    def GetNormalizationFactor(self,mirror_power,flux_map):
        return mirror_power / sum(flux_map*self.receiver.surface_area.flatten())


    def flux_by_section(self, section_id):
        """
        Given section ID returns dict of helios-flux_map pairing 

        Parameters
        ----------
        section_id :
            index of section

        Returns
        -------
        flux : dict
            key: Helio_index, Value: Flux_map

        """
        flux = {}
        if self.use_sp_flux:
            for h in self.field.helios_by_section[section_id]:
                flux_map = self.sp_flux.get_single_helio_flux(h,self.weather_data,self.hour_id,self.dni)
                flux[h] = numpy.array(flux_map)
        elif self.read_flux_from_file:
            import inputs
            num_cols = self.receiver.params["pts_per_len_dim"]
            num_rows = self.receiver.params["pts_per_ht_dim"]
            for h in self.field.helios_by_section[section_id]:
                fname = self.filenames["heliostat_file_dir"]+"heliostat"+str(h+1)+".csv"
                flux[h] = inputs.readFluxMapFromCSV(fname, num_rows, num_cols)
        else:
            for h in self.field.helios_by_section[section_id]:
                if self.receiver.params["receiver_type"] == "Flat plate":
                    center_idx = len(self.receiver.aimpoints) // 2
                else:
                    center_idx = len(self.receiver.aimpoints[h]) // 2
                solar_vector = self.getSolarVector(self.weather_data["solar_azimuth"][self.hour_id], self.weather_data["solar_zenith"][self.hour_id])
                flux_map = self.GenerateSingleFluxMap(h, center_idx, solar_vector, True)
                flux[h] = numpy.array(flux_map)
        return flux

    def get_halos_native_flux(self):
        """
        Get flux map per heliostat in parallel Using HALOS

        Returns
        -------
        None. Populates self.parallel_flux_maps

        """
        self.sp_flux = sp_module.SP_Flux(self.filenames, self.field)
        print(len(self.field.x))
        section_id = []
        for s in range(self.field.num_sections):
            section_id.append(s)
        self.parallel_flux_maps = {}
        for s in section_id:
            sect = self.flux_by_section(s)
            self.parallel_flux_maps.update(sect)
        for h in self.parallel_flux_maps.keys():
            self.parallel_flux_maps[h] = numpy.matrix(self.parallel_flux_maps[h]).round(3)  # round to nearest Watt
        print("HALOS Native Flux Calculation - Done!")

    def get_flux_maps_from_files(self):
        """
        Get flux maps by reading CSV files directly.

        Returns
        -------
        None. Populates self.parallel_flux_maps
        """
        section_id = []
        for s in range(self.field.num_sections):
            section_id.append(s)
        import multiprocessing
        p = multiprocessing.Pool(multiprocessing.cpu_count())
        results = p.map(self.flux_by_section, section_id)
        self.parallel_flux_maps = {}
        for sect in results:
            self.parallel_flux_maps.update(sect)
        for h in self.parallel_flux_maps.keys():
            self.parallel_flux_maps[h] = numpy.matrix(self.parallel_flux_maps[h]).round(3)  # round to nearest Watt
        print("Parallel Flux Calculation - Done!")
        del (p)

    def get_sp_flux_parallel(self):
        """
        Get flux map per heliostat in parallel Using Solarpilot

        Returns
        -------
        None. Populates self.parallel_flux_maps

        """
        self.sp_flux = sp_module.SP_Flux(self.filenames, self.field)
        section_id = []
        for s in range(self.field.num_sections):
            section_id.append(s)
        import multiprocessing
        p = multiprocessing.Pool(multiprocessing.cpu_count())
        results = p.map(self.flux_by_section, section_id)
        self.parallel_flux_maps = {}
        for sect in results:
            self.parallel_flux_maps.update(sect)
        for h in self.parallel_flux_maps.keys():
            self.parallel_flux_maps[h] = numpy.matrix(self.parallel_flux_maps[h]).round(3)  #round to nearest Watt
        print("Parallel Flux Calculation - Done!")
        del(p)


    def ShiftImage_GenerateSingleHeliostatFluxMaps(self,helio_idx,solar_vector,approx=True):
        """
        Given the heliostat, a solar vector, dni, and an efficiency map, 
        generates a flux map of the receiver for each aimpoint-measurement 
        point pair.  
        FLux map at center aimpoint is shifted on other aimpoints as well. 
        Here, the flux maps are flatteend into 1-dimensional 
        arrays for the optimization model.
        
        
        Returns
        -------
        maps:  Dictionary 
            Key is the index of the aimpoint whereas value has the repected 
            flux map generated through shifting center aim map.

        """
        def shift_2d_array(data, dx, dy, constant=False):
            """
            Shifts the array in two dimensions while setting rolled values to constant
    
            Parameters
            ----------
            data : 2D Array
                The 2d numpy array to be shifted
            dx : int
                The shift in x
            dy : int
                The shift in y
            constant :  optional
                The constant to replace rolled values with
    
            Returns
            -------
            shifted_array : 2D Array
                The shifted array with "constant" where roll occurs
    
            """
            
            shifted_data = numpy.roll(data, dx, axis=1)
            if dx < 0:
                shifted_data[:, dx:] = constant
            elif dx > 0:
                shifted_data[:, 0:dx] = constant
        
            shifted_data = numpy.roll(shifted_data, dy, axis=0)
            if dy < 0:
                shifted_data[dy:, :] = constant
            elif dy > 0:
                shifted_data[0:dy, :] = constant
            return shifted_data
        #Getting center map using same code of previous function 
        maps = {}
        aim_rows = int(self.receiver.params["aim_rows"])
        #sp_flux_dict = self.flux_dict_SP()
        #map_center = numpy.array(self.sp_flux_dict[helio_idx + 1])
        if self.receiver.params["receiver_type"] == "External cylindrical":
            center_idx = len(self.receiver.aimpoints[helio_idx]) // 2
            aim_cols = 1
            #PARALLEL FLUX
            map_center = self.parallel_flux_maps[helio_idx]
            x_shift_size = 0
            y_shift_size = round(len(map_center)/aim_rows)

        elif self.receiver.params["receiver_type"] == "Flat plate":
            center_idx = len(self.receiver.aimpoints) // 2
            map_center = self.parallel_flux_maps[helio_idx]
            aim_cols = int(self.receiver.params["aim_cols"])
            x_shift_size = round(len(map_center)/aim_cols)
            y_shift_size = round(len(map_center)/aim_rows)
        map_center = numpy.array(map_center)   
        center_map = map_center.flatten()
        if not (self.use_sp_flux or self.read_flux_from_file):
            mirror_power = self.dni * self.field.mirror_area
            factor = self.GetNormalizationFactor(mirror_power,center_map)
        else:
            factor = 1.0
        maps[center_idx] = factor * center_map # Center map stored in dict after flatting 
        # Image Shifting uses Center MAP before flatting
        for y_shift in range(0,(aim_rows//2)+1): #Aimpoint are 5 by 5 so we have to 2 steps up and 2 steps down. y_shift zero is center row
            #y_shift +ive goes downwards 
            for x_shift in range(0,(aim_cols//2)+1): #Similarly x_shift 0 is center column where as +2 is ahead and -2 is behind
                ## Center Line
                ## Positive in x moves ahead - Positive in y here moves down
                x_frwd = shift_2d_array(map_center, x_shift*x_shift_size, y_shift*y_shift_size, constant=0)
                ##TO PLOT
                # im = plt.imshow(x_frwd)
                # plt.colorbar(im)
                # plt.tight_layout()
                # plt.show()
                maps[center_idx + (aim_cols * y_shift) + x_shift] = factor * x_frwd.flatten()
                x_bcwd = shift_2d_array(map_center, - x_shift*x_shift_size, y_shift*y_shift_size, constant=0)
                ##TO PLOT
                # im = plt.imshow(x_bcwd)
                # plt.colorbar(im)
                # plt.tight_layout()
                # plt.show()
                maps[center_idx + (aim_cols * y_shift) - x_shift] = factor * x_bcwd.flatten()
                ## Negative in x moves behind - Negative in y here moves up from center
                x_frwd = shift_2d_array(map_center, x_shift*x_shift_size, - y_shift*y_shift_size, constant=0)
                maps[center_idx - (aim_cols * y_shift) + x_shift] = factor * x_frwd.flatten()
                x_bcwd = shift_2d_array(map_center, - x_shift*x_shift_size, - y_shift*y_shift_size, constant=0)
                maps[center_idx - (aim_cols * y_shift) - x_shift] = factor * x_bcwd.flatten()
        #print(maps.keys())   ##Indexes Tested - Image shift tested using plotting
        return maps


    def SectionFluxMap(self): 
        """
        Calculates section wise column sum of heliostat wise flux_maps

        Returns
        -------
        Dictionary with column sum of flux for each section
 
        """
        self.flux_sum_cols = {}
        for s in range(self.field.num_sections):
            flux = numpy.matrix(numpy.zeros_like(self.receiver.x))
            for h in self.field.helios_by_section[s]:
                map_center = numpy.matrix(self.parallel_flux_maps[h])
                flux = flux + map_center
            print("Section Number: ", s)
            sum_col = numpy.sum(numpy.array(flux), axis = 0)
            self.flux_sum_cols[s] = sum_col
        return self.flux_sum_cols 


    def getFractionMap(self): 
        """
        Using Section Column flux sum calculates the fraction for each column.
        The fraction is used to divide flux limits on the receiver in cylindrical
        case.

        Returns
        -------
        None, Populates self.fraction_maps

        """
        import pandas
        col_sum_each_map = self.SectionFluxMap()
        col_sums = []  # Stores sum as: 1st list is sum of 1st column of each section

        for ncol in range(self.receiver.params["pts_per_len_dim"]):
            col_sum_each_column = []   #Stores sum of one specific column from each section
            for s in range(self.field.num_sections):
                col_sum_each_column.append((col_sum_each_map[s])[ncol])
            col_sums.append(col_sum_each_column)
        
        self.fraction_maps = []
        for s in range(self.field.num_sections):
            fraction_map = pandas.DataFrame(numpy.zeros_like(self.receiver.x))
            for ncol in range(self.receiver.params["pts_per_len_dim"]):
                fraction_map[ncol] = fraction_map[ncol] + (col_sums[ncol])[s]/sum(col_sums[ncol])
            self.fraction_maps.append(numpy.array(fraction_map))


    def GenerateSingleHeliostatFluxMaps(self,helio_idx,solar_vector,dni,approx=True):
        """
        given the heliostat, a solar vector, dni, and an efficiency map, 
        generates a flux map of the receiver for each aimpoint-measurement 
        point pair.  Here, the flux maps are flatteend into 1-dimensional 
        arrays for the optimization model.
        """  
        center_idx = len(self.receiver.aimpoints) // 2
        maps = {}
        center_map = self.GenerateSingleFluxMap(helio_idx,center_idx,solar_vector,approx).flatten()
        mirror_power = dni * self.field.mirror_area
        factor = self.GetNormalizationFactor(mirror_power,center_map)
        maps[center_idx] = factor * center_map
        for aimpoint_idx in range(len(self.receiver.aimpoints)):
            if aimpoint_idx != center_idx:
                maps[aimpoint_idx] = factor * self.GenerateSingleFluxMap(helio_idx,aimpoint_idx,solar_vector,approx).flatten()
        return maps


    def GenerateFullFluxProblem(self,solar_vector,dni,outfile=None,approx=True):
        all_maps = {}
        for helio_idx in range(self.field.num_heliostats):
            all_maps[helio_idx] = self.GenerateSingleHeliostatFluxMaps(helio_idx,solar_vector,dni,approx)
        return all_maps
        
    
if __name__ == "__main__":
    import sun_shape
    import geometry
    import mirror_model
    import field

    sun = sun_shape.SinglePointSun(0)
#    params = {"length": 21, "height": 17, "pts_per_dim": 51}
    params = {"length": 21, "height": 17, "pts_per_dim": 10, "zenith_deg": 0, "azimuth_deg": 0,
              "rec_cent_offset": [-10, 5, 0]}
    h = 150
    receiver = geometry.FlatPlateReceiver(h, params)
    receiver.generateAimpointsGrid(5,5)
    mirror = mirror_model.SinglePointGaussianMirror(numpy.array([300, 300, 0]), 5.0)
    weather_file = "./../weather_files/USA NV Tonopah Airport (TMY3).csv"
    f = field.Field()
    f.GenerateCoords(numpy.array([-200,-100,0,100,200]),numpy.array([-300,-300,-300,-300,-300]),numpy.array([0,0,0,0,0]))
    fm = FluxModel(sun, mirror, receiver, None, weather_file, f)
    

