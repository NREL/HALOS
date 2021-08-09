
"""
HALOS, SolarPilot Integration Module 

Contains seprate subclasses for Solarpilot field and flux. Module for generating
flux images and layout (field) using the Solarpilot API. Can populate an 
optimization model or can independently run several solarpilot cases using the
co-pylot API.  

"""
import flux_model
import os, csv
import matplotlib.pyplot as plt
from copylot import CoPylot
import numpy
import pandas
import flux_model

cp = CoPylot()

class SolarPilot:
    def __init__(self, filenames):
        self.filenames = filenames
        self.get_input_data()
        
    def get_input_data(self):
        """
        Reads the input receiver csv file and assigns input parameters as 
        key/value to dictionary receiver_data

        Returns
        -------
        None.

        """
        self.receiver_data = {}
        reader = csv.reader(open(self.filenames["receiver_filename"],'r'))
        for line in reader:    
            if len(line)>1:
               self.receiver_data[line[0]] = line[1]
               
    def assign_inputs(self, weather_data = None, hour_id = None, dni = None, 
                      not_filter_helio = False, center_aimpoint = True, read_weather = False):
        """
        Creates solarpilot data pointer and assigns the inputs through data pointer 

        Parameters
        ----------
        weather_data : Dataframe, optional
            DESCRIPTION. The default is None.Provide weather file if 
            simulating for a specific hour.
        hour_id : 
            DESCRIPTION. The default is None.
            Specific hour number to simulate
        dni : 
            DESCRIPTION. The default is None.
            Flux Simulation DNI
        not_filter_helio : 
            DESCRIPTION. The default is False.
            If True: Doesn't filter heliostats based on efficiency and power
            delivery
        center_aimpoint : 
            DESCRIPTION. The default is True.
            Aimpoint Method 
        read_weather : 
            DESCRIPTION. Default is false for running HALOS optimization. If true: read weather file, for running SP optimization.

        Returns
        -------
        None.

        """
        self.r = cp.data_create()
        cp.reset_vars(self.r)
        ##Field Parameters:
        cp.data_set_number(self.r, "solarfield.0.dni_des", float(self.receiver_data["field_dni"]))   
        cp.data_set_number(self.r, "solarfield.0.q_des", ((float(self.receiver_data["power_rating"]))/10**6)) 
        cp.data_set_number(self.r, "solarfield.0.tht", float(self.receiver_data["tow_height"]))  
        ##Mirror Parameters - (Heliostats) 
        cp.data_set_number(self.r, "heliostat.0.height", float(self.receiver_data["mirror_ht"]))
        cp.data_set_number(self.r, "heliostat.0.width", float(self.receiver_data["mirror_len"]))
        ##Receiver Parameters:
        ## ZENITH DEF: Zenith Angle is measured relative to earth's normal (Z-Axis). (0 to 180) deg
        ## Elevation DEF: Elevation is measured relative to the XY plane. Has positive value towards upwards. (90 to -90) deg
        cp.data_set_number(self.r, "receiver.0.rec_height", float(self.receiver_data["height"]))  
        cp.data_set_number(self.r, "receiver.0.rec_azimuth", 180 - float(self.receiver_data["azimuth_deg"]))   
        cp.data_set_number(self.r, "receiver.0.rec_elevation",90 - float(self.receiver_data["zenith_deg"]))    
        cp.data_set_number(self.r, "receiver.0.rec_offset_x",float(self.receiver_data["rec_cent_offset_x"]))   
        cp.data_set_number(self.r, "receiver.0.rec_offset_y",float(self.receiver_data["rec_cent_offset_y"]))   
        cp.data_set_number(self.r, "receiver.0.rec_offset_z",float(self.receiver_data["rec_cent_offset_z"])) 
        if self.receiver_data["receiver_type"] == 'Flat plate':
            cp.data_set_number(self.r, "receiver.0.rec_width", float(self.receiver_data["length"]))   
        if self.receiver_data["receiver_type"] == 'External cylindrical':
            cp.data_set_number(self.r, "receiver.0.rec_diameter", float(self.receiver_data["diameter"]))
        #Receiver Type: External cylindrical or Flat plate
        # rec_type = bytes(self.receiver_data["receiver_type"], encoding="utf-8")
        cp.data_set_string(self.r, "receiver.0.rec_type",self.receiver_data["receiver_type"]) 
        ##Flux Parameters
        if dni is not None:
            cp.data_set_number(self.r, "fluxsim.0.flux_dni", dni)
        elif self.receiver_data.get('flux_dni') is not None: 
            cp.data_set_number(self.r, "fluxsim.0.flux_dni", float(self.receiver_data["flux_dni"]))
        try:
            cp.data_set_number(self.r, 'fluxsim.0.x_res', float(self.receiver_data["pts_per_len_dim"]))
            cp.data_set_number(self.r, 'fluxsim.0.y_res', float(self.receiver_data["pts_per_ht_dim"]))
        except KeyError:
            cp.data_set_number(self.r, 'fluxsim.0.x_res', float(self.receiver_data["pts_per_dim"]))
            cp.data_set_number(self.r, 'fluxsim.0.y_res', float(self.receiver_data["pts_per_dim"]))
        if hour_id is not None:
            if read_weather == True:
                weather_data = flux_model.ReadWeatherFile(weather_data)
            cp.data_set_number(self.r, "fluxsim.0.flux_day", weather_data['day'][hour_id])
            cp.data_set_number(self.r, "fluxsim.0.flux_hour", weather_data['hour'][hour_id])
            cp.data_set_number(self.r, "fluxsim.0.flux_month", weather_data['month'][hour_id])
            cp.data_set_number(self.r, "fluxsim.0.flux_dni", weather_data['dni'][hour_id])
        #Weather File 
        # w_file = bytes(self.filenames["weather_filename"], encoding="utf-8")
        cp.data_set_string(self.r, "ambient.0.weather_file", self.filenames["weather_filename"]) 
        if not_filter_helio is True:
            cp.data_set_string( self.r, "solarfield.0.des_sim_detail", "Do not filter heliostats" )
        if center_aimpoint is True:
            cp.data_set_string(self.r, "fluxsim.0.aim_method", "Simple aim points")
        
    def plot_field(self, field, name = None):
        x = field["x_location"].values.flatten()
        y = field["y_location"].values.flatten()
        plt.scatter(x, y, s=1.5)
        plt.tight_layout()
        if name is not None:
            plt.savefig(name) 
        
    def plot_flux_map(self, flux, name = None):
        im = plt.imshow(flux)
        plt.colorbar(im)
        plt.tight_layout()
        if name is not None:
            plt.savefig(name)
        
class SP_Field(SolarPilot):
    def __init__(self, filenames):
        super().__init__(filenames)
    
    def generate_field(self, not_filter_helio = False, plotting = False):
        """
        Generates Field Using SolarPilot

        Parameters
        ----------
        not_filter_helio : 
            DESCRIPTION. The default is False.
            If True: Doesn't filter heliostats based on efficiency and power
            delivery
        Returns
        -------
        solar_field : Dataframe with Heliostat location and performance parameters

        """
        self.assign_inputs(not_filter_helio = not_filter_helio)
        cp.generate_layout(self.r, nthreads = 8)
        print("Field Generated by SolarPilot")                                          
        field = cp.get_layout_info(self.r)   #Headers are [index, position-x, position-y, position-z, template, ranking metric value]
        solar_field = cp.detail_results(self.r) #Dataframe 
        solar_field.to_csv("current_field.csv",index=False)
        if plotting is True:
            self.plot_field(field)
        cp.data_free(self.r)
        return  solar_field
    
    def get_sp_layouts(self, weather_data, out_filename):
        """
        get annual power delivery of an Heliostat based on specific number of 
        hours. This method is used to optimize the layout. The field is obtained
        without filtering Heliostats. 

        Parameters
        ----------
        out_filename : CSV, Dataframe
            Contains Heliostat ID, location and annual power delivery to the receiver

        Returns
        -------
        None.

        """
        hour_id = []
        for i in range(len(weather_data["dni"])):
            if (weather_data["month"][i], weather_data["day"][i]) in [(2,3), (5,5), (8,4), (11,4)] and weather_data["dni"][i] >= 500:
                hour_id.append(i)      
        field = self.generate_field(not_filter_helio = True)
        sp_flux = SP_Flux(self.filenames, field = field)
        results = {}
        for i in range(len(hour_id)):
            print("Hour Number: ", hour_id[i])
            sp_flux.get_full_field_flux(weather_data,None,hour_id[i], dni =None, not_filter_helio = True, center_aimpoint = True)
            field_df = sp_flux.df_field
            results[hour_id[i]] = field_df
        #aggregate into annual energy estimate
        power_sum = []   
        for h in range(len(field)):
            sum_h = 0
            for hour in hour_id:
                sum_h = sum_h + results[hour]["power_to_receiver"].values[h]
            power_sum.append(sum_h)  
        print("results aggregated")
        #send to new dataframe
        annual_results = pandas.DataFrame(columns = ["id","X-pos","Y-pos","Z-pos","annual_power"])
        annual_results["id"] = field["id"]
        annual_results["X-pos"] = field["x_location"]
        annual_results["Y-pos"] = field["y_location"]
        annual_results["Z-pos"] = field["z_location"]
        annual_results["annual_power"] = power_sum
        annual_results.set_index("id")
        #get used field, add indicator to what's in the actual field
        annual_results["in_field"] = 0
        sp_field = self.generate_field(not_filter_helio = False)
        ID = sp_field["id"].values
        for i in range(len(ID)):
            annual_results.at[ID[i], "in_field"] = 1
        annual_results.to_csv(out_filename)
    
class SP_Flux(SolarPilot):
    def __init__(self, filenames, field = None, use_sp_field = False):
        super().__init__(filenames)
        self.field = field
        self.use_sp_field = use_sp_field
        self.get_helio_coord()
        
    def get_helio_coord(self):
        """
        Converts the Coordinates from field generated/assigned for simulation 
        into solarpilot template for flux simulation 

        Returns
        -------
        None.

        """
        if self.field is not None:
            try:
                x = self.field.x
                y = self.field.y
                z = self.field.z
            except AttributeError:
                x = self.field["x_location"].values
                y = self.field["y_location"].values
                z = self.field["z_location"].values
        else:
            if self.use_sp_field == True:
                spf = SP_Field(self.filenames)
                self.solar_field = spf.generate_field()
                x = self.solar_field["x_location"].values
                y = self.solar_field["y_location"].values
                z = self.solar_field["z_location"].values
            else:
                print("Reading Field From CSV file")
                df = pandas.read_csv(self.filenames["field_filename"])
                x = df["Pos-x"].values
                y = df["Pos-y"].values
                z = df["Pos-z"].values
        helio_template = [0]*len(x)                        #adding zeroes for SP API (Template)
        coords = numpy.array([helio_template,x,y,z]).transpose()  
        self.num_heliostats = len(x)
        ## simulate whole field to get SolarPILOTS aimpoints
        self.helio_data = []
        for h in range(self.num_heliostats):
            helio = coords[h].tolist()
            self.helio_data.append(helio)
                
    def get_single_helio_flux(self, helio_index, weather_data = None, hour_id = None, dni = None):
        """
        Returns Flux image/map of a single Heliostat 

        Parameters
        ----------
        helio_index : Heliostat ID
        weather_data : Dataframe, Required if simulating specific hour
             
        hour_id : 
            DESCRIPTION. Hour Number
        dni : 
            DESCRIPTION. Flux simulation dni

        Returns
        -------
        flux : Array
            Flux image of a single Heliostat

        """
        self.assign_inputs(weather_data, hour_id, dni)
        #Assiging layout and starting flux simulation
        h = []
        h.append(self.helio_data[helio_index])
        cp.assign_layout(self.r, h)                                                 
        field = cp.get_layout_info(self.r)                                         
        cp.simulate(self.r, nthreads = 8)                                                     
        flux = cp.get_fluxmap(self.r)
        cp.data_free(self.r)
        return flux
    
    def get_single_helio_flux_dict(self,aim_method = 'Simple aim points', weather_data = None, hour_id = None, dni = None):
        """
        Get a dictionary with Heliostat ID as key and flux image per heliostat
        as value. 

        Parameters
        ----------
        aim_method : 
            DESCRIPTION. Aimpoint Method

        Returns
        -------
        flux_dict : Dictionary 

        """
        self.assign_inputs(weather_data, hour_id, dni, read_weather=True)
        cp.data_set_string(self.r, "fluxsim.0.aim_method", aim_method)
        cp.data_set_string(self.r, "solarfield.0.des_sim_detail", "Do not filter heliostats")
        cp.assign_layout(self.r, self.helio_data)
        cp.simulate(self.r)
        ## Build modify field object
        res = cp.detail_results(self.r) 
        helio_dict = {}
        helio_dict['id'] = res['id'].tolist()
        helio_dict['aimpoint-x'] = res['x_aimpoint'].tolist()
        helio_dict['aimpoint-y'] = res['y_aimpoint'].tolist()
        helio_dict['aimpoint-z'] = res['z_aimpoint'].tolist()
        helio_dict['enabled'] = [0]*len(res['id']) # Start with all heliostats off
        flux_dict = {}
        for i,helio in enumerate(helio_dict['id']):
            helio_dict['enabled'][i] = 1 # Turn on a heliostat
            if i > 0:
                helio_dict['enabled'][i-1] = 0  # Turn off the previous
            # Modifying layout and starting flux simulation
            cp.modify_heliostats(self.r, helio_dict)
            cp.simulate(self.r,update_aimpoints=False)
            flux = cp.get_fluxmap(self.r)  
            flux_dict[helio] = flux   
            if i%10 == 0:
                print("Heliostat Field Flux Characterization {:5.2f} % Complete".format(100.*i/self.num_heliostats)) 
                if False:
                    im = plt.imshow(flux)
                    plt.colorbar(im)
                    plt.tight_layout()
                    plt.show()
        cp.data_free(self.r)
        return flux_dict
    
    def get_full_field_flux(self, weather_data = None, case_name = None, hour_id = None, dni = None, 
                            not_filter_helio = False, center_aimpoint = False):
        """
        Returns full field flux map
        
        Returns
        -------
        flux : 
            Flux image of the full field
        num_defocus : 
            Number of Heliostat Removed from field by Solarpilot 
        num_heliostats : 
            Total number Heliostats in the field

        """
        self.assign_inputs(weather_data, hour_id, dni, not_filter_helio, center_aimpoint, read_weather=True)
        cp.assign_layout(self.r, self.helio_data)
        cp.simulate(self.r, nthreads = 8)                                                             
        flux = cp.get_fluxmap(self.r)                                                
        field = cp.get_layout_info(self.r)  
        num_defocus = self.num_heliostats - len(field) - 1
        num_heliostats = len(field)
        res1 = cp.summary_results(self.r, save_dict=True)
        self.df_field = cp.detail_results(self.r) 
        self.full_field_flux = flux
        cp.data_free(self.r)
        return flux, num_defocus, num_heliostats
    
    def sp_aimpoint_defocus(self, weather_data = None, hour_id = None, dni = None):
        """
        Aimpoint heuristics using SolarPILOT flux calculations. Defocuses 
        Heliostats to avoid flux violations. 

        Returns
        -------
        flux_before_defocus : 
            Flux image before heuristic (defocusing)
        flux_after_defocus : 
            Flux image after heuristic (defocusing)
        defocused_helios : 
            Number of Heliostats defocused

        """
        import sp_aimpoint_heuristic as sp_aim
        import pickle
        
        field_flux_dict_file = "field_flux_image_priority.pkl"
        load_data = False
        if not load_data:
            flux_dict = self.get_single_helio_flux_dict(aim_method='Image size priority', weather_data=weather_data, hour_id=hour_id, dni=dni)
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
        comb_flux = numpy.zeros((len(flux_dict[0]), len(flux_dict[0][0])))
        for key in flux_dict.keys():
            comb_flux += numpy.array(flux_dict[key])
        flux_before_defocus = numpy.array(comb_flux)
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
    
        comb_flux, defocused_helios = sp_aim.flux_limit_heuristic(flux_dict, flux_limit)
        # print("Number of Heliostats Defocused: " + str(len(defocused_helios)))
        defocused_helios = len(defocused_helios)
        flux_after_defocus = numpy.array(comb_flux)
        return flux_before_defocus, flux_after_defocus, defocused_helios       
    
    def run_sp_case(self, weather_data = None, case_name = None, hour_id = None, dni = None, 
                    sp_aimpoint_heur = False, saveCSV = True):
        """
        Run full field SolarPilot Case with or without aimpoint heuristic

        Parameters
        ----------

        sp_aimpoint_heur : 
            DESCRIPTION. If True, uses sp_aimpoint_defocus to avoid flux
            violations
        saveCSV : The default is True.
            DESCRIPTION. Save results as CSV

        Returns
        -------
        results : Dict
            DESCRIPTION. Contains power delivery to receiver (obj_value), 
            maximum flux, number of helios, and number of helios filtered/defocus
            by solarpilot. If sp_aimpoint_heur is True also returns before and 
            after defocused obj_value. 

        """
        
        if self.receiver_data["receiver_type"] == 'Flat plate':
            area_m_point = (float(self.receiver_data["height"])/float (self.receiver_data["pts_per_dim"])) * (float(self.receiver_data["length"])/float(self.receiver_data["pts_per_dim"]))
        if self.receiver_data["receiver_type"] == 'External cylindrical':
            area_m_point = ((((float(self.receiver_data["diameter"])) * (numpy.pi))/(float (self.receiver_data["pts_per_dim"]))) * ((float(self.receiver_data["height"]) / (float(self.receiver_data["pts_per_dim"])))))
        flux, num_defocus, num_heliostats = self.get_full_field_flux(weather_data=weather_data, case_name=case_name, hour_id=hour_id)
        flux_with_area = numpy.array(flux) * area_m_point
        obj_value = numpy.sum(flux_with_area)
        max_flux = numpy.max(flux)
        results = {}
        results["num_defocus"] = num_defocus
        results["num_heliostats"] = num_heliostats
        results["max_flux"] = max_flux
        results["obj_value"] = obj_value
        if sp_aimpoint_heur:
            print("HALOS & SP Done! - SP Aimpoint Heuristic Started.......")
            before_flux, after_flux, defocused_helios = self.sp_aimpoint_defocus(weather_data = weather_data, hour_id = hour_id, dni = dni)
            before_obj = numpy.sum(before_flux*area_m_point)
            post_heur_obj = numpy.sum(after_flux*area_m_point)
            results["before_flux"] = before_obj
            results["post_heur_obj"] = post_heur_obj
            results["post_num_defocus"] = defocused_helios
        if saveCSV:
            results_filename = "SolarPilot_" + case_name +"_summary.csv"
            with open(results_filename, "w") as f:  
                w = csv.DictWriter(f, results.keys())
                w.writeheader()
                w.writerow(results) 
        
        if case_name is not None:
            plt.figure(1)
            self.plot_flux_map(flux, name = ("SolarPilot_Flux_map_" + case_name))
            if sp_aimpoint_heur:
                plt.figure(2)
                self.plot_flux_map(before_flux, name = (case_name+"_before_defocus_SP"))
                plt.figure(3)
                self.plot_flux_map(after_flux, name = (case_name+"_After_defocus_SP"))
        return results
    
                        
if __name__ == "__main__":
    import inputs
    import field

    case_filename = "./../case_inputs/radial_50_ca_case.csv"
    #case_filename = "./../case_inputs/radial_tonopah.csv"
    filenames = inputs.readCaseFile(case_filename)
    sp_field = SP_Field(filenames)
    field = sp_field.generate_field(plotting= True)
    sp_flux = SP_Flux(filenames,field)
    case_name = "50_MW_Radial_dagget"
    flux, num_defocus, num_heliostats = sp_flux.get_full_field_flux()
    # results = sp_flux.run_sp_case(case_name = case_name, sp_aimpoint_heur = True)
    
    # flux = sp.get_single_helio_flux_dict()
    sp_flux.plot_flux_map(flux)
