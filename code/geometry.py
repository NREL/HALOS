# -*- coding: utf-8 -*-
"""
Geometry Module
Builds different receiver modules with x, y, z coordinates
"""

import numpy as np
import scipy
import math
import csv
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pandas


class Receiver:
    def __init__(self, shape_type, tow_height, params,solar_field= None):
        """
                shape_type: string identifier
                tow_height: optical height of the tower
                params: dictionary of additional characteristics, such as length and
                        height for a flat plate
                solar_field: x,y,z location of Heliostats
        """
        self.solar_field = solar_field
        self.shape_type = shape_type
        self.tow_height = tow_height  # distance from ground to center of receiver
        self.params = params
        self.x = np.array([])
        self.y = np.array([])
        self.z = np.array([])
        self.normal = np.array([])
        # sets discretization for length and height dimension if uniform is given
        if "pts_per_len_dim" not in self.params or "pts_per_ht_dim" not in self.params:
            self.params["pts_per_len_dim"] = self.params["pts_per_dim"]  # length or circumferential dimension
            self.params["pts_per_ht_dim"] = self.params["pts_per_dim"]  # height dimension
        self.dni = params.get('flux_dni')


    def check_inputs(self, shape_type, tow_height, params):
        if not isinstance(params, dict):
            raise Warning("'Params' should be a dictionary")


    def getOpticalHeight(self, Heliostat):
        return self.tow_height - Heliostat._pylon_height


    def getCoords(self):
        """
        Converts x,y,z coordinates into a single array

        Returns
        -------
        None. Populates self.coords

        """
        self.coords = np.array([self.x.reshape(self.x.size),
                                   self.y.reshape(self.y.size),
                                   self.z.reshape(self.z.size)]).transpose()
        self.getNormals()


    def getAimpoints(self):
        """
        Converts aimpoint x,y,z into single array. 

        Returns
        -------
        None. Assigns the single array to self.aimpoints and the number of 
        aimpoints to self.num_aimpoints

        """
        self.aimpoints = np.array([self.aim_x.reshape(self.aim_x.size),
                                      self.aim_y.reshape(self.aim_y.size),
                                      self.aim_z.reshape(self.aim_z.size)]).transpose()
        self.num_aimpoints = self.aim_x.size


    def generateDynamicFluxLimits(self, lb, ub, Np):
        """
        Creates Dynamic flux limits on the receiver. The flux limit is 
        maximum when the fluid temperature is lowest. 
        
        Assumes flow enters from Top

        Parameters
        ----------
        lb : Lower Bound for Flux Limit - Minimum 
        ub : Upper bound for Flux Limit - Maximum 
        Np : Number of Fluid Circulations

        Returns
        -------
        None. Saves flux limits as CSV and Populates self.flux_upper_limits

        """
        def createposMat(Np, Nx, Ny):
            """
            Creates Flux Measurement point matrix

            Parameters
            ----------
            Np : Number of Fluid Circulations
            Nx : Number of Measurement point in the x direction - Horizontal
            Ny : Number of Measurement point in the y direction - Vertical

            Returns
            -------
            TYPE
                Position Matrix

            """
            ## creating half the matrix
            if (Nx % 2 != 0):
                print("ERROR: Nx must be divisible by 2")
                return False
            if (Nx % Np != 0):
                print("ERROR: Nx must be divisible by Np")
                return False
        
            NrepeatCols = int((Nx/2)/Np)                    # Number of repeat columns 
            NpointsPath = int((Nx/2)/NrepeatCols*Ny)        # number of points on receiver path
        
            # This method assume flow enters the top of the receiver
            posM = []
            for y in range(Ny):
                posM.append([])
                for x in range(int((Nx/2)/NrepeatCols)):
                    if (x % 2 == 0):    # this is assuming at serpentine flow
                        count = y + x*Ny
                    else:
                        count = (Ny - y) + x*Ny
                    for j in range(NrepeatCols):
                        posM[y].append(count/NpointsPath)
            
            RposM = np.array(posM)
            LposM = np.flip(RposM,1)
            posM = np.concatenate((LposM,RposM),axis=1)
            return posM   
        
        from scipy.interpolate import interp1d
        Nx = self.params["pts_per_len_dim"]               #36 # Number of measurement points in x-dir
        Ny = self.params["pts_per_ht_dim"]               #24 # Number of measurement points in y-dir

        pos = []
        AFD = []
        with open('AFD_flowpath.csv') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=',')
            line_count = 0 
            for row in csv_reader:
                if len(row) == 2:
                    pos.append(float(row[0]))
                    AFD.append(float(row[1]))
        #End effects
        pos.insert(0, 0)
        AFD.insert(0, 1)
        pos.append(1)
        AFD.append(AFD[-1])
        ## create a position maxtrix (normalized)
        posM = createposMat(Np, Nx, Ny)
        f = interp1d(pos, AFD)
        AFD_M = []
        for j in range(len(posM)):
            AFD_M.append([])
            for i in range(len(posM[0])):
                AFD_M[j].append(f(posM[j][i])*ub)
        AFD_M = np.array(AFD_M)
        self.flux_lower_limits = np.ones(self.x.size, dtype=float) * lb
        self.flux_upper_limits = AFD_M.flatten()
        AFD = pandas.DataFrame(AFD_M)
        AFD.to_csv('flux_limits.csv', index = False)


    def generateFixedFluxLimits(self, lb, ub):
        """
        Generate Constant Flux Limits 

        Parameters
        ----------
        lb : Lower Bound for Flux Limit - Minimum 
        ub : Upper bound for Flux Limit - Maximum 

        Returns
        -------
        None. populates self.flux_upper_limits

        """
        self.flux_lower_limits = np.ones(self.x.size, dtype=float) * lb
        self.flux_upper_limits = np.ones(self.x.size, dtype=float) * ub


    def getSurfaceArea(self):
        raise ValueError("Inherited classes must overwrite 'getSurfaceArea'")


    def getNormals(self):
        """temporary: everything faces north"""
        self.normals = self.coords * 1.0
        self.normals[:, 0] = 0.
        self.normals[:, 1] = -1.
        self.normals[:, 2] = 0.


    def generateAimpointsGrid(
            self, num_rows, num_cols, h_margin=0.,
            v_margin=0, epsilon=1e-6
    ):
        """
        Generates candidate aimpoints on the receiver surface. uses the
        measurement points on the surface.  The aimpoints are obtained
        via linear interpolation if the aimpoints are not already existing
        measurement points.

        ===== Parameters =====
        num_horizontal -- number of horizontal aimpoints
        num_cols -- number of vertical aimpoints
        h_margin -- number of horizontal-axis measurement points to
            skip (on both sides)
        v_margin -- number of z-axis measurement points to skip on each
            side (on both sides)

        ===== Returns =====
        None (assigns output to self.aimpoints)

        """
        shape = self.x.shape
        self.aim_x = np.zeros([num_rows, num_cols], dtype=float)
        self.aim_y = np.zeros_like(self.aim_x)
        self.aim_z = np.zeros_like(self.aim_x)
        for h in range(num_rows):
            eff_row = h_margin + ((h + 0.5) * (shape[0] - 1 - 2 * h_margin) / num_rows)
            row_start = int(eff_row)
            for v in range(num_cols):
                eff_col = v_margin + ((v + 0.5) * (shape[1] - 1 - 2 * v_margin) / num_cols)
                col_start = int(eff_col)
                self.aim_x[h, v] = self.x[row_start, col_start]
                self.aim_y[h, v] = self.y[row_start, col_start]
                self.aim_z[h, v] = self.z[row_start, col_start]
                if eff_col - col_start > 1e-4:
                    if eff_row - row_start > 1e-4:
                        # interpolate on row and column
                        p_row = eff_row - row_start
                        p_col = eff_col - col_start
                        self.aim_x[h, v] += p_row * (self.x[row_start + 1, col_start] - self.x[row_start, col_start])
                        self.aim_y[h, v] += p_row * (self.y[row_start + 1, col_start] - self.y[row_start, col_start])
                        self.aim_z[h, v] += p_row * (self.z[row_start + 1, col_start] - self.z[row_start, col_start])
                        self.aim_x[h, v] += p_col * (self.x[row_start, col_start + 1] - self.x[row_start, col_start])
                        self.aim_y[h, v] += p_col * (self.y[row_start, col_start + 1] - self.y[row_start, col_start])
                        self.aim_z[h, v] += p_col * (self.z[row_start, col_start + 1] - self.z[row_start, col_start])
                    else:
                        # interpolate on column but not row
                        p_col = eff_col - col_start
                        self.aim_x[h, v] += p_col * (self.x[row_start, col_start + 1] - self.x[row_start, col_start])
                        self.aim_y[h, v] += p_col * (self.y[row_start, col_start + 1] - self.y[row_start, col_start])
                        self.aim_z[h, v] += p_col * (self.z[row_start, col_start + 1] - self.z[row_start, col_start])
                elif eff_row - row_start > 1e-4:
                    # interpolate on row but not column
                    p_row = eff_row - row_start
                    self.aim_x[h, v] += p_row * (self.x[row_start + 1, col_start] - self.x[row_start, col_start])
                    self.aim_y[h, v] += p_row * (self.y[row_start + 1, col_start] - self.y[row_start, col_start])
                    self.aim_z[h, v] += p_row * (self.z[row_start + 1, col_start] - self.z[row_start, col_start])
        self.getAimpoints()


    def transReceiver(self, trans):
        # moves receiver location provided by user
        # expects trans an array of length 3 (x,y,z)
        self.x += trans[0]
        self.y += trans[1]
        self.z += trans[2]


    def build_zero_arrays(self):
        self.x = np.zeros([self.params["pts_per_ht_dim"], self.params["pts_per_len_dim"]], dtype=float)
        self.y = np.zeros([self.params["pts_per_ht_dim"], self.params["pts_per_len_dim"]], dtype=float)
        self.z = np.zeros([self.params["pts_per_ht_dim"], self.params["pts_per_len_dim"]], dtype=float)


    def build_measurement_points(self):
        raise ValueError("Inherited classes must overwrite 'build_measurement_points")


class FlatPlateReceiver(Receiver):
    def __init__(self, tow_height, params):
        super().__init__("flat_plate", tow_height, params)
        self.build_measurement_points()
        self.getSurfaceArea()
        try:
            self.generateAimpointsGrid(int(params["aim_rows"]),int(params["aim_cols"]),int(params["aim_h_margin"]),int(params["aim_v_margin"]))
        except KeyError:
            pass
        self.dni = params.get('flux_dni')


    def getSurfaceArea(self):
        """
        Calculates area per measurement point

        """
        """This is for the flat plate only"""
        area_per_point = self.params["length"] * self.params["height"] / (
                    self.params["pts_per_len_dim"] * self.params["pts_per_ht_dim"])
        self.surface_area = np.ones(self.x.size, dtype=float) * area_per_point


    def build_measurement_points(self):
        """
        Builds Measurement Points (x,y,z) on the Receiver Surface

        Returns
        -------
        None. Assigns coordinates of points to self.x, self.y, self.z and calls 
        get_coords

        """
        self.shape_type = "flat_plate"
        self.build_zero_arrays()
        # This assume receiver normal is [0, 1, 0] i.e., points directly north
        xs = (np.arange(self.params["pts_per_len_dim"], dtype=float) + 0.5) / (self.params["pts_per_len_dim"]) * \
             self.params["length"] - self.params["length"] / 2
        zs = (np.arange(self.params["pts_per_ht_dim"], dtype=float) + 0.5) / (self.params["pts_per_ht_dim"]) * \
             self.params["height"] - self.params["height"] / 2 + self.tow_height
        self.x += xs

        self.z = self.z.transpose()
        self.z += zs

        self.z = self.z.transpose()
        self.z = self.z[::-1]  # flips matrix for plotting to be consistent with z direction
        self.getCoords()
        # Determines zenith and azimuth angles if normal vector is given
        if "normal" in self.params:
            # note: normal is inward normal
            # calculates zenith and azimuth angles [radians] from normal provided by user
            norm = self.params["normal"]
            # normalize normal vector if not already done so by user
            norm_len = np.sqrt(np.sum(norm * norm))
            if not norm_len == 1:
                norm = norm / norm_len
                self.params["normal"] = norm  # returning normalized normal

            # calculate zenith angle (range [0, pi]), 0 = pointing up
            z = np.array([0, 0, 1])  # zenith axis
            self.params["zenith"] = math.acos(np.sum(norm * z))

            # calculate azimuth angle (defined as clockwise from south-facing (+ y-axis)) (range [0, 2*pi])
            north = np.array([0, 1, 0])  # y-axis or (north)
            norm_proj_xy = norm * np.array([1, 1, 0])
            if not sum(norm_proj_xy) == 0:
                self.params["azimuth"] = math.acos(
                    np.sum(norm_proj_xy * north) / np.sqrt(np.sum(norm_proj_xy * norm_proj_xy)))
            else:
                # assume azimuth equal to zero, this is only the case when receiver points directly down or up (norm [0,0,+-z])
                self.params["azimuth"] = 0.0

            # determine the appropriate quadrant to correct azimuth angle
            if norm[0] < 0:
                self.params["azimuth"] = 2 * np.pi - self.params["azimuth"]

        # Determines normal vector if zenith and azimuth angles are given (radians or degrees)
        elif ("zenith" in self.params and "azimuth" in self.params) or (
                "zenith_deg" in self.params and "azimuth_deg" in self.params):

            # converts degrees to radians
            if "zenith_deg" in self.params:
                # if given in degrees convert to radians
                self.params["zenith"] = self.params["zenith_deg"] * np.pi / 180
                self.params["azimuth"] = self.params["azimuth_deg"] * np.pi / 180

            norm = np.zeros([3, ], dtype=float)
            # calculate normal components
            norm[0] = math.sin(self.params["zenith"]) * math.sin(self.params["azimuth"])
            norm[1] = math.sin(self.params["zenith"]) * math.cos(self.params["azimuth"])
            norm[2] = math.cos(self.params["zenith"])
            self.params["normal"] = norm

        # default surface normal, zenith, and azimuth angles
        else:
            self.params["zenith"] = np.pi / 2
            self.params["azimuth"] = 0.0
            self.params["normal"] = np.array([0, 1, 0])

        for i in range(len(self.coords)):
            self.normals[i] = self.params["normal"]

        # Provides degrees of zenith and azimuth if not given
        if not "zenith_deg" in self.params:
            self.params["zenith_deg"] = self.params["zenith"] * 180 / np.pi
            self.params["azimuth_deg"] = self.params["azimuth"] * 180 / np.pi

        # rotates receiver to provided zenith and azimuth angles
        if not self.params["zenith_deg"] == 90 or not self.params["azimuth_deg"] == 0:
            # translate to the origin
            self.z -= self.tow_height

            zen_rot = self.params["zenith"] - np.pi / 2  # rotation in the zenith angle

            xprime = (math.cos(self.params["azimuth"]) * self.x
                      + math.sin(self.params["azimuth"]) * math.cos(zen_rot) * self.y
                      + math.sin(self.params["azimuth"]) * math.sin(zen_rot) * self.z)

            yprime = (-math.sin(self.params["azimuth"]) * self.x
                      + math.cos(self.params["azimuth"]) * math.cos(zen_rot) * self.y
                      + math.cos(self.params["azimuth"]) * math.sin(zen_rot) * self.z)

            zprime = -math.sin(zen_rot) * self.y + math.cos(zen_rot) * self.z

            # new receiver coordinates
            self.x = xprime
            self.y = yprime
            self.z = zprime + self.tow_height

        # translate receiver if offset provided
        if "rec_cent_offset" in self.params:
            self.transReceiver(self.params["rec_cent_offset"])
            
        self.getCoords()


class CylindricalPlateReceiver(Receiver):
    def __init__(self, tow_height, params,solar_field):
        super().__init__("cylindrical_plate", tow_height, params,solar_field)
        self.check_inputs("cylindrical_plate", tow_height, params)
        self.build_measurement_points()
        self.getSurfaceArea()
        try:
            #self.generateAimpointsGrid_rows_only(int(params["aim_rows"]))
            self.generateAimpointsGrid_rows_only(int(params["aim_rows"]),self.solar_field)
        except KeyError:
            pass
        #self.getSectionFluxColumn()
    def check_inputs(self, shape_type, tow_height, params):
        if not isinstance(params, dict):
            raise Warning("'Params' should be a dictionary")
        if 'diameter' not in params:
            raise ValueError("'diameter' needs to be included")
            
    def getSurfaceArea(self):
        """
        Calculates area per measurement point

        """
        """This is for the Cylindrical Plate Only"""
        #_area = _height * _radius * PI * 2.   => d = 2*r
        print("Cylindrical Receiver")
        area_per_point = self.params["diameter"] * np.pi * self.params["height"]/ (
                    self.params["pts_per_len_dim"] * self.params["pts_per_ht_dim"])
        self.surface_area = np.ones(self.x.size, dtype=float) * area_per_point
        
    def build_measurement_points(self):
        """
        Builds Measurement points for the cylindrical Receiver Case

        Returns
        -------
        None. Calls get_Coords to create single array of coordinates 

        """
        # Defaults if not provided
        if not "hoz_acpt_angle" in self.params:
            self.params["hoz_acpt_angle"] = 360

        if not "azimuth" in self.params and not "azimuth_deg" in self.params:
            self.params["azimuth"] = 0

        if "azimuth" in self.params:
            self.params["azimuth_deg"] = self.params["azimuth"] * 180 / np.pi  # convert to degrees
        elif "azimuth_deg" in self.params:
            self.params["azimuth"] = self.params["azimuth_deg"] * np.pi / 180  # convert to radians

        r = self.params["diameter"] / 2  # receiver radius

        # Builds simple cylindrical receiver (assumes a perfect cylinder)
        self.build_zero_arrays()
        if self.shape_type == "cylindrical_plate":

            # azimuth angle to each point (calculation in degree then converts to radians)
            self.azimuth = (np.arange(self.params["pts_per_len_dim"], dtype=float) + 0.5) / (
                self.params["pts_per_len_dim"]) * self.params["hoz_acpt_angle"]
            self.azimuth -= self.params["hoz_acpt_angle"] / 2  # rotates to center on north

            self.azimuth += self.params["azimuth_deg"]  # rotates to center on users azimuth angle provided
            np.place(self.azimuth, self.azimuth < 0, self.azimuth + 360)  # removes negative angles
            np.place(self.azimuth, self.azimuth > 360, self.azimuth - 360)  # removes angles greater than 360
            self.azimuth *= np.pi / 180  # convert to radians

            xs = r * np.sin(self.azimuth)
            ys = r * np.cos(self.azimuth)
            zs = (np.arange(self.params["pts_per_ht_dim"], dtype=float) + 0.5) / (self.params["pts_per_ht_dim"]) * \
                 self.params["height"] - self.params["height"] / 2 + self.tow_height

            self.x += xs
            self.y += ys

            self.z = self.z.transpose()
            self.z += zs

            self.z = self.z.transpose()
            self.z = self.z[::-1]  # flips matrix for plotting to be consistent with z direction

            # translate receiver if offset provided
            if "rec_cent_offset" in self.params:
                self.transReceiver(self.params["rec_cent_offset"])
            self.getCoords()
        ## BUILD PANELS 
        # Builds cylindrical receiver using panels
        if self.shape_type == "cylindrical_plate_panels":
            print("Number of Panels: ", self.params["n_panels"])
            azi_step = (self.params["hoz_acpt_angle"] / self.params["n_panels"]) * np.pi / 180  # convert to radians

            dis_pc = r * math.cos(azi_step / 2)  # distance to center of panel (for offset)
            panel_length = 2 * r * math.sin(azi_step / 2)  # panel length

            # "pts_per_len_dim" will be used for each panel
            constp_params = self.params.copy()
            constp_params["length"] = panel_length
            constp_params.pop("zenith_deg", None)
            constp_params["zenith"] = np.pi / 2  # 90 degrees

            # self.normal = np.zeros([self.params["n_panels"],3],dtype=float)   # for storing normals
            if 'normal' in self.params:
                self.normal = self.params['normal']
            else:
                self.NormalVectorCylinder()
            # For future reference, this was uncommented by Alex Mikulich because when these classes were all one class,
            # 'self.normal' was initialized in the flat plate function and then used by the cylindrical function. Now,
            # they are in different classes, so self.normal has to be initialized here.

            azi = azi_step / 2 - self.params["hoz_acpt_angle"] * (np.pi / 180) / 2 + self.params[
                "azimuth"]  # starting azimuth angle
            for x in range(self.params["n_panels"]):
                panel_params = constp_params.copy()
                panel_params["azimuth"] = azi
                panel_params["rec_cent_offset"] = [dis_pc * math.sin(azi), dis_pc * math.cos(azi), 0]

                # creating a flat plate receiver per panel
                #[KL] changed tow_height from zero
                panel = FlatPlateReceiver(self.params["tow_height"], panel_params)
                panel.build_measurement_points()
                #panel.getSurfaceArea()
                if not hasattr(self, 'x'):
                    self.x = panel.x.copy()
                    self.y = panel.y.copy()
                    self.z = panel.z.copy()
                    self.normal = np.array([panel.params["normal"]])
                else:
                    self.x = np.concatenate((self.x, panel.x), axis=1)
                    self.y = np.concatenate((self.y, panel.y), axis=1)
                    self.z = np.concatenate((self.z, panel.z), axis=1)
                    self.normal = np.append(self.normal, [panel.params["normal"]], axis=0)

                azi += azi_step
                del panel, x
                panel_params.clear()

            np.place(self.normal, abs(self.normal) < 1e-15, 0)  # makes small numbers equal to zero

            self.z += self.tow_height
            # translate receiver if offset provided
            if "rec_cent_offset" in self.params:
                self.transReceiver(self.params["rec_cent_offset"])
            self.getCoords()
            
    def generateAimpointsGrid_rows_only(self, num_rows,solar_field, epsilon=1e-6):
        """
        Dynamic Aimpoint Building Strategy for Cylindrical Receiver 
        
        Assumes aimpoints to be distributed in a vertical column at the center 
        of Receiver. 
        
        Dynamically moves aimpoint column from center to the surface of the 
        receiver using Heliostat Location. 

        Parameters
        ----------
        num_rows : number of aimpoints rows 
        solar_field : location of Heliostats

        Returns
        -------
        Assigns coordinate array to self.aimpoint by calling getAimpointsForField

        """
        num_cols = 1
        shape = self.x.shape
        self.aim_x = np.zeros([num_rows, num_cols], dtype=float)
        self.aim_y = np.zeros_like(self.aim_x)
        self.aim_z = np.zeros_like(self.aim_x)
        self.aimpoints = []
        for h_idx in range(len(solar_field.x)):
            for h in range(num_rows):
                eff_row = ((h + 0.5) * (shape[0] - 1) / num_rows)
                row_start = int(eff_row)
                for v in range(num_cols):
                    eff_col = ((v + 0.5) * (shape[1] - 1) / num_cols)
                    col_start = int(eff_col)
                    self.aim_z[h, v] = self.z[row_start, col_start]
                    #print(range(len(solar_field.x)))
                    if solar_field.x[h_idx] >= 0.0:
                        self.aim_x[h, v] = self.params["diameter"]//2
                    elif solar_field.x[h_idx] < 0.0:
                        self.aim_x[h, v] = -(self.params["diameter"]//2)
                    if solar_field.y[h_idx] >= 0.0:
                        self.aim_y[h, v] = self.params["diameter"]//2
                    elif solar_field.y[h_idx] < 0.0:
                        self.aim_y[h, v] = -(self.params["diameter"]//2)
                    if eff_col - col_start > 1e-4:
                        if eff_row - row_start > 1e-4:
                            p_row = eff_row - row_start
                            p_col = eff_col - col_start
                            self.aim_z[h, v] += p_row * (self.z[row_start + 1, col_start] - self.z[row_start, col_start])
                            self.aim_z[h, v] += p_col * (self.z[row_start, col_start + 1] - self.z[row_start, col_start])
                        else:
                        # interpolate on column but not row
                            p_col = eff_col - col_start
                            self.aim_z[h, v] += p_col * (self.z[row_start, col_start + 1] - self.z[row_start, col_start])
                    elif eff_row - row_start > 1e-4:
                    # interpolate on row but not column
                        p_row = eff_row - row_start
                        self.aim_z[h, v] += p_row * (self.z[row_start + 1, col_start] - self.z[row_start, col_start])
            #self.getAimpoints()
            aimpoints_helio = np.array([self.aim_x.reshape(self.aim_x.size),
                                      self.aim_y.reshape(self.aim_y.size),
                                      self.aim_z.reshape(self.aim_z.size)]).transpose() 
            self.aimpoints.append(aimpoints_helio)
        self.getAimpointsForField()


    def getAimpointsForField(self):
        """
        Get Aimpoints' location as an array for Cylindrical Case

        """
        self.aimpoints = self.aimpoints
        self.num_aimpoints = self.aim_x.size
      
    def NormalVectorCylinder(self):
        self.params['n_panels'] = int(self.params['n_panels'])
        angles = np.linspace(0, 2 * math.pi, self.params['n_panels'], endpoint=False)
        self.normal = [self.CalculateNormalVector([np.cos(angle), np.sin(angle), 0]) for angle in angles]

    def CalculateNormalVector(self, Hloc):
        """
        This subroutine should be used to calculate the normal vector to the receiver for a given heliostat location.
        Ultimately, the optical calculations should not use this method to calculate the normal vector. Instead, use
        the normal vector that is assigned to the receiver surface during setup.
        In the case of continuous cylindrical surfaces, this method can be called during optical calculations.
        Given a heliostat at point Hloc{x,y,z}, return a normal vector to the receiver absorber surface.
        """
        # TODO: is this necessary?
        try:
            magnitude = math.sqrt(Hloc[0]**2 + Hloc[1]**2 + Hloc[2]**2)
            return [Hloc[0]/magnitude, Hloc[1]/magnitude, Hloc[2]/magnitude]
        except ZeroDivisionError:
            self.normal = np.zeros(self.params['n_panels'], 3)
            raise Warning("'hloc' has zero magnitude, setting normal vector to [0,0,0].")


    def GetFractionUsingEff(self,section_id):
        """
        Optional Method: Not used by default  
        Builds Fraction Map for Cylindrical Case using Heliostat Efficiency
        
        The model currently uses fraction Maps generated using Flux per Heliostat
        in the flux_model.py

        Parameters
        ----------
        section_id : Section Number 

        Returns
        -------
        Array with fraction map of the section 

        """
        self.section_flux = [
            sum([self.solar_field.eff[idx] for idx in self.solar_field.helios_by_section[s]])
            / sum(self.solar_field.eff) for s in range(self.solar_field.num_sections)]
        
        self.fraction_maps = []
        sum_cols = []           #Stores lists - one list contains sum of 1 column of all sections
        for ncol in range(self.params["pts_per_len_dim"]):
            sums = []           #Stores sum of one specific column of each section
            for s in range(self.solar_field.num_sections):
                fraction_map = np.zeros_like(self.x)+self.section_flux[s]
                fraction_map = pandas.DataFrame(fraction_map)
                sum_col = sum(fraction_map[ncol])
                sums.append(sum_col)
            sum_cols.append(sums)       
        
        for s in range(self.solar_field.num_sections):
            map_fraction = pandas.DataFrame(np.zeros_like(self.x))
            for ncol in range(self.params["pts_per_len_dim"]): 
                map_fraction[ncol] = map_fraction[ncol] + (sum_cols[ncol][s]/sum(sum_cols[ncol]))
            self.map_fraction = np.array(map_fraction)
            self.fraction_maps.append(self.map_fraction.flatten())   
        return self.fraction_maps[section_id]
       

if __name__ == "__main__":
    # params = {"length": 40, "height": 25, "pts_per_dim": 50}

    # params = {"length": 40, "height": 25, "pts_per_len_dim": 50,
    #           "pts_per_ht_dim": 30}  # testing dimensional discretization

    # norm = np.array(
    #     [0, np.sqrt(2) / 2, -np.sqrt(2) / 2])  # 45 degrees tilt down (zenith = 135 deg) (normalized)
    # # norm = np.array([-1, 1, -np.sqrt(2)])  # 45 deg tilt down 45 deg azimuth (not normalized)
    # # norm = np.array([1, 1, 0])  # 45 deg azimuth (not normalized)
    # # norm = np.array([0, 0, 1])  # 45 deg azimuth (not normalized)
    # params = {"length": 40, "height": 25, "pts_per_dim": 50, "normal": norm}  # testing normal calculations

    # params = {"length": 40, "height": 25, "pts_per_dim": 50, "zenith_deg": 135, "azimuth_deg": 90,
    #           "rec_cent_offset": [-20, 10, 0]}  # testing zenith and azimuth calculations
    # h = 185
    # r = FlatPlateReceiver(h, params)
    # r.build_measurement_points()

    params = {"diameter": 6, "height": 17, "pts_per_dim": 15,
              "rec_cent_offset": [0, 0, 0]}  # testing simple external receiver with offset
    #params = {"diameter": 40, "height": 25, "pts_per_dim": 50, "hoz_acpt_angle": 180, "azimuth_deg": 15,
     #         "rec_cent_offset": [0, 0, 0]}  # testing horizontal acceptance angle
    h = 150
    #Test Helios (x,y,z)
    h1 = [10,15,0]
    h2 = [-12,-18,0]
    h3 = [14,-17,0]
    h4 = [-170,87,0]
    #Test field
    field = [h1,h2,h3,h4]
    #shape_type = "simple_ext_cylinder"
    r = CylindricalPlateReceiver(h, params,field)
    r.build_measurement_points()
    r.generateAimpointsGrid_rows_only(3,field)
    print("Number of aimpoints: ",r.num_aimpoints)
    # params = {"n_panels": 6, "diameter": 40, "height": 25, "pts_per_len_dim": 10,
    #           "pts_per_ht_dim": 30}  # testing panel external receiver with offset
    # params = {"n_panels": 3, "diameter": 40, "height": 25, "pts_per_len_dim": 10,
    #           "pts_per_ht_dim": 30}  # triangle receiver
    # params = {"n_panels": 4, "diameter": 40, "height": 25, "pts_per_len_dim": 10, "pts_per_ht_dim": 30,
    #           "hoz_acpt_angle": 90, "azimuth_deg": 45,
    #           "rec_cent_offset": [-20, 10, 0]}  # testing horizontal acceptance angle
    # # params = {"n_panels": 4,"diameter":40, "height":25, "pts_per_len_dim":10, "pts_per_ht_dim": 30, "azimuth_deg":45}    #square receiver

    # h = 150
    # r = CylindricalPlateReceiver(h, params)
    # r.build_measurement_points()
    # r.generateAimpointsGrid(5, 5, 0, 2)

    # print("Measurement point x coordinates:\n", r.x)
    # print("\n\nMeasurement point y coordinates:\n", r.y)
    # print("\n\nMeasurement point z coordinates:\n", r.z)
    print("\n\nAimpoint x coordinates:\n", r.aim_x)
    print("\n\nAimpoint y coordinates:\n", r.aim_y)
    print("\n\nAimpoint z coordinates:\n", r.aim_z)

    # set to true to plot receiver points (3D and axes projections)
    isplotting = True

    if isplotting:
        k = 4  # normal scalar multiplier

        # 3D plotting to test code
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        # ax.set_aspect('equal')  #this functionality doesn't work with (projection='3d')
        # ax.plot_surface(r.x, r.y, r.z)
        ax.scatter3D(r.x, r.y - 20, r.z)
        ax.axes.get_xaxis().set_ticklabels([])
        ax.axes.get_yaxis().set_ticklabels([])
        ax.axes.w_zaxis.set_ticklabels([])
        #        ax.set_xlabel("north-south coordinate")
        #        ax.set_ylabel("east-west coordinate")
        #        ax.set_zlabel("elevation")
        if "normal" in params:
            ax.plot3D([0, k * r.params["normal"][0]], [0, k * r.params["normal"][1]],
                      [h, k * r.params["normal"][2] + h], 'r')
            # normal looks off in the plot because aspect of axes are not equal (there are ways to get a 3D plot to have equal aspect axes)
        #       https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to

        # y-z plane
        fig = plt.figure()
        ax = fig.gca()
        ax.axis('equal')
        ax.scatter(r.y, r.z)
        if "normal" in params:
            ax.plot([0, k * r.params["normal"][1]], [h, k * r.params["normal"][2] + h], 'r')
        plt.title("Y-Z plane projection")

        # x-z plane
        fig = plt.figure()
        ax = fig.gca()
        ax.axis('equal')
        ax.scatter(r.x, r.z)
        if "normal" in params:
            ax.plot([0, k * r.params["normal"][0]], [h, k * r.params["normal"][2] + h], 'r')
        plt.title("X-Z plane projection")

        # x-y plane
        fig = plt.figure()
        ax = fig.gca()
        ax.axis('equal')
        ax.scatter(r.x, r.y)
        if "normal" in params:
            ax.plot([0, k * r.params["normal"][0]], [0, k * r.params["normal"][1]], 'r')
        plt.title("X-Y plane projection")
