# -*- coding: utf-8 -*-
"""
HALOS model
Heliostat class and methods
"""

import Toolbox as tb
import numpy as np

class Heliostat(object):
    def __init__(self, var_dict):
        """
        Inputs:
            var_dict:       A dictionary with all user inputs using the variables names below (without leading underscore)
        """
        # Default values
        self._helio_name = ""                    # string: [-] name of template heliostat uses
        self._id = None                          # int: [-] heliostat id number
        self._tht = 150.0                        # float: [m] Tower height -> Should come from Receiver object (need for panel canting)
        self._width = 0.0                        # float: [m] heliostat width
        self._height = 0.0                       # float: [m] heliostat height glass
        self._pylon_height = 0.0                 # float: [m] distance from ground to heliostat pivot point 
        self._is_multi_panels = True            # bool: [-] is heliostat made of multiple panels
        self._n_cant_x = 2                       # int: [-] number of canting panels in horizontal direction (x) 
        self._n_cant_y = 3                       # int: [-] number of canting panels in vertical direction (y) 
        self._x_gap = 0.                         # float: [m] distance between panels in x 
        self._y_gap = 0.                         # float: [m] distance between panels in y 
        self._cant_method = 0                    # int: [-] heliostat canting method {0: no canting, -1: on-axis at slant}
        self._err_elevation = 0                  # float: [rad] elevation pointing error 
        self._err_azimuth = 0                    # float: [rad] azimuth pointing error 
        self._err_surface_x = 0.00153            # float: [rad] surface slope error in x 
        self._err_surface_y = 0.00153            # float: [rad] surface slope error in y 
        self._err_reflect_x = 0.0002             # float: [rad] reflected beam error in x
        self._err_reflect_y = 0.0002             # float: [rad] reflected beam error in y 
        self._reflect_ratio = 0.97               # float: [-] the ratio of active reflective area to total structural area [-]
        self._reflectivity = 0.95                # float: [-] base surface reflectivity in clean state
        self._soiling = 0.95                     # float: [-] the fraction of light that is reflected after accounting for surface soiling

        self._focus_method = 0                   # int: [-] heliostat focusing type {0: flat, 1: at slant}
        self._xfocal = 1.e9                      # float: [m] focal length in the i^ direction (assume flat with no focusing)
        self._yfocal = 1.e9                      # float: [m] focal length in the j^ direction (assume flat with no focusing - TODO: update these base on focus method)

        self._is_enabled = False                 # bool: [-] is heliostat enabled (tracking)
        self._is_round = False                   # bool: [-] is heliostat round
        self._zenith = 0.0                       # float: [rad] zenith angle of heliostat
        self._azimuth = 0.0                      # float: [rad] azimuth angle of heliostat

        self._image_size_xy = np.zeros(2)           # float: [m/m] Image size on the receiver plane in {x,y}, normalized by tower height

        self._mu_MN = np.zeros((3,3))            # float: [-] Normalized moments of the Mirror Shape
        self._mu_M = np.zeros((3,3))             # float: [-] Moments of the Mirror Shape
        self._mu_S = np.zeros((3,3))             # float: [-] Moments of Sunshape
        self._mu_G = np.zeros((3,3))             # float: [-] Moments of the error distribution 
        self._mu_F = np.zeros((3,3))             # float: [-] Flux moment distrubtion - result 
        self._hcoef = np.zeros(3)                # float: [-] Hermite coefficients

        self._location = tb.Point(np.array([0,0,0]))       # numpy array (float): [m] location of heliostat (x, y, z coordinates) 
        self._aim_point = tb.Point(np.array([0,0,1]))      # numpy array (float): [-] location of aim point (x, y, z coordinates)
        self._track = tb.Vector(np.array([0,0,1]))          # numpy array (float): [-] the tracking vector for the heliostat
        self._tower_vect = tb.Vector(np.array([0,0,1]))     # numpy array (float): [-] heliostat-to-tower unit vector

        self._corners = list()                   # List of Point objects containing heliostat points {0: Upper right, 1: Upper left, 2: lower left, 3: lower right}

        # set class properties based on user input dictionary
        for key in var_dict:
            if isinstance(getattr(self,'_' + key), tb.Vector):
                # creating vector classes
                setattr(self, '_' + key, tb.Vector(var_dict[key]))
            elif isinstance(getattr(self,'_' + key), tb.Point):
                # creating point calsses
                setattr(self, '_' + key, tb.Point(var_dict[key]))
            else:
                setattr(self, '_' + key, var_dict[key])

        # is the heliostat made of multiple panels?
        if self._n_cant_x > 1. or self._n_cant_y > 1.:
            self._is_multi_panels = True

        self.updateCalculatedParameters()

        # set-up panels
        self._panels = [[Reflector() for j in range(self._n_cant_x)] for i in range(self._n_cant_y)]
        # Now define panel geometry
        self.installPanels()


    #Accessors
    def getFocalX(self): return self._xfocal
    def getFocalY(self): return self._yfocal

    def getRadialPos(self): return np.sqrt(sum(self._location.point[:2]**2))            # this assume the tower is at origin
    def getAzimuthalPos(self): return np.arctan2( self._location.x, self._location.y)   # returns radians

    def getLocation(self): return self._location
    def getTrackVector(self): return self._track
    def getAimPoint(self): return self._aim_point 
    def getTowerVector(self): return self._tower_vect
    def getSlantRange(self): return self._slant
    def getPanels(self): return self._panels

    def getMirrorShapeNormCoefObject(self): return self._mu_MN
    def getMirrorShapeCoefObject(self): return self._mu_M
    def getSunShapeCoefObject(self): return self._mu_S
    def getErrorDistCoefObject(self): return self._mu_G
    def getFluxMomentsObject(self): return self._mu_F
    def getHermiteCoefObject(self): return self._hcoef

    # Setters
    def setMirrorShapeNormCoefObject(self, mu_MN): self._mu_MN = mu_MN
    def setMirrorShapeCoefObject(self, mu_M): self._mu_M = mu_M
    def setSunShapeCoefObject(self, mu_S): self._mu_S = mu_S
    def setErrorDistCoefObject(self, mu_G): self._mu_G = mu_G
    def setFluxMomentsObject(self, mu_F): self._mu_F = mu_F
    def setHermiteCoefObject(self, hcoef): self._hcoef = hcoef

    def setSlantRange(self, L): self._slant = L
    def setTrackingAngles(self, azimuth, zenith): self._azimuth = azimuth; self._zenith = zenith

    def setLocation(self, nparray): self._location.setPoint(nparray)
    def setTrackVector(self,nparray): self._track.setVector(nparray)
    def setAimPoint(self, nparray): self._aim_point.setPoint(nparray)
    def setTowerVector(self,nparray): self._tower_vect.setVector(nparray)

    def setImageSize(self, sigx_n, sigy_n): self._image_size_xy[0] = sigx_n; self._image_size_xy[1] = sigy_n


    def updateCalculatedParameters(self):
        #calculate the collision radius and area
        self._r_collision = np.sqrt(self._height**2 + self._width**2)/2.

        self._area = (self._width * self._height * self._reflect_ratio   # width*height*structural density is the base area
                        - self._x_gap * self._height * (self._n_cant_x - 1)
                        - self._y_gap * self._width * (self._n_cant_y - 1) #subtract off gap areas
                        + (self._n_cant_y - 1) * (self._n_cant_x - 1) * self._y_gap * self._x_gap   #  don't double-count the little squares in both gaps
        )

        self._err_tot = np.sqrt( 4* (self._err_elevation**2 
                                    + self._err_azimuth**2 
                                    + self._err_surface_x**2
                                    + self._err_surface_y**2) 
                                + self._err_reflect_x**2
                                + self._err_reflect_y**2
                                ) * 1/np.sqrt(2)        # accounts for translation from 2 dimensions to 1

        self._slant = np.sqrt(self.getRadialPos()**2 + self._tht**2)
        
        return

    def installPanels(self):
        """
        This method uses the inputs to define the location and pointing vector of each
        panel on the heliostat. 

        DELSOL3 lines 6494-6520

        Note that in DELSOL3, this originally is part of the flux algorithm. The panel
        arrangement is more conveniently conceptualized as an attribute of the heliostat
        rather than as part of the flux algorithm, so it is placed here instead.
        """
        # Initialize the image plane image size for this heliostat to zero until it's calculated in the Flux methods
        self.setImageSize(0.,0.)

        if self._is_round:
            #TODO: add round heliostat (lines 410-426 in Heliostat.cpp)
            pass
        else:   #Rectangular heliostats
            dx = ( self._width - self._x_gap *( self._n_cant_x - 1.) )/ float(self._n_cant_x)   # [m] width of each canting panel
            dy = ( self._height - self._y_gap *( self._n_cant_y - 1.) )/ float(self._n_cant_y)  # [m] height of each panel

            
            # back-calculate the aim point
                    # need for Cant_method == OFFAXIS_DAY_AND_HOUR 
            #paim = tb.Point(self.getLocation + self.getSlantRange*self.getTowerVector)
            
            #initialize X and Y location of the facet
            y = - self._height*0.5 + dy*0.5

            idn = int(0)
            for j in range(self._n_cant_y):
                # initialize the X location of the facet
                x = - self._width*0.5 + dx*0.5

                for i in range(self._n_cant_x):
                    # Assign an ID
                    self._panels[j][i].setID(idn); idn += 1
                    self._panels[j][i].setType(1)    # Type=1, rectangular panel
                    self._panels[j][i].setWidth(dx)
                    self._panels[j][i].setHeight(dy)
                    # Set the position in the reflector plane. Assume the centroid is in the plane (z=0)
                    self._panels[j][i].setPosition(x, y, 0.0)

                    if self._cant_method == -1:            # {0: no canting, -1: on-axis at slant}
                        hyp = np.sqrt( self._slant**2 + x**2 + y**2)
                        self._panels[j][i].setAim(-x/hyp, -y/hyp, 2.*self._slant/hyp)
                    elif self._cant_method == 0:           # no canting
                        self._panels[j][i].setAim(0., 0., 1.)
                    else:
                        # other canting methods exist in SolarPilot Heliostat.cpp -> lines 461-539
                        print("Error: Cant Method is not supported!  Please add to installPanels method in Heliostat Class.")
                    
                    #increment the x panel position
                    x += dx + self._x_gap
                #increment the y panel position
                y += dy + self._y_gap
        return

    def updateTrackVector(self, sunVector):
        """
        Calculates the tracking vector given a solar position in "Ambient"
        and a receiver in "Receiver".

        Updates the coordinates of the heliostat corners for shadowing/blocking calculations.

        Do not update the aim point. This method uses the currently assigned aim point.

        This also updates:
        _track			| setTrackVector()		| The tracking vector for the heliostat
        _azimuth		| setTrackAngles()		| The tracking azimuth angle
        _zenith			| setTrackAngles()		| The tracking zenith angle
        _corners		| none					| The location of the heliostat corners in global coordinates (for shadowing and blocking)

        Store the new tracking vector in _track

        From Snell's law, n_hat = (s_hat + t_hat) / mag(s_hat + t_hat)
            where:
            n_hat is the normal tracking vector
            s_hat is the heliostat to sun vector
            t_hat is the helostat to receiver vector
        """

        # initialize vectors
        n_hat = tb.Vector(np.array([0.,0.,1.]))    # Tracking vector
        t_hat = tb.Vector(np.array([0.,0.,1.]))    # Heliostat to receiver
        s_hat = tb.Vector(sunVector)               # Heliostat to Sun


        if (self._is_enabled):
            t_hat.setVector(self.getAimPoint() - self.getLocation())

            ts = t_hat.vec + s_hat.vec
            ts_mag = np.sqrt(sum(ts**2))
            
            n_hat.setVector(ts/ts_mag)

            # set the tracking angles
            self.setTrackingAngles(np.arctan2(n_hat.i,n_hat.j), np.arccos(n_hat.k))
        else:
            # aimpoint vector is reflection of sun vector
            t_hat.setVector(s_hat*np.array([-1,-1,0]))

            # normal vector is zenith [0,0,1] 
            # (don't need to change due to default above if statement)

            # make tracking angles so that heliostat "faces" tower position when in stow
            self.setTrackingAngles(self.getAzimuthalPos(), 0.0)

        # Set the heliostat object tracking vector
        self.setTrackVector(n_hat)
        # Set the heliostat to tower vector
        self.setTowerVector(t_hat)

        """
        Calculate the location in global coordinates of the top two heliostat corners. Note that 
        by the azimuth convention where North is 0deg, the upper edges of the heliostat will begin on
        the southernmost edge of the heliostat.
            
        Assume that the heliostat is starting out facing upward in the z direction with the 
        upper and lower edges parallel to the global x axis (i.e. zenth=0, azimuth=0)
        """
        if not self._is_round:
            wm2 = self._width/2.
            hm2 = self._height/2.
            self._corners = list()  #Reset corner locations

            baseP = np.array([wm2,hm2,0.])
            self._corners.append(tb.Point(baseP*np.array([-1,-1,0])))       # southwest (upper right)
            self._corners.append(tb.Point(baseP*np.array([ 1,-1,0])))       # southeast (upper left)
            self._corners.append(tb.Point(baseP*np.array([ 1, 1,0])))       # northeast (lower left)
            self._corners.append(tb.Point(baseP*np.array([-1, 1,0])))       # northwest (lower right)

            for cornerPoint in self._corners:
                # Rotate first about the x axis (zenith)
                tb.rotation(self._zenith, 0, cornerPoint)
                # Now rotate about the z-axis (azimuth)
                tb.rotation(self._azimuth, 2, cornerPoint)
                # Move from heliostat coordinates to global coordinates
                cornerPoint.Add(self.getLocation())

        else:
            # no corner geometry to consider for round heliostats
            pass


        return

    def calcAndSetAimPointFluxPlane(self, aimPoint, Receiver):
        """
        Given a particular aim point in space, translate the position to an aimpoint on the actual
        flux plane of the receiver. The original aimpoint may not necessarily be on the plane of the 
        receiver, but the final aim point will be.
        """
        ## TODO: Fill out this function when needed
        pass

# --------------Reflector class methods ----------------------
class Reflector(object):
    def __init__(self):
        self._width = 0.
        self._height = 0.
        self._diameter = 0.
        self._focal_length = 0.
        self._id = -1
        self._type = 1
        
        self._locate_vector = tb.PointVect()
        self.setOrientation(0., 0., 0., 0., 0., 0.)

    # Get-Set methods
    def getID(self): return self._id
    def getWidth(self): return self._width
    def getHeight(self): return self._height
    def getDiameter(self): return self._diameter
    def getFocalLength(self): return self._focal_length
    def getType(self): return self._type
    def getOrientation(self): return self._locate_vector

    def setID(self, idnum): self._id = idnum
    def setType(self, Rtype): self._type = Rtype
    def setWidth(self, width): self._width = width
    def setHeight(self, height): self._height = height
    def setDiameter(self, diam): self._diameter = diam

    def setPosition(self, x, y, z):
        self._locate_vector.x = x
        self._locate_vector.y = y
        self._locate_vector.z = z

    def setAim(self, i, j, k):
        self._locate_vector.i = i
        self._locate_vector.j = j
        self._locate_vector.k = k


    def setOrientation(self, x, y, z, i, j, k):
        self.setPosition(x,y,z)
        self.setAim(i,j,k)
        return






if __name__ == "__main__":

    ## I think that we should have a default template accessable somewhere
    var_dict = {"helio_name": "template_1",
                "id": 1,
                "location" : np.array([300, 400, 0]),
                "width": 12.2,
                "height": 12.2,
                "is_multi_panels": False,
                "n_cant_x": 2,
                "n_cant_y": 3,
                "x_gap": 0,
                "y_gap": 0,
                "cant_method": -1,
                "focus_method": 0,
                "err_elevation": 0,
                "err_azimuth": 0,
                "err_surface_x": 0.00153,
                "err_surface_y": 0.00153,
                "err_reflect_x": 0.0002,
                "err_reflect_y": 0.0002,
                "reflect_ratio": 0.97,
                "reflectivity": 0.95,
                "soiling": 0.95,
                "tht": 200.0
                }

    h = Heliostat(var_dict)

    pass
    

