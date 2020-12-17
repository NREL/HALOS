# -*- coding: utf-8 -*-
"""
OHLASTS model
Mirror and slope error model
"""
import scipy.stats
import numpy


# import math

class Mirror(object):
    def __init__(self,locs,shape,error,error_type):
        self.mirror_locs = locs   #3-column array(float): x, y, z coordinates for each mirror in field
        self.shape = shape      #2-element tuple(float): length, width
        self.error = error      #data structure is a tuple for gaussian, not sure about pillbox yet
        self.error_type = error_type #string
        
    def getType(self):
        return self.error_type

    def getAtmAtt(self, helio_idx, measure_loc, atm_method=0, user_defined=None):
        """
        Parameters
        ----------
        helio_idx: int
            index of heliostat in field to use as input
        
        self.mirror_locs : numpy.ndarray (float)
            x, y, z coordinates of mirror (heliostat) location [m]
            
        measure_loc : numpy.ndarray (float)
            x, y, z coordinates of aimpoint [m]
            
        atm_method : INT
            Dictates the atmospheric attenuation method to be used for calculations
            0 := DELSOL3 clear day (Barstow 25km visibility)
            1 := DELSOL3 hazy day (Barstow 5km visibility)
            2 := user defined coefficients (polynomials) -> must have user_defined parameter specified
            
        user_defined : numpy.ndarray (float)
            Array of polynomial coefficients as a function of Slant range [km] starting with the zeroth degree
        
        Returns
        -------
        att : float
            Atmospheric attenuation losses for the reflected radiation between the heliostat and receiver [0,1]
    
        """
        atm_models = {
            0: numpy.array([0.006789, 0.1046, -0.017, 0.002845]),
            1: numpy.array([0.01293, 0.2748, -0.03394]),
            2: user_defined
        }

        atm_coefs = atm_models.get(atm_method)

        # Calculating slant range from the heliostat to the receiver
        slant_range = numpy.sqrt(numpy.sum(measure_loc - self.mirror_locs[helio_idx]) ** 2)
        slant_range *= 0.001  # convert from meters to kilometers

        att = 0.
        for i, c in enumerate(atm_coefs):
            att += c * (slant_range ** i)

        return 1. - att


class GaussianMirror(Mirror):
    def __init__(self, loc, shape, error):
        Mirror.__init__(self, loc, shape, error, "Gaussian")

    def getFlux(self,loc, aim_loc):
        pass


class PillboxMirror(Mirror):
    def __init__(self, loc, shape, error):
        Mirror.__init__(self, loc, shape, error, "Gaussian")


class SinglePointGaussianMirror(Mirror):
    """
    This mirror represents a single-point image, with azimuth and elevation
    aiming errors that are Gaussian in nature.
    """

    def __init__(self, loc, error):
        Mirror.__init__(self, loc, (0., 0.), error, "SinglePointGaussian")

    def getFlux(self, helio_idx, measure_loc, aim_loc):
        """
        Inputs
        loc -- x, y, z coordinates of measurement point in a numpy array
        aim_loc -- x, y, z coordinates of aimpoint 
        
        Returns
        flux: estimate of intensity at point (pdf of gaussian)
        """
        # a = distance from mirror to measurement
        a = numpy.sqrt(numpy.sum(measure_loc - self.mirror_locs[helio_idx]) ** 2)
        # b = distance from mirror to aimpoint
        b = numpy.sqrt(numpy.sum(aim_loc - self.mirror_locs[helio_idx]) ** 2)
        # c = distance from aimpoint to measurement
        c = numpy.sqrt(numpy.sum(measure_loc - aim_loc) ** 2)
        # Get angular error using law of cosines, convert to mrad
        angle = 1000 * numpy.arccos((a * a + b * b - c * c) / (2 * a * b))
        # return gaussian pdf at z = (angle / error)
        return scipy.stats.norm.pdf(angle / self.error)
