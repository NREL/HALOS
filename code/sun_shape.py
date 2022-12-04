# -*- coding: utf-8 -*-
"""
HALOS
sun_shape module
This module contains source code for building sun shape objects and their 
related error models.
"""

import scipy.stats
import numpy

class Sunshape(object):
    def __init__(self, model_type, model_params):
        self.type = model_type
        self.model_params = model_params


class PillBoxSunshape(Sunshape):
    def __init__(self, model_params):
        Sunshape.__init__(self, "Pillbox", model_params)

    def getIntensity(self, error):
        """ 
        Returns normalized intensity for a location, given angular error
        as input.  For the pillbox error type, is a simple binary distribution 
        in which flux is assumed to be equally distributed within a specified disc. 
        The region outside the disc is modeled with zero flux intensity.
        
        Inputs
        error -- angular error from center, [mrad] (float)
        
        Returns
        normalized intensity -- flux at the location [-] (float)
        """
        if error <= self.model_params["error"]:
            return 1.0
        else:
            return 0.0


class GaussianSunshape(Sunshape):
    def __init__(self, model_params):
        Sunshape.__init__(self, "Gaussian", model_params)

    def getIntensity(self, error):
        """
        Returns normalized intensity for a location, given angular error
        as input.  For the gaussian error type, this is equivalent to the pdf 
        of the standard normal distribution at z=angle_error/model_error. 
        
        Inputs
        error -- angular error from center, [mrad] (float)
        
        Returns
        normalized intensity -- flux at the location,  [-] (float)
        """
        return scipy.stats.norm.pdf(error / self.model_params["error"])


class LimbDarkenedSunshape(Sunshape):
    def __init__(self, model_params):
        Sunshape.__init__(self, "Limb-Darkened", model_params)

    def getIntensity(self, error):
        """
        Returns normalized intensity for a location, given angular error
        as input.  For the Limb-Darkened error type, this model specifies sun intensity 
        as a function of angular distance from the centroid of the sun disc.
        
        Inputs
        error -- angular error from center, [mrad] (float)
        
        Returns
        normalized intensity -- flux at the location,  [-] (float)
        """
        if error <= 5.4923:
            return 1 - 0.5138 * (error / 4.65) ** 4
        else:
            return 0.0


class SinglePointSun(Sunshape):
    def __init__(self, model_params):
        Sunshape.__init__(self, "SinglePoint", model_params)

    def getIntensity(self, error):
        return 1.0


class BuieCSRSunshape(Sunshape):
    def __init__(self, model_params):
        Sunshape.__init__(self, "BuieCSR", model_params)

    def getIntensity(self, error):
        """
        Returns normalized intensity for a location, given angular error
        as input.  For the Buie CSR error type, this model specifies sun intensity 
        as a function of angular distance from the centroid of the sun disc.
        
        Creates the Buie (2003) sun shape based on CSR
        [1] Buie, D., Dey, C., & Bosi, S. (2003). The effective size of the solar cone for solar concentrating systems.
        Solar energy, 74(2003), 417–427.
        [2] Buie, D., Monger, A., & Dey, C. (2003). Sunshape distributions for terrestrial solar simulations. Solar
        Energy, 74(March 2003), 113–122.

        Inputs
        error -- angular error from center, [mrad] (float)
        model_params:
            chi -- Circumsolar ratio, [-] (float)

            chi_corr -- corrected circumsolar ratio, [-] (float)
                        is calculated if not in model_params
        
        Returns
        normalized intensity -- flux at the location,  [-] (float)
        """
        if not 'chi_corr' in self.model_params:
            self.updateBuieCalcParameters()

        if error > 4.65:
            return numpy.exp(self.model_params['buie_kappa']) * (error ** self.model_params['buie_gamma'])
        else:
            return numpy.cos(0.326 * error) / numpy.cos(0.308 * error)

    def updateBuieCalcParameters(self):
        """
        Updates Buie Sunshape calculated parameters.

        CSR correction is from Tonatiuh:
            https://github.com/pablo-benito/tonatiuh    Main repository
            https://github.com/pablo-benito/tonatiuh/blob/master/TonatiuhProject/plugins/SunshapeBuie/src/SunshapeBuie.cpp  starting at line 175
        
        Inputs
        model_params:
            chi -- Circumsolar ratio, [-] (float)
        
        Returns
        model_params:
            chi_corr -- Corrected Circumsolar ratio, [-] (float)
            buie_kappa -- Sunshape variable, [-] (float)
            buie_gamma -- Sunshape variable, [-] (float)
        """
        csr = self.model_params['chi']
        if csr > 0.145:
            self.model_params['chi_corr'] = -0.04419909985804843 + csr * (1.401323894233574 + csr * (
                        -0.3639746714505299 + csr * (-0.9579768560161194 + 1.1550475450828657 * csr)))
        elif csr > 0.035:
            self.model_params['chi_corr'] = 0.022652077593662934 + csr * (
                        0.5252380349996234 + (2.5484334534423887 - 0.8763755326550412 * csr) * csr)
        else:
            self.model_params['chi_corr'] = 0.004733749294807862 + csr * (4.716738065192151 + csr * (
                        -463.506669149804 + csr * (
                            24745.88727411664 + csr * (-606122.7511711778 + 5521693.445014727 * csr))))

        self.model_params['buie_kappa'] = 0.9 * numpy.log(13.5 * self.model_params['chi_corr']) * (
                    self.model_params['chi_corr'] ** (-0.3))
        self.model_params['buie_gamma'] = 2.2 * numpy.log(0.52 * self.model_params['chi_corr']) * (
                    self.model_params['chi_corr'] ** 0.43) - 0.1
        return None


if __name__ == "__main__":
    test = BuieCSRSunshape(model_params={"chi": 0.1})
    I = test.getIntensity(5.7)
    print(I)
