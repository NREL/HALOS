# -*- coding: utf-8 -*-
"""

OHLASTS miscellaneous scripts

This module calculates a Hermite expansion 

"""

import Toolbox as tb

import math
import numpy as np
import scipy
import scipy.stats



class Flux(object):
    def __init__(self):
        """
            This class was taken form SolarPILOT (source code in Flux.cpp)
        """
        self._n_order = 6                                                   # int: [-] Order of the Hermite expansion
        self._n_terms = self._n_order + 1                                   # int: [-] Number of terms (n_order+1), to include constant term

        #self._hermitePoly = np.zeros((1, 3))                               # Float: [?] Contains a vector of values of the evaluated hermite polynomial coefficients (I dont think this is needed)
        self._fact_odds = np.zeros((self._n_terms*2))                       # Float: [-] Contains the factorial values for the odd terms in the hermite expansion
        self._fact_d = np.zeros((self._n_terms*2))                          # float: [-] Contains factorial values
        self._binomials = np.zeros((self._n_terms,self._n_terms))           # float: [-] Contains binomial coefficients
        self._binomials_hxn = np.zeros((self._n_terms,self._n_terms))       # float: [-] Contains binomial coefficients for the HxN array

        self._mu_SN = np.zeros((self._n_terms,self._n_terms))               # float: [-] Normalized moment of sunshape
        self._mu_GN = np.zeros((self._n_terms,self._n_terms,4))             # float: [-] Normalized moments of the error distribution

        #Calculate the array of factorial terms in the Hermite expansion
        self.factOdds()
        self._fact_d = np.array([math.factorial(i) for i in range(self._n_terms*2)])
        self.Binomials()
        self.Binomials_hxn()

        #Set up coefficient weighting arrays for Hermite integral
        self._ci = np.array([.196584, .115194, .000344, .019527])
        self._ag = np.array([.02715246, .06225352, .09515851, .12462897,
					         .14959599, .16915652, .18260341, .18945061,
					         .14959599, .16915652, .18260341, .18945061,
					         .02715246, .06225352, .09515851, .12462897])
        self._xg = np.array([.98940093, .94457502, .86563120, .75540441,
					         .61787624, .45801678, .28160355, .09501251,
					        -.61787624,-.45801678,-.28160355,-.09501251,
					        -.98940093,-.94457502,-.86563120,-.75540441])

        #Create jmin and jmax arrays
        self._jmin = np.zeros((self._n_terms),dtype=int)
        self._jmax = np.zeros((self._n_terms),dtype=int)
        for i in range(self._n_terms):
            self._jmin[i] = i%2+1
            self._jmax[i] = self._n_terms - i
        return

    def factOdds(self):
        #Calculate the factorial values for the odd terms in the hermite expansion
        #factorial of odds.. 1*3*5... etc
        self._fact_odds[1] = 1.       #first value is 1
        for i in range(3,self._n_terms*2,2):
            self._fact_odds[i] = self._fact_odds[i-2]*np.double(i)
        return

    def Binomials(self):
        #Calculate the binomial coefficients
        for i in range(1,self._n_terms+1):
            for j in range(1, self._n_terms+1):
                self._binomials[i-1,j-1] = self._fact_d[i-1]/self._fact_d[j-1]/self._fact_d[i-j]
        return

    def Binomials_hxn(self):
        """
            Calculation of "binomial coefficients" (different from previously calculated array) assigned
	        to the array HXN, DELSOL lines 741-754.
	    	
            I don't know where this comes from mathematically.
	    """
        self._binomials_hxn[0,0] = 1.
        self._binomials_hxn[1,1] = 1.

        for i in range(3,self._n_terms+1):
            fi = np.float(i-2)
            self._binomials_hxn[i-1,0] = -fi*self._binomials_hxn[i-3,0]

            for j in range(2,self._n_terms+1):
                self._binomials_hxn[i-1,j-1] = self._binomials_hxn[i-2][j-2] - fi*self._binomials_hxn[i-3][j-1]
        return

    def JMN(self, i):
        """
      	//Frequently used bounds calculation
    	//Treat as an array, should return [1,2,1,2,1,2]...
	    //return i%2+1;
        """
        return self._jmin[i]
    
    def JMX(self, i):
        """
    	//Frequently used bounds calculation
	    //Treat JMX like an array for conventions to match DELSOL. i.e. where in DELSOL we would call JMX(1) and 
    	//expect to get back 7, here we call JMX(0) and get back 7. ???? this should return 6? - WTH
	    //return _n_terms - i;
        """
        return self._jmax[i]
    
    def IMN(self, i):
        """
        //Frequently used bounds calculation
	    //return JMN(i);
        """
        return self._jmin[i]

    def hermitePoly(self, x):
        # SolarPILOT says this is not used - Flux.cpp (lines 288-320)
        pass

    def initHermiteCoefs(self, SunShape):
        # Fills out the constant coefficients that don't change during the simulation

        # Called in SolarField.cpp
        
        # Sun shape
        self.hermiteSunCoefs(SunShape)

        # Error distribution coefficients
        self.hermiteErrDistCoefs()
        return

    def hermiteSunCoefs(self, SunShape): 
        """
        Directly from SolarPILOT Flux.cpp (lines 337 - 552)
        ###############################################################################################
        -------WHEN TO CALL------
        Call this subroutine once to determine the sunshape coefficients. The coefficients do not 
        currently depend on solar position, time of day, weather, etc. These coefficients will be used
        in the subroutine "imagePlaneIntercept" to determine the moments of sunshape.

        ---INFORMATION REQUIRED--
        * Requires the sunshape model type -> Ambient
        * the _fact_odds array must have been calculated
        * for a user-specified sunshape, the sunshape array must be filled out -> Ambient.getUserSun()

        ---------OUTPUT----------
        Fills out the "mSun" (NxN) array
        ###############################################################################################
        
        This function calculates the coefficients of sunshape distribution. The moments are used
        in the analytical Hermite polynomial formulation. Each sunshape moment
        is of form [Dellin, 1979]:
        mu_S_(i,j) = (1/(i+j+2) - 0.5138/(i+j+6))/(1/2 - .5138/6) * (i!!*j!!/((i+j)/2)!)*2^((i+j)/2)r_0^(i+j)
        where r_0 = (4.65e-3 mrad) * (slant range)
            
        The sun-shape can take one of several forms. The limb-darkening expres-
        sion is given by [Dellin, 1979]:
        S(r) = 1 - .5138 * (r/r_0)^4
        where 'r' is the radius about the aim poing in the image plane,
        and r_0 is given above. The r_0 value is not actually applied in this
        algorithm.
            
        Other options for sunshape include point-source (model=0), square-wave
        sun (model=2), or user-defined sunshape (model=3).
            
        For the user-defined sunshape, the user must provide a 2-D array of
        sun intensity and corresponding angular deviation from the center of the
        solar disc.
        suntype = [[angle0, intens0], [angle1, intens1] ...]
            
        Otherwise, sunshape can be declared using:
        suntype = <option number>
            
        Default sunshape is limb-darkened (option 1) 
        """
        # get the sun type (current options: "Pillbox","Gaussian","Limb-Darkened","SinglePoint","BuieCSR")
        # Select a suntype case using if statement
        if SunShape.type == "Pillbox":
            #---Square-wave sunshape--- see DELSOL3 lines 6416-6425
            for i in range(1,self._n_terms+1,2):                # iterate 'i' from 1 to N_order by 2
                factdum1 = 1.
                if (i>1):
                    factdum1 = self._fact_odds[i-2]             # Hold on to the factorial i-1 value to avoid constant recalculation. If i==1, set to 1.
                
                for j in range(1,self._n_terms+1,2):
                    factdum2 = 1.
                    if (j>1):                                   # Hold on to the factorial j-1 value to avoid recalc. if j==1, set to 1.
                        factdum2 = self._fact_odds[j-2]         # Hold on to the i+j value to avoid multiple recalculations
                    ij = int(i+j)
                    self._mu_SN[i-1,j-1] = (2.*factdum1*factdum2/math.factorial(ij/2-1)/(2**((ij-2)/2))/np.double(ij))*(SunShape.model_params["error"]/1000.)**(ij-2) #4.645e-3**(ij-2)

        elif SunShape.type == "Limb-Darkened":
            #---Limb-darkened sunshape--- see DELSOL3 lines 6399-6414
            for i in range(1,self._n_terms+1,2):                # iterate 'i' from 1 to N_order by 2
                jmax = self._n_terms-i+1                        # Set the upper bound on the iteration limit for j
                factdum1 = 1.
                if (i>1):
                    factdum1 = self._fact_odds[i-2]             # Hold on to the factorial i-1 value to avoid constant recalculation. If i==1, set to 1.
                
                for j in range(1,jmax+1,2):
                    factdum2 = 1.
                    if (j>1):                                   # Hold on to the factorial j-1 value to avoid recalc. if j==1, set to 1.
                        factdum2 = self._fact_odds[j-2]         # Hold on to the i+j value to avoid multiple recalculations
                    ij = int(i+j)
                    # Calculate the moment for this i-j combination. Algorithm taken from DELSOL3, lines 6399-6414
                    dfact = np.double(math.factorial((ij)/2-1))
                    self._mu_SN[i-1,j-1] = ((1./np.double(ij) - 0.5138/np.double(ij+4))/(0.5 - 0.5138/6.) * factdum1 * factdum2 / dfact / 2**((ij-2)/2)) * 4.65e-3**(ij-2)
            
        elif SunShape.type == "SinglePoint":
            #---Point of sun unit intensity---
            self._mu_SN[0,0] = 1.0

            # removed for loops because _mu_SN is initialized as an zero array
            """
            for i in range(1,self._n_terms+1,2):
                jmax = self._n_terms-i+1
                for j in range(1,jmax+1,2):
                    self._mu_SN[i-1,j-1] = 0.
            """
        elif SunShape.type == "Gaussian" or SunShape.type == "BuieCSR" or SunShape.type == "User_Sun":
            if SunShape.type == "Gaussian" or SunShape.type == "BuieCSR":
                npt = 50
                temp_sun = np.zeros((npt,2))
                theta = np.array([x*25./npt for x in range(npt)])  # equally spaced to 25 mrad
                for i in range(npt):
                    temp_sun[i, 0] = theta[i]                      # mrad
                    temp_sun[i, 1] = SunShape.getIntensity(theta[i])       # Calls sunshape to get Intensity at theta

                user_sun = temp_sun                                # Assign 
            else:
                user_sun = SunShape.user_sun                       # Assign
            
            azmin = np.zeros((12))
            nn = int(len(user_sun) - 1)
            for n in range(1, nn+1):
                # The disc angle and corresponding intensity
                disc_angle = user_sun[n-1,0]/1000.      #converted mrad to radians
                intens = user_sun[n-1,1]
                # The next disc angle and intensity pair in the array
                disc_angle_next = user_sun[n,0]/1000.
                intens_next = user_sun[n,1]/1000.

                # fractional step size
                rel_step = 1./(disc_angle_next - disc_angle)
                # Relative to the step size, how far away are we from the centroid?
                    # Is this making the assumptions that the steps are equal?
                r_steps = disc_angle*rel_step

                for m in range(1, 8, 2):
                    temp1 = (disc_angle_next**(m+1) + disc_angle**(m+1))/np.double(m+1)
                    temp2 = (disc_angle_next**(m+2) + disc_angle**(m+2))/np.double(m+2)
                    azmin[m-1] += intens*(temp1*(1+r_steps) - temp2*rel_step) + intens_next*(-temp1*r_steps + temp2*rel_step)
            
            xnorm = 1.       # Initialize the normalizing variable.. it will be reset once the new value is calculated below
            #Also initialize an array that's needed for this calculation - see DELSOL3 lines 6238-6244
            RSPA = np.zeros((self._n_terms,self._n_terms))
            RSPA[0] = [2.,0.,1.,0.,.75,0.,.625]
            RSPA[2] = [1.,0.,.25,0.,.125,0.,0.]
            RSPA[4] = [.75,0.,.125,0.,0.,0.,0.]
            RSPA[6] = [.625,0.,0.,0.,0.,0., 0.]

            for i in range(1, self._n_terms+1, 2):
                jmax = self._n_terms-i+1
                for j in range(1, jmax+1, 2):
                    ij = int(i+j)
                    self._mu_SN[i-1,j-1] = azmin[ij-2]*RSPA[i-1,j-1]/xnorm*np.pi
                    xnorm = self._mu_SN[0,0]

            self._mu_SN[0,0] = 1.
        else:
            print("ERROR: Unsupported sunshape model")

        return

    def hermiteErrDistCoefs(self):
        """
        ###############################################################################################
        -------WHEN TO CALL------
        Call once. Coefficients will be used to calculate the moments of the error distribution in the
        subroutine "imagePlaneIntercept" below.

        ---INFORMATION REQUIRED--
        No dependencies.

        ---------OUTPUT----------
        Fills in the "errDM" (N x N x 4) array for the coefficients of the error distribution moments.
        ###############################################################################################

        This method calculates the moments of error distribution "G". 

        The array is N_terms x N_terms x 4 (3D array)

        From Dellin (1979), pp14:
            G represents the normalized probability distribution that the reflected
            vector t_hat will be displaced from its nominal value by and amount dt in the image
            plane (i_hat_t, j_hat_t) due to errors in the system. These displacements result 
            from the cumulative effect of many individual error sources and the magnitude
            depends on the detailed design of the system.
        
        DELSOL3 lines 6461-6483
        """

        """
        # REMOVED: Not needed
        temp1 = 1.
        # I'm not sure the purpose of this for loop -> setting temp1 to _fact_odds[7] = 105
        for i in range(1,self._n_terms+1,2):
            if (i>1):
                temp1 = self._fact_odds[i]
        """

        # Calculate each moment
        for i in range(1,self._n_terms+1):
            jmax = self.JMX(i-1)
            jmin = self.JMN(i-1)
            for j in range(jmin,jmax+1,2):
                ii = int((i-1)/2+1)
                for k in range(1,ii+1):
                    temp1 = 1.              # but then resets it to 1 here???
                    if (i+j > 2*k):
                        temp1 = self._fact_odds[i+j-2*k-1]
                    self._mu_GN[i-1,j-1,k-1] = temp1 * self._fact_d[i-1]/self._fact_d[i-2*k+1]*self._fact_d[k-1]

        return

    def hermiteMirrorCoefs(self, Heliostat, tht):
        """
        ###############################################################################################
        -------WHEN TO CALL------
        This method should be called once for each heliostat template (NOT for each heliostat!).

        ---INFORMATION REQUIRED--
        Heliostat geometry, tower height (to normalize information). The heliostat templates should be
        previously initialized, then when flux initialization occurs, set the mirror coefficients

        ---------OUTPUT----------
        Fills out the mirror moment coefficient array "errMM", which is normalized by the tower height
        and applies to each unique heliostat geometry.
        ###############################################################################################

        This method calculates the moments of mirror shape, 'M'. Moments are based on heliostat
        dimensions that are normalized by the tower height.

        This method references an established heliostat geometry and requires a receiver 
        object with an associated tower height.
        -> Heliostat
        -> Receiver

        From Dellin (1979), pp 17:
            M represents the flux produced by a perfect heliostat reflecting a point sun.
            For a flat heliostat M is given by the geometrical projection of all point
            on the mirror that are neither shaded nor blocked. WLV showed that the moments
            of a flat rectangular heliostat (including the effects of shading and blocking)
            can be evaluated analytically.

        DELSOL3 lines 6494-6525

        Round heliostats: 
        DELSOL3 lines 6529-6540
        """
        # binding Heliostat to H
        H = Heliostat   # Same ID (changes to H impact Heliostat object)

        # Assign some values from the heliostat instance
        wm = H._width
        hm = H._height
        ncantx = H._n_cant_x
        ncanty = H._n_cant_y

        # Calculate the effective mirror width to use for 
        if H._is_multi_panels:
            #Use the smaller cant panel dimensions for image calculations
            wm2s = 0.
            hm2s = 0.
            #Use average panel width/height
            ff = 1./(np.double(ncantx*ncanty)*2.*tht)
            for r in range(ncanty):
                for c in range(ncantx):
                    wm2s += H._panels[r][c]._width*ff
                    hm2s += H._panels[r][c]._height*ff
                    pass
        else:
            # Use the large structure width (or no canting is specified) for image calculations
            wm2s = wm/(tht*2.)
            hm2s = hm/(tht*2.)

        errMM = np.zeros((self._n_terms,self._n_terms))

        # Calculate the moments depending on whether the heliostats are circular or rectangular
        if H._is_round:
            #----Round heliostats----
            pass
        else:
            #----Rectangular heliostats ----
            wm2sk = wm2s                # is pow(wm2s,k)
            wm2s2 = wm2s*wm2s
            for k in range(1,self._n_terms+1, 2):
                wm2sk *= wm2s2

                hm2sl = hm2s            # is pos(hm2s,1)
                hm2s2 = hm2s*hm2s
                for l in range(1,self._n_terms+1, 2):
                    kl = k*l
                    # Calculate the moments
                    hm2sl*= hm2s2
                    errMM[k-1,l-1] = 4./np.double(kl)*wm2sk*hm2sl

        # Updating heliostat object
        Heliostat.setMirrorShapeNormCoefObject(errMM)   # sets _mu_MN in heliostat object    # this also updates object H

        #Potentially add other geometries here. No other geometries are defined in DELSOL3.
        return

    def imagePlaneIntercept(self, Heliostat, Receiver, SunVector):
        """
        ###############################################################################################
        -------WHEN TO CALL------
        Call this method to evaluate the sunshape, heliostat error distribution, and heliostat mirror
        shape distribution for EACH HELIOSTAT IN THE FIELD. Output will change when any input except
        tower height changes. 

        *This is the main optical intercept algorithm that should be called during simulation. This 
        calls the other Hermite flux characterization algorithms.*

        ---INFORMATION REQUIRED--
        * Heliostat geometry
        * Tower height
        * Receiver geometry
        * Sun position
        * Evaluated Hermite flux coefficients (mirror geometry, error distribution, sunshape)
        
        ---------OUTPUT----------
        Calculates:
            * moments of sunshape "_mu_S"
            * moments of mirror shape "_mu_M"
            * moments of error distribution "_mu_G"
            * convolved moments "_mu_F"
            --> All of these are normalized by the tower height
        Calls the subroutine chain:
            hermiteIntEval() -> hermiteIntegralSetup() -> hermiteIntegral()
        ...and returns the result (double) of the chain.
        ###############################################################################################	

        This subroutine calculates the moments of the sunshape, heliostat error
        distribution, and heliostat shape on the image plane and then convolves
        these to find the moments of the heliostat image. 

        This method is derived from the DELSOL3 subroutine "MOMENT", beginning on 
        line 8249.
        """ 
        H = Heliostat       #local objects
        R = Receiver

        #First check to make sure the heliostat can see the aperture plane. If not, don't bother.
        # TODO: create a calculate or get normal vector function for receiver Receiver.cpp line 622 (depending on type)
        # TODO: this needs to be cleaned up 
        rnorm = r.params['normal']
        # how is TowerVector set? -> solarField.cpp
        """
        view_ok = False
        if (np.dot(rnorm,H.getTowerVector()) < 0.):    # We actually want the opposite of the tower vector, so take negative dot products to be ok.
            view_ok = True
        if not view_ok:
            return 0.
        """

        h_rad = H.getRadialPos()        #[m] Get the heliostat radius from the tower

        #--------Calculate the moments of the sun-----------
        tht = R.tow_height

        # vector to store calculated sun moments
        mu_S = H.getSunShapeCoefObject()
        mu_S = np.zeros((self._n_terms,self._n_terms)) # reshape and fill

        # slant range could be calculated outside of this function
        
        if R.shape_type == 'simple_ext_cylinder' or R.shape_type == 'ext_cylinder':
            pass
        else:
            rec_opt_ht = R.getOpticalHeight(H)
            slant = np.sqrt(h_rad**2 + rec_opt_ht**2)
        
        srange = np.zeros((self._n_terms))
        srange[0] = 1.
        for i in range(1,self._n_terms):
            srange[i] = srange[i-1]*slant

        for i in range(1,self._n_terms+1,2):
            jmax = self.JMX(i-1)
            #jmin = int(2-i+2*(i/2))    # original code -> [1,2,1,2,1,2,1]
            jmin = self.JMN(i-1)        # does the same thing

            for j in range(jmin,jmax+1,2):
                # RUS -> hmSun
                mu_S[i-1,j-1] = srange[i+j-2] * self._mu_SN[i-1,j-1]

        H.setSunShapeCoefObject(mu_S)       # update heliostat object

        # ---------Calculate moments of the error distribution -------------
	    # DELSOL 8304-8339
        """
        The error distribution moments are evaluated using the formulation from Dellin (1979)
        equations 3-6 and 3-7 on page 16. The angles for this equation are defined according
        to the figures on page 38-39 of Dellin, and also on page 26-27 of the DELSOL manual.

        Equation 3-6 contains terms for sigma_x^2, sigma_y^2, and rho, each of which are evaluated
        by their formulations described in Eq. 3-7 of Dellin. These three relationships contain 
        several phrases that are identical, so to save computational effort, they are assigned
        unique variable names.

        The expansion coefficients are calculated in the "B(:,:)" array.

        Description of variables and angles from DELSOL:
        ---Nomenclature system	---
        -Char 1	
        C	Cosine
        S	Sine
        -Char 2	
        S	Sun vector
        T	Helio-to-receiver vector
        N	Heliostat tracking vector
        -Char 3	
        P	with reference to zenith (polar)
        A	with reference to azimuth
        Suffix	
        ANG	Array
        Variable	Primary description---
        AZMANG		Heliostat aiming azimuth angle
        CBNS		CSA*CNA+SNA*SSA
        CBNT		j_hat component of the normal tracking vector
        CDEC		Cosine of the declination angle
        CLAT		cosine of latitude
        CNA			check CTA
        CNP			Cosine of the tracking normal zenith angle
        CSA			Equals CSAANG(I,J)
        CSAANG		Cosine of solar azimuth angle
        CSP			Equals CSPANG(I,J)
        CSPANG		Cosine of solar zenith angle
        CT			Cosine of hour angle
        CTA			Equals CTAANG(L)
        CTAANG		cosine of heliostat location azimuth angle
        CTP			Equals CTPANG(K)
        CTPANG		cosine of TPANG
        ELVANG		Heliostat aiming elevation (?) angle
        HCOS		Heliostat cosine loss
        HINSOL	
        I			Iterator, 1 -> Number of azimuth angles in table
        J			Iterator, 1-> Number of zenith angles in table
        K			Iterator, 1-> Number of azimuthal zones
        L			Iterator, 1-> Number of radial zones
        SBNS		SNA*CSA-CNA*SSA
        SBNT		i_hat component of the normal tracking vector
        SDEC		Sine of the declination angle
        SLAT		sine of latitude
        SNA			Sine of the tracking azimuth angle
        SNP			Sine of the tracking zenith angle
        SSA			Equals SSAANG(I,J)
        SSAANG		Sine of solar azimuth angle
        SSP			Equals SSPANG(I,J)
        SSPANG		Sine of solar zenith angle
        STA			Equals STAANG(L)
        STAANG		sine of heliostat location azimuth angle
        STP			Equals STPANG(K)
        STPANG		sin of TPANG
        T			Hour angle
        TPANG		Heliostat to receiver zenith angle
        UAZ			Array of user-defined solar azimuth angles
        UEL			Array of user-defined solar ZENITH angles
        """

        # calculate relevant angles
        # SolarPILOT Flux.cpp starting at line 899

        z_hat = tb.Vector([0.,0.,1.])        # zenith unit vector
        n_hat = H.getTrackVector()          # tracking vector
        t_hat = H.getTowerVector()          # heliostat to tower vector

        s_hat = SunVector                   # Sun poistion vector

        cos_s_zen = s_hat.k
        sin_s_zen = np.sqrt(1. - s_hat.k**2)
        if sin_s_zen == 0.: sin_s_zen = 1.e-6
        cos_s_az = s_hat.j/sin_s_zen
        sin_s_az = s_hat.i/sin_s_zen

        #---------------------------------------------------

        theta_n_zen = np.arccos(n_hat.k)
        #cos_n_zen = n_hat.k
        sin_n_zen = np.sin(theta_n_zen)
        if sin_n_zen == 0.: sin_n_zen = 1.e-6

        #theta_n_az = np.arctan2(n_hat.i,n_hat.j)
        #sin_n_az = np.sin(theta_n_az)
        #cos_n_az = np.cos(theta_n_az)

        theta_t_zen = np.arccos(np.dot(t_hat.vec, z_hat.vec))
        cos_t_zen = np.dot(t_hat.vec, z_hat.vec)            # cos of zenith of helio-tower vector
        sin_t_zen = np.sin(theta_t_zen)
        
        if sin_t_zen == 0.: sin_t_zen =  1.e-6 
        theta_t_az = np.arctan2(t_hat.i,t_hat.j)            # azimuth angle of the heliostat-to-receiver vector
        sin_t_az = np.sin(theta_t_az)
        cos_t_az = np.cos(theta_t_az)

        #---------------------------------------------------	
	    # Calculate the heliostat cosine loss
        eta_consine = (np.sqrt(2.)/2.)*np.sqrt(1.+cos_s_zen*cos_t_zen+sin_s_zen*sin_t_zen*(cos_t_az*cos_s_az + sin_t_az*sin_s_az))

        # get error terms - See Kistler pp. 184 for definition
        err_angular = np.zeros(2)
        err_surface = np.zeros(2)
        err_reflected = np.zeros(2)

        err_angular[0] = H._err_azimuth
        err_angular[1] = H._err_elevation

        err_surface[0] = H._err_surface_x
        err_surface[1] = H._err_surface_y

        err_reflected[0] = H._err_reflect_x
        err_reflected[1] = H._err_reflect_y

        # Depending on the canting method, calculate the A[], B[] arrays differently.
        A11 = A12 = A21 = A22 = B11 = B12 = B21 = B22 = 0.0
        cant_method = H._cant_method        # {0=none, -1=on-axis at slant, 1=on-axis at user def., 3=off-axis at hour-day}

        # reused terms:
        # SAVE=SIGAZ2*SNP**2+SIGSX2 | 8304
        term1 = err_angular[0] * sin_n_zen
        term1 *= term1
        term1 += err_surface[0]*err_surface[0]

        # SAVE2=SIGEL2+SIGSY2
        term2 = err_angular[1]*err_angular[1] + err_surface[1] * err_surface[1]  # second reused term

        if cant_method == 1:
            #TODO: add onaxis userdefined method
            pass
        elif cant_method == 3:
            #TODO: add offaxis day and hour canting
            pass
        elif cant_method == 4:
            #TODO: add user-vector cant method
            # Currently not supported in SolarPILOT
            pass

        # if the focal distance of the heliostat is not equal to the slant range, do additional 
	    # calculations here
        
        #DELTA1=SIGTX2+SAVE*B(1,1)**2+SAVE2*B(1,2)**2
        delta_1 = err_reflected[0]*err_reflected[0] + term1 * B11 * B11 + term2 * B12 * B12
        #DELTA2=SIGTY2+SAVE*B(2,1)**2+SAVE2*B(2,2)**2
        delta_2 = err_reflected[1]*err_reflected[1] + term1 * B21 * B21 + term2 * B22 * B22
	    #DELTAA=SAVE*B(1,1)*B(2,1)+SAVE2*B(2,2)*B(1,2)
        delta_a = term1 * B11 * B21 + term2 * B22 * B12
	    #SIGX=SRANGE(2)*SQRT(DELTA1)
        sigma_x = srange[1] * np.sqrt(delta_1)      #Standard deviation of the image error in the X-direction
	    #SIGY=SRANGE(2)*SQRT(DELTA2)
        sigma_y = srange[1] * np.sqrt(delta_2)      #Standard deviation of the image error in the Y-direction

        # The argument to part of the RHO equation (sigma_a) determines the remaining calculations. 
	    # If the value of sigma_a^2 is zero, handle with a separate method.
        mu_G = H.getErrorDistCoefObject()           # this is unnecessary
        mu_G = np.zeros((self._n_terms, self._n_terms))   # reshape and fill

        if delta_a**2 < 1.e-20:
            # The argument delta_a^2 is very close to zero, handle separately
            fact1 = fact2 = 0.0
            # Calculate the factorial terms
            factarr = np.zeros((self._n_terms, self._n_terms))
            fact1 = 1. # use for calculating the factorial terms || 6461-6470
            for i in range(1, self._n_terms+1, 2):
                if i>1:
                    fact1 = self._fact_odds[i-2]
                fact2 = 1.
                for j in range(1, self._n_terms+1, 2):
                    if j>1:
                        fact2 = self._fact_odds[j-2]
                    factarr[i-1,j-1] = fact1*fact2

            # Calculate the moments
            for i in range(1, self._n_terms+1):
                #term2 = float(i - 2*(i/2))     # old version (C++)
                term2 = i%2                     # produces (1,0,1,0,1,0,1)
                term1 = sigma_x**(i-1)
                jmax = self.JMX(i-1)
                jmin = self.JMN(i-1)

                for j in range(jmin, jmax+1, 2):
                    mu_G[i-1,j-1] = factarr[i-1,j-1] * term1 * sigma_y**(j-1) * term2
        else:
            rho_tmp = delta_a / np.sqrt(delta_1 * delta_2) + 1.e-10
            rho_sq = rho_tmp**2
            rho = (1. - rho_sq)/2.

            for i in range(1, self._n_terms+1):
                # set iteration bounds
                jmax = self.JMX(i-1)
                jmin = self.JMN(i-1)
                for j in range(jmin, jmax+1, 2):
                    rs2 = rho_sq    # initialize variable, used in moment calculation
                    rs = 1./rho
                    term1 = 0.      # Another temp variable, initialize here
                    k_ct = int((i-1)/2 + 1) # Iteration limit for moment calculation

                    for k in range(1, k_ct+1):
                        rs = rs*rho
                        rs2 = rs2/rho_sq
                        term1 += self._mu_GN[i-1,j-1,k-1] * rs * rs2       #_mu_GN are the constant error distribution coefficients calculated once (see hermiteErrDistCoefs() )

                    # Calculate the moments mu_G
                    mu_G[i-1,j-1] = (sigma_x * rho_tmp)**(i-1) * sigma_y**(j-1) * term1
        
        H.setErrorDistCoefObject(mu_G)          # update heliostat object
        # -----end moments of error distribution ---------------

        # ------ begin moments of the mirror projection "M"-----
        """
        8343-8377	Flat, Focused, and Canted Heliostats

        For a flat heliostat, projection along t_hat of the differential mirror element (dxm, dym)
        at (xm,ym) onto the image plane located at the receiver is:
            xt = A11*xm + A12*ym
            yt = A21*xm + A22*ym
        Aij are defined by Eq (A-2) (Dellin) and are functions of time and position with respect to 
        the tower. By definition, the center of the heliostat (xm=0, ym=0) is projected onto the 
        origin of the image plane.

        Focusing is induced by displacements delta_x and delta_y in the mirror along the i_n j_n 
        directions. With displacements, the size of the image becomes:
            x = A11*xm + A12*ym + B11*R*delta_x + B12*R*delta_y
            y = A21*xm + A22*ym + B21*R*delta_x + B22*R*delta_y
        where R is the slant range from heliostat and Bij are defined in Eq A-6 (Dellin).
        """
        mu_M = H.getMirrorShapeCoefObject()
        mu_M = np.zeros((self._n_terms, self._n_terms))   # reshape and fill
        mu_Msave = np.zeros((self._n_terms, self._n_terms))

        # hermiteMirrorCoefs(errm_M, H, tht);	//Get the moments
        errm_M = H.getMirrorShapeNormCoefObject()
        xfocal = -tht/H.getFocalX()
        yfocal = -tht/H.getFocalY()

        E11 = A11+B11*xfocal*srange[1]/2.+1.e-10
        E12 = A12+B12*yfocal*srange[1]/2.+1.e-10
        E21 = A21+B21*xfocal*srange[1]/2.+1.e-10
        E22 = A22+B22*yfocal*srange[1]/2.+1.e-10

        step = E21*E21/(E22*E22)
        temp_res = e_ratio = binoms = S11 = S22 = 0
        start = np.zeros(2)

        # Loop over each term in the polynomial
        for i in range(1, self._n_terms+1):
            jmax = self.JMX(i-1) +1
            jmin = self.JMN(i-1)
            for j in range(jmin, jmax, 2):
                temp_res = 0.
                ij = i+j+1
                if j>1:
                    mstep = 1
                else:
                    mstep = 2
                
                S11 = E12**(i-1)
                e_ratio = (E11/E12)**mstep
                start[0] = E22**(j-1)
                start[1] = start[0]*E21/E22

                l_min = self.IMN(i-1)

                for k in range(1, i+1, mstep):
                    ijk = ij - k
                    binoms = self._binomials[i-1,k-1]*S11   # e_ratio;  moved the s11 calculation after this to reduce number of operations

                    S11 *= e_ratio
                    S22 = start[l_min - 1]
                    for l in range(l_min, j+1, 2):
                        S22 *= step
                        temp_res += binoms*S22*self._binomials[j-1,l-1]*errm_M[k+l-2,ijk-l-1]/step

                # calculate the moments to save
                mu_M[i-1,j-1] = temp_res * eta_consine
                mu_Msave[i-1,j-1] = mu_M[i-1,j-1]       # save a copy for other calcs

        panels = H.getPanels()                          # Get the matrix of all of the panels

        ncantx = H._n_cant_x
        ncanty = H._n_cant_y

        if H._is_multi_panels:
            # No canting=0; On-axis at slant=-1
            if H._cant_method == 0:
                # no canting
                gcanta = 0.0
                gcantx = 0.0
                gcantb = 0.0
                gcanty = 0.0
            elif H._cant_method == -1:
                # method -1 for on-axis at default slant range (6 tht)
                gcanta = 0.0
                gcanty = 0.0
                
                gcantx = -0.5*tht/H.getSlantRange()*srange[1]
                gcantb = gcantx
            else:
                print("ERROR: This canting method is not supported.")

            xcent = np.zeros((ncanty,ncantx))    
            ycent = np.zeros((ncanty,ncantx)) 

            xc = np.zeros((ncanty,ncantx, self._n_terms))
            yc = np.zeros((ncanty,ncantx, self._n_terms))

            # calculate optical terms for each canted panel
            thtinv = 1./tht
            for i in range(ncanty):         # over the rows of panels
                for j in range(ncantx):     # over the columns of panels
                    ploc = panels[i][j].getOrientation()
                    ploc_x = ploc.x * thtinv
                    ploc_y = ploc.y * thtinv
                    
                    xcent[i,j] = (A11 + B11*gcantx + B12*gcanty)*ploc_x + (A12 + B11*gcanta + B12*gcantb)*ploc_y
                    ycent[i,j] = (A21 + B21*gcantx + B22*gcanty)*ploc_x + (A22 + B21*gcanta + B22*gcantb)*ploc_y
                    xc[i,j,0] = 1.
                    yc[i,j,0] = 1.

                    # loop over each term in the Hermite expansion
                    xcentij = xcent[i,j]
                    ycentij = ycent[i,j]
                    for k in range(1,self._n_terms):
                        xc[i,j,k] = xc[i,j,k-1]*xcentij
                        yc[i,j,k] = yc[i,j,k-1]*ycentij

            # Adjust the image for multiple facets
            # DELSOL 8393-8407
            for i in range(1,self._n_terms):
                jmin = self.JMN(i-1)
                jmax = self.JMX(i-1)
                im1 = i-1
                ip1 = i+1
                for j in range(jmin, jmax+1, 2):
                    jp1 = j+1
                    jm1 = j-1
                    
                    temp_res = 0.
                    if j>1:
                        mstep = 1
                    else:
                        mstep = 2
                    
                    for m in range(ncanty):
                        for n in range(ncantx):
                            for ii in range(1, ip1, mstep):
                                iim1 = ii-1
                                S11 = self._binomials[im1, iim1]*xc[m,n,i-ii]
                                for jj in range(iim1%2+1, jp1, 2):
                                    jjm1 = jj-1
                                    temp_res += S11*self._binomials[jm1,jjm1]*mu_Msave[iim1,jjm1]*yc[m,n,j-jj]
                    # moments of M
                    mu_M[im1,jm1] = temp_res
     
        H.setMirrorShapeCoefObject(mu_M)    # update heliostat object

        # ------------end of mirror projection moments--------------

        # -----Combine moments of S, G, and M to get moments of F---------
	    # DELSOL3 8412-8450
        mu_F = H.getFluxMomentsObject()                 # unnecessary
        mu_F = np.zeros((self._n_terms,self._n_terms))  # reshape and fill

        comb_ord = mu_S[0,0] * mu_G[0,0] * mu_M[0,0]    # Combined moments at the ordinate (using local varables)
        comb_ord_inv = 1./comb_ord

        for m in range(1, self._n_terms+1):
            nmin = self.JMN(m-1)
            nmax = self.JMX(m-1)

            for n in range(nmin, nmax+1, 2):
                #ipak += 1
                temp_res = 0.0

                if (n>1):
                    mstep = 1
                else:
                    mstep = 2
                
                for k in range(1, m+1, mstep):
                    lmin = self.JMN(k-1)
                    binom_temp0 = self._binomials[m-1,k-1]
                    mk = m-k+1
                    kp1 = k+1
                    km1 = k-1
                    for l in range(lmin, n+1, 2):
                        ugs = mu_G[mk-1,n-l]*binom_temp0*self._binomials[n-1,l-1]

                        for i in range(1, kp1, 2):
                            ki = kp1 - i
                            binom_temp1 = self._binomials[km1,i-1]

                            for j in range(1, l+1, 2):
                                muS = mu_S[i-1, j-1]
                                muM = mu_M[ki-1, l-j]
                                term1 = self._binomials[l-1,j-1]*binom_temp1*muS*muM*ugs
                                temp_res += term1

                mu_F[m-1,n-1] = temp_res*comb_ord_inv

        H.setFluxMomentsObject(mu_F)        # update heliostat object

        # -------end combination-------------

        # Now evaluate the hermite coefficients and assign the spillage intercept value
        heval = self.hermiteIntEval(H, R)

        # TODO: Translate lines below
        """
        if(Rv->rec_type.mapval() == var_receiver::REC_TYPE::EXTERNAL_CYLINDRICAL 
            && h_rad < Receiver::getReceiverWidth( *Rv )/2.)
            return 0.0;
        else
            return heval;
        """
        return heval
    
    def hermiteIntEval(self, Heliostat, Receiver):
        """
        ###############################################################################################
        -------WHEN TO CALL------
        Automatically called by imagePlaneIntercept. Should be called for each flux intercept calculation.

        ---INFORMATION REQUIRED--
        Managed by imagePlaneIntercept. Requires all heliostat, flux, and receiver geometry; ambient

        ---------OUTPUT----------
        Returns spillage efficiency (1=no spillage) {type = double}
        ###############################################################################################

        "HERMIT"

        This subroutine evaluates the Hermite coefficients.

        Returns spillage efficiency (1=no spillage)
        """
        H = Heliostat
        R = Receiver

        mu_F = H.getFluxMomentsObject()

        nrf = 0
        for i in range(1, self._n_terms+1):
            for j in range(self.JMN(i-1), self.JMX(i-1)+1, 2):
                nrf += 1

        sig_x2 = np.sqrt(mu_F[2,0])             # x-direction image standard deviation
        sig_y2 = np.sqrt(mu_F[0,2])             # y-direction image standard deviation
        SigXY = np.array([sig_x2, sig_y2])      # put in an array
        # Set the standard deviation of the image (scaled by tower height) for the heliostat
        H.setImageSize(sig_x2, sig_y2)

        npak = int(nrf*4)

        # for optimization runs, don't do spillage calculations
        h_spill = np.zeros(nrf)
        axi = np.zeros(self._n_terms)
        ayi = np.zeros(self._n_terms)

        hcoef = H.getHermiteCoefObject()
        hcoef = np.zeros(npak/4)

        # detailed spillage calculations
        self.hermiteIntegralSetup(SigXY, H, h_spill, R)
        eta_spill = 0.

        axi[0] = 1.
        ayi[0] = 1.

        for i in range(1, self._n_terms):
            axi[i] = axi[i-1]/sig_x2
            ayi[i] = ayi[i-1]/sig_y2

        # Evaluate the hermite coefs
        ipak = 0
        save1 = axi[1]*ayi[1]/6.2832
        for i in range(self._n_terms+1):
            jmin = self.JMN(i-1)
            jmax = self.JMX(i-1)
            kmin = jmin     # small difference in calculation method, but the hard-coded values appear to always give the same results. DELSOL 1393

            save2 = save1 / self._fact_d[i-1]

            for j in range(jmin, jmax+1, 2):
                lmin =self.JMN(j-1)
                ipak += 1
                temp_res = 0.

                for k in range(kmin, i+1, 2):
                    for l in range(lmin, j+1, 2):
                        temp_res += self._binomials_hxn[i-1,k-1] * mu_F[k-1,l-1] * axi[k-1] * ayi[l-1] * self._binomials_hxn[j-1,l-1]

                temp_res /= self._fact_d[j-1] * save2
                eta_spill += h_spill[ipak-1]*temp_res
                hcoef[ipak-1] = temp_res            # This array is used to evaluate the flux density

        H.setHermiteCoefObject(hcoef)       # update heliostat object

        return eta_spill

    def hermiteIntegralSetup(self, SigXY, Heliostat, hspill, Receiver):
        """
        ###############################################################################################
        -------WHEN TO CALL------
        Automatically called by the hermiteIntEval method. Should be called each time the flux 
        intercept is calculated.

        ---INFORMATION REQUIRED--
        Automatically managed by hermiteIntEval

        ---------OUTPUT----------
        Passes data from hermiteIntegral through 'hspill' matrix
        ###############################################################################################	
        
        Calculate the amount of flux intercepted by the receiver. The receiver number is assigned 
        using the variable "rec", by default = 0.

        This subroutine is called by the hermiteIntEval method in evaluating the Hermite coefficients.

        DELSOL code 2041-2506

        The code manipulates and returns the hspill array of length(npak) ~typically x16
        """
        # Flux.cpp line 1564 to 2007

        pass

class AshleyGaussianCalculation(object):
    def __init__(self,sigma):
        """
        this object is used to generate calculations of flux using a 
        gaussian method.  This implementation uses the notation from the paper
        in Energy: Ashley et al. (2017),
        "Optimisation of aiming strategies in Solar Power Tower plants"
        
        sigma -- standard error of convolution of mirror, sunshape and geometry
            (mrad)
        """
        self.method_type = "Gaussian"
        self.solar_vec = scipy.array([0,0,1])  #high noon
        self.normal_vec = scipy.array([0,1,0])  #direct to y direction
        self.sigma = sigma  
    
    #Accessors
    def GetFluxMethod(self): return self.method_type
    def GetError(self): return self.sigma
    def GetSolarVec(self): return self.solar_vec
    def GetNormalVec(self): return self.normal_vec
    
    #Mutators
    def SetError(self,sigma):
        self.sigma = sigma
    def SetSolarVec(self,dx,dy,dz): 
        self.solar_vec = scipy.array([dx,dy,dz])   
    def SetNormalVec(self,dx,dy,dz): 
        self.normal_vec = scipy.array([dx,dy,dz])   
        
    
    def GetFlux(self,u,v,w,x,y,mu1=0,mu2=0):
        """
        Obtains the flux calculation from heliostat location (x,y) in the 
        direction of the vector w given receiver aimpoint (u,v) and error
        level theta.
        
        =============Parameters===========
        x -- x coordinate of heliostat location relative to receiver midpoint 
            (float)
        y -- y coordinate of heliostat location relative to receiver aimpoint 
            (float)
        u -- x coordinate of heliostat aimpoint (float)
        v -- y coordinate of heliostat aimpoint (float)
        w -- vector indicating direction of heliostat to measured target 
                (3-element scipy array)
        mu1 -- mean bias on x-axis?  stdev?
        mu2 -- mean bias on y-axis?  stdev?
        =============Returns==============
        flux: amount of thermal flux delivered from heliostat
        """
        return self.f1(w,x,y,mu1,mu2) * scipy.exp(
                -1.0*self.f2(w,u,v,x,y) / 
                2*(self.f3(w,x,y,mu1,mu2)**2)
                )
        
    def f1(self,w,x,y,mu1,mu2): 
        return self.f4(w) / (2 * scipy.pi 
                      * self.f3(w,x,y,mu1,mu2)**2 
                      * scipy.sum(w*w)
                      )
    
    def f2(self,w,u,v,x,y):
        return ( (u*u + v*v) / (2*scipy.sum(w*w)) ) * (
                (1 + self.f4(w)**2) + (scipy.absolute(u) * (1 - self.f4(w)**2))
                 / scipy.sqrt(u*u + v*v)
                )
    
    def f3(self,w,x,y,mu1,mu2):
        return scipy.sqrt(mu1**2 + ((mu2*(1-self.fcos(w)))/(4*scipy.sqrt(scipy.sum(w*w))))**2 )
    
    def f4(self,w):
        return (
                scipy.absolute(scipy.sum(w*self.GetNormalVec())) 
                / scipy.sqrt(scipy.sum(w*w))
                )
        
    def fcos(self,w):
        return scipy.sqrt(0.5 + (scipy.sum(w*self.GetSolarVec()))/(scipy.sqrt(scipy.sum(w*w))))

class HermiteCalculation(object):
    def __init__(self):
        """
        This class calculates a Hermite expansion for a collection of 
        measurement-based locations.  
        """
        self.hermite_terms = scipy.array([
            [1,0,0,0,0,0,0],
            [0,1,0,0,0,0,0],
            [-1,0,1,0,0,0,0],
            [0,-3,0,1,0,0,0],
            [3,0,-6,0,1,0,0],
            [0,15,0,-10,0,1,0],
            [-15,0,45,0,-15,0,1]
            ],dtype=float)
        self.gaussian_moments = scipy.array([1,0,1,0,3,0,15],dtype=float)
        
    def get_moments(self,func):
        pass
        
    def calculate_hermite_const(self,moments,i,j):
        return (moments*self.hermite_terms[i]).sum() * (moments*self.hermite_terms[j]).sum()

    def orthonormal_const(self,n):
        return (math.factorial(n))**-0.5

    def evaluate_hermite(self,x,herm_idx):
    #    print(x,herm_idx)
        Hx = 0
        for i in range(7):
            Hx += self.hermite_terms[herm_idx,i]*(x**i)
    #        print(i, hermite_terms[herm_idx,i], Hx)
        return Hx

    def full_hermite_eval(self,x,y,mx,my,sx,sy,moments):
        dx = (x-mx)/sx
        dy = (y-my)/sy
        pdf = (1./scipy.sqrt(2*scipy.pi*sx*sy))*scipy.exp(-0.5*(dx**2+dy**2))
        
        herm_sum = 0.
        for i in range(7):
            for j in range(7):
                herm_sum +=  self.calculate_hermite_const(moments,i,j) * self.evaluate_hermite(dx,i) * self.evaluate_hermite(dy,j) / (math.factorial(i) * math.factorial(j))
        return pdf * herm_sum

    def calculate_moment_xy(self,power_x,power_y,func,pts_per_dim,lb1,ub1,lb2,ub2,fname=None):
        """
        calcuate \mu_{xy} across a known support.
        """
        xs = lb1 + (0.5 + scipy.arange(pts_per_dim,dtype=float))*(ub1-lb1)/pts_per_dim
        ys = (lb2 + (0.5 + scipy.arange(pts_per_dim,dtype=float))*(ub2-lb2)).reshape([pts_per_dim,1])/pts_per_dim
        evals = scipy.zeros([pts_per_dim,pts_per_dim],dtype=float)
        for i in range(pts_per_dim):
            evals[i] = func(xs,ys[i])
        if fname != None: 
    #        normal_tests.plot_obj_heatmap(evals,fname)
            pass
        if power_x == 0: 
            xs = 1.
        if power_y == 0:
            ys = 1.
        return ((xs**power_x) * evals * (ys**power_y)).sum() * ((ub1-lb1)*(ub2-lb2)) / ((pts_per_dim)*(pts_per_dim))
    
    
    def calculate_1d_moment(self,xs,mx,sx,lp,power):
        p = scipy.stats.norm.pdf((xs-mx)/sx)
        px = (xs**power) * p * lp
        return px.sum()
    
    
    def gaussian_2d_func(self,x,y,mx=0,my=0,sx=1,sy=1):
        return scipy.stats.norm.pdf((x-mx)/sx)*scipy.stats.norm.pdf((y-my)/sy)/(sx*sy)
    
class SimpleNormalFluxCalc(object):
    def __init__(self,error=0, mirror=None):
        if mirror == None:
            self.sigma = error #in mrad
        else:
            self.sigma = mirror.error
        
    def GetFlux(self,helio,aim,measurement,normal,solar_vec,dni,approx=False,center_aim=None):
        """
        obtains the flux intensity at the measurement point given the heliostat
        location and aimpoint, and the normal vector to the heliostat as input.
        ======Parameters======
        helio - x,y,z coordinates of heliostat (3-element array of floats)
        aim - x,y,z coordinates of heliostat (3-element array of floats)
        measurement - x,y,z coordinates of heliostat 
            (3-element array of floats)
        solar_vec - solar vector (for cosine efficiency). i.e., vector from
            ground to the sun
        normal -- x,y,z vector of inward normal vector to receiver 
            at the measurement point
        dni -- direct normal irradiance from sun (W)
        =============Returns==============
        flux: amount of thermal flux delivered from heliostat at measurement
            point (W)
        """
        if approx:
            return self.GetFlux(helio,center_aim,(center_aim+measurement-aim),normal,solar_vec,dni)
        #generate three vectors (heliostat to aimpoint (h2a); heliostat
        # to measurement point (h2m); aimpoint to measurement point (a2m))
        h2a = aim-helio
        h2m = measurement-helio
        a2m = measurement-aim
        # generate the length of each vector; together these are the sides of
        # a triangle
        a = scipy.sqrt(scipy.sum(h2a*h2a))
        b = scipy.sqrt(scipy.sum(h2m*h2m))
        c = scipy.sqrt(scipy.sum(a2m*a2m))
        # Use the law of cosines to get the angle between h2a and h2m; convert
        # to milliradians (math.acos returns the angle in radians). This is the
        # angular 'error' which we divide by self.sigma to obtain the z-statistic
        # (and, in turn, the flux delivery to the measurement point.)
        err = math.acos((a*a+b*b-c*c)/(2*a*b)) * 1000  
        # Calculate cosine image efficiency; this is taken from Ashley et al. (2017)
        image_eff = scipy.sqrt(0.5 + 0.5*scipy.sum(h2m*solar_vec)/scipy.sqrt(scipy.sum(h2m*h2m)))  #taken from Ashley et al.
        # Cosine efficiency due to angle between mirror and inward normal
        #cosine_eff = scipy.sum(h2m*normal)/scipy.sqrt(scipy.sum(h2m*h2m))
        # There are three things that we still need to consider: 
        # (1) a normalization factor, since the integration of the normal 
        # will vary by the angle.  
        # (2) Atmospheric attenuation.
        # (3) relative blocking and 
        # But, for a very quick case study, this should be sufficient to 
        # get some quick insights via the optimization model.
#        reflective_eff = getReceiverReflectiveLosses(
#                getIncidenceAngle(helio, measurement, normal)
#                )
        #print(dni * scipy.stats.norm.pdf(err / self.sigma) * cosine_eff)
        return dni * scipy.stats.norm.pdf(err / self.sigma) * image_eff# * reflective_eff
    
def getIncidenceAngle(mirror_loc, measure_loc, measure_normal):
    """
    Parameters
    ----------
    mirror_loc : numpy.ndarray (float)
        x, y, z coordinates of mirror (heliostat) location [m]
        
    measure_loc : numpy.ndarray (float)
        x, y, z coordinates of aimpoint [m]
        
    measure_normal : numpy.ndarray (float)
        normal vector for receiver flux measurement point (unit vector) [m] 

    Returns
    -------
    theta : float
        Angle of Incidence with respect to receiver normal [degrees]
    """
        # vector from the receiver measurement point to heliostat
    r2h = mirror_loc - measure_loc
        # angle of incidence 
    theta = np.arccos(scipy.sum(measure_normal*r2h)/(scipy.sqrt(scipy.sum(r2h*r2h))))
    theta *= 180/np.pi   #convert to degrees
    
    if theta > 90:
        print("ERROR: Heliostat is behind surface", theta)

    return theta
     
def getReceiverReflectiveLosses(incidence_ang):
    """
    Parameters
    ----------
    incidence_ang : float
        Angle of Incidence with respect to receiver normal [degrees]
        MUST BE WITHIN 0 AND 90 DEGREES

    Returns
    -------
    f_abs : float
        Ratio of solar absorptance and solar absorptance at normal incidence
    """
    if incidence_ang > 90 or incidence_ang < 0:
        print("ERROR: Incidence angle is out of range!")
        
    # polynominal coefficients from equation 4.11.1 page 196 of Duffie and Beckman
    # Solar Engineering of Thermal Processes
    a = [x for x in range(8)]
    a[0] = 1.0
    a[1] = -1.5870e-3
    a[2] = 2.7314e-4
    a[3] = -2.3026e-5
    a[4] = 9.0244e-7
    a[5] = -1.8e-8
    a[6] = 1.7734e-10
    a[7] = -6.9937e-13
    
    f_abs = 0.0
    for i,coeff in enumerate(a):
        f_abs += coeff*(incidence_ang**i)
        
    if f_abs < 0.:  #for extreme angles, change small negative value to zero
        return 0.
    return f_abs

if __name__ == "__main__":
    #build geometry
    import geometry as rec
    import plotting

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    import sun_shape
    import heliostat

    # build a receiver
    norm = np.array( [0, 1, 0])
    params = {"length": 40, "height": 25, "pts_per_dim": 50, "normal": norm}  # testing normal calculations
    tht = 185
    shape_type = "flat_plate"
    r = rec.FlatPlateReceiver(shape_type, tht, params)
    r.build_measurement_points()

    var_dict = {"helio_name": "template_1",
                "id": 1,
                "location" : np.array([300, 400, 0]),
                "width": 12.2,
                "height": 12.2,
                "is_multi_panels": False,
                "cant_method": 0,
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
                "tht": tht
                }

    h = heliostat.Heliostat(var_dict)

    # Making Sunshape
    model_params = {"error": 2.78}
    ## sunshape models
    #sun = sun_shape.PillBoxSunshape(model_params)
    #sun = sun_shape.LimbDarkenedSunshape(model_params)
    #sun = sun_shape.SinglePointSun(model_params)
    #sun = sun_shape.GaussianSunshape(model_params)
    model_params = {"error": 2.78, "chi": 0.1}
    sun = sun_shape.BuieCSRSunshape(model_params)
    sunVector = tb.Vector(np.array([0,0,1]))

    #TESTING FLUX CLASS
    f = Flux()
    f.initHermiteCoefs(sun)
    f.hermiteMirrorCoefs(h,tht)
    f.imagePlaneIntercept(h, r, sunVector)

    ## Checking moment calculations
   

    # mu_MN
    # fig, ax = plt.subplots()
    # cs = ax.imshow(h._mu_MN,norm=LogNorm())
    # cbar = fig.colorbar(cs)
    # plt.show()

    # mu_S
    moments = ['_mu_S', '_mu_M', '_mu_G', '_mu_F' ]
    fig, axs = plt.subplots(2,2)
    cs0 = axs[0, 0].imshow(h._mu_S,norm=LogNorm())
    axs[0, 0].set_title(moments[0])
    plt.colorbar(cs0, ax = axs[0, 0])

    cs1 = axs[0, 1].imshow(h._mu_M,norm=LogNorm())
    axs[0, 1].set_title(moments[1])
    plt.colorbar(cs1, ax = axs[0, 1])

    cs2 = axs[1, 0].imshow(h._mu_G,norm=LogNorm())
    axs[1, 0].set_title(moments[2])
    plt.colorbar(cs2, ax = axs[1, 0])

    cs3 = axs[1, 1].imshow(h._mu_F,norm=LogNorm())
    axs[1, 1].set_title(moments[3])
    plt.colorbar(cs3, ax = axs[1, 1])
    plt.show()


    
    """
    #TESTING ASHLEY MODEL
    pts_per_dim = 21
    params = {"length":40, "height":25, "pts_per_dim":pts_per_dim}
    h = 185
    shape_type = "flat_plate"
    r = geometry.Receiver(shape_type,h,params)
    r.buildFlatPlate()
    
    c = GaussianCalculation(2)
    
    #create flux map
    flux = scipy.zeros([pts_per_dim,pts_per_dim])
    
    u = 1
    v = 1
    helio_x = 0
    helio_y = 500
    helio_z = 0
    aim_x = 0
    aim_y = 0
    aim_z = h
    for xidx in range(pts_per_dim):
        for zidx in range(pts_per_dim):
            x = r.x[xidx,zidx]
            y = 0
            z = r.z[xidx,zidx]  
            w = scipy.array([
                        helio_x-x,
                        helio_y-y,
                        helio_z-z
                    ])
#            print(w)
            flux[xidx,zidx] = c.GetFlux(u,v,w,helio_x,helio_y,5,25)
    
    plotting.plot_obj_heatmap(flux,"ashley_map.pdf")
    """
    
    """
    #TESTING SIMPLE GAUSSIAN, building example map
    pts_per_dim = 51
    params = {"length":25, "height":25, "pts_per_dim":pts_per_dim}
    
    norm = scipy.array([0, scipy.sqrt(2)/2, -scipy.sqrt(2)/2])  # 45 degrees tilt down (zenith = 135 deg) (normalized)
    params = {"length":25, "height":25, "pts_per_dim":pts_per_dim, "normal":norm}
    #params = {"length":25, "height":25, "pts_per_len_dim":51, "pts_per_ht_dim":31, "normal":norm}
    
    h = 185
    shape_type = "flat_plate"
    r = geometry.Receiver(shape_type,h,params)
    r.buildFlatPlate()
    
    c = SimpleNormalFluxCalc(5)
    
    #create flux map
    flux = scipy.zeros_like(r.x)
    f1 = scipy.zeros_like(r.x)
    
    solar_vec = scipy.array([0,0,-1])
    
    helio = scipy.array([500,500,0])
    aim = scipy.array([4,0,h+7.5])
    for i in range(pts_per_dim):
        for j in range(pts_per_dim):
            measurement = scipy.array([r.x[i,j],r.y[i,j],r.z[i,j]])
            flux[i,j] += c.GetFlux(helio,aim,measurement,solar_vec,1000/scipy.sqrt(2))
            f1[i,j] = c.GetFlux(helio,aim,measurement,solar_vec,600)

#    print("helio1:",f1.sum()/(pts_per_dim*pts_per_dim))
    
    helio = scipy.array([250,500,0])
    aim = scipy.array([4,0,h-4])
    for i in range(pts_per_dim):
        for j in range(pts_per_dim):
            measurement = scipy.array([r.x[i,j],r.y[i,j],r.z[i,j]])
            flux[i,j] += c.GetFlux(helio,aim,measurement,solar_vec,1000)
            f1[i,j] = c.GetFlux(helio,aim,measurement,solar_vec,1000)
#    print("helio2:",f1.sum()/(pts_per_dim*pts_per_dim))
            
    helio = scipy.array([-250,500,0])
    aim = scipy.array([-4,0,h+4])
    for i in range(pts_per_dim):
        for j in range(pts_per_dim):
            measurement = scipy.array([r.x[i,j],r.y[i,j],r.z[i,j]])
            flux[i,j] += c.GetFlux(helio,aim,measurement,solar_vec,1000)
            f1[i,j] = c.GetFlux(helio,aim,measurement,solar_vec,1000)
#    print("helio3:",f1.sum()/(pts_per_dim*pts_per_dim))
            
    helio = scipy.array([-500,500,0])
    aim = scipy.array([-7.5,0,h-4])
    for i in range(pts_per_dim):
        for j in range(pts_per_dim):
            measurement = scipy.array([r.x[i,j],r.y[i,j],r.z[i,j]])
            flux[i,j] += c.GetFlux(helio,aim,measurement,solar_vec,1000/scipy.sqrt(2))
            f1[i,j] = c.GetFlux(helio,aim,measurement,solar_vec,600)
#    print("helio4:",f1.sum()/(pts_per_dim*pts_per_dim))

#    helio = scipy.array([0,500,h])
#    aim = scipy.array([0,0,h])
#    for i in range(pts_per_dim):
#        for j in range(pts_per_dim):
#            measurement = scipy.array([r.x[i,j],r.y[i,j],r.z[i,j]])
#            flux[i,j] += c.GetFlux(helio,aim,measurement,solar_vec,1000)
#            f1[i,j] = c.GetFlux(helio,aim,measurement,solar_vec,1000)
#    print("helio5:",f1.sum()/(pts_per_dim*pts_per_dim))
    
    plotting.plot_obj_heatmap(flux,"gaussian_map_4_mirrors_diag_200.pdf")
    
    diag = flux.sum()
    diagmax = flux.max()
    
    flux = scipy.zeros_like(r.x)
    
    helio = scipy.array([500,500,0])
    aim = scipy.array([4,0,h])
    solar_vec = scipy.array([0,0,1]) #'high noon'
    for i in range(pts_per_dim):
        for j in range(pts_per_dim):
            measurement = scipy.array([r.x[i,j],r.y[i,j],r.z[i,j]])
            flux[i,j] += c.GetFlux(helio,aim,measurement,solar_vec,600)
    
    helio = scipy.array([250,500,0])
    aim = scipy.array([0,0,h-4])
    for i in range(pts_per_dim):
        for j in range(pts_per_dim):
            measurement = scipy.array([r.x[i,j],r.y[i,j],r.z[i,j]])
            flux[i,j] += c.GetFlux(helio,aim,measurement,solar_vec,1000)
            
    helio = scipy.array([-250,500,0])
    aim = scipy.array([0,0,h+4])
    for i in range(pts_per_dim):
        for j in range(pts_per_dim):
            measurement = scipy.array([r.x[i,j],r.y[i,j],r.z[i,j]])
            flux[i,j] += c.GetFlux(helio,aim,measurement,solar_vec,1000)
            
    helio = scipy.array([-500,500,0])
    aim = scipy.array([-4,0,h])
    for i in range(pts_per_dim):
        for j in range(pts_per_dim):
            measurement = scipy.array([r.x[i,j],r.y[i,j],r.z[i,j]])
            flux[i,j] += c.GetFlux(helio,aim,measurement,solar_vec,600)
            
    #plotting.plot_obj_heatmap(flux,"gaussian_map_4_mirrors_normal_200.pdf")
    
    
    cross = flux.sum()
    crossmax = flux.max()
    
    
    print("diag:",diag,diagmax,diag*(params["length"]*params["height"])/(pts_per_dim**2))
    print("cross:",cross, crossmax,cross*(params["length"]*params["height"])/(pts_per_dim**2))
    """
    
    
    """
    #TESTING HERMITE CALCULATION
    pts_per_dim = 50
    dist_per_pt = 10 / 50.
    xs = -5 + (0.5 + scipy.arange(pts_per_dim,dtype=float))*(10)/pts_per_dim
    print(xs)
    #lp = 
    
    map1 = scipy.zeros([pts_per_dim,pts_per_dim],dtype=float)
    map2 = scipy.zeros([pts_per_dim,pts_per_dim],dtype=float)
    
    for i in range(pts_per_dim):
        print (i,"working")
        for j in range(pts_per_dim):
            map1[i,j] = gaussian_2d_func(xs[i],xs[j],mx=0,my=0,sx=1,sy=1)
            map2[i,j] = full_hermite_eval(xs[i],xs[j],0,0,1,1,gaussian_moments)
    print(map1, map2)
    """
#normal_tests.plot_obj_heatmap(map1,"gaussian_map.pdf")
#normal_tests.plot_obj_heatmap(map2,"hermite_map.pdf")


