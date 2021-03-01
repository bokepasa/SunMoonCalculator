#MIT License

#Copyright (c) 2020 Leandro Castro Pérez

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.


# *** CREDIT ***
#This implementation is based on a Java class by T. Alonso Albi from OAN (Spain)
#published at http://conga.oan.es/~alonso/doku.php?id=blog:sun_moon_position

#################################################################################################################################
# A very simple yet accurate Sun/Moon calculator without using JPARSEC library.
# @author T. Alonso Albi - OAN (Spain), email t.alonso@oan.es
# @version August 7, 2020 (forgot to add equation of equinoxes to lst, affects 1s at most to rise/set/transit times)
# @version July 1, 2020 (improved aberration in getSun, illumination phase for planets fixed, refraction with ambient P/T and extensible to radio wavelengths)
# @version June 15, 2020 (more terms for Sun, equinoxes/solstices methods revised to improve accuracy and performance)
# @version June 12, 2020 (added planetary ephemerides from another class, better performance)
# *** CODE PARTIALLY REWRITTEN, PLEASE UPDATE EVERYTHING ***
# *** INCLUDING THE SECTIONS BELOW FOR A CORRECT VERSION ***
# @version June 3, 2020 (fix rise/set problem for high latitudes, return azimuth and elevation even in geocentric computations)
# @version November 26, 2018 (two new methods getCulminationTime and getAzimuthTime)
# @version November 6, 2018 (better accuracy for Moon, angular radius in ephemeris, cosmetic improvements)
# @version July 24, 2018 (new class to hold results, illumination phase, moon phases, equinoxes and solstices)
# @version May 25, 2017 (fixed nutation correction and moon age, better accuracy in Moon)
#################################################################################################################################

import sys
import enum
import math
import numpy as np


# Radians to degrees.
RAD_TO_DEG = 180.0 / math.pi

# Degrees to radians.
DEG_TO_RAD = 1.0 / RAD_TO_DEG


class SunMoonCalculator:

    # Arcseconds to radians
    ARCSEC_TO_RAD = (DEG_TO_RAD / 3600.0)

    # Astronomical Unit in km. As defined by JPL.
    AU = 149597870.691

    # Earth equatorial radius in km. IERS 2003 Conventions.
    EARTH_RADIUS = 6378.1366

    # Two times Pi.
    TWO_PI = 2.0 * math.pi

    # Pi divided by two.
    PI_OVER_TWO = math.pi / 2.0

    # Julian century conversion constant = 100 * days per year.
    JULIAN_DAYS_PER_CENTURY = 36525.0

    # Seconds in one day.
    SECONDS_PER_DAY = 86400

    # Light time in days for 1 AU. DE405 definition.
    LIGHT_TIME_DAYS_PER_AU = 0.00577551833109

    # Our default epoch. The Julian Day which represents noon on 2000-01-01.
    J2000 = 2451545.0


    # The set of twilights to calculate (types of rise/set events).
    class TWILIGHT(enum.Enum):
        #/**
        #* Event ID for calculation of rising and setting times for astronomical
        #* twilight. In this case, the calculated time will be the time when the
        #* center of the object is at -18 degrees of geometrical elevation below the
        #* astronomical horizon. At this time astronomical observations are possible
        #* because the sky is dark enough.
        #*/
        ASTRONOMICAL = 4
        #/**
		# * Event ID for calculation of rising and setting times for nautical
		# * twilight. In this case, the calculated time will be the time when the
		# * center of the object is at -12 degrees of geometric elevation below the
		# * astronomical horizon.
		# */
        NAUTICAL = 3
        #/**
		# * Event ID for calculation of rising and setting times for civil twilight.
		# * In this case, the calculated time will be the time when the center of the
		# * object is at -6 degrees of geometric elevation below the astronomical
		# * horizon.
		# */
        CIVIL = 2
        #/**
		# * The standard value of 34' for the refraction at the local horizon.
		# */
        HORIZON_34arcmin = 1
    
    # The set of events to calculate (rise/set/transit events).
    class EVENT(enum.Enum):
        RISE = 'RISE'
        SET = 'SET'
        TRANSIT = 'TRANSIT'
    
    # The set of phases to compute the moon phases.
    class MOONPHASE(enum.Enum):
        NEW_MOON = ("New Moon:        ", 0)
        CRESCENT_QUARTER = ("Crescent quarter:", 0.25)
        FULL_MOON = ("Full Moon:       ", 0.5)
        DESCENT_QUARTER = ("Descent quarter: ", 0.75)
    
    # The set of bodies to compute ephemerides.
    class BODY(enum.Enum):
        
        MERCURY = (0, 2439.7)
        VENUS = (1, 6051.8)
        MARS = (3, 3396.19)
        JUPITER = (4, 71492)
        SATURN = (5, 60268)
        URANUS = (6, 25559)
        NEPTUNE = (7, 24764)
        Moon = (-2, 1737.4)
        Sun = (-1, 696000)
        EMB = (2, 0)

        @property
        def eqRadius(self):
            return self.value[1]

    # Input values and nutation/obliquity parameters only calculated once.
    jd_UT = 0
    t = 0
    obsLon = 0 
    obsLat = 0
    TTminusUT = 0
    twilight = TWILIGHT.HORIZON_34arcmin
    nutLon = 0 
    nutObl = 0
    meanObliquity = 0

    lst = 0

    #/**
	# * Class to hold the results of ephemerides.
	# * @author T. Alonso Albi - OAN (Spain)
	# */
    class Ephemeris(object):

        azimuth = 0
        elevation = 0
        rise = 0
        set = 0
        transit = 0
        transitElevation = 0
        rightAscension = 0
        declination = 0
        distance = 0
        illuminationPhase=100
        eclipticLongitude = 0
        eclipticLatitude = 0
        angR = 0

        def __init__(self, azimuth, elevation, rise, set, transit, transitElevation, rightAscension, declination, distance,
			eclipticLongitude, eclipticLatitude, angR):
            self.azimuth = azimuth
            self.elevation = elevation
            self.rise = rise
            self.set = set
            self.transit = transit
            self.transitElevation = transitElevation
            self.distance = distance
            self.rightAscension = rightAscension
            self.declination = declination
            self.eclipticLongitude = eclipticLongitude
            self.eclipticLatitude = eclipticLatitude
            self.angularRadius = angR
                 
    # Ephemeris for the Sun and Moon bodies.
    sun = Ephemeris
    moon = Ephemeris

    # Moon's age in days as an independent variable.
    moonAge=None


    def __init__(self, year, month, day, hour, minute, second, obsLon, obsLat):
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.minute = minute
        self.second = second
        self.obsLon = obsLon
        self.obsLat = obsLat

        jd = self.toJulianDay(self.year,self.month,self.day, self.hour, self.minute,self.second)
        TTminusUT = 0.0
        if (year > -600 and year < 2200):
            x = year + (month - 1 + day / 30.0) / 12.0
            x2 = x * x 
            x3 = x2 * x 
            x4 = x3 * x
            if (year < 1600):
                self.TTminusUT = 10535.328003 - 9.9952386275 * x + 0.00306730763 * x2 - 7.7634069836E-6 * x3 + 3.1331045394E-9 * x4 + 8.2255308544E-12 * x2 * x3 - 7.4861647156E-15 * x4 * x2 + 1.936246155E-18 * x4 * x3 - 8.4892249378E-23 * x4 * x4
            else:
                self.TTminusUT = -1027175.3477559977 + 2523.256625418965 * x - 1.885686849058459 * x2 + 5.869246227888417E-5 * x3 + 3.3379295816475025E-7 * x4 + \
                    1.7758961671447929E-10 * x2 * x3 - 2.7889902806153024E-13 * x2 * x4 + 1.0224295822336825E-16 * x3 * x4 - 1.2528102370680435E-20 * x4 * x4
    
        self.setUTDate(jd)    
    
    # /**
    # Sets the rise/set times to return. Default is for the local horizon.
    # @param t The Twilight.
    # */
    def setTwilight(self,t):
        self.twilight = t

    # **
	# * Sets the UT date from the provided Julian day and computes the nutation, obliquity, and 
	# * sidereal time. TT minuts UT1 is not updated since it changes very slowly with time.
	# * Use this only to update the computation time for the same year as the one used when the 
	# * instance was created.
	# * @param jd The new Julian day in UT.
	# */
    def setUTDate(self,jd):
        self.jd_UT = jd
        self.t = (jd + self.TTminusUT / self.SECONDS_PER_DAY - self.J2000) / self.JULIAN_DAYS_PER_CENTURY

        # Compute nutation
        M1 = (124.90 - 1934.134 * self.t + 0.002063 * self.t * self.t) * DEG_TO_RAD
        M2 = (201.11 + 72001.5377 * self.t + 0.00057 * self.t * self.t) * DEG_TO_RAD
        self.nutLon = (-0.0047785 * math.sin(M1) - 0.0003667 * math.sin(M2)) * DEG_TO_RAD
        self.nutObl = (0.002558 * math.cos(M1) - 0.00015339 * math.cos(M2)) * DEG_TO_RAD

        # Compute mean obliquity
        t2 = self.t / 100.0
        tmp = t2 * (27.87 + t2 * (5.79 + t2 * 2.45))
        tmp = t2 * (-249.67 + t2 * (-39.05 + t2 * (7.12 + tmp)))
        tmp = t2 * (-1.55 + t2 * (1999.25 + t2 * (-51.38 + tmp)))
        tmp = (t2 * (-4680.93 + tmp)) / 3600.0
        self.meanObliquity = (23.4392911111111 + tmp) * DEG_TO_RAD

        # Obtain local apparent sidereal time
        jd0 = math.floor(self.jd_UT - 0.5) + 0.5
        T0 = (jd0 - self.J2000) / self.JULIAN_DAYS_PER_CENTURY
        secs = (self.jd_UT - jd0) * self.SECONDS_PER_DAY
        gmst = (((((-6.2e-6 * T0) + 9.3104e-2) * T0) + 8640184.812866) * T0) + 24110.54841
        msday = 1.0 + (((((-1.86e-5 * T0) + 0.186208) * T0) + 8640184.812866) / (self.SECONDS_PER_DAY * self.JULIAN_DAYS_PER_CENTURY))
        gmst = (gmst + msday * secs) * (15.0 / 3600.0) * DEG_TO_RAD
        self.lst = self.normalizeRadians(gmst + obsLon + self.nutLon * math.cos(self.meanObliquity + self.nutObl))

    # ** Calculates everything for the Sun and the Moon.
    def calcSunAndMoon(self):
        jd = self.jd_UT

        # First the Sun
        self.sun = self.doCalc(self.getSun(),False)

        niter = 15 # // Number of iterations to get accurate rise/set/transit times
        self.sun.rise = self.obtainAccurateRiseSetTransit(self.sun.rise, self.EVENT.RISE, niter, True)
        self.sun.set= self.obtainAccurateRiseSetTransit(self.sun.set, self.EVENT.SET, niter, True)
        self.sun.transit = self.obtainAccurateRiseSetTransit(self.sun.transit, self.EVENT.TRANSIT, niter, True)
        if (self.sun.transit == -1):
            self.sun.transitElevation = 0
        else:
            #// Update Sun's maximum elevation
            self.setUTDate(self.sun.transit)
            self.sun.transitElevation = self.doCalc(self.getSun(),False).transitElevation
        
        # Now Moon
        self.setUTDate(jd)
        self.moon = self.doCalc(self.getMoon(),False)
        ma = self.moonAge
        #//niter = 15 #// Number of iterations to get accurate rise/set/transit times
        self.moon.rise = self.obtainAccurateRiseSetTransit(self.moon.rise, self.EVENT.RISE, niter, False)
        self.moon.set = self.obtainAccurateRiseSetTransit(self.moon.set, self.EVENT.SET, niter, False)
        self.moon.transit = self.obtainAccurateRiseSetTransit(self.moon.transit, self.EVENT.TRANSIT, niter, False)
        if (self.moon.transit == -1):
	        self.moon.transitElevation = 0
        else: 
	        #// Update Moon's maximum elevation
	        self.setUTDate(self.moon.transit)
	        self.moon.transitElevation = self.doCalc(self.getMoon(),False).transitElevation

        self.setUTDate(jd)
        self.moonAge = ma
		
        #// Compute illumination phase percentage for the Moon
        self.getIlluminationPhase(self.moon)        
        
            
    # Sun data from the expansion "Planetary Programs 
	# and Tables" by Pierre Bretagnon and Jean-Louis 
	# Simon, Willman-Bell, 1986
    sun_elements = np.array(([403406.0, 0.0, 4.721964, 1.621043], 
    [195207.0, -97597.0, 5.937458, 62830.348067], 
    [119433.0, -59715.0, 1.115589, 62830.821524],
    [ 112392.0, -56188.0, 5.781616, 62829.634302 ], 
    [ 3891.0, -1556.0, 5.5474, 125660.5691 ], 
    [ 2819.0, -1126.0, 1.512, 125660.9845 ], 
    [ 1721.0, -861.0, 4.1897, 62832.4766 ], 
    [ 0.0, 941.0, 1.163, .813 ], 
    [ 660.0, -264.0, 5.415, 125659.31 ], 
    [ 350.0, -163.0, 4.315, 57533.85 ], 
    [ 334.0, 0.0, 4.553, -33.931 ], 
    [ 314.0, 309.0, 5.198, 777137.715 ], 
    [ 268.0, -158.0, 5.989, 78604.191 ], 
    [ 242.0, 0.0, 2.911, 5.412 ], 
    [ 234.0, -54.0, 1.423, 39302.098 ], 
    [ 158.0, 0.0, .061, -34.861 ], 
    [ 132.0, -93.0, 2.317, 115067.698 ], 
    [ 129.0, -20.0, 3.193, 15774.337 ], 
    [ 114.0, 0.0, 2.828, 5296.67 ], 
    [ 99.0, -47.0, .52, 58849.27 ], 
    [ 93.0, 0.0, 4.65, 5296.11 ], 
    [ 86.0, 0.0, 4.35, -3980.7 ], 
    [ 78.0, -33.0, 2.75, 52237.69 ], 
    [ 72.0, -32.0, 4.5, 55076.47 ], 
    [ 68.0, 0.0, 3.23, 261.08 ], 
    [ 64.0, -10.0, 1.22, 15773.85 ], 
    [ 46.0, -16.0, .14, 188491.03 ], 
    [ 38.0, 0.0, 3.44, -7756.55 ], 
    [ 37.0, 0.0, 4.37, 264.89 ], 
    [ 32.0, -24.0, 1.14, 117906.27 ], 
    [ 29.0, -13.0, 2.84, 55075.75 ], 
    [ 28.0, 0.0, 5.96, -7961.39 ], 
    [ 27.0, -9.0, 5.09, 188489.81 ], 
    [ 27.0, 0.0, 1.72, 2132.19 ], 
    [ 25.0, -17.0, 2.56, 109771.03 ], 
    [ 24.0, -11.0, 1.92, 54868.56 ], 
    [ 21.0, 0.0, .09, 25443.93 ], 
    [ 21.0, 31.0, 5.98, -55731.43 ], 
    [ 20.0, -10.0, 4.03, 60697.74 ], 
    [ 18.0, 0.0, 4.27, 2132.79 ], 
    [ 17.0, -12.0, .79, 109771.63 ], 
    [ 14.0, 0.0, 4.24, -7752.82 ], 
    [ 13.0, -5.0, 2.01, 188491.91 ], 
    [ 13.0, 0.0, 2.65, 207.81 ], 
    [ 13.0, 0.0, 4.98, 29424.63 ], 
    [ 12.0, 0.0, .93, -7.99 ], 
    [ 10.0, 0.0, 2.21, 46941.14 ], 
    [ 10.0, 0.0, 3.59, -68.29 ], 
    [ 10.0, 0.0, 1.5, 21463.25 ], 
    [ 10.0, -9.0, 2.55, 157208.4 ]))


    def getSun(self):
        L = 0.0
        R = 0.0
        t2 = self.t * 0.01
        Lp = 0.0
        deltat = 0.5
        t2p = (self.t + deltat / self.JULIAN_DAYS_PER_CENTURY) * 0.01

        for i in range(len(self.sun_elements)):
            v = self.sun_elements[i][2] + self.sun_elements[i][3] * t2
            u = self.normalizeRadians(v)
            L = L + self.sun_elements[i][0] * math.sin(u)
            R = R + self.sun_elements[i][1] * math.cos(u)
            vp = self.sun_elements[i][2] + self.sun_elements[i][3] * t2p
            up = self.normalizeRadians(vp)
            Lp = Lp + self.sun_elements[i][0] * math.sin(up)

        lon = self.normalizeRadians(4.9353929 + self.normalizeRadians(62833.196168 * t2) + L / 10000000.0) * RAD_TO_DEG
        sdistance = 1.0001026 + R / 10000000.0

        # // Now subtract aberration.
        dlon = ((Lp - L) / 10000000.0 + 62833.196168 * (t2p - t2)) / deltat
        aberration = dlon * sdistance * 0.00577551833109
        lon -= aberration * RAD_TO_DEG

        slongitude = lon * DEG_TO_RAD # apparent longitude (error<0.001 deg)
        slatitude = 0 # Sun's ecliptic latitude is always negligible

        array = [slongitude, slatitude, sdistance, math.atan(self.BODY.Sun.eqRadius / (self.AU * sdistance))]

        return array
    
    def getMoon(self):
        t = self.t
        #// MOON PARAMETERS (Formulae from "Calendrical Calculations")
        phase = self.normalizeRadians((297.8502042 + 445267.1115168 * t - 0.00163 * t * t + t * t * t / 538841 - t * t * t * t / 65194000) * DEG_TO_RAD)

		#// Anomalistic phase
        anomaly = (134.9634114 + 477198.8676313 * t + .008997 * t * t + t * t * t / 69699 - t * t * t * t / 14712000)
        anomaly = anomaly * DEG_TO_RAD

		#// Degrees from ascending node
        node = (93.2720993 + 483202.0175273 * t - 0.0034029 * t * t - t * t * t / 3526000 + t * t * t * t / 863310000)
        node = node * DEG_TO_RAD

        E = 1.0 - (.002495 + 7.52E-06 * (t + 1.0)) * (t + 1.0)

		#// Solar anomaly
        sanomaly = (357.5291 + 35999.0503 * t - .0001559 * t * t - 4.8E-07 * t * t * t) * DEG_TO_RAD

		#// Now longitude, with the three main correcting terms of evection,
		#// variation, and equation of year, plus other terms (error<0.01 deg)
		#// P. Duffet's MOON program taken as reference
        l = (218.31664563 + 481267.8811958 * t - .00146639 * t * t + t * t * t / 540135.03 - t * t * t * t / 65193770.4)
        l += 6.28875 * math.sin(anomaly) + 1.274018 * math.sin(2 * phase - anomaly) + .658309 * math.sin(2 * phase)
        l +=  0.213616 * math.sin(2 * anomaly) - E * .185596 * math.sin(sanomaly) - 0.114336 * math.sin(2 * node)
        l += .058793 * math.sin(2 * phase - 2 * anomaly) + .057212 * E * math.sin(2 * phase - anomaly - sanomaly) + .05332 * math.sin(2 * phase + anomaly)
        l += .045874 * E * math.sin(2 * phase - sanomaly) + .041024 * E * math.sin(anomaly - sanomaly) - .034718 * math.sin(phase) - E * .030465 * math.sin(sanomaly + anomaly)
        l += .015326 * math.sin(2 * (phase - node)) - .012528 * math.sin(2 * node + anomaly) - .01098 * math.sin(2 * node - anomaly) + .010674 * math.sin(4 * phase - anomaly)
        l += .010034 * math.sin(3 * anomaly) + .008548 * math.sin(4 * phase - 2 * anomaly)
        l += -E * .00791 * math.sin(sanomaly - anomaly + 2 * phase) - E * .006783 * math.sin(2 * phase + sanomaly) + .005162 * math.sin(anomaly - phase) + E * .005 * math.sin(sanomaly + phase)
        l += .003862 * math.sin(4 * phase) + E * .004049 * math.sin(anomaly - sanomaly + 2 * phase) + .003996 * math.sin(2 * (anomaly + phase)) + .003665 * math.sin(2 * phase - 3 * anomaly)
        l += E * 2.695E-3 * math.sin(2 * anomaly - sanomaly) + 2.602E-3 * math.sin(anomaly - 2*(node+phase))
        l += E * 2.396E-3 * math.sin(2*(phase - anomaly) - sanomaly) - 2.349E-3 * math.sin(anomaly+phase)
        l += E * E * 2.249E-3 * math.sin(2*(phase-sanomaly)) - E * 2.125E-3 * math.sin(2*anomaly+sanomaly)
        l += -E * E * 2.079E-3 * math.sin(2*sanomaly) + E * E * 2.059E-3 * math.sin(2*(phase-sanomaly)-anomaly)
        l += -1.773E-3 * math.sin(anomaly+2*(phase-node)) - 1.595E-3 * math.sin(2*(node+phase))
        l += E * 1.22E-3 * math.sin(4*phase-sanomaly-anomaly) - 1.11E-3 * math.sin(2*(anomaly+node))
        longitude = l * DEG_TO_RAD
				
        Psin = 29.530588853
        if (self.sun != None):
			#// Get Moon age, more accurate than 'phase' but we need the Sun position
            self.moonAge = self.normalizeRadians(longitude - self.sun.eclipticLongitude) * Psin / self.TWO_PI
        else:
			#// Use the phase variable as estimate, less accurate but this is used only when we don't need an accurate value
            self.moonAge = phase * Psin / self.TWO_PI			
		
		#// Now Moon parallax
        parallax = .950724 + .051818 * math.cos(anomaly) + .009531 * math.cos(2 * phase - anomaly)
        parallax += .007843 * math.cos(2 * phase) + .002824 * math.cos(2 * anomaly)
        parallax += 0.000857 * math.cos(2 * phase + anomaly) + E * .000533 * math.cos(2 * phase - sanomaly)
        parallax += E * .000401 * math.cos(2 * phase - anomaly - sanomaly) + E * .00032 * math.cos(anomaly - sanomaly) - .000271 * math.cos(phase)
        parallax += -E * .000264 * math.cos(sanomaly + anomaly) - .000198 * math.cos(2 * node - anomaly)
        parallax += 1.73E-4 * math.cos(3 * anomaly) + 1.67E-4 * math.cos(4*phase-anomaly)

		#// So Moon distance in Earth radii is, more or less,
        distance = 1.0 / math.sin(parallax * DEG_TO_RAD)

		#// Ecliptic latitude with nodal phase (error<0.01 deg)
        l = 5.128189 * math.sin(node) + 0.280606 * math.sin(node + anomaly) + 0.277693 * math.sin(anomaly - node)
        l += .173238 * math.sin(2 * phase - node) + .055413 * math.sin(2 * phase + node - anomaly)
        l += .046272 * math.sin(2 * phase - node - anomaly) + .032573 * math.sin(2 * phase + node)
        l += .017198 * math.sin(2 * anomaly + node) + .009267 * math.sin(2 * phase + anomaly - node)
        l += .008823 * math.sin(2 * anomaly - node) + E * .008247 * math.sin(2 * phase - sanomaly - node) + .004323 * math.sin(2 * (phase - anomaly) - node)
        l += .0042 * math.sin(2 * phase + node + anomaly) + E * .003372 * math.sin(node - sanomaly - 2 * phase)
        l += E * 2.472E-3 * math.sin(2 * phase + node - sanomaly - anomaly)
        l += E * 2.222E-3 * math.sin(2 * phase + node - sanomaly)
        l += E * 2.072E-3 * math.sin(2 * phase - node - sanomaly - anomaly)
        latitude = l * DEG_TO_RAD

        array = [longitude, latitude, distance * self.EARTH_RADIUS / self.AU, math.atan(self.BODY.Moon.eqRadius / (distance * self.EARTH_RADIUS))]
		
        return array 
    
    #/**
	# * Compute the position of the body.
	# * @param pos Values for the ecliptic longitude, latitude, distance and so on from previous methods for the specific body.
	# * @param geocentric True to return geocentric position. Set this to false generally.
	# * @return The ephemeris object with the output position
	# */
    def doCalc(self, pos, geocentric):
        #// Correct for nutation in longitude and obliquity
        pos[0] = pos[0] + self.nutLon
        pos[1] = pos[1] + self.nutObl
        
        #// Ecliptic to equatorial coordinates
        cl = math.cos(pos[1])
        x = pos[2] * math.cos(pos[0]) * cl
        y = pos[2] * math.sin(pos[0]) * cl
        z = pos[2] * math.sin(pos[1])
        sinEcl = math.sin(self.meanObliquity)
        cosEcl = math.cos(self.meanObliquity)
        tmp = y * cosEcl - z * sinEcl
        z = y * sinEcl + z * cosEcl
        y = tmp

        #// Obtain topocentric rectangular coordinates
        sinLat = math.sin(self.obsLat)
        cosLat = math.cos(self.obsLat)
        if (geocentric==True):
            radiusAU =  0 
        else:
            radiusAU = self.EARTH_RADIUS / self.AU

        correction = np.array([radiusAU * cosLat * math.cos(self.lst), 
                radiusAU * cosLat * math.sin(self.lst), radiusAU * sinLat])
        
            
        xtopo = x - correction[0]
        ytopo = y - correction[1]
        ztopo = z - correction[2]

        # Obtain topocentric equatorial coordinates
        ra = 0.0
        dec = self.PI_OVER_TWO
        if (ztopo < 0.0): 
            dec = - (dec)
        if (ytopo != 0.0 or xtopo != 0.0):
            ra = math.atan2(ytopo, xtopo)
            dec = math.atan2(ztopo / math.hypot(xtopo, ytopo), 1.0)
        dist = math.sqrt(xtopo * xtopo + ytopo * ytopo + ztopo * ztopo)

        # Hour angle
        angh = self.lst - ra
        
        # Obtain azimuth and geometric alt
        sinDec = math.sin(dec) 
        cosDec = math.cos(dec)
        h = sinLat * sinDec + cosLat * cosDec * math.cos(angh)
        alt = math.asin(h)
        azy = math.sin(angh)
        azx = math.cos(angh) * sinLat - sinDec * cosLat / cosDec
        azi = math.pi + math.atan2(azy, azx) #// 0 = north
		
        if (geocentric==True): 
            return self.Ephemeris(azi, alt, -1, -1, -1, -1, normalizeRadians(ra), dec, dist, pos[0], pos[1], pos[3])

		# Get apparent elevation
        alt = self.refraction(alt)

        if (self.twilight.value==self.TWILIGHT.HORIZON_34arcmin.value):
            tmp = -(34.0 / 60.0) * DEG_TO_RAD - pos[3]

        if (self.twilight.value==self.TWILIGHT.CIVIL.value):
            tmp = -6 * DEG_TO_RAD
        
        if (self.twilight.value==self.TWILIGHT.NAUTICAL.value):
            tmp = -12 * DEG_TO_RAD

        if (self.twilight.value==self.TWILIGHT.ASTRONOMICAL.value):
            tmp = -18 * DEG_TO_RAD
        
        # // Compute cosine of hour angle
        tmp = (math.sin(tmp) - sinLat * sinDec) / (cosLat * cosDec)
		# /** Length of a sidereal day in days according to IERS Conventions. */
        siderealDayLength = 1.00273781191135448
        celestialHoursToEarthTime = 1.0 / (siderealDayLength * self.TWO_PI)

        # // Make calculations for the meridian
        transit_time1 = celestialHoursToEarthTime * self.normalizeRadians(ra - self.lst)
        transit_time2 = celestialHoursToEarthTime * (self.normalizeRadians(ra - self.lst) - self.TWO_PI)
        transit_alt = math.asin(sinDec * sinLat + cosDec * cosLat)
        transit_alt = self.refraction(transit_alt)

        # // Obtain the current event in time
        transit_time = transit_time1
        jdToday = math.floor(self.jd_UT - 0.5) + 0.5
        transitToday2 = math.floor(self.jd_UT + transit_time2 - 0.5) + 0.5
        # // Obtain the transit time. Preference should be given to the closest event
        # // in time to the current calculation time
        if (jdToday == transitToday2 and abs(transit_time2) < abs(transit_time1)):
            transit_time = transit_time2
		
        transit = self.jd_UT + transit_time

        # // Make calculations for rise and set
        rise = -1
        set = -1
        if (abs(tmp) <= 1.0):
            ang_hor = abs(math.acos(tmp))
            rise_time1 = celestialHoursToEarthTime * self.normalizeRadians(ra - ang_hor - self.lst)
            set_time1 = celestialHoursToEarthTime * self.normalizeRadians(ra + ang_hor - self.lst)
            rise_time2 = celestialHoursToEarthTime * (self.normalizeRadians(ra - ang_hor - self.lst) - self.TWO_PI)
            set_time2 = celestialHoursToEarthTime * (self.normalizeRadians(ra + ang_hor - self.lst) - self.TWO_PI)

			# // Obtain the current events in time. Preference should be given to the closest event
			# // in time to the current calculation time (so that iteration in other method will converge)
            rise_time = rise_time1
            riseToday2 = math.floor(self.jd_UT + rise_time2 - 0.5) + 0.5
            if (jdToday == riseToday2 and abs(rise_time2) < abs(rise_time1)):
                rise_time = rise_time2
            set_time = set_time1
            setToday2 = math.floor(self.jd_UT + set_time2 - 0.5) + 0.5
            if (jdToday == setToday2 and abs(set_time2) < abs(set_time1)):
                set_time = set_time2
            rise = self.jd_UT + rise_time
            set = self.jd_UT + set_time

        out = self.Ephemeris(azi, alt, rise, set, transit, transit_alt, self.normalizeRadians(ra), dec, dist, pos[0], pos[1], pos[3])
        return out
        
    #/**
	# * Corrects geometric elevation for refraction if it is greater than -3 degrees.
	# * @param alt Geometric elevation in radians.
	# * @return Apparent elevation.
	# */
    def refraction(self,alt):
        if (alt <= -3 * DEG_TO_RAD):
            return alt
        Ps = 1010 # Pressure in mb
        Ts = 10 + 273.15; # Temperature in K
        altDeg = alt * RAD_TO_DEG

        # Bennet 1982 formulae for optical wavelengths, do the job but not accurate close to horizon
        r = DEG_TO_RAD * abs(math.tan(self.PI_OVER_TWO - (altDeg + 7.31 / (altDeg + 4.4)) * DEG_TO_RAD)) / 60.0
        refr = r * (0.28 * Ps / Ts)
        return min(alt + refr, self.PI_OVER_TWO)
    
    #/**
	# * Sets the illumination phase field for the provided body. 
	# * Sun position must be computed before calling this method.
	# * @param body The ephemeris object for this body.
	# */
    def getIlluminationPhase(self,body):
        dlon = body.rightAscension - self.sun.rightAscension
        cosElong = (math.sin(self.sun.declination) * math.sin(body.declination) + math.cos(self.sun.declination) * math.cos(body.declination) * math.cos(dlon))
        
        RE = self.sun.distance
        RO = body.distance
        #Use elongation cosine as trick to solve the rectangle and get RP (distance body - sun)
        RP = math.sqrt(-(cosElong * 2.0 * RE * RO - RE * RE - RO * RO))

        DPH = ((RP * RP + RO * RO - RE * RE) / (2.0 * RP * RO))
        body.illuminationPhase = 100 * (1.0 + DPH) * 0.5

    #/**
	# * Transforms a common date into a Julian day number (counting days from Jan 1, 4713 B.C. at noon).
	# * Dates before October, 15, 1582 are assumed to be in the Julian calendar, after that the Gregorian one is used.
	# * @param year Year.
	# * @param month Month.
	# * @param day Day.
	# * @param h Hour.
	# * @param m Minute.
	# * @param s Second.
	# * @return Julian day number.
	# * @throws Exception For an inexistent date.
	# */
    def toJulianDay(self,year,month,day,h,m,s):
        
        julian = False
        jd = 0
    
        if (year < 1582 or (year == 1582 and month <= 10) or (year == 1582 and month == 10 and day < 15)):
            julian=True
        D = day
        M = month
        Y = year
        if (M<3):
            Y = Y - 1
            M += 12
        A = int(Y / 100)
        if (julian==True):
            B = 0
        else:
            B = 2 - A + A / 4
        B = int(B)
        dayFraction = (h + (m + (s / 60.0)) / 60.0) / 24.0
        jd = dayFraction + (int) (365.25 * (Y + 4716)) + (int) (30.6001 * (M + 1)) + D + B - 1524.5

        if (jd < 2299160.0 and jd >= 2299150.0):
            print("invalid julian day " + jd + ". This date does not exist.")

        return jd
    
    #/**
	# * Transforms a Julian day (rise/set/transit fields) to a common date.
	# * @param jd The Julian day.
	# * @return A set of integers: year, month, day, hour, minute, second.
	# * @throws Exception If the input date does not exists.
	# */
    def getDate(self,jd):
        if (jd < 2299160.0 and jd >= 2299150.0):
            print("invalid julian day " + jd + ". This date does not exist.")
        
        #// The conversion formulas are from Meeus,
        # Chapter 7
        Z = math.floor(jd + 0.5)
        F = jd + 0.5 - Z
        A = Z
        if (Z >= 2299161):
            a = (int) ((Z - 1867216.25) / 36524.25)
            A += 1 + a - a / 4
        B = A + 1524
        C = (int) ((B - 122.1) / 365.25)
        D = (int) (C * 365.25)
        E = (int) ((B - D) / 30.6001)
        exactDay = F + B - D - (int) (30.6001 * E)
        day = int(exactDay)
        if (E<14):
            month = E-1
        else:
            month = E-13
        year = C - 4715
        if (month > 2):
            year=year-1
        h = ((exactDay - day) * self.SECONDS_PER_DAY) / 3600.0

        hour = int(h)
        m = (h - hour) * 60.0
        minute = int(m)
        second = (int) ((m - minute) * 60.0)
        
        array = np.array([year, month, day, hour, minute, second])

        return  array

    #/**
	# * Returns a date as a string.
	# * @param jd The Julian day.
	# * @return The String.
	# * @throws Exception If the date does not exists.
	# */
    def getDateAsString(self,jd):
        if (jd == -1): 
            return "NO RISE/SET/TRANSIT FOR THIS OBSERVER/DATE"
        date = np.array(self.getDate(jd))
        zyr = ""
        zmo = ""
        zh = ""
        zm = ""
        zs = ""

        if (date[1] < 10): zyr = "0"
        if (date[2] < 10): zmo = "0"
        if (date[3] < 10): zh = "0"
        if (date[4] < 10): zm = "0"
        if (date[5] < 10): zs = "0"

        return str(date[0])+"-"+zyr+str(date[1])+"-"+zmo+str(date[2])+" "+zh+str(date[3])+":"+zm+str(date[4])+":"+zs+str(date[5])+" UT"

    #/**
	# * Reduce an angle in radians to the range (0 - 2 Pi).
	# * @param r Value in radians.
	# * @return The reduced radians value.
	# */
    def normalizeRadians(self,r):
	
        if (r < 0 and r >= -self.TWO_PI): 
            return r + self.TWO_PI
        if (r >= self.TWO_PI and r < 2*self.TWO_PI):
            return r - self.TWO_PI
        if (r >= 0 and r < self.TWO_PI): 
            return r
        
        r -= self.TWO_PI * math.floor(r / self.TWO_PI)
        if (r < 0.): 
            r += self.TWO_PI
    
        return(r)
      
    #/**
	# * Computes an accurate rise/set/transit time for a moving object.
	# * @param riseSetJD Start date for the event.
	# * @param index Event identifier.
	# * @param niter Maximum number of iterations.
	# * @param sun True for the Sun.
	# * @return The Julian day in UT for the event, 1s accuracy.
	# */
    def obtainAccurateRiseSetTransit(self,riseSetJD, index, niter, sun):
        step = -1
        i=0
        while (i < niter):
            if (riseSetJD == -1): 
                return riseSetJD #// -1 means no rise/set from that location
            self.setUTDate(riseSetJD)
            out = self.Ephemeris
            out = None
            if (sun):
                out = self.doCalc(self.getSun(),False)
            else:
                out = self.doCalc(self.getMoon(),False)
            val = out.rise
            if (index == self.EVENT.SET): 
                val = out.set
            if (index == self.EVENT.TRANSIT): 
                val = out.transit
            step = abs(riseSetJD - val)
            riseSetJD = val
            if (step <= 1.0 / self.SECONDS_PER_DAY):
                break # // convergency reached
            i=i+1
        if (step > 1.0 / self.SECONDS_PER_DAY): 
            return -1 # // did not converge => without rise/set/transit in this date
        return riseSetJD

    #/**
    #* Returns the instant of a given moon phase.
    #* @param phase The phase.
    #* @return The instant of that phase, accuracy around 1 minute or better.
    #*/
    def getMoonPhaseTime(self, phase):
        accuracy = 10 / (30 * self.SECONDS_PER_DAY) #// 10s / lunar cycle length in s => 10s accuracy
        refPhase = phase
        oldJD = self.jd_UT 
        oldMoonAge = self.moonAge
        while (True):
            age = self.normalizeRadians((self.getMoon()[0] - self.getSun()[0])) / self.TWO_PI - refPhase
            if (age < 0): age += 1
            if (age < accuracy or age > 1 - accuracy): break
            if (age < 0.5):
                self.jd_UT -= age
            else:
                self.jd_UT += 1-age
            
            self.setUTDate(self.jd_UT)

        out = self.jd_UT
        self.setUTDate(oldJD)
        moonAge = oldMoonAge
        
        return out


############################################## MAIN PROGRAM ####################################################################

try:
    #Time in UT !!!      
    year = 2021
    month = 3
    day = 1
    h = 21
    m = 0
    s = 0

    #ADD LONGITUDE-LATITUDE
    obsLon = math.radians(LONGITUDE) #lon is negative to the west.
    obsLat = math.radians(LATITUDE)


    smc = SunMoonCalculator(year, month, day, h, m, s, obsLon, obsLat)
    smc.calcSunAndMoon()
    degSymbol = "\u00b0"
    print("###SUN###")
    print(" Az:       ",(math.degrees(smc.sun.azimuth)),degSymbol)
    print(" El:       ",(math.degrees(smc.sun.elevation)),degSymbol)
    print(" Dist:     ",(smc.sun.distance), " AU")
    print(" RA:       ",(math.degrees(smc.sun.rightAscension)),degSymbol)
    print(" DEC:      ",(math.degrees(smc.sun.declination)),degSymbol)
    print(" Ill:      ",(smc.sun.illuminationPhase),  "%")
    print(" Rise:     "+smc.getDateAsString(smc.sun.rise))
    print(" Set:      "+smc.getDateAsString(smc.sun.set))
    print(" Transit:  "+smc.getDateAsString(smc.sun.transit)," (elev. ", (math.degrees(smc.sun.transitElevation)),degSymbol+")")

    print("###MOON###")
    print(" Az:       ",(math.degrees(smc.moon.azimuth)),degSymbol)
    print(" El:       ",(math.degrees(smc.moon.elevation)),degSymbol)
    print(" Dist:     ",int(smc.moon.distance * smc.AU),"km")
    print(" RA:       ",(math.degrees(smc.moon.rightAscension)),degSymbol)
    print(" DEC:      ",(math.degrees(smc.moon.declination)),degSymbol)
    print(" Ill:      ",(smc.moon.illuminationPhase),  "%")
    print(" Age:      ",(smc.moonAge)," days")
    print(" Rise:     "+smc.getDateAsString(smc.moon.rise))
    print(" Set:      "+smc.getDateAsString(smc.moon.set))
    print(" Transit:  "+smc.getDateAsString(smc.moon.transit)," (elev. ", (math.degrees(smc.moon.transitElevation)),degSymbol+")")

    smc.setTwilight(SunMoonCalculator.TWILIGHT.ASTRONOMICAL)
    smc.calcSunAndMoon()
 
    print("")
    print("Astronomical twilights:")
    print("Sun")
    print(" Rise:     "+smc.getDateAsString(smc.sun.rise))
    print(" Set:      "+smc.getDateAsString(smc.sun.set))
    print("Moon")
    print(" Rise:     "+smc.getDateAsString(smc.moon.rise))
    print(" Set:      "+smc.getDateAsString(smc.moon.set))

    print("")
    print("Closest Moon phases:")
    for s in smc.MOONPHASE:
        print(" "+s.value[0]+"  "+smc.getDateAsString(smc.getMoonPhaseTime(s.value[1])))


except:
    print("Unexpected error:", sys.exc_info()[0])
    raise


	