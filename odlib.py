# odlib.py
# Fiona Luo, NMT 2
# 7/6/2020
# This library contains functions which will help in the OD of an asteroid

import numpy as np
import math
from math import radians, sin, cos, sqrt, pi, degrees, acos, asin
import ephem

# CONSTANTS
gauss = 58.132358929894
k = 0.01720209895  # Gaussian gravitational constant
G = 6.67428 * 10 ** -11  # Gravitational constant
Msun = 1.98892 * 10 ** 30
Mearth = 5.9722 * 10 ** 24
mu = G * (Msun + Mearth)
cAU = 173.144632674240  # speed of light in au/(mean solar)day 
eps = math.radians(23.4366)  # Earth's obliquity

# Returns an angle in the correct quadrant given its sine and cosine value (RADIANS)
def findQuadrant(sine, cosine):
    if cosine > 0 and sine > 0: #1
        return math.asin(sine)

    if cosine < 0 and sine > 0: #2
        return math.acos(cosine)

    if cosine < 0 and sine < 0: #3
        return math.pi - math.asin(sine)

    if cosine > 0 and sine < 0: #4
        return 2*math.pi + math.asin(sine)

# Returns PERCENTAGE of error between two values
def percenterror(value, truevalue):
    return abs( (truevalue - value) / truevalue ) * 100


# problem 2
# last 2 params for radian conversion, normalization
def convert_angle(degrees, minutes, seconds, radians):
    # handle a negative angle
    minutes = math.copysign(minutes, degrees)
    seconds = math.copysign(seconds, degrees)

    # perform angle conversion
    result = degrees + minutes/60 + seconds/3600
    result %= 360
    if radians:
        result *= math.pi/180
        
    # return result
    return result

# Takes input as hh:mm:ss.s or dd:mm:ss.s and converts to decimal RADIANS
# hours boolean if first term in hours
def sextodec(s, hours):
    x = 0
    h = float( s[0:s.index(":")] )
    m = float( s[s.index(":") + 1 : s.index(":") + 3])
    s = float( s[s.rfind(":") + 1:])

    m = math.copysign(m, h)
    s = math.copysign(s, h)

    x += h * math.pi / 180
    x += m * math.pi / 180 / 60
    x += s * math.pi / 180 / 60 / 60
    if (hours):
        x *= 15
    return x

# Convert a decimal DEGREE angle to HH:MM:SS.S OR DD:MM:SS.S
# (Seconds is ROUNDED)
def dectosex(n, hours):
    negative = False
    if (n < 0):
        negative = True
        n = -1 * n

    if(hours): 
        n /= 15  # there are 15 degrees in an hour
    d = int(n)
    n -= d
    m = int (n * 60)
    n -= m/60
    s = 0
    if (hours):
        s = int (n * 60 * 60 * 100) / 100
    else:
        s = int (n * 60 * 60 * 10) / 10

    # Add zeroes for single digit numbers
    if (d < 10):
        d = "0" + str(d)
    if (m < 10):
        m = "0" + str(m)
    if (s < 10):
        s = "0" + str(s)

    result = str(d) + ":" + str(m) + ":" + str(s)

    if hours == False:
        if negative:
            result = "-" + result
        else:
            result = "+" + result
    return result


# Function to get data dictionary (ddict) of specified time in data file formatted as:
# X     Y     Z
# VX    VY    VZ
# LT    RG    RR
# date: yyyy-Mon-dd
# time: hh:mm:ss.ssss
def getData (filename, date, time):
    # First read in the input line by line
    inputfile = open(filename)
    text = inputfile.readlines()

    # Initialize variables
    ddict = {}

    i = -1 # keep track of which type of line we're on; title, (x,y,z), (vx,vy,vz) etc
    correctday = False # if we're on the correct day this is true
    # In this loop, we use the direct indeces because they're the same for each line
    for line in text:
        # Check if we're on the correct day
        if date in line and time in line:
            i = 0
            correctday = True
        # read in positions
        if i == 1:
            ddict["x"] = float(line[4:26].strip() ) 
            ddict["y"] = float(line[30:52].strip() ) 
            ddict["z"] = float(line[56:78].strip())  
        # read in velocities
        if i == 2:
            # read/convert velocities to gaussian units
            ddict["vx"] = float(line[4:26].strip() ) / k 
            ddict["vy"] = float(line[30:52].strip() ) / k 
            ddict["vz"] = float(line[56:78].strip()) / k  
            break
        if correctday:
            i += 1
    return ddict

# function to get dictionary of orbit elements with copy pasted data from Horizons 
# (Make sure Ephemeris type is ELEMENTS)
def getElements(filename, date, time):
    # First read in the input line by line
    inputfile = open(filename)
    text = inputfile.readlines()

    # Element dictionary
    edict = {}

    i = -1 # keep track of which type of line we're on; title, (x,y,z), (vx,vy,vz) etc
    correctday = False # if we're on the correct day this is true
    # In this loop, we use the direct indeces because they're the same for each line
    for line in text:
        # Check if we're on the correct day
        if date in line and time in line:
            i = 0
            correctday = True
            edict["jd"] = float(line[0:18].strip() ) # Julian Date
        if i == 1:
            edict["e"] = float(line[4:26].strip() ) # eccentricity 
            edict["qr"] = float(line[30:52].strip() ) # periapsis distance (au)
            edict["i"] = float(line[56:78].strip())  # inclination (degrees)
        if i == 2:
            edict["Omega"] = float(line[4:26].strip() ) # longitude of ascending node (degrees)
            edict["w"] = float(line[30:52].strip() ) # argument of perihelion (degrees)
            edict["tp"] = float(line[56:78].strip() ) # time of periapsis (Julian Day Number)
        if i == 3:
            edict["n"] = float(line[4:26].strip() ) # mean motion (degrees/day)
            edict["ma"] = float(line[30:52].strip() ) # mean anomaly (degrees)
            edict["ta"] = float(line[56:78].strip() ) # true anomaly (degrees)
        if i == 4:
            edict["a"] = float(line[4:26].strip() ) # semi-major axis (au)
            edict["ad"] = float(line[30:52].strip() ) # apoapsis distance (au)
            edict["pr"] = float(line[56:78].strip() ) # sidereal orbit period (day)
            break
        if correctday:
            i += 1
    return edict


# r = position (x, y, z)
# v = velocity (vx, vy, vz)
def getangularmomentum (r, v):
    return np.cross(r, v)

# r: a numpy vector [x,y,z] for position
# v: a numpy vector [vx,vy,vz] for velocity
def semimajoraxis(r, v):
    term1 = 2/np.linalg.norm(r)
    term2 = np.dot(v, v.transpose())
    return 1/(term1 - term2)

# a = semimajor axis
def eccentricity(r, v):
    a = semimajoraxis(r, v)
    term = np.linalg.norm(np.cross(r, v)) ** 2
    ec = (1 - term/(a)) ** 0.5
    return ec

# i = inclination
def inclination(r, v, degrees = False):
    h = getangularmomentum(r, v)
    i = math.atan( (h[0]**2 + h[1]**2)**0.5 / h[2] )
    if degrees:
        return i * 180 / math.pi
    return i

# ongitude of ascending node = Omega
def lonascending(r, v, degrees = False):
    h = getangularmomentum(r, v)
    i = inclination(r, v)
    sine = h[0] / (np.linalg.norm(h) * math.sin(i))
    cosine = -1 * h[1] / (np.linalg.norm(h) * math.sin(i))
    angle = findQuadrant(sine, cosine)
    if degrees:
        return angle * 180 / pi
    return angle

def trueanomaly(r, v, degrees = False):
    a = semimajoraxis(r, v)
    e = eccentricity(r, v)
    h = getangularmomentum(r, v)
    cosinev = ( a * (1 - e**2) / np.linalg.norm(r) - 1) / e
    sinev = (a * (1 - e**2) / np.linalg.norm(h) * np.dot(r, v) / np.linalg.norm(r)) / e
    v = findQuadrant(sinev, cosinev)
    v = v % (math.pi * 2)
    if degrees:
        return v * 180 / math.pi
    return v


# omega
# r = position vector [x, y, z]
# i = inclination (RADIANS)
def argperihelion(r, rdot, degrees = False):
    Omega = lonascending(r, rdot, degrees = False)
    i = inclination(r, rdot, degrees = False)
    x = r[0]
    y = r[1]
    z = r[2]
    rmag = np.linalg.norm(r)
    e = eccentricity(r, rdot)
    a = semimajoraxis(r, rdot)
    hmag = np.linalg.norm(getangularmomentum(r, rdot))

    # eq 68 and 71 to find U
    cosU = (x*cos(Omega) + y*sin(Omega))/rmag
    sinU = z / (rmag * sin(i))
    U = findQuadrant(sinU, cosU)

    # eq 73 and 75 for v
    cosv = (a*(1-e**2)/rmag - 1)/e
    sinv = (a*(1-e**2)/hmag * np.dot(r, rdot)/rmag)/e
    v = findQuadrant(sinv, cosv)
    # eq 65 for w
    w = U - v 
    w %= 2*pi
    if(degrees):
        return w * 180 / pi
    return w


# ma ; mean anomaly is how much of the orbit the asteroid has passed through
def meananomaly(r, v, degrees = False):
    # Solve for E
    e = eccentricity(r, v)
    a = semimajoraxis(r, v)
    E = math.acos(1/e * (1 - np.linalg.norm(r)/a) )

    # Solve for mt2
    M = E - e * math.sin(E)
    M %= 2 * math.pi
    if degrees:
        return M * 180 / math.pi
    return M

# Find a NEW MEAN ANOMALY given an old jd and a new jd
def meananomalyjd(r, v, jdold, jdnew, k2 = 1, degrees = False): #, wrap180 = False):
    M = meananomaly(r, v)
    n = meanmotion(r, v, k1 = k2)
    Mnew = M + n * (jdnew - jdold)

    Mnew %= 2 * pi

    if degrees:
        Mnew = Mnew * 180 / pi
    return Mnew

# E = the initial E (eccentric anomaly) value
# M = The mean anomaly
# e = the eccentricity
# tol = the tolerance
# Get E, the eccentric anomaly, from M and e
def eccentricanomaly(E, M, e, tol):
    Mnew = E - e*sin(E)
    while abs(Mnew - M) > tol:
        Mnew = E - e*sin(E)
        #E = E - f(E, M, e) / (1 - e*cos(E))
        E = E - (E - e*sin(E) - M) / (1 - e*cos(E)) #replaced f(x)
        #print(Mnew)
    return E

# Returns the eccentric anomaly in RADIANS using r and v
def eccentricanomalyrv(r, v, degrees = False):
    a = semimajoraxis(r, v)
    e = eccentricity(r, v)
    E = acos((1 - np.linalg.norm(r)/a)/e)
    if degrees:
        return E * 180 / pi
    return E


# ma ; mean anomaly is how much of the orbit the asteroid has passed through
def meananomaly(r, v, degrees = False):
    # Solve for E
    e = eccentricity(r, v)
    a = semimajoraxis(r, v)
    E = math.acos(1/e * (1 - np.linalg.norm(r)/a) )
    if (np.dot(r, v) < 0):
        E = -1 * E

    # Solve for mt2
    M = E - e * math.sin(E)
    #M %= 2 * math.pi
    if degrees:
        return M * 180 / math.pi
    return M

# Find a NEW MEAN ANOMALY given an old jd and a new jd
def meananomalyjd(r, v, jdold, jdnew, k2 = 1, degrees = False): #, wrap180 = False):
    M = meananomaly(r, v)
    n = meanmotion(r, v, k1 = k2)
    Mnew = M + n * (jdnew - jdold)

    if degrees:
        Mnew = Mnew * 180 / pi
    return Mnew

# E = the initial E (eccentric anomaly) value
# M = The mean anomaly
# e = the eccentricity
# tol = the tolerance
# Get E, the eccentric anomaly, from M and e
def eccentricanomaly(E, M, e, tol):
    Mnew = E - e*sin(E)
    while abs(Mnew - M) > tol:
        Mnew = E - e*sin(E)
        #E = E - f(E, M, e) / (1 - e*cos(E))
        E = E - (E - e*sin(E) - M) / (1 - e*cos(E)) #replaced f(x)
        #print(Mnew)
    return E

# Returns the eccentric anomaly in RADIANS using r and v
def eccentricanomalyrv(r, v, degrees = False):
    a = semimajoraxis(r, v)
    e = eccentricity(r, v)
    E = acos((1 - np.linalg.norm(r)/a)/e)
    if np.dot(r, v) < 0:
        E = -1 * E
    if degrees:
        return E * 180 / pi
    return E


# jd
# date = "yyyy-Mon-dd"
        # 0123456789012
# time = "hh:mm:ss.ssss"
# source: https://quasar.as.utexas.edu/BillInfo/JulianDatesG.html
def juliandate(date, time):
    Y = int( date[0:4] )
    monthcodes = ["none","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    M = monthcodes.index(date[5:8])
    if M == 1 or M == 2:
        M += 12
        Y -= 1
    D = int( date[9:] )
    hours = float (time[0:2])
    minutes = float( time[3:5])
    seconds = float( time[6:])
    time = (hours + minutes/60 + seconds/3600)/24
    
    A = int(Y/100)
    B = int(A/4)    
    C = 2-A+B
    E = int(365.25 * (Y+4716))
    F = int(30.6001 * (M+1))
    JD= C+D+E+F-1524.5 + time

    return JD

# tp: time of periapsis
def timeofperiapsis(r, v, date, time):
    JD = juliandate(date, time)
    a = semimajoraxis(r, v)
    M = meananomaly(r, v)
    
    T = JD - M * (a**3 / k**2) ** 0.5
    return T

# tp: time of periapsis
def timeofperiapsisjd(r, v, JD):
    a = semimajoraxis(r, v)
    M = meananomaly(r, v)
    
    T = JD - M * (a**3 / k**2) ** 0.5
    return T

# Determine the initial value of r2 using the scalar equation of lagrange (SEL)
# taus = [tau1,tau3,tau]
# sun2 = [Sunx,Suny,Sunz]
# rhohat2 = [rhohatx,rhohat2y,rhohat2z]
# Ds = [D0,D21,D22,D23]
def SEL(taus, Sun2, rhohat2, Ds):
    roots = [0.,0.,0.] # for up to three real, positive roots 
    rhos = [0.,0.,0.] # range values for each real, positive root

    t1 = taus[0]
    t3 = taus[1]
    t = taus[2]

    A1 = t3/t
    B1 = A1/6 * (t**2 - t3**2)
    A3 = -1 * t1 / t
    B3 = A3/6 * (t**2 - t1**2)

    D0 = Ds[0] 
    D21 = Ds[1]
    D22 = Ds[2]
    D23 = Ds[3]

    A = ( A1 * D21 - D22 + A3 * D23 ) / (-1 * D0)
    B = ( B1 * D21 + B3 * D23) / (-1 * D0)

    R2 = np.array(Sun2)
    p2 = np.array(rhohat2)
    E = -2 * (np.dot(p2, R2))
    F = np.dot(R2, R2)

    a = -1 * (A**2 + A * E + F)
    b = -1 * (2*A*B + B*E)
    c = -1 * B**2

    r2 = np.roots([1,0,a,0,0,b,0,0,c])
    roots = []
    for num in r2:
        if np.isreal(num) == True and np.real(num) > 0 :
            roots.append(np.real(num))
    
    rhos = []
    for r in roots:
        rhos.append(A + B/(r**3) )
    
    roots.sort()
    rhos.sort()
    return(roots,rhos) # returns roots for SEL, and rhos (momentum)


# GENERIC NEWTON RALPHSON METHOD
# init = initial value
# f, fprime = functions for iterations
# tol = tolerance between subsequent value before value is returned
def newtonralphson(init, f, fprime, tol):
    x = init
    xnew = init - f(init) / fprime(init)
    while (abs(xnew - x) > tol):
        x = xnew
        xnew = xnew - f(xnew) / fprime(xnew)
    return xnew

# Return the mean motion / how many radians it travels each day
# k is the gaussian gravitational constant
# k = 0.01720209895 sometimes
def meanmotion(r, v, k1 = 1, degrees = False):
    a = semimajoraxis(r, v)
    n = (k1**2 / a**3) ** 0.5
    if degrees:
        return n * 180 / pi
    return n

# Period of orbit
# k2 is the gaussian gravitational constant, can be 0.01720209895
def period(r, v, k2 = 1):
    n = meanmotion(r, v, k1 = k2)
    return 2 * pi / n


# Reads input file to generate ephemeris in the format (LuoInput.txt)
# Julian Date, time (initial)
# X = , Y = , Z = (r on initial date)
# Orbital elements on correct day (Aug 3)
# Ra and Dec on correct day
# returns a dictionary?
def readinputephemeris(filename):
    # Element dictionary
    edict = {}
    # First read in the input line by line
    inputfile = open(filename)
    text = inputfile.readlines()

    # Read in initial r (x, y, z)

    i = 1 # keep track of which type of line we're on
    # In this loop, we use the direct indeces because they're the same for each line
    for line in text:
        if i == 1:
            edict["jd"] = float(line[0:18].strip() ) # Julian Date
        if i == 2:
            edict["x"] = float(line[4:26].strip() ) # x (of r from sun to asteroid)
            edict["y"] = float(line[30:52].strip() ) # y (of r from sun to asteroid)
            edict["z"] = float(line[56:78].strip() ) # z (of r from sun to asteroid)
        if i == 3:
            edict["e"] = float(line[4:26].strip() ) # eccentricity 
            edict["qr"] = float(line[30:52].strip() ) # periapsis distance (au)
            edict["i"] = float(line[56:78].strip())  # inclination (degrees)
        if i == 4:
            edict["Omega"] = float(line[4:26].strip() ) # longitude of ascending node (degrees)
            edict["w"] = float(line[30:52].strip() ) # argument of perihelion (degrees)
            edict["tp"] = float(line[56:78].strip() ) # time of periapsis (Julian Day Number)
        if i == 5:
            edict["n"] = float(line[4:26].strip() ) # mean motion (degrees/day)
            edict["ma"] = float(line[30:52].strip() ) # mean anomaly (degrees)
            edict["ta"] = float(line[56:78].strip() ) # true anomaly (degrees)
        if i == 6:
            edict["a"] = float(line[4:26].strip() ) # semi-major axis (au)
            edict["ad"] = float(line[30:52].strip() ) # apoapsis distance (au)
            edict["pr"] = float(line[56:78].strip() ) # sidereal orbit period (day)
        if i == 7: 
            n = line.split()
            edict["ra"] = ( float(n[0]) + float(n[1])/60 + float(n[2])/3600 ) * 15 
            edict["dec"] =  float(n[3]) + float(n[4])/60 + float(n[5])/3600 
            break
        i += 1
    return edict

# predict ra and dec for given date/time
# edict = dictionary of elements from the starting date; should have elements from readinputephemeris() method
# startdate = "yyyy-Mon-dd" (initial elements given)
# enddate =   "yyyy-Mon-dd" (ephemeris calculated for this date)
def ephemeris(edict, startdate, enddate, juliandate = True, degree = True): # Instructions in OD4_ephemeris.pdf
    # 6 orbital elements
    # r in orbital plane coords
    # r in ecliptic coords
    # r in equatorial coords
    # ra and dec

    # Get julian dates for start and end times
    t1 = startdate
    t2 = enddate
    if juliandate is False:
        t1 = juliandate(startdate, "00:00:00.0000")
        t2 = juliandate(enddate, "00:00:00.0000")

    # 6 orbital elements
    a = edict["a"] # semi major axis
    e = edict["e"] # eccentricity
    ma = radians(edict["M"]) # mean anomaly (sometimes i write "ma" instead)
    i = radians(edict["i"]) # inclination
    w = radians(edict["w"]) # argument of perihelion
    Omega = radians(edict["Omega"]) # longitude of ascending node
    n = (k**2/a**3) ** 0.5 # mean motion
    X = edict["x"] # sun vector
    Y = edict["y"]
    Z = edict["z"]

    # obtaining NEW M (Mean anomaly) and E (eccentric anomaly) at final time
    Mnew = n * (t2 - t1) + ma
    E = eccentricanomaly(Mnew, Mnew, e, 0.000001)

    # get r in orbital plane coords
    r = np.array([a*cos(E) - a*e, a * (1-e**2)**0.5 * sin(E), 0])

    # get r in ecliptic plane coords by doing 3 rotational matrix multiplications
    m1 = np.array([[cos(w), -sin(w), 0],[sin(w),cos(w),0],[0,0,1]])
    r = np.dot(m1, r)
    m2 = np.array([[1, 0, 0],[0, cos(i), -1*sin(i)],[0,sin(i),cos(i)]])
    r = np.dot(m2, r)
    m3 = np.array([[cos(Omega), -1*sin(Omega),0],[sin(Omega),cos(Omega),0],[0,0,1]])
    r = np.dot(m3, r)

    # get r in equatorial plane coordinates
    eps1 = radians( 23.4366 ) # obliquity of ecliptic plane
    #eps1 = radians (23.44)
    m4 = np.array([[1,0,0],[0,cos(eps1),-1*sin(eps1)],[0,sin(eps1),cos(eps1)]])
    r = np.dot(m4, r)

    # get ra and dec
    R = np.array([X, Y, Z]) # Sun vector from earth to sun (AU)
    p = np.add(R, r) # rho is vector from sun to asteroid (AU); add R and r from vector triangle
    phat = np.divide(p, float(np.linalg.norm(p)))

    dec = asin(phat[2])
    cosra = phat[0]/cos(dec)
    sinra = phat[1]/cos(dec)

    ra = findQuadrant(sinra, cosra)

    if degree:
        return degrees(ra), degrees(dec)
    return ra, dec


# Function used in eccentricanomaly2() in newton ralphson process
def f1(x, r2, r2dot, tau):
    a = semimajoraxis(r2, r2dot)
    n = meanmotion(r2, r2dot)
    return x - (1-np.linalg.norm(r2)/a)*sin(x) + \
        np.dot(r2, r2dot)/(n*a**2) * (1 - cos(x)) - n*tau

# Function used in eccentricanomaly2() in newton ralphson process
def fprime1(x, r2, r2dot, tau):
    a = semimajoraxis(r2, r2dot)
    n = meanmotion(r2, r2dot)
    return 1 - (1 - np.linalg.norm(r2)/a)*cos(x) + \
        np.dot(r2, r2dot) / (n * a**2) * sin(x)

# E = the initial E (eccentric anomaly) value
# M = The mean anomaly
# e = the eccentricity
# tol = the tolerance
# Get E, the eccentric anomaly, from M and e

# use newton's method to find dE
# newton's method: xnew = x - f(x)/f'(x)
def eccentricanomaly2(r2, r2dot, tau, tol):
    a = semimajoraxis(r2, r2dot)
    e = eccentricity(r2, r2dot)
    n = meanmotion(r2, r2dot)

    # initial value for delta E; depends on eccentricity
    init = n * tau
    # Add extra term to E if eccentricity is too high
    if e > 0.1: # (Info in Finding Delta E.pdf; sign of x0 depends on signtest)
        term = np.dot(r2, r2dot) / (n * a**2) # commonly ocurring term
        signtest = term * cos(n*tau - term) + \
            ( 1 - np.linalg.norm(r2)/a ) * sin(n*tau - term)
        if (signtest >= 0):
            init = init + 0.85 * e - term
        else:
            init = init - 0.85 * e - term
    
    # Use newton ralphson method with init value and f1/fprime1 function
    x = init
    xnew = x - f1(x, r2, r2dot, tau) / fprime1(x, r2, r2dot, tau)
    while (abs(xnew - x) >= tol):
        x = xnew
        xnew = xnew - f1(xnew, r2, r2dot, tau) / fprime1(xnew, r2, r2dot, tau)
    return xnew


# tau = time (gaussian i think? not sure)
# r2 = position vector
# r2dot = velocity vector
# flag = number (3, 4); or 0 for function directly
def fg(tau, r2, r2dot, flag):
    # cast position/velocity vectors to numpy arrays, initialize mu
    r2 = np.array(r2)
    r2dot = np.array(r2dot)
    mu = 1

    # Formula in OD2020 page 29
    # use r and rdot to find a and n
    # use newton's method to get delta E (start at n delta t)
    if flag == 2:
        u = mu / np.linalg.norm(r2) ** 3
        f = 1 - 1/2 * u * tau**2 
        g = tau
        return f, g
    if flag == 3 or flag == 4:
        u = mu / np.linalg.norm(r2) ** 3
        z = np.dot(r2, r2dot) / np.linalg.norm(r2) **2
        q = np.dot(r2dot, r2dot) / np.linalg.norm(r2) **2 - u
        # add the first 3 terms of the series
        f = 1 - 1/2 * u * tau**2 + 1/2 * u * z * tau**3
        g = tau - 1/6*u*tau**3 + 1/4 * u * z * tau**4

        # add the fourth term
        if (flag == 4):
            f += 1/24 * (3*u*q - 15*u*z**2 + u**2) * tau**4
            g += 1/4 * u * z * tau**4
        return f, g  # return, break out of method
    
    # if flag is 0 and we want the exact f and g functions
    tol = 10**-12
    n = meanmotion(r2, r2dot)
    a = semimajoraxis(r2, r2dot)
    dE = eccentricanomaly2(r2, r2dot, tau, tol)

    f = 1 - a/np.linalg.norm(r2) * (1 - cos(dE))
    g = tau + 1/n * (sin(dE) - dE)

    return f, g

# Return values for r1 and r3 using r2 and r2dot
# taus = [tau1, tau3, tau]
def r1r3(r2, r2dot, tau1, tau3, tau):
    # Convert input to numpy arrays
    r2 = np.array(r2)
    r2dot = np.array(r2dot)
    mu = 1

    u = mu / np.linalg.norm(r2) ** 3
    z = np.dot(r2, r2dot) / np.linalg.norm(r2) **2
    q = np.dot(r2dot, r2dot) / np.linalg.norm(r2) **2 - u

    r1 = (1 - u/(2*r2**3) + u*(np.dot(r2, r2dot))/(2*np.linalg.norm(r2)**5))*r2
    r3 = [0,0,0]
    print(r1, r3)


# Read in the od input file, returning values in an ARRAY of observation DICTIONARIES
def readodfile(filename):
    # Array (observations) of dictionaries (element values)
    odicts = []
    
    # First read in the input line by line
    inputfile = open(filename)
    text = inputfile.readlines()
    
    # for each line/observation, parse the values with indeces
    i = 0
    for line in text:
        values = line.split()
        odicts.append({})
        odicts[i]["jd"] = float(values[0]) # Julian Date of Observation (JD)
        odicts[i]["ra"] = sextodec(values[1], True) # Right ascension (hours->radians)
        odicts[i]["dec"] = sextodec(values[2].strip(), False) # Declination (degrees->radians)
        x = float(values[3]) # x of sun vector
        y = float(values[4]) # y of sun vector
        z = float(values[5]) # z of sun vector
        odicts[i]["R"] = np.array([x, y, z]) # Sun Vector R (AU)
        i += 1

    return odicts

def printfile(filename):
     # First read in the input line by line
    inputfile = open(filename)
    text = inputfile.readlines()
    for line in text:
        print(line.strip())

# return an array for rho components [x, y, z]
def getrhohat(ra, dec):
    x = cos(ra) * cos(dec)
    y = sin(ra) * cos(dec)
    z = sin(dec)
    return np.array([x, y, z])

# returns various orbital elements
# r2ec, r2dotec in ecliptic coordinates
# returns an element dictionary
# a, e, i, w, Omega, M, E, n, lastper, P
def orbitalelements(r2ec, r2dotec, jd):
    edict = {}
    edict["a"] = semimajoraxis(r2ec, r2dotec)
    edict["e"] = eccentricity(r2ec, r2dotec)
    edict["i"] = inclination(r2ec, r2dotec, degrees = True)
    edict["w"] = argperihelion(r2ec, r2dotec, degrees = True)
    edict["Omega"] = lonascending(r2ec, r2dotec, degrees = True)
    edict["M"] = meananomaly(r2ec, r2dotec, degrees = True)
    edict["E"] = eccentricanomalyrv(r2ec, r2dotec, degrees = True) 
    edict["n"] = meanmotion(r2ec, r2dotec, k1 = k)
    edict["lastper"] = timeofperiapsisjd(r2ec, r2dotec, jd) # Time of last perihelion passage
    edict["P"] = period(r2ec, r2dotec, k2 = k) / 365 # Period in days
    edict["jd"] = jd
    return edict

# a, e, i, w, Omega, M, E, n, lastper, P, jd
# prints orbital elements
def printorbitalelements(edict):
    a = edict["a"]
    e = edict["e"]
    i = edict["i"]
    w = edict["w"]
    Omega = edict["Omega"]
    E = edict["E"]
    M = edict["M"]
    n = edict["n"]
    lastper = edict["lastper"]
    P = edict["P"]
    jd = edict["jd"]
    print("\nORBITAL ELEMENTS")
    print("a =", round(a, 6), "AU")
    print("e =", round(e, 6))
    print("i =", round(i, 6), "degrees")
    print("omega =", round(w, 6), "degrees")
    print("Omega =", round(Omega, 6) , "degrees")
    print("E =", round(E, 6) , "degrees")
    print("M =", round(M, 6), "degrees at JD =", round(jd, 6))
    print("n =", round(n, 6), "rad/day")
    print("JD of last perihelion passage =", round(lastper, 6))
    print("P =", round(P, 6), "yrs")

def printuncertainties(avg_elements, std_elements, sdom_elements):
    #PRINT AVERAGES
    a = avg_elements["a"]
    e = avg_elements["e"]
    i = avg_elements["i"]
    w = avg_elements["w"]
    Omega = avg_elements["Omega"]
    E = avg_elements["E"]
    M = avg_elements["M"]
    n = avg_elements["n"]
    lastper = avg_elements["lastper"]
    P = avg_elements["P"]
    print("\nAVERAGE ORBITAL ELEMENTS")
    print("a =", a, "AU")
    print("e =", e)
    print("i =", i, "degrees")
    print("omega =", w, "degrees")
    print("Omega =", Omega , "degrees")
    print("E =", E, "degrees")
    print("M =", M)
    print("n =", n, "rad/day")
    print("JD of last perihelion passage =", lastper)
    print("P =", P, "yrs")

    # PRINT STANDARD DEVIATIONS
    a = std_elements["a"]
    e = std_elements["e"]
    i = std_elements["i"]
    w = std_elements["w"]
    Omega = std_elements["Omega"]
    E = std_elements["E"]
    M = std_elements["M"]
    n = std_elements["n"]
    lastper = std_elements["lastper"]
    P = std_elements["P"]
    print("\nSTANDARD DEVIATIONS")
    print("a =", a, "AU")
    print("e =", e)
    print("i =", i, "degrees")
    print("omega =", w, "degrees")
    print("Omega =", Omega , "degrees")
    print("E =", E, "degrees")
    print("M =", M)
    print("n =", n, "rad/day")
    print("JD of last perihelion passage =", lastper)
    print("P =", P, "yrs")

    # PRINT SDOM
    a = sdom_elements["a"]
    e = sdom_elements["e"]
    i = sdom_elements["i"]
    w = sdom_elements["w"]
    Omega = sdom_elements["Omega"]
    E = sdom_elements["E"]
    M = sdom_elements["M"]
    n = sdom_elements["n"]
    lastper = sdom_elements["lastper"]
    P = sdom_elements["P"]
    print("\nSDOM's")
    print("a =", a, "AU")
    print("e =", e)
    print("i =", i, "degrees")
    print("omega =", w, "degrees")
    print("Omega =", Omega , "degrees")
    print("E =", E, "degrees")
    print("M =", M)
    print("n =", n, "rad/day")
    print("JD of last perihelion passage =", lastper)
    print("P =", P, "yrs")

# Generate the ra and dec given the orbital elements of an asteroid
# date example: "2020/07/05 07:08:00"
def ephemerispyephem(edict, lon1, lat1, elev1, datetime):
    # Ephemeris Generator
    # Make an observer
    obs = ephem.Observer()
    '''
    obs.lon = "-81:24:52.9" #longitude of obs.
    obs.lat = "36:15:09.7" #lat of obs.
    obs.elev = 922. #elevation of obs, in meters (not a string!) 
    #obs.date = datetime
    obs.date = "2020/07/05 07:08:00" #(UTC date/time of observation)
    '''
    obs.lon = lon1
    obs.lat = lat1
    obs.elev = elev1
    obs.date = datetime
    
    a = edict["a"]
    e = edict["e"]
    i = edict["i"]
    w = edict["w"]
    w = 1.049186122879573E+02
    Omega = edict["Omega"]
    E = edict["E"]
    M = edict["M"]
    n = edict["n"]
    lastper = edict["lastper"]
    P = edict["P"]
    jd = edict["jd"]
    '''
    e = 1.305475148579241E-01
    i = 1.917247159927224E+01
    Omega = 1.761888605611246E+02
    w = 1.049186122879573E+02
    lastper =  2459041.857668430079
    n = 3.749324289784072E-01
    M = 1.093962878703165E+00
    a = 1.904732734778862E+00
    jd = edict["jd"]
    '''
    
    dateformat = datetime[5:7] + "/" + datetime[8:10] + "." + str(round(jd - int(jd), 5))[2:] + "/" + datetime[0:4]
    # make an xephem-style database entry for your asteroid 
    # It's a comma-delimited string with
    alist = [10737, "e", i, Omega, w, a, "",\
        e, M, dateformat, 2000,"",""]
    alist = [str(element) for element in alist]
    line = ",".join(alist)

    #create the asteroid
    asteroid = ephem.readdb(line)
    #compute its position for the observer 
    asteroid.compute(obs)
    return sextodec(str(asteroid.a_ra), hours = True), sextodec(str(asteroid.a_dec), hours = False)
