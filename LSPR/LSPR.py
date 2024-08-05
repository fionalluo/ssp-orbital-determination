# LSPR.py
# Fiona Luo, NMT2
# 07-03-2020

'''
Description of program:
Input:
- (x, y) centroids, (RA, dec) of 12 stars (INPUT FILE)
- (x, y) centroid of object (COMMAND LINE)

Output:
- 6 plate constants (degrees, degrees/pix)
- Uncertainty of fit in RA and dec (arsec)
- (RA, dec) of unknown object (hh:mm:ss.ss, dd:mm:ss.s)
'''

import math
import numpy as np

print("\n\n")

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

#---------------------------------------------------------------------------------------------

# Read in the user's input using command line
# (x, y) centroid of UNKNOWN OBJECT
filename = input("test input file = ")
testposition = input("test position (x,y) = ")

xn = float( testposition[1:testposition.index(",")] ) # x coord of new object
yn = float( testposition[testposition.index(",") + 1:len(testposition)-1] ) # y coord of new object

# Variables for the text file inputs of 12 stars
x = []
y = []
RA = []  # capital bc stars have fixed values; lowercase ra used for variable in later code
DEC = [] 

# Read in (x,y) and (RA, dec) of 12 stars in the INPUT FILE
file1 = open(filename, 'r') 
lines = file1.readlines()
for line in lines: 
    line = line.split()
    x.append(float(line[0]))
    y.append(float(line[1]))
    raraw = line[2]
    decraw = line[3]
    RA.append( sextodec(raraw, True) )
    DEC.append( sextodec(decraw, False))

# Matrix manipulations to solve for 6 plate values
# Calculate SUMMATIONS first
N = 12
Sx = 0
Sy = 0
Sx2 = 0
Sy2 = 0
Sxy = 0
Sra = 0
Srax = 0
Sray = 0
Sdec = 0
Sdecx = 0
Sdecy = 0

for i in range(N):
    Sx += x[i]
    Sy += y[i]
    Sx2 += (x[i])**2
    Sy2 += (y[i])**2
    Sxy += x[i] * y[i]
    Sra += RA[i]
    Srax += RA[i] * x[i]
    Sray += RA[i] * y[i]
    Sdec += DEC[i]
    Sdecx += DEC[i] * x[i]
    Sdecy += DEC[i] * y[i]

# Define 3x3 matrix from equation and invert it
m = np.array([[N,Sx,Sy], [Sx,Sx2,Sxy],[Sy,Sxy,Sy2]])
m = np.linalg.inv(m)
mra = np.array([Sra, Srax, Sray])
mdec = np.array([Sdec, Sdecx, Sdecy])

m1 = np.dot(m, mra)
m2 = np.dot(m, mdec)

b1 = m1[0]
a11 = m1[1]
a12 = m1[2]
b2 = m2[0]
a21 = m2[1]
a22 = m2[2]

# Find the error for RA and DEC (era, edec)
era = 0
edec = 0
for i in range(N):
    era += ( RA[i] - b1 - a11*x[i] - a12*y[i] ) ** 2
    edec += ( DEC[i] - b2 - a21*x[i] - a22*y[i]) **2
era = ( era/(N-3) ) ** 0.5
edec = ( edec/(N-3) ) ** 0.5
# convert era and edec from RADIANS to ARCSEC
era = era * 180 / math.pi * 3600
edec = edec * 180 / math.pi * 3600

# Use equations to find final ra and dec of unknown object
ra = b1 + a11 * xn + a12 * yn
dec = b2 + a21 * xn + a22 * yn
ra = ra * 180 / math.pi # convert to degrees
dec = dec * 180 / math.pi
ra = dectosex(ra, True)
dec = dectosex(dec, False)

# Convert everything to degrees
b1 = math.degrees(b1)
a11 = math.degrees(a11)
a12 = math.degrees(a12)
b2 = math.degrees(b2)
a21 = math.degrees(a21)
a22 = math.degrees(a22)


# Print out results with formatting
print("***************")
print("plate constants")
print("***************")
print("b1:", "{0:.12g}".format(b1), "deg")
print("b2:", "{0:.12g}".format(b2), "deg")
print("a11:", "{0:.12g}".format(a11), "deg/pix")
print("a12:", "{0:.12g}".format(a12), "deg/pix")
print("a21:", "{0:.12g}".format(a21), "deg/pix")
print("a22:", "{0:.12g}".format(a22), "deg/pix")
print("***********")
print("uncertainty")
print("***********")
print(" RA:", '%.2f' % era, "arsec")
print("Dec:", '%.2f' % edec, "arsec")
print("*********************************")
print("astrometry for")
print("(x,y)=",testposition)
print("*********************************")
print(" RA: ", ra)
print("Dec:", dec)

print("\n\n")