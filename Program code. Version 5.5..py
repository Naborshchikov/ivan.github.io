#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import numpy
from numpy import *
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib.pyplot import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy.optimize import curve_fit
from scipy.interpolate import PchipInterpolator
from scipy.signal import savgol_filter
from prettytable import PrettyTable
from collections import namedtuple
import csv

"""A table with the initial knowledge with which the calculation begins"""

Initial_conditions = [[1, 'Project implemented in Python', 'ver. 3.7.6 \n'],
                      [2, 'ID - Anaconda', 'ver. 2020 02 \n'],
                      [3, 'All data presented in the SI system', 'nist.gov/\n'],
                      [4, 'All constants are taken from the\n'
                       'data the US NIST.\n', 'nist.gov/ \n'],
                      [5, 'Particle structure -\n'                       
                       'published works of Nobel laureates.\n', 'published scientific works \n'],
                      [6, 'Protons, neutrons\n'
                       'have a core and two shells.\n', 'Robert Hofstadter \n'],
                      [7, 'Speed of light in a vacuum, c = 299792458\n', 'nist.gov/ \n'],
                      [8, 'Electrical constant, ε0 = 8.8541878128E−12\n', 'nist.gov/ \n'],
                      [9, 'Gravitational constant, G = 6.67430E-11\n', 'nist.gov/ \n'],
                      [10, 'Electric charge of an electron\n'
                       '-1.602176634e-19 \n', 'nist.gov/ \n'],
                      [11, 'π = 3.14159265358979', "Scientific American\n"],
                      [12, "Planck's constant, h = 6.62607015E−34", 'nist.gov/\n'],
                      [13, 'Electron diameter 10e−22,\n', 'Hans D. Dehmelt Experiments\n'],
                      [14, 'The proton consists of two quarks \n', 'Murray Gell-Mann\n'],
                      [15, 'The newneutron consists of two quarks \n', 'Murray Gell-Mann \n'],
                      [16, 'Quark radius − (0.47 · 10E−16 cm)E2\n'
                       '< RE2 < (0.43 · 10E−16 cm)E2 \n', 'arxiv.org/pdf/1604.01280.pdf \n'],
                      [17, 'Additional information\n', 'Data from available sources. \n'],
                     [18, "Quark condensate provides about 9\n"
                      "percent of the proton's mass\n", 'Physical Review Letters, 2018\n,'
                      ' website arXiv.org\n'],
                     [19, 'Electron diameter: 10e−22 \n', 'Nobel lecture, December, 8, 1989,\n'
                      ' Hans D. Dehmelt Experiments with \n'
                      'an isolated subatomic particle at rest\n'],
                     [20, 'proton mass: 1.67262192369E-27\n', 'nist.gov/\n'],
                     [21, 'neutron mass: 1.67492749804E-27\n', 'nist.gov/\n'],
                     [22, 'The magnitude of the charge\n'
                      'of the core, shells in the proton\n'
                     'respectively: 0.35; 0.5; 0.15\n', 'Robert Hofstadter the\n'
                      'Nobel laureate\n'], 
                     [23, 'The magnitude of the charge of the core,\n'
                      ' shells in the neutron\n'
                     'respectively: 0.35; - 0.5; 0.15\n', 'Robert Hofstadter the\n'],
                     [24, 'The proton radius: 0.84 fm\n', 'aps.org/publications/apsnews/201806/proton.cfm\n'],
                     [25, 'The neutron radius: 0.8e−15\n', 'Povh, B.; Rith, K.(2002).\n'],
                     [26, 'Rradius of the proton core: 0.23 ± 0.03 F\n', 'https://doi.org/10.1103/PhysRevD.18.2484\n'],
                     [27, 'Rradius of the neutron core:\n'
                      ' from 0.3 to 0.36 fm\n', 'arxiv.org/pdf/1810.00486.pdf\n'],
                     [28, 'The radius of the inner shell of the neutron\n'
                      'is approximately 0.6 fm.\n', 'actaphys.uj.edu.pl/fulltext?series=Reg&vol=30&page=119\n'],
                     [29, 'Proton, neutron consists of a nucleus and two\n'
                      'shells, or three quarks, ... \n', 'https://cerncourier.com/a/the-proton-laid-bare/ \n'
                     'https://www.nature.com/articles/s41586-019-0925-9'],
                     [30, 'This calculation starts with the fact that \n'
                      'a quark consists of a nucleus and two shells \n', 'This is a conditional division']] 
table1 = PrettyTable(['#', 'Description', 'Link to source/ comments'])
for rec in Initial_conditions:
    table1.add_row(rec)
    
"""Program class for obtaining data for analysis"""

class Algorithm():
# Assigning values to constants.
# Constants with more characters than constants according to US NIST data are index two.

    constantε0 = 8.8541878128e-12
    constantε02 = 8.85418781762039e-12
    
    constantc = 299792458
        
    constantg = 6.67430E-11
    constantg2 = 6.67448478E-11
    
    constanth = 6.62607015e-34     

# Data from different research groups may differ from each other.

    π = 3.14159265358979
    
# electron mass
    me = 9.1093837015e-31
    
# electron diameter
    de = 10e-22
    
# Electric charge of an electron
    qe = 1.602176634e-19
    qe2 = 1.602176620898e-19
    
# proton mass  
    mp = 1.67262192369E-27
# radius of a proton estimated by electric charge
    rp = 0.84e-15
# Rradius of the proton core
    rpc = 0.23e-15
# The radius of the inner layer (assumption).
    rpi = 0.6e-15 
# neutron mass
    mn = 1.67492749804E-27
# radius of a neutron
    rn = 0.8e-15
# Radius of the neutron core, 
# following from the hadronic and nuclear matter properties
    rnc = 0.33e-15
# The radius of the inner layer.
    rni = 0.6e-15
    
# quark radius
# The third sign "n" - for negative radius
# The third sign "p" - for positive radius
    qrn = - 0.47 * 10e-18
    qrp = 0.43 * 10e-18
    
# The magnitude of the charge of the core, shells respectively

# proton
    SHELLP0 = 0.35
    SHELLP1 = 0.5
    SHELLP2 = 0.15
    
# neutron
    SHELLN0 = 0.35
    SHELLN1 = -0.5
    SHELLN2 = 0.15
    
    def __init__ (self, xq02, xq13, xv02, xv13, xm02, xm13):
        
        
# The first symbol is the name of the quark, the second symbol is:
# 0 - core, 1 - inner shell, 2 - outer shell.

# RULE 1:
# The calculation takes into account that the quarks of the nucleus
# can not fall on a single line, as it will mean the synthesis of quarks
# and the loss of their identity.

# RULE 2:
# Quarks are connected if there is their intersection is at least one shell.

# RULE 3:
# The combination of quarks is obliged to provide the densest arrangement.

        self.xq02 = xq02
        self.xq13 = xq13
                
        self.xv02 = xv02
        self.xv13 = xv13
                
        self.xm02 = xm02
        self.xm13 = xm13
                
# Enter the number of different quarks in a particle.
n = 2
# Enter the total number of quarks in the particle.
k = 3

i = 1
aa = []
bb = []
cc = []
for i in range(n):
    quark_3n_2 = 3 * (i+1) - 2
    aa.append(quark_3n_2)
    quark_3n_1 = 3 * (i+1) - 1
    bb.append(quark_3n_1)
    quark_3n = 3 * (i+1) 
    cc.append(quark_3n)
A = np.array([aa, bb, cc])

arr = [k for k in range(k+1)]
def myfun(k):
    if k>n:
        return n
    else:
        return k
arr.pop(0)
list(map(myfun, arr))

def permutation(num_array):
    res=[]
    if len(num_array) <= 1:
        return [num_array]
    for num in set(num_array):
        temp_array = num_array.copy()
        temp_array.remove(num)
        res += [[num] + perm for perm in permutation(temp_array)]
    return res
my_list = list(map(myfun, arr)) 

C = []
for i in range(n):
    C += (permutation(my_list)[0][i], A[:, i].tolist())

if C[0] == 1:
    N1 = np.diag(C[1])
if C[2] == 2:
    N2 = np.diag(C[3])
O1 = np.zeros((k))
O2 = np.zeros((2 * k))
if np.array_equal(permutation(my_list)[0], np.array([1, 2, 2])) == True: 
    P1 = np.vstack((N1,O1)) 
    P2 = np.vstack((O1,N2))    
    P3 = np.hstack((P1,P2))
    P4 = np.vstack((P3,O2))    
    P5 = np.vstack((O1,P2))
    P = np.hstack((P4,P5))    
if np.array_equal(permutation(my_list)[1], np.array([2, 1, 2])) == True:    
    P21 = np.vstack((O1,N2))
    P51 = np.vstack((O1,P21))
    P7 = np.vstack((N2,O1)) 
    P8 = np.vstack((O1,N1))    
    P9 = np.hstack((P7,P8))    
    P10 = np.vstack((P9,O2))
    H = np.hstack((P10,P51))
if np.array_equal(permutation(my_list)[2], np.array([2, 2, 1])) == True: 
    P22 = np.vstack((O1,N2))
    P71 = np.vstack((N2,O1))
    P81 = np.vstack((O1,N1))
    P12 = np.hstack((P71,P22))    
    P13 = np.vstack((P12,O2))
    P14 = np.vstack((O1,P81))    
    L = np.hstack((P13,P14))
unique, counts = numpy.unique(P, return_counts=True)
x1 = dict(zip(unique, counts))
unique, counts = numpy.unique(H, return_counts=True)
x2 = dict(zip(unique, counts))
unique, counts = numpy.unique(L, return_counts=True)
x3 = dict(zip(unique, counts))
if x1 == x2 and x2 == x3:    
    print(f'')
else:
    print(f'Check the code')
    
K = np.split(P, [3])
K1 = np.split(K[1], [1])
unique, counts = numpy.unique(K[0], return_counts=True)
x4 = dict(zip(unique, counts))

unique, counts = numpy.unique(K1[0], return_counts=True)
x5 = dict(zip(unique, counts))

unique, counts = numpy.unique(K1[1], return_counts=True)
x6 = dict(zip(unique, counts))

H = np.split(H, [3])
H1 = np.split(H[1], [1])
unique, counts = numpy.unique(H[0], return_counts=True)
x7 = dict(zip(unique, counts))

unique, counts = numpy.unique(H1[0], return_counts=True)
x8 = dict(zip(unique, counts))

unique, counts = numpy.unique(H1[1], return_counts=True)
x9 = dict(zip(unique, counts))

L = np.split(L, [3])
L1 = np.split(L[1], [1])
unique, counts = numpy.unique(L[0], return_counts=True)
x10 = dict(zip(unique, counts))

unique, counts = numpy.unique(L1[0], return_counts=True)
x11 = dict(zip(unique, counts))

unique, counts = numpy.unique(L1[1], return_counts=True)
x12 = dict(zip(unique, counts))

# Let's compose the matrices for the neutron.
# Form for calculating quarks.
N1 = np.zeros((k, k * 2))
N2 = np.zeros((k, k * 2))
N3 = np.zeros((k, k * 2))

N1[0][0] = x4[1]
N1[0][1] = x4[2]
N1[0][2] = x4[3]
N1[0][3] = x4[4]
N1[0][4] = x4[5]
N1[1][4] = x5[5]
N1[1][5] = x5[6]
N1[2][5] = x6[6]

N2[0][0] = x7[1]
N2[0][1] = x7[2]
N2[0][3] = x7[4]
N2[0][4] = x7[5]
N2[0][5] = x7[6]
N2[1][2] = x8[3]
N2[1][4] = x8[5]
N2[2][5] = x9[6]

N3[0][0] = x10[1]
N3[0][3] = x10[4]
N3[0][4] = x10[5]
N3[0][5] = x10[6]
N3[1][1] = x11[2]
N3[1][5] = x11[6]
N3[2][2] = x12[3]

# Let's compose the matrices for the proton.
# Form for calculating quarks.
S = np.hsplit(N1, 2)
P1 = np.concatenate((S[1], S[0]), axis=1)
S1 = np.hsplit(N2, 2)
P2 = np.concatenate((S1[1], S1[0]), axis=1)
S3 = np.hsplit(N3, 2)
P3 = np.concatenate((S3[1], S3[0]), axis=1)

# Prepare matrices for calculating quarks.
Q1 = np.vstack((P1, N1))
Q2 = np.vstack((P1, N2))
Q3 = np.vstack((P1, N3))
Q4 = np.vstack((P2, N1))
Q5 = np.vstack((P2, N2))
Q6 = np.vstack((P2, N3))
Q7 = np.vstack((P3, N1))
Q8 = np.vstack((P3, N2))
Q9 = np.vstack((P3, N3))

# Let us determine the characteristics of the quarks "u" and "d" 
# for each matrix.

bq = array([Algorithm.SHELLP0, Algorithm.SHELLP1, Algorithm.SHELLP2, 
            Algorithm.SHELLN0, Algorithm.SHELLN1, Algorithm.SHELLN2])

try:    
    xq13 = linalg.solve(Q3, bq)
    xq16 = linalg.solve(Q6, bq)
    xq17 = linalg.solve(Q7, bq)
    xq18 = linalg.solve(Q8, bq)
    print('')
    B13 = np.linalg.inv(xq13)
    B16 = np.linalg.inv(xq16)
    B17 = np.linalg.inv(xq17)
    B18 = np.linalg.inv(xq18)
except:
    print('')
else:
    print("")
xq11 = linalg.solve(Q1, bq)
xq12 = linalg.solve(Q2, bq)
xq14 = linalg.solve(Q4, bq)
xq15 = linalg.solve(Q5, bq)
xq19 = linalg.solve(Q9, bq) 

# Error estimation.
Calcul = linalg.solve(Q5, bq)
Calcu = linalg.solve(Q9, bq)
Calculation_error_u2 = (Calcul[0] + Calcul[1] + Calcul[2] - 2/3)/(2/3) * 100
Calculation_error_d2 = (- 1/3 - Calcul[3] - Calcul[4] - Calcul[5])/(-1/3) * 100
Calculation_error_u = (2/3 - Calcu[0] - Calcu[1] - Calcu[2])/(2/3) * 100
Calculation_error_d = (- 1/3 - Calcu[3] - Calcu[4] - Calcu[5])/(1/3) * 100

# Comparison of matrices           
if np.array_equal(xq11, xq12) == True:
    xq12
if np.array_equal(xq12, xq14) == True:
    xq14
if np.array_equal(xq14, xq15) == True:
    xq15
    
# This formula is used at the INITIAL stage of calculations.
# The results obtained allow the interested person to obtain an 
#accurate volumetric model.
# V = 4/3πR**3
# In accordance with the existing representation
# v - volume; u - quark "u"; d - quark "d"; n - neutron; p - proton
# i - inner; o - outer
    
vu = 4/3 * Algorithm.π * (abs(Algorithm.qrn)**3 + Algorithm.qrp**3)/2
vd = vu
vn = 4/3 * Algorithm.π * Algorithm.rn**3
vp = 4/3 * Algorithm.π * Algorithm.rp**3
vrpc = 4/3 * Algorithm.π * Algorithm.rpc**3
vrnc = 4/3 * Algorithm.π * Algorithm.rnc**3
vrpi = 4/3 * Algorithm.π * Algorithm.rpi**3 - vrpc
vrpo = vp - vrpc - vrpi
vrni = 4/3 * Algorithm.π * Algorithm.rni**3 - vrnc
vrno = vn - vrnc - vrni
ve = 4/3 * Algorithm.π * (Algorithm.de/2)**3
mpc = Algorithm.mp * 0.91
mnc = Algorithm.mn * 0.91
mpi = Algorithm.mp * 0.09 * (vrpi/vrpo)
mni = Algorithm.mn * 0.09 * (vrni/vrno)
mpo = Algorithm.mp - mpc - mpi
mno = Algorithm.mn - mnc - mni
   
bv = array ([vrpc, vrpi, vrpo, vrnc, vrni, vrno])
    
bm = array ([mpc, mpi, mpo, mnc, mni, mno])

# Calculation of the charge in the electric charges of an electron for the
# core and shells of the "u" and "d" quarks.
# The numbers from [0] to [2] refer to the "u" quark.
# The numbers from [3] to [5] refer to the "d" quark.
xq02 = xq19
xq13 = xq15

# Calculation of volume for core and shells of the "u" and "d" quarks.
# The numbers from [0] to [2] refer to the "u" quark.
# The numbers from [3] to [5] refer to the "d" quark.
xv02 = linalg.solve(Q9, bv)
xv13 = linalg.solve(Q5, bv)

# Calculation of mass for core and shells of the "u" and "d" quarks.
# The numbers from [0] to [2] refer to the "u" quark.
# The numbers from [3] to [5] refer to the "d" quark.
xm02 = linalg.solve(Q9, bm)
xm13 = linalg.solve(Q5, bm)

# Calculation of the charge for the core and shells of the "u" and "d" 
#quarks.
# The numbers from [0] to [2] refer to the "u" quark.
# The numbers from [3] to [5] refer to the "d" quark.
for i, item in enumerate(xq02):
    xq02[i] *= Algorithm.qe
    
for i, item in enumerate(xq13):
    xq13[i] *= Algorithm.qe

Tuxq02 = xq02[0] + xq02[1] + xq02[2]
Tdxq02 = xq02[3] + xq02[4] + xq02[5]

Tuxq13 = xq13[0] + xq13[1] + xq13[2]
Tdxq13 = xq13[3] + xq13[4] + xq13[5]


unit = Algorithm(xq02, xq13, xv02, xv13, xm02, xm13)

Quarks = [[1, "Part of the matrices for calculating quarks \n"
           "xq13, xq16, xq17, xq18 \n", "Matrix inverse metrics do \n"  
           "not exist\n"],
          [2, "S to P conversions", "Mathematical conversion of \n"
           "neutron \n to proton"], 
          [3, "The proof of the algorithm: \n \n"
           "At the input, the magnitude of the \n" 
           "electric charge of the proton, \n" 
           "neutron in the shells.\n \n"
           "Output: The magnitude of \n"
           "the electric charge of quarks \n"
           "in the shells. The discrepancy with generally \n"
           "accepted experimental data: \n \n", "\n\n %"],
          [4, "Quark 'u' \n \n", {Calculation_error_u}],
          [5, "Quark 'd' \n", {Calculation_error_d}],
          [6, "Quark 'u2' \n \n", {Calculation_error_u2}],
          [7, "Quark 'd2' \n", {Calculation_error_d2}], 
          [8, "The electric charge of a quark: \n", "Coulomb"], 
          [9,"quark core 'u': \n \n", {unit.xq13[0]}], 
          [10,"inner quark shell 'u': \n\n", {unit.xq13[1]}], 
          [11,"outer quark shell 'u': \n \n", {unit.xq13[2]}], 
          [12,"total quark charge 'u': \n \n", {Tuxq13}], 
          [13,"quark core 'd': \n \n", {unit.xq13[3]}], 
          [14,"inner quark shell 'd': \n\n", {unit.xq13[4]}], 
          [15,"outer quark shell 'd': \n \n", {unit.xq13[5]}], 
          [12,"total quark charge 'd': \n \n", {Tdxq13}], 
          [16,"quark core 'u2': \n \n", {unit.xq02[0]}], 
          [17,"inner quark shell 'u2': \n\n", {unit.xq02[1]}], 
          [18,"outer quark shell 'u2': \n \n", {unit.xq02[2]}], 
          [19,"total quark charge 'u2': \n \n", {Tuxq02}], 
          [20,"quark core 'd2': \n \n", {unit.xq02[3]}], 
          [21,"inner quark shell 'd2': \n\n", {unit.xq02[4]}], 
          [22,"outer quark shell 'd2': \n \n", {unit.xq02[5]}], 
          [23,"total quark charge 'd2': \n \n", {Tdxq02}]]
                      
table2 = PrettyTable(["#", "Name", "Comments"])
for rec in Quarks:
    table2.add_row(rec)

"""Preparation of a data set for a protons, and neutrons"""
"""Calculating tachyon values"""

class Particles():
    def __init__ (self, proton0, proton1, neutron0, neutron1, tachyon_charge,
                 tachyon_mass, tachyon_volume):
        self.proton0 = proton0
        self.proton1 = proton1
        self.neutron0 = neutron0
        self.neutron1 = neutron1
        self.tachyon_charge = tachyon_charge
        self.tachyon_mass = tachyon_mass
        self.tachyon_volume = tachyon_volume        
        
# Matrices from a0 to a3 from the Algorithm class are used to form a data set
# for protons, and neutrons. 
# The x...02 values are used for the matrices a0 and a2.
# The x...13 values are used for the matrices a1 and a3.
Proton2 = namedtuple('Proton2', 'name1 charge name2 mass name3 volume')

proton0 = [[1, 'pq1 \n', unit.xq02[0], 'pm1', unit.xm02[0], 'pv1', unit.xv02[0]],
           [2, 'pq2 \n', unit.xq02[1], 'pm2', unit.xm02[1], 'pv2', unit.xv02[1]],
           [3, 'pq3 \n', unit.xq02[0], 'pm3', unit.xm02[0], 'pv3', unit.xv02[0]],
           [4, 'pq4 \n', unit.xq02[2], 'pm4', unit.xm02[2], 'pv4', unit.xv02[2]],
           [5, 'pq5 \n', unit.xq02[1], 'pm5', unit.xm02[1], 'pv5', unit.xv02[1]],           
           [6, 'pq6 \n', unit.xq02[3], 'pm6', unit.xm02[3], 'pv6', unit.xv02[3]],           
           [7, 'pq7 \n', unit.xq02[2], 'pm7', unit.xm02[2], 'pv7', unit.xv02[2]],           
           [8, 'pq8 \n', unit.xq02[4], 'pm8', unit.xm02[4], 'pv8', unit.xv02[4]],
           [9, 'pq9 \n', unit.xq02[5], 'pm9', unit.xm02[5], 'pv9', unit.xv02[5]]] 

table3 = PrettyTable(['#', 'Q symbol', 'Charge in Cl', 'M symbol',
                      'Mass in kg.', 'V symbol', 'Volume in cbm'])

for rec in proton0:
    table3.add_row(rec) 

Proton = namedtuple('Proton', 'name1 charge name2 mass name3 volume')
proton1 = [[1, 'pq1 \n', unit.xq13[0], 'pm1', unit.xm13[0], 'pv1', unit.xv13[0]], 
           [2, 'pq2 \n', unit.xq13[1], 'pm2', unit.xm13[1], 'pv2', unit.xv13[1]], 
           [3, 'pq3 \n', unit.xq13[3], 'pm3', unit.xm13[3], 'pv3', unit.xv13[3]],
           [4, 'pq4 \n', unit.xq13[2], 'pm4', unit.xm13[2], 'pv4', unit.xv13[2]],
           [5, 'pq5 \n', unit.xq13[4], 'pm5', unit.xm13[4], 'pv5', unit.xv13[4]],           
           [6, 'pq6 \n', unit.xq13[0], 'pm6', unit.xm13[0], 'pv6', unit.xv13[0]],
           [7, 'pq7 \n', unit.xq13[5], 'pm7', unit.xm13[5], 'pv7', unit.xv13[5]],           
           [8, 'pq8 \n', unit.xq13[1], 'pm8', unit.xm13[1], 'pv8', unit.xv13[1]],
           [9, 'pq9 \n', unit.xq13[2], 'pm9', unit.xm13[2], 'pv9', unit.xv13[2]]] 

table4 = PrettyTable(['#', 'Q symbol', 'Charge in Cl', 'M symbol',
                      'Mass in kg.', 'V symbol', 'Volume in cbm'])
for rec in proton1:
    table4.add_row(rec)

Neutron2 = namedtuple('Neutron2', 'name1 charge name2 mass name3 volume')
neutron0 = [[1, 'nq1 \n', unit.xq02[3], 'nm1', unit.xm02[3], 'nv1', unit.xv02[3]], 
            [2, 'nq2 \n', unit.xq02[4], 'nm2', unit.xm02[4], 'nv2', unit.xv02[4]],
            [3, 'nq3 \n', unit.xq02[3], 'nm3', unit.xm02[3], 'nv3', unit.xv02[3]],
            [4, 'nq4 \n', unit.xq02[5], 'nm4', unit.xm02[5], 'nv4', unit.xv02[5]],
            [5, 'nq5 \n', unit.xq02[4], 'nm5', unit.xm02[4], 'nv5', unit.xv02[4]],            
            [6, 'nq6 \n', unit.xq02[0], 'nm6', unit.xm02[0], 'nv6', unit.xv02[0]],
            [7, 'nq7 \n', unit.xq02[4], 'nm7', unit.xm02[4], 'nv7', unit.xv02[4]],            
            [8, 'nq8 \n', unit.xq02[1], 'nm8', unit.xm02[1], 'nv8', unit.xv02[1]],
            [9, 'nq9 \n', unit.xq02[2], 'nm9', unit.xm02[2], 'nv9', unit.xv02[2]]]

table5 = PrettyTable(['#', 'Q symbol', 'Charge in Cl', 'M symbol',
                      'Mass in kg.', 'V symbol', 'Volume in cbm'])
for rec in neutron0:
    table5.add_row(rec)

Neutron = namedtuple('Neutron', 'name1 charge name2 mass name3 volume')    
neutron1 = [[1, 'nq1 \n', unit.xq13[3], 'nm1', unit.xm13[3], 'nv1', unit.xv13[3]], 
            [2, 'nq2 \n', unit.xq13[4], 'nm2', unit.xm13[4], 'nv2', unit.xv13[4]], 
            [3, 'nq3 \n', unit.xq13[0], 'nm3', unit.xm13[0], 'nv3', unit.xv13[0]],
            [4, 'nq4 \n', unit.xq13[5], 'nm4', unit.xm13[5], 'nv4', unit.xv13[5]],
            [5, 'nq5 \n', unit.xq13[1], 'nm5', unit.xm13[1], 'nv5', unit.xv13[1]],
            [6, 'nq6 \n', unit.xq13[3], 'nm6', unit.xm13[3], 'nv6', unit.xv13[3]],            
            [7, 'nq7 \n', unit.xq13[2], 'nm7', unit.xm13[2], 'nv7', unit.xv13[2]],
            [8, 'nq8 \n', unit.xq13[4], 'nm8', unit.xm13[4], 'nv8', unit.xv13[4]],
            [9, 'nq9 \n', unit.xq13[5], 'nm9', unit.xm13[5], 'nv9', unit.xv13[5]]]

table6 = PrettyTable(['#', 'Q symbol', 'Charge in Cl', 'M symbol',
                      'Mass in kg.', 'V symbol', 'Volume in cbm'])
for rec in neutron1:
    table6.add_row(rec)

proton0 = list(zip(* proton0))
proton1 = list(zip(* proton1))
neutron0 = list(zip(* neutron0))
neutron1 = list(zip(* neutron1))

"""Algorithm for finding tachyon""" 

proton0_min_charge = min((proton0)[2], key=abs)
proton1_min_charge = min((proton1)[2], key=abs)
neutron0_min_charge = min((neutron0)[2], key=abs)
neutron1_min_charge = min((neutron1)[2], key=abs)

# Let's compare the minimum values of charges in protons, and neutrons, 
# and find the value of a tachyon

if (proton0_min_charge == neutron0_min_charge and  
    proton1_min_charge == neutron1_min_charge and 
    proton0_min_charge == proton1_min_charge):    

# ATTENTION! THE CYCLE DOES NOT CONTAIN A FORCED INTERRUPTION.
   
    a = proton0_min_charge
    b = unit.qe2
    while a != b:
        if a > b:
            a = a - b
        else:
            b = b - a

# Electric charge of the tachyon.
tachyon_charge = a
        
# Find the mass of a tachyon.
tachyon_mass = unit.me/(unit.qe2/tachyon_charge)
 
# Minimum_volume tachyon. 
tachyon_volume = ve/(unit.qe2/tachyon_charge)

# Let's define protons, and neutrons, through the definition
# of electric charge.

if (sum(neutron0[2]) > sum(neutron1[2]) and sum(proton0[2]) == sum(proton1[2])): 
    neutron = neutron1 
    new_neutron = neutron0 
    proton = proton1 
    new_proton = proton0
else:
    print('Algorithm requires verification.')     

unit1 = Particles(proton0, proton1, neutron0, neutron1, tachyon_charge,
                 tachyon_mass, tachyon_volume)

tachyon_electrical_charge = [[1, 'Electric charge \n', unit1.tachyon_charge],
                              [2, 'Mass \n', unit1.tachyon_mass],
                              [3, 'Volume', unit1.tachyon_volume]]
table7 = PrettyTable(['#', 'Description', 'Design data'])

for rec in tachyon_electrical_charge:
                              table7.add_row(rec)       

# We divide matrices of protons, neutrons by time segments
# NEUTRON
# We divide matrices of protons, neutrons by time segments
# NEUTRON
NEUTRON_Present_matrix = []
for i in range(9): 
    if unit1.neutron1[2][i] > 0 and unit1.neutron1[6][i] > 0:
        NEUTRON_Present_matrix += [[unit1.neutron1[2][i], unit1.neutron1[4][i], unit1.neutron1[6][i]]]    
    
NEUTRON_Past_matrix = []
for j in range(9): 
    if unit1.neutron1[2][j] < 0 and unit1.neutron1[6][j] < 0:
        NEUTRON_Past_matrix += [[unit1.neutron1[2][j], unit1.neutron1[4][j], unit1.neutron1[6][j]]]    
    
NEUTRON_Future_matrix = []
for j in range(9): 
    if unit1.neutron1[2][j] > 0 and unit1.neutron1[6][j] < 0:
        NEUTRON_Future_matrix += [[unit1.neutron1[2][j], unit1.neutron1[4][j], unit1.neutron1[6][j]]]    
    
# NEUTRON 2
NEUTRON2_Present_matrix = []
for i in range(9): 
    if unit1.neutron0[2][i] > 0 and unit1.neutron0[6][i] > 0:
        NEUTRON2_Present_matrix += [[unit1.neutron0[2][i], unit1.neutron0[4][i], unit1.neutron0[6][i]]]    
    
NEUTRON2_Past_matrix = []
for j in range(9): 
    if unit1.neutron0[2][j] < 0 and unit1.neutron0[6][j] < 0:
        NEUTRON2_Past_matrix += [[unit1.neutron0[2][j], unit1.neutron0[4][j], unit1.neutron0[6][j]]]    
    
NEUTRON2_Future_matrix = []
for j in range(9): 
    if unit1.neutron0[2][j] > 0 and unit1.neutron0[6][j] < 0:
        NEUTRON2_Future_matrix += [[unit1.neutron0[2][j], unit1.neutron0[4][j], unit1.neutron0[6][j]]]    
    
# PROTON
PROTON_Present_matrix = []
for i in range(9): 
    if unit1.proton1[2][i] > 0 and unit1.proton1[6][i] > 0:
        PROTON_Present_matrix += [[unit1.proton1[2][i], unit1.proton1[4][i], unit1.proton1[6][i]]]    
    
PROTON_Past_matrix = []
for j in range(9): 
    if unit1.proton1[2][j] < 0 and unit1.proton1[6][j] < 0:
        PROTON_Past_matrix += [[unit1.proton1[2][j], unit1.proton1[4][j], unit1.proton1[6][j]]]    
    
PROTON_Future_matrix = []
for j in range(9): 
    if unit1.proton1[2][j] > 0 and unit1.proton1[6][j] < 0:
        PROTON_Future_matrix += [[unit1.proton1[2][j], unit1.proton1[4][j], unit1.proton1[6][j]]]    
    
# PROTON 2
PROTON2_Present_matrix = []
for i in range(9): 
    if unit1.proton0[2][i] > 0 and unit1.proton0[6][i] > 0:
        PROTON2_Present_matrix += [[unit1.proton0[2][i], unit1.proton0[4][i], unit1.proton0[6][i]]]    
    
PROTON2_Past_matrix = []
for j in range(9): 
    if unit1.proton0[2][j] < 0 and unit1.proton0[6][j] < 0:
        PROTON2_Past_matrix += [[unit1.proton0[2][j], unit1.proton0[4][j], unit1.proton0[6][j]]]    
    
PROTON2_Future_matrix = []
for j in range(9): 
    if unit1.proton0[2][j] > 0 and unit1.proton0[6][j] < 0:
        PROTON2_Future_matrix += [[unit1.proton0[2][j], unit1.proton0[4][j], unit1.proton0[6][j]]]

# We find the value of each time segment
# Let's calculate the value for each column
NEUTRON_Present = numpy.sum(NEUTRON_Present_matrix, axis=0, dtype=None, out=None)
NEUTRON_Past = numpy.sum(NEUTRON_Past_matrix, axis=0, dtype=None, out=None)
NEUTRON_Future = numpy.sum(NEUTRON_Future_matrix, axis=0, dtype=None, out=None)
NEUTRON2_Present = numpy.sum(NEUTRON2_Present_matrix, axis=0, dtype=None, out=None)
NEUTRON2_Past = numpy.sum(NEUTRON2_Past_matrix, axis=0, dtype=None, out=None)
NEUTRON2_Future = numpy.sum(NEUTRON2_Future_matrix, axis=0, dtype=None, out=None)

PROTON_Present = numpy.sum(PROTON_Present_matrix, axis=0, dtype=None, out=None)
PROTON_Past = numpy.sum(PROTON_Past_matrix, axis=0, dtype=None, out=None)
PROTON_Future = numpy.sum(PROTON_Future_matrix, axis=0, dtype=None, out=None)
PROTON2_Present = numpy.sum(PROTON2_Present_matrix, axis=0, dtype=None, out=None)
PROTON2_Past = numpy.sum(PROTON2_Past_matrix, axis=0, dtype=None, out=None)
PROTON2_Future = numpy.sum(PROTON2_Future_matrix, axis=0, dtype=None, out=None)

# Justification of neutron decay and proton stability
singularly = np.array(NEUTRON_Past_matrix[0]) - np.array(NEUTRON_Past_matrix[1])

# Neutron 2 

NEUTRON2_Present_V = NEUTRON2_Present[2]
NEUTRON2_Past_V = NEUTRON2_Past[2]
NEUTRON2_Future_V = NEUTRON2_Future[2]
    
NEUTRON2_Present_Q = NEUTRON2_Present[0]
NEUTRON2_Past_Q = NEUTRON2_Past[0]
NEUTRON2_Future_Q = NEUTRON2_Future[0]
    
NEUTRON2_Present_M = NEUTRON2_Present[1]
NEUTRON2_Past_M = NEUTRON2_Past[1]
NEUTRON2_Future_M = NEUTRON2_Future[1]

# Proton 2
 
PROTON2_Present_V = PROTON2_Present[2]
PROTON2_Past_V = PROTON2_Past[2]
PROTON2_Future_V = PROTON2_Future[2]
    
PROTON2_Present_Q = PROTON2_Present[0]
PROTON2_Past_Q = PROTON2_Past[0]
PROTON2_Future_Q = PROTON2_Future[0]
    
PROTON2_Present_M = PROTON2_Present[1]
PROTON2_Past_M = PROTON2_Past[1]
PROTON2_Future_M = PROTON2_Future[1]

# Neutron 

NEUTRON_Present_V = NEUTRON_Present[2]
NEUTRON_Past_V = NEUTRON_Past[2]
NEUTRON_Future_V = NEUTRON_Future[2]
    
NEUTRON_Present_Q = NEUTRON_Present[0]
NEUTRON_Past_Q = NEUTRON_Past[0]
NEUTRON_Future_Q = NEUTRON_Future[0]
    
NEUTRON_Present_M = NEUTRON_Present[1]
NEUTRON_Past_M = NEUTRON_Past[1]
NEUTRON_Future_M = NEUTRON_Future[1]

# Proton
 
PROTON_Present_V = PROTON_Present[2]
PROTON_Past_V = PROTON_Past[2] 
PROTON_Future_V = PROTON_Future[2]
    
PROTON_Present_Q = PROTON_Present[0]
PROTON_Past_Q = PROTON_Past[0]
PROTON_Future_Q = PROTON_Future[0]
    
PROTON_Present_M = PROTON_Present[1]
PROTON_Past_M = PROTON_Past[1]
PROTON_Future_M = PROTON_Future[1]                      
                        
class Tachion():   
# The class makes calculations for protons, neutrons, 
# taking into account tachyons

# The distribution for protons, neutrons over time segments 

# Neutron     
    
    def __init__ (self, neutron_present_tachyon_quantity, neutron_future_tachyon_quantity, 
                  neutron_present_tachyon_mass, neutron_future_tachyon_mass,
                  neutron_present_withouttachyon_mass, neutron_future_withouttachyon_mass,
                  neutron_present_withouttachyon_energy, neutron_future_withouttachyon_energy,
                  neutron_present_at1tachyon_energy, neutron_future_at1tachyon_energy,
                  
                  proton_present_tachyon_quantity, proton_future_tachyon_quantity, 
                  proton_present_tachyon_mass, proton_future_tachyon_mass,
                  proton_present_withouttachyon_mass, proton_future_withouttachyon_mass,
                  proton_present_withouttachyon_energy, proton_future_withouttachyon_energy,
                  proton_present_at1tachyon_energy, proton_future_at1tachyon_energy,
                  
                  neutron2_present_tachyon_quantity, neutron2_future_tachyon_quantity, 
                  neutron2_present_tachyon_mass, neutron2_future_tachyon_mass,
                  neutron2_present_withouttachyon_mass, neutron2_future_withouttachyon_mass,
                  neutron2_present_withouttachyon_energy, neutron2_future_withouttachyon_energy,
                  neutron2_present_at1tachyon_energy, neutron2_future_at1tachyon_energy,
                  
                  proton2_present_tachyon_quantity, proton2_future_tachyon_quantity, 
                  proton2_present_tachyon_mass, proton2_future_tachyon_mass,
                  proton2_present_withouttachyon_mass, proton2_future_withouttachyon_mass,
                  proton2_present_withouttachyon_energy, proton2_future_withouttachyon_energy,
                  proton2_present_at1tachyon_energy, proton2_future_at1tachyon_energy):
        
        self.neutron_present_tachyon_quantity = neutron_present_tachyon_quantity
        self.neutron_future_tachyon_quantity = neutron_future_tachyon_quantity
        self.neutron_present_tachyon_mass = neutron_present_tachyon_mass
        self.neutron_future_tachyon_mass = neutron_future_tachyon_mass
        self.neutron_present_withouttachyon_mass = neutron_present_withouttachyon_mass
        self.neutron_future_withouttachyon_mass = neutron_future_withouttachyon_mass
        self.neutron_present_withouttachyon_energy = neutron_present_withouttachyon_energy
        self.neutron_future_withouttachyon_energy = neutron_future_withouttachyon_energy
        self.neutron_present_at1tachyon_energy = neutron_present_at1tachyon_energy
        self.neutron_future_at1tachyon_energy = neutron_future_at1tachyon_energy
        
        self.proton_present_tachyon_quantity = proton_present_tachyon_quantity
        self.proton_future_tachyon_quantity = proton_future_tachyon_quantity
        self.proton_present_tachyon_mass = proton_present_tachyon_mass
        self.proton_future_tachyon_mass = proton_future_tachyon_mass
        self.proton_present_withouttachyon_mass = proton_present_withouttachyon_mass
        self.proton_future_withouttachyon_mass = proton_future_withouttachyon_mass
        self.proton_present_withouttachyon_energy = proton_present_withouttachyon_energy
        self.proton_future_withouttachyon_energy = proton_future_withouttachyon_energy
        self.proton_present_at1tachyon_energy = proton_present_at1tachyon_energy
        self.proton_future_at1tachyon_energy = proton_future_at1tachyon_energy
        
        self.neutron2_present_tachyon_quantity = neutron2_present_tachyon_quantity
        self.neutron2_future_tachyon_quantity = neutron2_future_tachyon_quantity
        self.neutron2_present_tachyon_mass = neutron2_present_tachyon_mass
        self.neutron2_future_tachyon_mass = neutron2_future_tachyon_mass
        self.neutron2_present_withouttachyon_mass = neutron2_present_withouttachyon_mass
        self.neutron2_future_withouttachyon_mass = neutron2_future_withouttachyon_mass
        self.neutron2_present_withouttachyon_energy = neutron2_present_withouttachyon_energy
        self.neutron2_future_withouttachyon_energy = neutron2_future_withouttachyon_energy
        self.neutron2_present_at1tachyon_energy = neutron2_present_at1tachyon_energy
        self.neutron2_future_at1tachyon_energy = neutron2_future_at1tachyon_energy
        
        self.proton2_present_tachyon_quantity = proton2_present_tachyon_quantity
        self.proton2_future_tachyon_quantity = proton2_future_tachyon_quantity
        self.proton2_present_tachyon_mass = proton2_present_tachyon_mass
        self.proton2_future_tachyon_mass = proton2_future_tachyon_mass
        self.proton2_present_withouttachyon_mass = proton2_present_withouttachyon_mass
        self.proton2_future_withouttachyon_mass = proton2_future_withouttachyon_mass
        self.proton2_present_withouttachyon_energy = proton2_present_withouttachyon_energy
        self.proton2_future_withouttachyon_energy = proton2_future_withouttachyon_energy
        self.proton2_present_at1tachyon_energy = proton2_present_at1tachyon_energy
        self.proton2_future_at1tachyon_energy = proton2_future_at1tachyon_energy
        
# We calculate the possible number of tachyon
neutron_present_tachyon_quantity = NEUTRON_Present_Q/unit1.tachyon_charge
neutron_future_tachyon_quantity = NEUTRON_Future_Q/unit1.tachyon_charge

proton_present_tachyon_quantity = PROTON_Present_Q/unit1.tachyon_charge
proton_future_tachyon_quantity = PROTON_Future_Q/unit1.tachyon_charge

neutron2_present_tachyon_quantity = NEUTRON2_Present_Q/unit1.tachyon_charge
neutron2_future_tachyon_quantity = NEUTRON2_Future_Q/unit1.tachyon_charge

proton2_present_tachyon_quantity = PROTON2_Present_Q/unit1.tachyon_charge
proton2_future_tachyon_quantity = PROTON2_Future_Q/unit1.tachyon_charge

# Determine the mass of tachyons in the segment
neutron_present_tachyon_mass = neutron_present_tachyon_quantity * unit1.tachyon_mass
neutron_future_tachyon_mass = neutron_future_tachyon_quantity * unit1.tachyon_mass

proton_present_tachyon_mass = proton_present_tachyon_quantity * unit1.tachyon_mass
proton_future_tachyon_mass = proton_future_tachyon_quantity * unit1.tachyon_mass

neutron2_present_tachyon_mass = neutron2_present_tachyon_quantity * unit1.tachyon_mass
neutron2_future_tachyon_mass = neutron2_future_tachyon_quantity * unit1.tachyon_mass

proton2_present_tachyon_mass = proton2_present_tachyon_quantity * unit1.tachyon_mass
proton2_future_tachyon_mass = proton2_future_tachyon_quantity * unit1.tachyon_mass

# Determine the mass of the segment without tachyons
neutron_present_withouttachyon_mass = NEUTRON_Present_M - neutron_present_tachyon_mass
neutron_future_withouttachyon_mass = NEUTRON_Future_M - neutron_future_tachyon_mass

proton_present_withouttachyon_mass = PROTON_Present_M - proton_present_tachyon_mass
proton_future_withouttachyon_mass = PROTON_Future_M - proton_future_tachyon_mass

neutron2_present_withouttachyon_mass = NEUTRON2_Present_M - neutron_present_tachyon_mass
neutron2_future_withouttachyon_mass = NEUTRON2_Future_M - neutron_future_tachyon_mass

proton2_present_withouttachyon_mass = PROTON2_Present_M - proton_present_tachyon_mass
proton2_future_withouttachyon_mass = PROTON2_Future_M - proton_future_tachyon_mass

# We calculate the amount of energy in a segment without taking into account tachyons
neutron_present_withouttachyon_energy = neutron_present_withouttachyon_mass * Algorithm.constantc ** 2 
neutron_future_withouttachyon_energy = neutron_future_withouttachyon_mass * Algorithm.constantc ** 2 

proton_present_withouttachyon_energy = proton_present_withouttachyon_mass * Algorithm.constantc ** 2 
proton_future_withouttachyon_energy = proton_future_withouttachyon_mass * Algorithm.constantc ** 2

neutron2_present_withouttachyon_energy = neutron2_present_withouttachyon_mass * Algorithm.constantc ** 2 
neutron2_future_withouttachyon_energy = neutron2_future_withouttachyon_mass * Algorithm.constantc ** 2 

proton2_present_withouttachyon_energy = proton2_present_withouttachyon_mass * Algorithm.constantc ** 2 
proton2_future_withouttachyon_energy = proton2_future_withouttachyon_mass * Algorithm.constantc ** 2

# We calculate the amount of energy in a segment of one tachyon
neutron_present_at1tachyon_energy = neutron_present_withouttachyon_energy/neutron_present_tachyon_quantity
neutron_future_at1tachyon_energy = neutron_future_withouttachyon_energy/neutron_future_tachyon_quantity

proton_present_at1tachyon_energy = proton_present_withouttachyon_energy/proton_present_tachyon_quantity
proton_future_at1tachyon_energy = proton_future_withouttachyon_energy/proton_future_tachyon_quantity

neutron2_present_at1tachyon_energy = neutron2_present_withouttachyon_energy/neutron2_present_tachyon_quantity
neutron2_future_at1tachyon_energy = neutron2_future_withouttachyon_energy/neutron2_future_tachyon_quantity

proton2_present_at1tachyon_energy = proton2_present_withouttachyon_energy/proton2_present_tachyon_quantity
proton2_future_at1tachyon_energy = proton2_future_withouttachyon_energy/proton2_future_tachyon_quantity

unit2 = Tachion(neutron_present_tachyon_quantity, neutron_future_tachyon_quantity, 
                neutron_present_tachyon_mass, neutron_future_tachyon_mass,
                neutron_present_withouttachyon_mass, neutron_future_withouttachyon_mass,
                neutron_present_withouttachyon_energy, neutron_future_withouttachyon_energy,
                neutron_present_at1tachyon_energy, neutron_future_at1tachyon_energy,
                
                proton_present_tachyon_quantity, proton_future_tachyon_quantity, 
                proton_present_tachyon_mass, proton_future_tachyon_mass,
                proton_present_withouttachyon_mass, proton_future_withouttachyon_mass,
                proton_present_withouttachyon_energy, proton_future_withouttachyon_energy,
                proton_present_at1tachyon_energy, proton_future_at1tachyon_energy,
                
                neutron2_present_tachyon_quantity, neutron2_future_tachyon_quantity, 
                neutron2_present_tachyon_mass, neutron2_future_tachyon_mass,
                neutron2_present_withouttachyon_mass, neutron2_future_withouttachyon_mass,
                neutron2_present_withouttachyon_energy, neutron2_future_withouttachyon_energy,
                neutron2_present_at1tachyon_energy, neutron2_future_at1tachyon_energy,
                
                proton2_present_tachyon_quantity, proton2_future_tachyon_quantity, 
                proton2_present_tachyon_mass, proton2_future_tachyon_mass,
                proton2_present_withouttachyon_mass, proton2_future_withouttachyon_mass,
                proton2_present_withouttachyon_energy, proton2_future_withouttachyon_energy,
                proton2_present_at1tachyon_energy, proton2_future_at1tachyon_energy)

# Table for Neutron_time_segments
Neutron_time_segments = namedtuple('Particle', 'name1 name2 name3')

neutron_time_segments = [[1, 'Present \n', NEUTRON_Present_Q, NEUTRON_Present_M, NEUTRON_Present_V],
                          [2, 'Past \n', NEUTRON_Past_Q, NEUTRON_Past_M, NEUTRON_Past_V],
                          [3, 'Future \n', NEUTRON_Future_Q, NEUTRON_Future_M, NEUTRON_Future_V]]
                      
table8 = PrettyTable(['#', 'Distribution of parts of a particle in time', 'Charge in Cl', 
                      'Mass in kg.', 'Volume in cbm'])

for rec in neutron_time_segments:
    table8.add_row(rec)  
    
# Table for Proton_time_segments
Proton_time_segments = namedtuple('Particle', 'name1 name2 name3')

proton_time_segments = [[1, 'Present \n', PROTON_Present_Q, PROTON_Present_M, PROTON_Present_V],
                          [2, 'Past \n', PROTON_Past_Q, PROTON_Past_M, PROTON_Past_V],
                          [3, 'Future \n', PROTON_Future_Q, PROTON_Future_M, PROTON_Future_V]]
                      
table9 = PrettyTable(['#', 'Distribution of parts of a particle in time', 'Charge in Cl', 
                      'Mass in kg.', 'Volume in cbm'])

for rec in proton_time_segments:
    table9.add_row(rec) 
    
# Table for Neutron2_time_segments    
Neutron2_time_segments = namedtuple('Particle', 'name1 name2 name3')

neutron2_time_segments = [[1, 'Present \n', NEUTRON2_Present_Q, NEUTRON2_Present_M, NEUTRON2_Present_V],
                          [2, 'Past \n', NEUTRON2_Past_Q, NEUTRON2_Past_M, NEUTRON2_Past_V],
                          [3, 'Future \n', NEUTRON2_Future_Q, NEUTRON2_Future_M, NEUTRON2_Future_V]]
                      
table10 = PrettyTable(['#', 'Distribution of parts of a particle in time', 'Charge in Cl', 
                       'Mass in kg.', 'Volume in cbm'])

for rec in neutron2_time_segments:
    table10.add_row(rec) 
    
# Table for Proton2_time_segments    
Proton2_time_segments = namedtuple('Particle', 'name1 name2 name3')

proton2_time_segments = [[1, 'Present \n', PROTON2_Present_Q, PROTON2_Present_M, PROTON2_Present_V],
                          [2, 'Past \n', PROTON2_Past_Q, PROTON2_Past_M, PROTON2_Past_V],
                          [3, 'Future \n', PROTON2_Future_Q, PROTON2_Future_M, PROTON2_Future_V]]
                      
table11 = PrettyTable(['#', 'Distribution of parts of a particle in time', 'Charge in Cl', 
                       'Mass in kg.', 'Volume in cbm'])

for rec in proton2_time_segments:
    table11.add_row(rec)   

"""Analysis of the data obtained for protons, neutrons, tachyon."""

# Histograms for protons, neutrons
# The histogram of the volume of shells of protons, neutrons
np_par = [f"P{i}" for i in range(9)]
width = 0.2
x = np.arange(len(np_par))
fig, ax = plt.subplots(figsize=(14,5))
rects1 = ax.bar(x - width/4 -0.2, ([unit.xv02[0], unit.xv02[1], unit.xv02[0], unit.xv02[2], 
                                    unit.xv02[1], unit.xv02[3], unit.xv02[2], unit.xv02[4], 
                                    unit.xv02[5]]), width, label='proton2')
rects2 = ax.bar(x + width/4 -0.2, ([unit.xv13[0], unit.xv13[1], unit.xv13[3], unit.xv13[2], 
                                    unit.xv13[4], unit.xv13[0], unit.xv13[5], unit.xv13[1],
                                    unit.xv13[2]]), width, label='proton')
rects3 = ax.bar(x - width/4 +0.2, ([unit.xv02[3], unit.xv02[4], unit.xv02[3], unit.xv02[5], 
                                    unit.xv02[4], unit.xv02[0], unit.xv02[4], unit.xv02[1],
                                    unit.xv02[2]]), width, label='neutron2')
rects4 = ax.bar(x + width/4 +0.2, ([unit.xv13[3], unit.xv13[4], unit.xv13[0], unit.xv13[5], 
                                    unit.xv13[1], unit.xv13[3], unit.xv13[2], unit.xv13[4], 
                                    unit.xv13[5]]), width, label='neutron')
ax.set_title('The histogram of the volume of shells of protons, neutrons\n'
             'Graph#1', fontsize = 20)
ax.set_xticks(x)
ax.set_xticklabels(np_par, fontsize = 14)
ax.legend(fontsize = 14)

# The histogram of the distribution of electric charge over the shells\n' 
# for protons, neutrons.
np_par = [f"P{i}" for i in range(9)]
width = 0.2
x = np.arange(len(np_par))
fig, ax = plt.subplots(figsize=(14,5))
rects1 = ax.bar(x - width/4 -0.2, ([unit.xq02[0], unit.xq02[1], unit.xq02[0], unit.xq02[2], 
                                    unit.xq02[1], unit.xq02[3], unit.xq02[2], unit.xq02[4], 
                                    unit.xq02[5]]), width, label='proton2')
rects2 = ax.bar(x + width/4 -0.2, ([unit.xq13[0], unit.xq13[1], unit.xq13[3], unit.xq13[2], 
                                    unit.xq13[4], unit.xq13[0], unit.xq13[5], unit.xq13[1],
                                    unit.xq13[2]]), width, label='proton')
rects3 = ax.bar(x - width/4 +0.2, ([unit.xq02[3], unit.xq02[4], unit.xq02[3], unit.xq02[5], 
                                    unit.xq02[4], unit.xq02[0], unit.xq02[4], unit.xq02[1],
                                    unit.xq02[2]]), width, label='neutron2')
rects4 = ax.bar(x + width/4 +0.2, ([unit.xq13[3], unit.xq13[4], unit.xq13[0], unit.xq13[5], 
                                    unit.xq13[1], unit.xq13[3], unit.xq13[2], unit.xq13[4], 
                                    unit.xq13[5]]), width, label='neutron')
ax.set_title('The histogram of the distribution of electric charge over the shells\n' 
             'for protons, neutrons. Graph#2\n', fontsize = 20)
ax.set_xticks(x)
ax.set_xticklabels(np_par, fontsize = 14)
ax.legend(fontsize = 14)

# The histogram of the mass distribution over the shells\n'
# for protons, neutrons.
np_par = [f"P{i}" for i in range(9)]
width = 0.2
x = np.arange(len(np_par))
fig, ax = plt.subplots(figsize=(14,5))
rects1 = ax.bar(x - width/4 -0.2, ([unit.xm02[0], unit.xm02[1], unit.xm02[0], unit.xm02[2], 
                                    unit.xm02[1], unit.xm02[3], unit.xm02[2], unit.xm02[4], 
                                    unit.xm02[5]]), width, label='proton2')
rects2 = ax.bar(x + width/4 -0.2, ([unit.xm13[0], unit.xm13[1], unit.xm13[3], unit.xm13[2], 
                                    unit.xm13[4], unit.xm13[0], unit.xm13[5], unit.xm13[1],
                                    unit.xm13[2]]), width, label='proton')
rects3 = ax.bar(x - width/4 +0.2, ([unit.xm02[3], unit.xm02[4], unit.xm02[3], unit.xm02[5], 
                                    unit.xm02[4], unit.xm02[0], unit.xm02[4], unit.xm02[1],
                                    unit.xm02[2]]), width, label='neutron2')
rects4 = ax.bar(x + width/4 +0.2, ([unit.xm13[3], unit.xm13[4], unit.xm13[0], unit.xm13[5], 
                                    unit.xm13[1], unit.xm13[3], unit.xm13[2], unit.xm13[4], 
                                    unit.xm13[5]]), width, label='neutron')
ax.set_title('The histogram of the mass distribution over the shells\n'
             'for protons, neutrons. Graph#3\n', fontsize = 20)
ax.set_xticks(x)
ax.set_xticklabels(np_par, fontsize = 14)
ax.legend(fontsize = 14)

plt.show()

# We will build pie charts based on the obtained data visualized using histograms.
# The values of the protons, neutrons are the positive and negative volume.
# Neutron 2 

data_names = ['In the 3D space \n', 
              'In the invariant space \n']
data_names2 = ['In the 3D space, Electric charge is positive \n', 
               'In the invariant space, Electric charge is negative \n']
data_names3 = ['In the 3D space, Electric charge is positive \n', 
               'In the invariant space, Electric charge is positive \n']

data_values = [(PROTON2_Present_V) * 10e45, 
               -(PROTON2_Past_V + PROTON2_Future_V) * 10e45]
data_values2 = [(PROTON2_Present_Q) * 10e19, 
                -(PROTON2_Past_Q + PROTON2_Future_Q) * 10e19]
data_values3 = [(PROTON2_Present_M) * 10e28, 
                (PROTON2_Past_M + PROTON2_Future_M) * 10e28]
data_values4 = [(PROTON_Present_V) * 10e45, 
                -(PROTON_Past_V + PROTON_Future_V) * 10e45]
data_values5 = [(PROTON_Present_Q) * 10e19, 
                (PROTON_Past_Q + PROTON_Future_Q) * 10e19]
data_values6 = [(PROTON_Present_M) * 10e28, 
                (PROTON_Past_M + PROTON_Future_M) * 10e28]
data_values7 = [(NEUTRON2_Present_V) * 10e45, 
                -(NEUTRON2_Past_V + NEUTRON2_Future_V) * 10e45]
data_values8 = [(NEUTRON2_Present_Q) * 10e19, 
                -(NEUTRON2_Past_Q + NEUTRON2_Future_Q) * 10e19]
data_values9 = [(NEUTRON2_Present_M) * 10e28, 
                (NEUTRON2_Past_M + NEUTRON2_Future_M) * 10e28]
data_values10 = [(NEUTRON_Present_V) * 10e45, 
                 -(NEUTRON_Past_V + NEUTRON_Future_V) * 10e45]
data_values11 = [(NEUTRON_Present_Q) * 10e19, 
                 -(NEUTRON_Past_Q + NEUTRON_Future_Q) * 10e19]
data_values12 = [(NEUTRON_Present_M) * 10e28, 
                 (NEUTRON_Past_M + NEUTRON_Future_M) * 10e28]

fig = plt.figure(figsize=plt.figaspect(0.3))
ax = fig.add_subplot(4, 3, 1)

mpl.rcParams.update({'font.size': 10})

plt.title("Distribution to \n placement in spaces in (%) \n"
              "for proton 2 for volume \n Graph#4")

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1, 
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.0, - 0.3, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

ax = fig.add_subplot(4, 3, 2)

mpl.rcParams.update({'font.size': 10})

plt.title("Distribution to \n placement in spaces in (%) \n"
              "for proton 2 for electric charge \n Graph#5")

xs = range(len(data_names))

plt.pie(data_values2, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names2) - 1)])
plt.legend(bbox_to_anchor = (-0.5, - 0.3,  0.25, 0.25),
           loc = 'lower left', labels = data_names2)

ax = fig.add_subplot(4, 3, 3)

mpl.rcParams.update({'font.size': 10})

plt.title("Distribution to \n placement in spaces in (%) \n"
              "for proton 2 for mass \n Graph#6")

xs = range(len(data_names))

plt.pie(data_values3, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.0, - 0.3, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

ax = fig.add_subplot(4, 3, 4)

mpl.rcParams.update({'font.size': 10})

plt.title("Distribution to \n placement in spaces in (%) \n"
              "for proton for volume \n Graph#7")

xs = range(len(data_names))

plt.pie(data_values4, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.0, - 0.3, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

ax = fig.add_subplot(4, 3, 5)

mpl.rcParams.update({'font.size': 10})

plt.title("Distribution to \n placement in spaces in (%) \n"
              "for proton for electric charge \n Graph#8")

xs = range(len(data_names3))

plt.pie(data_values5, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names3) - 1)])
plt.legend(bbox_to_anchor = (-0.5, - 0.3, 0.25, 0.25),
           loc = 'lower left', labels = data_names3)

ax = fig.add_subplot(4, 3, 6)

mpl.rcParams.update({'font.size': 10})

plt.title("Distribution to \n placement in spaces in (%) \n"
              "for proton for mass \n Graph#9")

xs = range(len(data_names))

plt.pie(data_values6, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.0, - 0.3, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

ax = fig.add_subplot(4, 3, 7)

mpl.rcParams.update({'font.size': 10})

plt.title("Distribution to \n placement in spaces in (%) \n"
              "for neutron 2 for volume \n Graph#10")

xs = range(len(data_names))

plt.pie(data_values7, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.0, - 0.3, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

ax = fig.add_subplot(4, 3, 8)

mpl.rcParams.update({'font.size': 10})

plt.title("Distribution to \n placement in spaces in (%) \n"
              "for neutron 2 for electric charge \n Graph#11")

xs = range(len(data_names2))

plt.pie(data_values8, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names2) - 1)])
plt.legend(bbox_to_anchor = (-0.5, - 0.3, 0.25, 0.25),
           loc = 'lower left', labels = data_names2)

ax = fig.add_subplot(4, 3, 9)

mpl.rcParams.update({'font.size': 10})

plt.title("Distribution to \n placement in spaces in (%) \n"
              "for neutron 2 for mass \n Graph#12")

xs = range(len(data_names))

plt.pie(data_values9, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.0, - 0.3, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

ax = fig.add_subplot(4, 3, 10)

mpl.rcParams.update({'font.size': 10})

plt.title("Distribution to \n placement in spaces in (%) \n"
              "for neutron for volume \n Graph#13")

xs = range(len(data_names))

plt.pie(data_values10, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.0, - 0.3, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

ax = fig.add_subplot(4, 3, 11)

mpl.rcParams.update({'font.size': 10})

plt.title("Distribution to \n placement in spaces in (%) \n"
              "for neutron for electric charge \n Graph#14")

xs = range(len(data_names2))

plt.pie(data_values11, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names2) - 1)])
plt.legend(bbox_to_anchor = (-0.5, - 0.3, 0.25, 0.25),
           loc = 'lower left', labels = data_names2)

ax = fig.add_subplot(4, 3, 12)

mpl.rcParams.update({'font.size': 10})

plt.title("Distribution to \n placement in spaces in (%) \n"
              "for neutron for mass \n Graph#15")

xs = range(len(data_names))

plt.pie(data_values12, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.0, - 0.3, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

plt.subplots_adjust(left=0.3,
                    bottom=0.001, 
                    right=1.0, 
                    top=2.9, 
                    wspace=0.5, 
                    hspace=0.8)

# Time distribution of proton, neutron

# Time distribution of proton volume in (%)
fig = plt.figure(figsize=plt.figaspect(0.3))
data_names = ['Present time \n', 
              'Past time: \n',
              'Future']
data_values = [(PROTON_Present_V) * 10e46,
               -(PROTON_Past_V) * 10e46, 
               -(PROTON_Future_V) * 10e46]
ax = fig.add_subplot(2, 3, 1)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n proton volume in (%)\n Graph#16')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Time distribution of proton electric charge in (%)
data_names = ['Present time \n', 
              'Past time: \n',
              'Future']
data_values = [(PROTON_Present_Q) * 10e46,
               -(PROTON_Past_Q) * 10e46, 
               (PROTON_Future_Q) * 10e46]
ax = fig.add_subplot(2, 3, 2)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n proton electric charge in (%)\n Graph#17')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Time distribution of proton mass in (%)
data_names = ['Present time \n', 
              'Past time: \n',
              'Future']
data_values = [(PROTON_Present_M) * 10e46,
               (PROTON_Past_M) * 10e46, 
               (PROTON_Future_M) * 10e46]
ax = fig.add_subplot(2, 3, 3)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n proton mass in (%)\n Graph#18')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Time distribution of proton2 volume in (%)
data_names = ['Present time \n', 
              'Past time: \n',
              'Future']
data_values = [(PROTON2_Present_V) * 10e46,
               -(PROTON2_Past_V) * 10e46, 
               -(PROTON2_Future_V) * 10e46]
ax = fig.add_subplot(2, 3, 4)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n proton2 volume in (%)\n Graph#19')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Time distribution of proton2 electric charge in (%)
data_names = ['Present time \n', 
              'Past time: \n',
              'Future']
data_values = [(PROTON2_Present_Q) * 10e20,
               -(PROTON2_Past_Q) * 10e20, 
               (PROTON2_Future_Q) * 10e20]
ax = fig.add_subplot(2, 3, 5)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n proton2 electric charge in (%)\n Graph#20')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)] )
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Time distribution of proton2 mass in (%)
data_names = ['Present time \n', 
              'Past time: \n',
              'Future']
data_values = [(PROTON2_Present_M) * 10e29,
               (PROTON2_Past_M) * 10e29, 
               (PROTON2_Future_M) * 10e29]

ax = fig.add_subplot(2, 3, 6)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n proton2 mass in (%)\n Graph#21')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)] )
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

plt.subplots_adjust(left=0.3,
                    bottom=0.001, 
                    right=1.0, 
                    top=2.5, 
                    wspace=0.1, 
                    hspace=0.1)

# Time distribution of neutron volume in (%)
fig = plt.figure(figsize=plt.figaspect(0.3))
data_names = ['Present time \n', 
              'Past time: \n',
              'Future']
data_values = [(NEUTRON_Present_V) * 10e46,
               -(NEUTRON_Past_V) * 10e46, 
               -(NEUTRON_Future_V) * 10e46]
ax = fig.add_subplot(2, 3, 1)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n neutron volume in (%)\n Graph#22')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Time distribution of neutron electric charge in (%)
data_names = ['Present time \n', 
              'Past time: \n',
              'Future']
data_values = [(NEUTRON_Present_Q) * 10e20,
               -(NEUTRON_Past_Q) * 10e20, 
               (NEUTRON_Future_Q) * 10e20]
ax = fig.add_subplot(2, 3, 2)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n neutron electric charge in (%)\n Graph#23')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Time distribution of neutron mass in (%)
data_names = ['Present time \n', 
              'Past time: \n',
              'Future']
data_values = [(NEUTRON_Present_M) * 10e29,
               (NEUTRON_Past_M) * 10e29, 
               (NEUTRON_Future_M) * 10e29]
ax = fig.add_subplot(2, 3, 3)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n neutron mass in (%)\n Graph#24')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Time distribution of neutron2 volume in (%)
data_names = ['Present time \n', 
              'Past time: \n',
              'Future']
data_values = [(NEUTRON2_Present_V) * 10e46,
               -(NEUTRON2_Past_V) * 10e46, 
               -(NEUTRON2_Future_V) * 10e46]
ax = fig.add_subplot(2, 3, 4)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n neutron2 volume in (%)\n Graph#25')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)] )
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Time distribution of neutron2 electric charge in (%)
data_names = ['Present time \n', 
              'Past time: \n',
              'Future']
data_values = [(NEUTRON2_Present_Q) * 10e20,
               -(NEUTRON2_Past_Q) * 10e20, 
               (NEUTRON2_Future_Q) * 10e20]
ax = fig.add_subplot(2, 3, 5)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n neutron2 electric charge in (%)\n Graph#26')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)] )
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Time distribution of neutron2 mass in (%)
data_names = ['Present time \n', 
              'Past time: \n',
              'Future']
data_values = [(NEUTRON2_Present_M) * 10e29,
               (NEUTRON2_Past_M) * 10e29, 
               (NEUTRON2_Future_M) * 10e29]
ax = fig.add_subplot(2, 3, 6)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n neutron2 mass in (%)\n Graph#27')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)] )
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

plt.subplots_adjust(left=0.3,
                    bottom=0.001, 
                    right=1.0, 
                    top=2.5, 
                    wspace=0.1, 
                    hspace=0.1)

# Demo graph comparing charge and mass from tachyon to neutron
x = np.array([tachyon_mass, unit.me, unit.xm02[3] + unit.xm02[4] + 
              unit.xm02[3] + unit.xm02[5] + unit.xm02[4] + 
              unit.xm02[0] + unit.xm02[4] + unit.xm02[1] + 
              unit.xm02[2], unit.xm02[0] + unit.xm02[1] + 
              unit.xm02[0] + unit.xm02[2] + unit.xm02[1] + 
              unit.xm02[3] + unit.xm02[2] + unit.xm02[4] + 
              unit.xm02[5], unit.xm13[0] + unit.xm13[1] + 
              unit.xm13[3] + unit.xm13[2] + unit.xm13[4] + 
              unit.xm13[0] + unit.xm13[5] + unit.xm13[1] + 
              unit.xm13[2], unit.xm13[3] + unit.xm13[4] + 
              unit.xm13[0] + unit.xm13[5] + unit.xm13[1] + 
              unit.xm13[3] + unit.xm13[2] + unit.xm13[4] + 
              unit.xm13[5]])

znp = np.array([tachyon_charge, -unit.qe, unit.xq02[3] + 
                unit.xq02[4] + unit.xq02[3] + unit.xq02[5] + 
                unit.xq02[4] + unit.xq02[0] + unit.xq02[4] + 
                unit.xq02[1] + unit.xq02[2] , unit.xq02[0] + 
                unit.xq02[1] + unit.xq02[0] + unit.xq02[2] + 
                unit.xq02[1] + unit.xq02[3] + unit.xq02[2] + 
                unit.xq02[4] + unit.xq02[5], unit.xq13[0] + 
                unit.xq13[1] + unit.xq13[3] + unit.xq13[2] + 
                unit.xq13[4] + unit.xq13[0] + unit.xq13[5] +
                unit.xq13[1] + unit.xq13[2], unit.xq13[3] + 
                unit.xq13[4] + unit.xq13[0] + unit.xq13[5] + 
                unit.xq13[1] + unit.xq13[3] + unit.xq13[2] + 
                unit.xq13[4] + unit.xq13[5]])  

fig, axs = plt.subplots(1, 1, figsize=(14, 11))

axs.plot(x, znp, 'bs', label= 'Weight in kg and charge \n in coulombs, respectively')

plt.ylabel('The amount of charge \n \n in Cl х Е-19', fontsize=15)
plt.xlabel('The amount of mass', fontsize=15)

plt.text(0, 0.01e-18, "Tachyon")
plt.text(0.02e-27, 0, tachyon_mass)
plt.text(0.02e-27, -0.01e-18, tachyon_charge)

plt.text(0.05e-27, -unit.qe + 0.015e-18, "Electron")
plt.text(0.05e-27, -unit.qe + 0.005e-18, unit.me)
plt.text(0.05e-27, -unit.qe - 0.025e-19, -unit.qe)

plt.text(1.5e-27 - 0.05e-27, unit.qe + 0.05e-19, "Proton2&Proton")

plt.text(1.5e-27 - 0.5e-27, unit.qe, unit.xm02[0] + unit.xm02[1] + 
         unit.xm02[0] + unit.xm02[2] + unit.xm02[1] + 
         unit.xm02[3] + unit.xm02[2] + unit.xm02[4] + 
         unit.xm02[5])

plt.text(1.5e-27 - 0.5e-27, unit.qe - 0.1e-19, unit.xq02[0] + unit.xq02[1] + 
         unit.xq02[0] + unit.xq02[2] + unit.xq02[1] + 
         unit.xq02[3] + unit.xq02[2] + unit.xq02[4] + 
         unit.xq02[5])


plt.text(1.5e-27 - 0.2e-27, unit.qe - 0.2e-19, unit.xm13[0] + unit.xm13[1] + 
         unit.xm13[3] + unit.xm13[2] + unit.xm13[4] + 
         unit.xm13[0] + unit.xm13[5] + 
         unit.xm13[1] + unit.xm13[2])

plt.text(1.5e-27 - 0.2e-27, unit.qe - 0.3e-19, unit.xq13[0] + unit.xq13[1] + 
         unit.xq13[3] + unit.xq13[2] + unit.xq13[4] + 
         unit.xq13[0] + unit.xq13[5] +
         unit.xq13[1] + unit.xq13[2])

plt.text(1.55e-27, 0.4e-19, "Neutron2")

plt.text(1.15e-27, 0.25e-19, unit.xm02[3] + unit.xm02[4] + 
         unit.xm02[3] + unit.xm02[5] + unit.xm02[4] + 
         unit.xm02[0] + unit.xm02[4] + unit.xm02[1] + 
         unit.xm02[2])

plt.text(1.15e-27, 0.15e-19, unit.xq02[3] + unit.xq02[4] + 
         unit.xq02[3] + unit.xq02[5] + unit.xq02[4] + 
         unit.xq02[0] + unit.xq02[4] + 
         unit.xq02[1] + unit.xq02[2])

plt.text(1.5e-27, 0 -0.05e-19, "Neutron")
plt.text(1.3e-27, 0 -0.15e-19, unit.xm13[3] + unit.xm13[4] + 
              unit.xm13[0] + unit.xm13[5] + unit.xm13[1] + 
              unit.xm13[3] + unit.xm13[2] + unit.xm13[4] + 
              unit.xm13[5])

plt.text(1.3e-27, 0-0.25e-19, unit.xq13[3] + 
                unit.xq13[4] + unit.xq13[0] + unit.xq13[5] + 
                unit.xq13[1] + unit.xq13[3] + unit.xq13[2] + 
                unit.xq13[4] + unit.xq13[5])


yticks(fontsize=12)
plt.legend(loc='upper left', fontsize=18)
grid()         
plt.title('Mass and charge from tachyon to neutron\n' 
          'Graph#28\n', fontsize=20)


# THE DISTRIBUTION ELECTRIC CHARGE FOR NEUTRON
# FOR SEGMENT AT PRESENT`S TIME

x = np.array([0, 1, 2, 3, 4])

# neutronv free state
znp = np.array([unit.xq13[3], unit.xq13[5], unit.xq13[3], unit.xq13[2], unit.xq13[5]])  

xx = np.linspace(x.min(),x.max(), 1000)
fig, axs = plt.subplots(1, 1, figsize=(14, 11))

itp2 = PchipInterpolator(x,znp)
window_size, poly_order = 3, 1

znpznp_sg = savgol_filter(itp2(xx), window_size, poly_order)

axs.plot(x, znp, 'bs', label= 'The neutron')
axs.plot(xx, znpznp_sg, 'b', label= "Smoothed curve")

# or fit to a global function for neutron2 
def func(x, A, B, x0, sigma):
    return A+B*np.tanh((x-x0)/sigma)
    
fit, _ = curve_fit(func, x, znp)
znpznp_fit = func(xx, *fit)

axs.plot(xx, znpznp_fit, 'b--', 
         label=r"$f(xnn) = |A| + B \tanh\left(\frac{x-x_0}{\sigma}\right)$")

plt.ylabel('The amount of charge \n \n in Cl х Е-20', fontsize=15)
plt.xlabel('Shells & present time', fontsize=15)

yticks(fontsize=12)
plt.legend(loc='upper left', fontsize=16)
grid()         
plt.title('THE DISTRIBUTION ELECTRIC CHARGE FOR NEUTRON \n' 
          'FOR SEGMENT AT PRESENT`S TIME\n Graph#29.\n', fontsize=17)

# Interrelation of electric charge, volume and mass,
# segment for the present time, 3D graf

fig = plt.figure(figsize=plt.figaspect(0.3))

ax = fig.add_subplot(1, 2, 1, projection='3d')

Xnn = ([unit.xq02[0], unit.xq02[0], unit.xq02[2], unit.xq02[2], unit.xq02[5]])
Ynn = ([unit.xm02[0], unit.xm02[0], unit.xm02[2], unit.xm02[2], unit.xm02[5]])
Znn = ([unit.xv02[0], unit.xv02[0], unit.xv02[2], unit.xv02[2], unit.xv02[5]])

ax.plot(Xnn,Ynn,Znn)

ax.set_xlabel('\n \n \n Electric charge \n ', fontsize = 15)
ax.set_zlabel('\n \n \n \n \n Volume \n ', fontsize = 15)
ax.set_ylabel('\n \n \n \n Mass\n ', fontsize = 15)

ax.text2D(0.2, 0.95,         
          "Proton 2 segment & present time \n" 
          "Graph # 30", 
          transform=ax.transAxes, fontsize = 16)

# Interrelation of Q, M and V in neytron, 3D graf

ax = fig.add_subplot(1, 2, 2, projection='3d')

Xnn = ([unit.xq13[3], unit.xq13[5], unit.xq13[3], unit.xq13[2], unit.xq13[5]])
Ynn = ([unit.xm13[3], unit.xm13[5], unit.xm13[3], unit.xm13[2], unit.xm13[5]])
Znn = ([unit.xv13[3], unit.xv13[5], unit.xv13[3], unit.xv13[2], unit.xv13[5]])

ax.plot(Xnn,Ynn,Znn)

ax.set_xlabel('\n \n \n Electric charge \n ', 
              fontsize = 15)
ax.set_zlabel('\n \n \n \n \n Mass \n ', fontsize = 15)
ax.set_ylabel('\n \n \n \n Volume\n ', fontsize = 15)

ax.text2D(0.2, 0.95, 
          
          "Neytron segment & present time \n"
          "Graph # 31", 
          transform=ax.transAxes, fontsize = 16)

# NEUTRON 2
fig = plt.figure(figsize=plt.figaspect(0.2))
ax = fig.add_subplot(2, 2, 1, projection='3d')
n = 500
# Present
X = np.array([unit.xq02[5], unit.xq02[0], unit.xq02[2]])
Y = np.array([unit.xv02[5], unit.xv02[0], unit.xv02[2]])
Z = np.array([unit.xm02[5], unit.xm02[0], unit.xm02[2]])
# Past
X1 = np.array([unit.xq02[3], unit.xq02[3],unit.xq02[1]])
Y1 = np.array([unit.xv02[3], unit.xv02[3],unit.xv02[1]])
Z1 = np.array([unit.xm02[3], unit.xm02[3],unit.xm02[1]])
# Future
X2 = np.array([unit.xq02[4], unit.xq02[4], unit.xq02[4]])
Y2 = np.array([unit.xv02[4], unit.xv02[4], unit.xv02[4]])
Z2 = np.array([unit.xm02[4], unit.xm02[4], unit.xm02[4]])

def randrange(n, vmin, vmax):
    
    """ Helper function to make an array of random numbers having 
    shape (n) with each number distributed Uniform(vmin, vmax)."""

    return (vmax - vmin)*np.random.rand(n) + vmin

plt.scatter(X,Y,Z)

for m, zlow, zhigh in [('s', -2e-30, 0), ('^', 0, 10e-30)]:

    xs = randrange(n, -2e-19, 2e-19)
    ys = randrange(n, -1e-45,2e-45)
    zs = randrange(n, 0,5e-29)
    
    ax.scatter(xs, ys, zs, marker=m,
              s = 5)    
    ax.scatter(X, Y, Z,
           c = 'r',
           s = 200)
    ax.scatter(X1, Y1, Z1,
           c = [[0.1, 0.63, 0.55]],
           s = 200)
    ax.scatter(X2, Y2, Z2,
           c = '#ad09a3',
           s = 200)
    
ax.text2D(0.2, 0.95,         
          "NEUTRON 2 by time segments&shells \n" 
          "Graph # 32", 
          transform=ax.transAxes, fontsize = 16)
ax.set_xlabel('\n \n \n The magnitude of the \n charge in coulombs')
ax.set_ylabel('\n \n \n The volume of the shell \n in cubic meters')
ax.set_zlabel(' \n \n Weight in kg.')

# PROTON 2
ax = fig.add_subplot(2, 2, 2, projection='3d')
n = 500
# Present
X = np.array([unit.xq02[0], unit.xq02[0], unit.xq02[2], 
                    unit.xq02[2], unit.xq02[5]])
Y = np.array([unit.xv02[0], unit.xv02[0], unit.xv02[2], 
                    unit.xv02[2], unit.xv02[5]])
Z = np.array([unit.xm02[0], unit.xm02[0], unit.xm02[2], 
                    unit.xm02[2], unit.xm02[5]])
# Past
X1 = np.array([unit.xq02[1], unit.xq02[1], unit.xq02[3]])
Y1 = np.array([unit.xv02[1], unit.xv02[1], unit.xv02[3]])
Z1 = np.array([unit.xm02[1], unit.xm02[1], unit.xm02[3]])
# Future
X2 = np.array([unit.xq02[4]])
Y2 = np.array([unit.xv02[4]])
Z2 = np.array([unit.xm02[4]])

def randrange(n, vmin, vmax):
    return (vmax - vmin)*np.random.rand(n) + vmin

plt.scatter(X,Y,Z)

for m, zlow, zhigh in [('s', -2e-30, 0), ('^', 0, 10e-30)]:

    xs = randrange(n, -2e-19, 2e-19)
    ys = randrange(n, -1e-45,2e-45)
    zs = randrange(n, 0,5e-29)
    
    ax.scatter(xs, ys, zs, marker=m,
              s = 5)
    
    ax.scatter(X, Y, Z,
           c = 'r',
           s = 200)
    ax.scatter(X1, Y1, Z1,
           c = [[0.1, 0.63, 0.55]],
           s = 200)
    ax.scatter(X2, Y2, Z2,
           c = '#ad09a3',
           s = 200)    
ax.text2D(0.2, 0.95,         
          "PROTON 2 by time segments&shells \n" 
          "Graph # 33", 
          transform=ax.transAxes, fontsize = 16)
ax.set_xlabel('\n \n \n The magnitude of the \n charge in coulombs')
ax.set_ylabel('\n \n \n The volume of the shell \n in cubic meters')
ax.set_zlabel('\n \n Weight in kg.')

ax.text2D(- 1.0, - 0.15,         
          "Red - Present  \n\n Aquamarine - Past \n\n Lilac - Future \n\n" 
          "Small dots - 500 \n generated numbers", 
          transform=ax.transAxes, fontsize = 16)

# NEUTRON
ax = fig.add_subplot(2, 2, 3, projection='3d')
n = 500
# Present
X = np.array([unit.xq13[3], unit.xq13[5], unit.xq13[3], 
                      unit.xq13[2], unit.xq13[5]])
Y = np.array([unit.xv13[3], unit.xv13[5], unit.xv13[3], 
                      unit.xv13[2], unit.xv13[5]])
Z = np.array([unit.xm13[3], unit.xm13[5], unit.xm13[3], 
                      unit.xm13[2], unit.xm13[5]])
# Past
X1 = np.array([unit.xq13[4], unit.xq13[4]])
Y1 = np.array([unit.xv13[4], unit.xv13[4]])
Z1 = np.array([unit.xm13[4], unit.xm13[4]])
# Future
X2 = np.array([unit.xq13[0], unit.xq13[1]])
Y2 = np.array([unit.xv13[0], unit.xv13[1]])
Z2 = np.array([unit.xm13[0], unit.xm13[1]])

def randrange(n, vmin, vmax):
    return (vmax - vmin)*np.random.rand(n) + vmin

plt.scatter(X,Y,Z)

for m, zlow, zhigh in [('s', -2e-30, 0), ('^', 0, 10e-30)]:

    xs = randrange(n, -2e-19, 2e-19)
    ys = randrange(n, -1e-45,2e-45)
    zs = randrange(n, 0,5e-29)

    
    ax.scatter(xs, ys, zs, marker=m,
              s = 5)
    
    ax.scatter(X, Y, Z,
           c = 'r',
           s = 200)
    ax.scatter(X1, Y1, Z1,
           c = [[0.1, 0.63, 0.55]],
           s = 200)
    ax.scatter(X2, Y2, Z2,
           c = '#ad09a3',
           s = 200)
    
ax.text2D(0.2, 0.95,         
          "NEUTRON by time segments&shells \n" 
          "Graph # 34", 
          transform=ax.transAxes, fontsize = 16)
ax.set_xlabel('\n \n \n The magnitude of the \n charge in coulombs')
ax.set_ylabel('\n \n \n The volume of the shell \n in cubic meters')
ax.set_zlabel('\n \n Weight in kg.')

# PROTON
ax = fig.add_subplot(2, 2, 4, projection='3d')
n = 500
# Present
X = np.array([unit.xq13[3], unit.xq13[2], unit.xq13[5], unit.xq13[2]])
Y = np.array([unit.xv13[3], unit.xv13[2], unit.xv13[5], unit.xv13[2]])
Z = np.array([unit.xm13[3], unit.xm13[2], unit.xm13[5], unit.xm13[2]])
# Past
X1 = np.array([unit.xq13[4]])
Y1 = np.array([unit.xv13[4]])
Z1 = np.array([unit.xm13[4]])
# Future
X2 = np.array([unit.xq13[0], unit.xq13[1], unit.xq13[0], unit.xq13[1]])
Y2 = np.array([unit.xv13[0], unit.xv13[1], unit.xv13[0], unit.xv13[1]])
Z2 = np.array([unit.xm13[0], unit.xm13[1], unit.xm13[0], unit.xm13[1]])

def randrange(n, vmin, vmax):
    return (vmax - vmin)*np.random.rand(n) + vmin

plt.scatter(X,Y,Z)

for m, zlow, zhigh in [('s', -2e-30, 0), ('^', 0, 10e-30)]:

    xs = randrange(n, -2e-19, 2e-19)
    ys = randrange(n, -1e-45,2e-45)
    zs = randrange(n, 0,5e-29)

    
    ax.scatter(xs, ys, zs, marker=m,
              s = 5)
    
    ax.scatter(X, Y, Z,
           c = 'r',
           s = 200)
    ax.scatter(X1, Y1, Z1,
           c = [[0.1, 0.63, 0.55]],
           s = 200)
    ax.scatter(X2, Y2, Z2,
           c = '#ad09a3',
           s = 200)
ax.text2D(0.2, 0.95,         
          "PROTON by time segments&shells \n" 
          "Graph # 35", 
          transform=ax.transAxes, fontsize = 16)    
ax.set_xlabel('\n \n \n The magnitude of the \n charge in coulombs')
ax.set_ylabel('\n \n \n The volume of the shell \n in cubic meters')
ax.set_zlabel('\n \n Weight in kg.')
plt.subplots_adjust(left=0.3,
                    bottom=0.001, 
                    right=1.0, 
                    top=2.9, 
                    wspace=1.4, 
                    hspace=0.2)
plt.show()

# Analyzing density

# NEUTRON

NEUTRON_Present_p = NEUTRON_Present_M/NEUTRON_Present_V
NEUTRON_Past_p = NEUTRON_Past_M/NEUTRON_Past_V
NEUTRON_Future_p = NEUTRON_Future_M/NEUTRON_Future_V

# Bulk density of electric charge

NEUTRON_Present_pq = NEUTRON_Present_Q/NEUTRON_Present_V
NEUTRON_Past_pq = NEUTRON_Past_Q/NEUTRON_Past_V
NEUTRON_Future_pq = NEUTRON_Future_Q/NEUTRON_Future_V

# PROTON

PROTON_Present_p = PROTON_Present_M/PROTON_Present_V
PROTON_Past_p = PROTON_Past_M/PROTON_Past_V
PROTON_Future_p = PROTON_Future_M/PROTON_Future_V

# Bulk density of electric charge

PROTON_Present_pq = PROTON_Present_Q/PROTON_Present_V
PROTON_Past_pq = PROTON_Past_Q/PROTON_Past_V
PROTON_Future_pq = PROTON_Future_Q/PROTON_Future_V

# NEUTRON2

NEUTRON2_Present_p = NEUTRON2_Present_M/NEUTRON2_Present_V
NEUTRON2_Past_p = NEUTRON2_Past_M/NEUTRON2_Past_V
NEUTRON2_Future_p = NEUTRON2_Future_M/NEUTRON2_Future_V

# Bulk density of electric charge

NEUTRON2_Present_pq = NEUTRON2_Present_Q/NEUTRON2_Present_V
NEUTRON2_Past_pq = NEUTRON2_Past_Q/NEUTRON2_Past_V
NEUTRON2_Future_pq = NEUTRON2_Future_Q/NEUTRON2_Future_V

# PROTON2

PROTON2_Present_p = PROTON2_Present_M/PROTON2_Present_V
PROTON2_Past_p = PROTON2_Past_M/PROTON2_Past_V
PROTON2_Future_p = PROTON2_Future_M/PROTON2_Future_V

# Bulk density of electric charge

PROTON2_Present_pq = PROTON2_Present_Q/PROTON2_Present_V
PROTON2_Past_pq = PROTON2_Past_Q/PROTON2_Past_V
PROTON2_Future_pq = PROTON2_Future_Q/PROTON2_Future_V

# Comparison of the values of density by mass, density 
# by volumetric electric charge in the time segments
# between a proton and a neutron, a proton2 and a neutron2.

# Comparison of density for proton and neutron by time segments
P_N_Present_p = NEUTRON_Present_p == PROTON_Present_p
P_N_Past_p = NEUTRON_Past_p == PROTON_Past_p
P_N_Future_p = NEUTRON_Future_p == PROTON_Future_p

# Comparison of the volumetric density of the electric 
# charge for a proton and a neutron by time segments
P_N_Present_pq = NEUTRON_Present_pq == PROTON_Present_pq
P_N_Past_pq = NEUTRON_Past_pq == PROTON_Past_pq
P_N_Future_pq = NEUTRON_Future_pq == PROTON_Future_pq

# Comparison of density for proton2 and neutron2 by time segments
P2_N_Present_p = NEUTRON2_Present_p == PROTON2_Present_p
P2_N_Past_p = NEUTRON2_Past_p == PROTON2_Past_p
P2_N_Future_p = NEUTRON2_Future_p == PROTON2_Future_p

# Comparison of the volumetric density of the electric 
# charge for a proton2 and a neutron2 by time segments
P2_N_Present_pq = NEUTRON2_Present_pq == PROTON2_Present_pq
P2_N_Past_pq = NEUTRON2_Past_pq == PROTON2_Past_pq
P2_N_Future_pq = NEUTRON2_Future_pq == PROTON2_Future_pq

# The results obtained are presented in the table
Density_comparison = namedtuple('Density_comparison', 'Present Past Future')
density = [[1, 'The density of the proton \n'
            'and neutron are equal \n', P_N_Present_p, P_N_Past_p, P_N_Future_p],
           [2, 'The bulk density of electric charge of \n'
            'the proton and neutron are equal \n', P_N_Present_pq, P_N_Past_pq, P_N_Future_pq],
           [3, 'The density of the proton2 \n'
            'and neutron2 are equal \n', P2_N_Present_p, P2_N_Past_p,  P2_N_Future_p],
           [4, 'The bulk density of electric charge of\n'
            'the proton2 and neutron2 are equal', P2_N_Present_pq, P2_N_Past_pq, P2_N_Future_pq]]

table12 = PrettyTable(['#', 'Description','Present', 'Past', 'Future'])
for rec in density:
    table12.add_row(rec)

print('\n \n \n Initial conditions. Table 1.\n')
print(table1)

print("\n \n Quarks 'u' and 'd', proof of the algorithm. Table2")
print(table2)

print('\nValues of electric charge, mass, volume by shells for the proton2.\n'
     'Table 3.')
print(table3)

print('\nValues of electric charge, mass, volume by shells for the proton.\n'
     'Table 4.')
print(table4)

print('\nValues of electric charge, mass, volume by shells for the neutron2.\n'
     'Table 5.')
print(table5)

print('\nValues of electric charge, mass, volume by shells for the neutron.\n'
     'Table 6.')
print(table6)

print('\n Mass, electric charge and volume of the tachyon.\n'
      'Table 7.')
print(table7)

print('\n The distribution of characteristics for neutron over time segments. Table 8.\n')
print(table8)

print('\n The distribution of characteristics for proton over time segments. Table 9.\n')
print(table9)

print('\n The distribution of characteristics for neutron2 over time segments. Table 10.\n')
print(table10)

print('\n The distribution of characteristics for proton2 over time segments. Table 11.\n')
print(table11)

print('\n Comparison of the density and volumetric electric \n density of'
      'protons and neutrons by time segments. Table 12')
print(table12)

# Time distribution of neutron & proton volume in (%)
# for past time
fig = plt.figure(figsize=plt.figaspect(0.3))
data_names = ['Protont past time \n', 
              'Neutron past time: \n']
data_values = [-(PROTON_Past_V) * 10e46,
               -(NEUTRON_Past_V) * 10e46]
ax = fig.add_subplot(2, 3, 1)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n neutron&proton \n'
          'volume in (%)\n'
          'for past time. Graph#36')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Time distribution of neutron&proton electric charge in (%)
# for past time
data_names = ['Protont past time \n', 
              'Neutron past time: \n']
data_values = [-(PROTON_Past_Q) * 10e20,
               -(NEUTRON_Past_Q) * 10e20]
ax = fig.add_subplot(2, 3, 2)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n neutron&proton \n'
          'electric charge in (%)\n'
          'for past time. Graph#37')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Time distribution of neutron&proton mass in (%)
# for past time
data_names = ['Proton past time \n', 
              'Neutron past time: \n']
data_values = [(PROTON_Past_M) * 10e29,
               (NEUTRON_Past_M) * 10e29]
ax = fig.add_subplot(2, 3, 3)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n neutron&proton \n'
          'mass in (%)\n'
          'for past time. Graph#38')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Time distribution of neutron&proton volume in (%)
# for future
data_names = ['Proton future time \n', 
              'Neutron future time: \n']
data_values = [-(PROTON_Future_V) * 10e46,
               -(NEUTRON_Future_V) * 10e46]
ax = fig.add_subplot(2, 3, 4)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n neutron&proton \n'
          'volume in (%)\n'
          'for future time. Graph#39')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)] )
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Time distribution of neutron&proton electric charge in (%)
# for future
data_names = ['Proton future time \n', 
              'Neutron future time: \n']
data_values = [(PROTON_Future_Q) * 10e20,
               (NEUTRON_Future_Q) * 10e20]
ax = fig.add_subplot(2, 3, 5)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n neutron&proton \n'
          'electric charge in (%)\n'
          'for future time. Graph#40')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)] )
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Time distribution of neutron&proton mass in (%)
# for future
data_names = ['Proton future time \n', 
              'Neutron future time: \n']
data_values = [(PROTON_Future_M) * 10e29,
               (NEUTRON_Future_M) * 10e29]
ax = fig.add_subplot(2, 3, 6)
mpl.rcParams.update({'font.size': 14})

plt.title('Time distribution of \n neutron&proton \n'
          'mass in (%)\n'
          'for future time. Graph#41')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)] )
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

plt.subplots_adjust(left=0.3,
                    bottom=0.001, 
                    right=1.0, 
                    top=2.5, 
                    wspace=0.1, 
                    hspace=0.1)

