#!/usr/bin/env python
# coding: utf-8

# In[10]:


import numpy as np
from numpy import array
from scipy import linalg

from collections import namedtuple
from prettytable import PrettyTable

class Algorithm(): 
    
    """ The 'class Algorithm()' class is used to calculate nuclei, shells of quarks"""
    """'u', 'd', and contains the data of constants used in the program.""" 

# Combinatorics, matrix calculus is used in this class.
    
# Assigning values to constants.
# Constants with more characters than constants according to US NIST data are index two.
# Assigning values to data that has become common knowledge.

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

# The combination of "u" and "d" quarks allows you to get several options 
# for matrices for the analysis of the proton, neutron.

# The resulting matrix may consist of the union of the variations of 
# the quark matrices "u" and "d":

# for the "u" quark:

#u00 = ['u0', 0, 0, 0, 0]
#      [ 0, 'u1', 0, 0, 0]
#      [ 0, 0, 'u2', 0, 0]

#u01 = [0, 'u0', 0, 0, 0]
#      [0, 0, 'u1', 0, 0]
#      [0, 0, 0, 'u2', 0]

#u02 = [0, 0, 'u0', 0, 0]
#      [0, 0, 0, 'u1', 0]
#      [0, 0, 0, 0, 'u2']
        
# for the "d" quark:

#d00 = ['d0', 0, 0, 0, 0]
#      [ 0, 'd1', 0, 0, 0]
#      [ 0, 0, 'd2', 0, 0]

#d01 = [0, 'd0', 0, 0, 0]
#      [0, 0, 'd1', 0, 0]
#      [0, 0, 0, 'd2', 0]

#d02 = [0, 0, 'd0', 0, 0]
#      [0, 0, 0, 'd1', 0]
#      [0, 0, 0, 0, 'd2']

# Let's calculate the values for quarks.

"""One of the matrix options for calculation in general form will look like this.

[['u0' '0' '0' '0' '0' '0' '0' '0' '0']
 ['0' 'u1' '0' 'd0' '0' '0' '0' '0' '0']
 ['0' '0' 'u2' '0' 'd1' '0' 'u0' '0' '0']
 ['0' '0' '0' '0' '0' 'd2' '0' 'u1' '0']
 ['0' '0' '0' '0' '0' '0' '0' '0' 'u2']]
[['d0' '0' '0' '0' '0' '0' '0' '0' '0']
 ['0' 'd1' '0' 'd0' '0' '0' '0' '0' '0']
 ['0' '0' 'd2' '0' 'd1' '0' 'u0' '0' '0']
 ['0' '0' '0' '0' '0' 'd2' '0' 'u1' '0']
 ['0' '0' '0' '0' '0' '0' '0' '0' 'u2']]"""

"""The same version of the matrix in the form of an array for calculation.

array([[2.0, 1.0, 1.0, 1.0, 1.0, 0.0], 
       [0.0, 1.0, 0.0, 0.0, 0.0, 1.0], 
       [0.0, 0.0, 1.0, 0.0, 0.0, 0.0], 
       [1.0, 0.0, 0.0, 2.0, 2.0, 1.0], 
       [0.0, 1.0, 0.0, 0.0, 0.0, 1.0], 
       [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]])"""

# The number of variants is determined by combinatorics.

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
unique, counts = np.unique(P, return_counts=True)
x1 = dict(zip(unique, counts))
unique, counts = np.unique(H, return_counts=True)
x2 = dict(zip(unique, counts))
unique, counts = np.unique(L, return_counts=True)
x3 = dict(zip(unique, counts))
if x1 == x2 and x2 == x3:    
    print(f'')
else:
    print(f'Check the code')
    
K = np.split(P, [3])
K1 = np.split(K[1], [1])
unique, counts = np.unique(K[0], return_counts=True)
x4 = dict(zip(unique, counts))

unique, counts = np.unique(K1[0], return_counts=True)
x5 = dict(zip(unique, counts))

unique, counts = np.unique(K1[1], return_counts=True)
x6 = dict(zip(unique, counts))

H = np.split(H, [3])
H1 = np.split(H[1], [1])
unique, counts = np.unique(H[0], return_counts=True)
x7 = dict(zip(unique, counts))

unique, counts = np.unique(H1[0], return_counts=True)
x8 = dict(zip(unique, counts))

unique, counts = np.unique(H1[1], return_counts=True)
x9 = dict(zip(unique, counts))

L = np.split(L, [3])

L1 = np.split(L[1], [1])
unique, counts = np.unique(L[0], return_counts=True)
x10 = dict(zip(unique, counts))

unique, counts = np.unique(L1[0], return_counts=True)
x11 = dict(zip(unique, counts))

unique, counts = np.unique(L1[1], return_counts=True)
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


class Particles():
    
    """This class forms the date set for protons, neutrons, calculates"""
    """the characteristics of the tachyon."""
    
    
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

# ATTENTION! If you change the values of the constants.
# THE CYCLE DOES NOT CONTAIN A FORCED INTERRUPTION.
   
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

class Segments():
    
    """This class defines segments of protons, neutrons conditionally located"""
    """in the past, present and future time."""
    
# Imagine a train. Part of the train is behind the railway crossing, part is at 
# the railway crossing, and the last part of the train did not reach the railway crossing.
# At this moment you are in the car in front of the railway crossing.
    
    def __init__ (self, NEUTRON_Present_matrix, NEUTRON_Past_matrix, NEUTRON_Future_matrix, 
                  NEUTRON2_Present_matrix, NEUTRON2_Past_matrix, NEUTRON2_Future_matrix,
                  
                  PROTON_Present_matrix, PROTON_Past_matrix, PROTON_Future_matrix,                   
                  PROTON2_Present_matrix, PROTON2_Past_matrix, PROTON2_Future_matrix, 
                  
                  NEUTRON_Present, NEUTRON_Past, NEUTRON_Future, 
                  NEUTRON2_Present, NEUTRON2_Past, NEUTRON2_Future,
                  
                  PROTON_Present, PROTON_Past, PROTON_Future, 
                  PROTON2_Present, PROTON2_Past, PROTON2_Future,
                  
                  NEUTRON2_Present_V, NEUTRON2_Past_V, NEUTRON2_Future_V,
                  NEUTRON2_Present_Q, NEUTRON2_Past_Q, NEUTRON2_Future_Q,
                  NEUTRON2_Present_M, NEUTRON2_Past_M, NEUTRON2_Future_M,
                  
                  PROTON2_Present_V, PROTON2_Past_V, PROTON2_Future_V,
                  PROTON2_Present_Q, PROTON2_Past_Q, PROTON2_Future_Q,
                  PROTON2_Present_M, PROTON2_Past_M, PROTON2_Future_M,
                  
                  NEUTRON_Present_V, NEUTRON_Past_V, NEUTRON_Future_V,
                  NEUTRON_Present_Q, NEUTRON_Past_Q, NEUTRON_Future_Q,
                  NEUTRON_Present_M, NEUTRON_Past_M, NEUTRON_Future_M,
                  
                  PROTON_Present_V, PROTON_Past_V, PROTON_Future_V,
                  PROTON_Present_Q, PROTON_Past_Q, PROTON_Future_Q,
                  PROTON_Present_M, PROTON_Past_M, PROTON_Future_M):
        
        self.NEUTRON_Present_matrix = NEUTRON_Present_matrix
        self.NEUTRON_Past_matrix = NEUTRON_Past_matrix
        self.NEUTRON_Future_matrix = NEUTRON_Future_matrix
        
        self.NEUTRON2_Present_matrix = NEUTRON2_Present_matrix
        self.NEUTRON2_Past_matrix = NEUTRON2_Past_matrix
        self.NEUTRON2_Future_matrix = NEUTRON2_Future_matrix 
        
        self.PROTON_Present_matrix = PROTON_Present_matrix
        self.PROTON_Past_matrix = PROTON_Past_matrix
        self.PROTON_Future_matrix = PROTON_Future_matrix
        
        self.PROTON2_Present_matrix = PROTON2_Present_matrix
        self.PROTON2_Past_matrix = PROTON2_Past_matrix
        self.PROTON2_Future_matrix = PROTON2_Future_matrix
        
        self.NEUTRON_Present = NEUTRON_Present
        self.NEUTRON_Past = NEUTRON_Past
        self.NEUTRON_Future = NEUTRON_Future
        
        self.NEUTRON2_Present = NEUTRON2_Present
        self.NEUTRON2_Past = NEUTRON2_Past
        self.NEUTRON2_Future = NEUTRON2_Future 
        
        self.PROTON_Present = PROTON_Present
        self.PROTON_Past = PROTON_Past
        self.PROTON_Future = PROTON_Future
        
        self.PROTON2_Present = PROTON2_Present
        self.PROTON2_Past = PROTON2_Past
        self.PROTON2_Future = PROTON2_Future
        
        self.NEUTRON2_Present_V = NEUTRON2_Present_V
        self.NEUTRON2_Past_V = NEUTRON2_Past_V
        self.NEUTRON2_Future_V = NEUTRON2_Future_V
        
        self.NEUTRON2_Present_Q = NEUTRON2_Present_Q
        self.NEUTRON2_Past_Q = NEUTRON2_Past_Q
        self.NEUTRON2_Future_Q = NEUTRON2_Future_Q
        
        self.NEUTRON2_Present_M = NEUTRON2_Present_M
        self.NEUTRON2_Past_M = NEUTRON2_Past_M
        self.NEUTRON2_Future_M = NEUTRON2_Future_M
        
        self.PROTON2_Present_V = PROTON2_Present_V
        self.PROTON2_Past_V = PROTON2_Past_V
        self.PROTON2_Future_V = PROTON2_Future_V
        
        self.PROTON2_Present_Q = PROTON2_Present_Q
        self.PROTON2_Past_Q = PROTON2_Past_Q
        self.PROTON2_Future_Q = PROTON2_Future_Q
        
        self.PROTON2_Present_M = PROTON2_Present_M
        self.PROTON2_Past_M = PROTON2_Past_M
        self.PROTON2_Future_M = PROTON2_Future_M
        
        self.NEUTRON_Present_V = NEUTRON_Present_V
        self.NEUTRON_Past_V = NEUTRON_Past_V
        self.NEUTRON_Future_V = NEUTRON_Future_V
        
        self.NEUTRON_Present_Q = NEUTRON_Present_Q
        self.NEUTRON_Past_Q = NEUTRON_Past_Q
        self.NEUTRON_Future_Q = NEUTRON_Future_Q
        
        self.NEUTRON_Present_M = NEUTRON_Present_M
        self.NEUTRON_Past_M = NEUTRON_Past_M
        self.NEUTRON_Future_M = NEUTRON_Future_M
        
        self.PROTON_Present_V = PROTON_Present_V
        self.PROTON_Past_V = PROTON_Past_V
        self.PROTON_Future_V = PROTON_Future_V
        
        self.PROTON_Present_Q = PROTON_Present_Q
        self.PROTON_Past_Q = PROTON_Past_Q
        self.PROTON_Future_Q = PROTON_Future_Q
        
        self.PROTON_Present_M = PROTON_Present_M
        self.PROTON_Past_M = PROTON_Past_M
        self.PROTON_Future_M = PROTON_Future_M
        
        
        
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
NEUTRON_Present = np.sum(NEUTRON_Present_matrix, axis=0, dtype=None, out=None)
NEUTRON_Past = np.sum(NEUTRON_Past_matrix, axis=0, dtype=None, out=None)
NEUTRON_Future = np.sum(NEUTRON_Future_matrix, axis=0, dtype=None, out=None)
NEUTRON2_Present = np.sum(NEUTRON2_Present_matrix, axis=0, dtype=None, out=None)
NEUTRON2_Past = np.sum(NEUTRON2_Past_matrix, axis=0, dtype=None, out=None)
NEUTRON2_Future = np.sum(NEUTRON2_Future_matrix, axis=0, dtype=None, out=None)

PROTON_Present = np.sum(PROTON_Present_matrix, axis=0, dtype=None, out=None)
PROTON_Past = np.sum(PROTON_Past_matrix, axis=0, dtype=None, out=None)
PROTON_Future = np.sum(PROTON_Future_matrix, axis=0, dtype=None, out=None)
PROTON2_Present = np.sum(PROTON2_Present_matrix, axis=0, dtype=None, out=None)
PROTON2_Past = np.sum(PROTON2_Past_matrix, axis=0, dtype=None, out=None)
PROTON2_Future = np.sum(PROTON2_Future_matrix, axis=0, dtype=None, out=None)

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


unit01 = Segments(NEUTRON_Present_matrix, NEUTRON_Past_matrix, NEUTRON_Future_matrix, 
                  NEUTRON2_Present_matrix, NEUTRON2_Past_matrix, NEUTRON2_Future_matrix,
                  
                  PROTON_Present_matrix, PROTON_Past_matrix, PROTON_Future_matrix,                   
                  PROTON2_Present_matrix, PROTON2_Past_matrix, PROTON2_Future_matrix, 
                  
                  NEUTRON_Present, NEUTRON_Past, NEUTRON_Future, 
                  NEUTRON2_Present, NEUTRON2_Past, NEUTRON2_Future,
                  
                  PROTON_Present, PROTON_Past, PROTON_Future, 
                  PROTON2_Present, PROTON2_Past, PROTON2_Future,
                  
                  NEUTRON2_Present_V, NEUTRON2_Past_V, NEUTRON2_Future_V,
                  NEUTRON2_Present_Q, NEUTRON2_Past_Q, NEUTRON2_Future_Q,
                  NEUTRON2_Present_M, NEUTRON2_Past_M, NEUTRON2_Future_M,
                  
                  PROTON2_Present_V, PROTON2_Past_V, PROTON2_Future_V,
                  PROTON2_Present_Q, PROTON2_Past_Q, PROTON2_Future_Q,
                  PROTON2_Present_M, PROTON2_Past_M, PROTON2_Future_M,
                  
                  NEUTRON_Present_V, NEUTRON_Past_V, NEUTRON_Future_V,
                  NEUTRON_Present_Q, NEUTRON_Past_Q, NEUTRON_Future_Q,
                  NEUTRON_Present_M, NEUTRON_Past_M, NEUTRON_Future_M,
                  
                  PROTON_Present_V, PROTON_Past_V, PROTON_Future_V,
                  PROTON_Present_Q, PROTON_Past_Q, PROTON_Future_Q,
                  PROTON_Present_M, PROTON_Past_M, PROTON_Future_M)

class Molecularhydrogen():
    
# This class calculates the characteristics of the common area for protons 
# in a hydrogen molecule.

    def __init__ (self, Molecularhydrogen1, Molecularhydrogen2, 
                  Molecularhydrogen3, Molecularhydrogen4, 
                  Molecularhydrogen5, Molecularhydrogen6, 
                  Molecularhydrogen7, Molecularhydrogen8):
        
        self.Molecularhydrogen1 = Molecularhydrogen1
        self.Molecularhydrogen2 = Molecularhydrogen2
        self.Molecularhydrogen3 = Molecularhydrogen3
        self.Molecularhydrogen4 = Molecularhydrogen4
    
        self.Molecularhydrogen5 = Molecularhydrogen5
        self.Molecularhydrogen6 = Molecularhydrogen6
        self.Molecularhydrogen7 = Molecularhydrogen7
        self.Molecularhydrogen8 = Molecularhydrogen8
        
# Head-to-tail connection. 
if unit01.PROTON2_Future_Q > 0 and unit01.PROTON_Past_Q < 0: 
    Molecularhydrogen1 = unit01.PROTON2_Future + unit01.PROTON_Past
else:
    print('Molecular_hydrogen1 does not exist')
            
if unit01.PROTON_Future_Q > 0 and unit01.PROTON2_Past_Q < 0: 
    Molecularhydrogen2 = unit01.PROTON_Future + unit01.PROTON2_Past 
else:
    print('Molecular_hydrogen2 does not exist')
        
if unit01.PROTON2_Future_Q > 0 and unit01.PROTON2_Past_Q < 0: 
    Molecularhydrogen3 = unit01.PROTON2_Future + unit01.PROTON2_Past 
else:
    print('Molecular_hydrogen3 does not exist')
            
if unit01.PROTON_Future_Q > 0 and unit01.PROTON_Past_Q < 0: 
    Molecularhydrogen4 = unit01.PROTON_Future + unit01.PROTON_Past 
else:
    print('Molecular_hydrogen4 does not exist')
        
# Core-tail connection. 
if unit01.PROTON2_Present_Q > 0 and unit01.PROTON_Past_Q < 0: 
    Molecularhydrogen5 = unit01.PROTON2_Present + unit01.PROTON_Past 
else:
    print('Molecular_hydrogen5 does not exist')
            
if unit01.PROTON_Present_Q > 0 and unit01.PROTON2_Past_Q < 0: 
    Molecularhydrogen6 = unit01.PROTON_Present + unit01.PROTON2_Past 
else:
    print('Molecular_hydrogen6 does not exist')
        
if unit01.PROTON2_Present_Q > 0 and unit01.PROTON2_Past_Q < 0: 
    Molecularhydrogen7 = unit01.PROTON2_Present + unit01.PROTON2_Past 
else:
    print('Molecular_hydrogen7 does not exist')
            
if unit01.PROTON_Present_Q > 0 and unit01.PROTON_Past_Q < 0: 
    Molecularhydrogen8 = unit01.PROTON_Present + unit01.PROTON_Past
else:
    print('Molecular_hydrogen8 does not exist')

unit02 = Molecularhydrogen(Molecularhydrogen1, Molecularhydrogen2, 
                            Molecularhydrogen3, Molecularhydrogen4, 
                            Molecularhydrogen5, Molecularhydrogen6, 
                            Molecularhydrogen7, Molecularhydrogen8)

