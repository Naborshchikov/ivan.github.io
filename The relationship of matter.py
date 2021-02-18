#!/usr/bin/env python
# coding: utf-8

# In[6]:


"""The relationship of matter.""" 
# Program ver 5.1
import pandas as pd
import numpy as np
from numpy import *
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib.pyplot import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import PchipInterpolator
from scipy.signal import savgol_filter
from prettytable import PrettyTable
from collections import namedtuple
import pandas_profiling
import cufflinks as cf
import plotly.offline

# Uncomment the line below if you plan to view 3D graphics from different angles.
# %matplotlib notebook

# Error elimination, since it does not affect the values obtained
# The graph is rendered in 3D, the package for this type of charts uses the square root
# The presence of a pair of negative numbers excludes their visualization
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

Initial_conditions = [[1, 'Project implemented in Python', 'ver. 3.7.6 \n'],
                      [2, 'ID - Anaconda', 'ver. 2020 02 \n'],
                      [3, 'All data presented in the SI system', 'nist.gov/\n'],
                      [4, 'All constants are taken from the\n'
                       'data the US NIST.\n', 'nist.gov/ \n'],
                      [5, 'Particle structure -\n'                       
                       'published works of Nobel laureates.\n', 
                       'published scientific works \n'],
                      [6, 'Protons, neutrons\n'
                       'have a core and two shells.\n', 'Robert Hofstadter \n'],
                      [7, 'Speed of light in a vacuum, c = 299792458\n', 'nist.gov/ \n'],
                      [8, 'Electrical constant, ε0 = 8.8541878128E−12\n', 'nist.gov/ \n'],
                      [9, 'Gravitational constant, G = 6.67430E-11\n', 'nist.gov/ \n'],
                      [10, 'Electric charge of an electron\n'
                       '-1.602176634e-19 \n', 'nist.gov/ \n'],
                      [11, 'π = 3.14159265358979', "Scientific American\n"],
                      [12, "Planck's constant, h = 6.62607015E−34", 'nist.gov/\n'],
                      [13, 'Electron diameter 10e−22,\n', 
                       'Hans D. Dehmelt Experiments\n'],
                      [14, 'The proton consists of two quarks \n', 
                       'Murray Gell-Mann\n'],
                      [15, 'The newneutron consists of two quarks \n', 
                       'Murray Gell-Mann \n'],
                      [16, 'Quark radius − (0.47 · 10E−16 cm)E2\n'
                       '< RE2 < (0.43 · 10E−16 cm)E2 \n', 
                       'arxiv.org/pdf/1604.01280.pdf \n'],
                      [17, 'Additional information\n', 
                       'Data from available sources. \n'],
                     [18, "Quark condensate provides about 9\n"
                      "percent of the proton's mass\n",
                     'Physical Review Letters, 2018\n,'
                      ' website arXiv.org\n'],
                     [19, 'Electron diameter: 10e−22 \n', 
                      'Nobel lecture, December, 8, 1989,\n'
                      ' Hans D. Dehmelt Experiments with \n'
                      'an isolated subatomic particle at rest\n'],
                     [20, 'proton mass: 1.67262192369E-27\n', 
                      'nist.gov/\n'],
                     [21, 'neutron mass: 1.67492749804E-27\n',
                      'nist.gov/\n'],
                     [22, 'The magnitude of the charge\n'
                      'of the core, shells in the proton\n'
                     'respectively: 0.35; 0.5; 0.15\n', 
                      'Robert Hofstadter the\n'
                      'Nobel laureate\n'], 
                     [23, 'The magnitude of the charge of the core,\n'
                      ' shells in the neutron\n'
                     'respectively: 0.35; - 0.5; 0.15\n', 
                      'Robert Hofstadter the\n'],
                     [24, 'The proton radius: 0.84 fm\n', 
                      'aps.org/publications/apsnews/201806/proton.cfm\n'],
                     [25, 'The neutron radius: 0.8e−15\n', 
                      'Povh, B.; Rith, K.(2002).\n'],
                     [26, 'Rradius of the proton core: 0.23 ± 0.03 F\n', 
                      'https://doi.org/10.1103/PhysRevD.18.2484\n'],
                     [27, 'Rradius of the neutron core:\n'
                      ' from 0.3 to 0.36 fm\n', 
                      'arxiv.org/pdf/1810.00486.pdf\n'],
                     [28, 'The radius of the inner shell of the neutron\n'
                      'is approximately 0.6 fm.\n', 
                      'actaphys.uj.edu.pl/fulltext?series=Reg&vol=30&page=119\n']] 

table1 = PrettyTable(['#', 'Description', 'Link to source/ comments'])

for rec in Initial_conditions:
    table1.add_row(rec)

# Proton, neutron consist of a nucleus and two shells, or three quarks.
# Therefore, quarks must consist of a nucleus and shells.
# A quark can have several shells, like layers of air near the Earth.
# This calculation is limited to two shells.
# Further study of quarks by scientists will possibly increase the number
# shells for the calculation, make it more accurate.

class Algorithm():
# Assigning values to constants.
# Constants with more characters than constants according to US NIST data are index two.

    constantε0 = 8.8541878128e-12
    constantε02 = 8.85418781762039e-12
    
    constantc = 299792458
        
    constantg = 6.67430E-11
    constantg2 = 6.67448478E-11
    
    constanth = 6.62607015e-34
    
    
# Assigning values to data that has become common knowledge.
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
# Rradius of the neutron core, 
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

# Rule 1:
# The calculation takes into account that the quarks of the nucleus
# can not fall on a single line, as it will mean the synthesis of quarks
# and the loss of their identity.

# Rule 2:
# Quarks are connected if there is their intersection is at least one shell.

# Rule 3:
# The combination of quarks is obliged to provide the most dense arrangement.

        self.xq02 = xq02
        self.xq13 = xq13
                
        self.xv02 = xv02
        self.xv13 = xv13
                
        self.xm02 = xm02
        self.xm13 = xm13
                
# The combination of "u" and "d" quarks makes it possible to obtain several 
# variants of matrices for the proton, neutron.

# The matrixes for the proton.
a000 = ['u0', 0,   0,   0,  0]
a001 = [ 0,  'u1', 0,   0,  0]
a002 = [ 0,   0,  'u2', 0,  0]

a003 = [0, 'u0',  0,   0,   0]
a004 = [0,  0,   'u1', 0,   0]
a005 = [0,  0,    0,  'u2', 0]

a006 = [0,  0, 'd0', 0,    0]
a007 = [0,  0,  0,  'd1',  0]
a008 = [0,  0,  0,   0,   'd2']

#a0 = list(zip(a000, a001, a002, a003, a004, a005, a006, a007, a008))

"""The result is a matrix.
a0 = [('u0', 0,   0,   0,   0,   0,   0,   0,   0), 
      (0,   'u1', 0,  'u0', 0,   0,   0,   0,   0), 
      (0,    0,  'u2', 0,  'u1', 0,  'd0', 0,   0),
      (0,    0,   0,   0,   0,  'u2', 0,  'd1', 0), 
      (0,    0,   0,   0,   0,   0,   0,   0,  'd2')]"""

      
# The matrixes for the neutron.
a020 = ['d0', 0,   0,   0,   0]
a021 = [ 0,  'd1', 0,   0,   0]
a022 = [ 0,   0,  'd2', 0,   0]

a023 = [0, 'd0',  0,   0,   0]
a024 = [0,  0,   'd1', 0,   0]
a025 = [0,  0,    0,  'd2', 0]

a026 = [0,  0, 'u0', 0,   0]
a027 = [0,  0,  0,  'u1', 0]
a028 = [0,  0,  0,   0,  'u2']

#a2 = list(zip(a020, a021, a022, a023, a024, a025, a026, a027, a028))

"""The result is a matrix.
a2 = [['d0',  0,  0,   0,   0,   0,   0,   0,   0], 
      [0,   'd1', 0,  'd0', 0,   0,   0,   0,   0], 
      [0,    0,  'd2', 0,  'd1', 0,  'u0', 0,   0],
      [0,    0,   0,   0,   0,  'd2', 0,  'u1', 0], 
      [0,    0,   0,   0,   0,   0,   0,   0,  'u2']]"""

# Since we know the values for the nuclei and shells of the proton, neutron, 
# for the calculation we use the matrices a0 with a2.

"""It looks visually.
[['u0' '0' '0' '0' '0' '0' '0' '0' '0']
 ['0' 'u1' '0' 'u0' '0' '0' '0' '0' '0']
 ['0' '0' 'u2' '0' 'u1' '0' 'd0' '0' '0']
 ['0' '0' '0' '0' '0' 'u2' '0' 'd1' '0']
 ['0' '0' '0' '0' '0' '0' '0' '0' 'd2']]
[['d0' '0' '0' '0' '0' '0' '0' '0' '0']
 ['0' 'd1' '0' 'd0' '0' '0' '0' '0' '0']
 ['0' '0' 'd2' '0' 'd1' '0' 'u0' '0' '0']
 ['0' '0' '0' '0' '0' 'd2' '0' 'u1' '0']
 ['0' '0' '0' '0' '0' '0' '0' '0' 'u2']]"""

# All lines with 0 in the second character form a core.
# The remaining two lines form the inner and outer shell.

# Matrices are converted into an array, taking into account the available
# data for the calculation. 
# The array represents the equations for the proton and neutron.
      
# The top three lines of the array are proton (coefficients for the array)
# (u0+u0 = 2; u1+u1 = 2; u2 = 1; d0 = 1) - core for a0; (d1 = 1; u2 = 1) - 
# inner shell for a0; d2 = 1 - outer shell for a0

# The bottom three lines of the array are a neutron (coefficients for the array)
# (d0+d0 = 2; d1+d1 = 2; d2 = 1; u0 = 1) - core for a2; (d2 = 1; u1 = 1) - 
# inner shell for a2; u2 = 1 - outer shell for a2

x00 = (a000.count('u0') + a001.count('u0') + a002.count('u0') +
       a003.count('u0') + a004.count('u0') + a006.count('u0')) 

x01 = (a000.count('u1') + a001.count('u1') + a002.count('u1') +
       a003.count('u1') + a004.count('u1') + a006.count('u1'))

x02 = (a000.count('u2') + a001.count('u2') + a002.count('u2') +
       a003.count('u2') + a004.count('u2') + a006.count('u2'))

x03 = (a000.count('d0') + a001.count('d0') + a002.count('d0') +
       a003.count('d0') + a004.count('d0') + a006.count('d0'))

x04 = (a000.count('d1') + a001.count('d1') + a002.count('d1') +
       a003.count('d1') + a004.count('d1') + a006.count('d1'))

x05 = (a000.count('d2') + a001.count('d2') + a002.count('d2') +
       a003.count('d2') + a004.count('d2') + a006.count('d2'))

an20 = [0]
an20.insert(0, x00)
an20.insert(1, x01)
an20.insert(2, x02)
an20.insert(3, x03)
an20.insert(4, x04)
an20.insert(5, x05)
an20.pop(6)

x10 = (a005.count('u0') + a007.count('u0')) 
x11 = (a005.count('u1') + a007.count('u1'))
x12 = (a005.count('u2') + a007.count('u2'))
x13 = (a005.count('d0') + a007.count('d0'))
x14 = (a005.count('d1') + a007.count('d1'))
x15 = (a005.count('d2') + a007.count('d2'))

an21 = [0]
an21.insert(0, x10)
an21.insert(1, x11)
an21.insert(2, x12)
an21.insert(3, x13)
an21.insert(4, x14)
an21.insert(5, x15)
an21.pop(6)

x20 = a008.count('u0') 
x21 = a008.count('u1')
x22 = a008.count('u2')
x23 = a008.count('d0')
x24 = a008.count('d1')
x25 = a008.count('d2')

an22 = [0]
an22.insert(0, x20)
an22.insert(1, x21)
an22.insert(2, x22)
an22.insert(3, x23)
an22.insert(4, x24)
an22.insert(5, x25)
an22.pop(6)

x30 = (a020.count('u0') + a021.count('u0') + a022.count('u0') +
       a023.count('u0') + a024.count('u0') + a026.count('u0')) 

x31 = (a020.count('u1') + a021.count('u1') + a022.count('u1') +
       a023.count('u1') + a024.count('u1') + a026.count('u1'))

x32 = (a020.count('u2') + a021.count('u2') + a022.count('u2') +
       a023.count('u2') + a024.count('u2') + a026.count('u2'))

x33 = (a020.count('d0') + a021.count('d0') + a022.count('d0') +
       a023.count('d0') + a024.count('d0') + a026.count('d0'))

x34 = (a020.count('d1') + a021.count('d1') + a022.count('d1') +
       a023.count('d1') + a024.count('d1') + a026.count('d1'))

x35 = (a020.count('d2') + a021.count('d2') + a022.count('d2') +
       a023.count('d2') + a024.count('d2') + a026.count('d2'))

a120 = [0]
a120.insert(0, x30)
a120.insert(1, x31)
a120.insert(2, x32)
a120.insert(3, x33)
a120.insert(4, x34)
a120.insert(5, x35)
a120.pop(6)

x40 = (a025.count('u0') + a027.count('u0')) 
x41 = (a025.count('u1') + a027.count('u1'))
x42 = (a025.count('u2') + a027.count('u2'))
x43 = (a025.count('d0') + a027.count('d0'))
x44 = (a025.count('d1') + a027.count('d1'))
x45 = (a025.count('d2') + a027.count('d2'))

a121 = [0]
a121.insert(0, x40)
a121.insert(1, x41)
a121.insert(2, x42)
a121.insert(3, x43)
a121.insert(4, x44)
a121.insert(5, x45)
a121.pop(6)

x50 = a028.count('u0') 
x51 = a028.count('u1')
x52 = a028.count('u2')
x53 = a028.count('d0')
x54 = a028.count('d1')
x55 = a028.count('d2')

a122 = [0]
a122.insert(0, x50)
a122.insert(1, x51)
a122.insert(2, x52)
a122.insert(3, x53)
a122.insert(4, x54)
a122.insert(5, x55)
a122.pop(6)

a02 = [an20, an21, an22, a120, a121, a122]
a02 = array(a02)

"""It looks visually.
a02 = array([[2.0 , 2.0, 1.0, 1.0, 0.0, 0.0],
             [0.0, 0.0, 1.0, 0.0, 1.0, 0.0], 
             [0.0, 0.0, 0.0, 0.0, 0.0, 1.0], 
             [1.0, 0.0, 0.0, 2.0, 2.0, 1.0], 
             [0.0, 1.0, 0.0, 0.0, 0.0, 1.0], 
             [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]])"""

a010 = ['u0', 0,   0,   0,   0]
a011 = [ 0,  'u1', 0,   0,   0]
a012 = [ 0,   0,  'u2', 0,   0]

a013 = [0, 'd0',  0,   0,   0]
a014 = [0,  0,   'd1', 0,   0]
a015 = [0,  0,    0,  'd2', 0]

a016 = [0,  0, 'u0', 0,   0]
a017 = [0,  0,  0,  'u1', 0]
a018 = [0,  0,  0,   0,  'u2']

#a1 = list(zip(a010, a011, a012, a013, a014, a015, a016, a017, a018))

"""The result is a matrix.
a1 = [('u0',  0,  0,   0,   0,   0,   0,   0,   0), 
      (0,   'u1', 0,  'd0', 0,   0,   0,   0,   0), 
      (0,    0,  'u2', 0,  'd1', 0,  'u0', 0,   0),
      (0,    0,   0,   0,   0,  'd2', 0,  'u1', 0), 
      (0,    0,   0,   0,   0,   0,   0,   0,  'u2')]"""

a030 = ['d0', 0,   0,   0,   0]
a031 = [ 0,  'd1', 0,   0,   0]
a032 = [ 0,   0,  'd2', 0,   0]

a033 = [0, 'u0',  0,   0,   0]
a034 = [0,  0,   'u1', 0,   0]
a035 = [0,  0,    0,  'u2', 0]

a036 = [0,  0, 'd0', 0,   0]
a037 = [0,  0,  0,  'd1', 0]
a038 = [0,  0,  0,   0,  'd2']

#a3 = list(zip(a030, a031, a032, a033, a034, a035, a036, a037, a038))

"""The result is a matrix.
a3 = [('d0', 0,   0,   0,   0,   0,   0,   0,   0), 
      (0,   'd1', 0,  'u0', 0,   0,   0,   0,   0), 
      (0,    0,  'd2', 0,  'u1', 0,  'd0', 0,   0),
      (0,    0,   0,   0,   0,  'u2', 0,  'd1', 0), 
      (0,    0,   0,   0,   0,   0,   0,   0,  'd2')]"""

# Since we know the values for the nuclei and shells of the proton, neutron, 
# for the calculation we use the matrices a1 with a3.
"""It looks visually.
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

# The top three lines of the array are a proton (coefficients for the array)
# (u0+u0 = 2; u1 = 1; u2 = 1; d0 = 1; d1 = 1) - core for a1; (d2 = 1; u1 = 1) - 
# inner shell for a1; u2 = 1 - outer shell for a1

# The bottom three lines of the array are neutron (coefficients for the array)
# (d0+d0 = 2; d1 = 1; d2 = 1; u0 = 1; u1 = 1) - core for a3; (u2 = 1; d1 = 1) - 
# inner shell for a3; d2 = 1 - outer shell for a3

x60 = (a010.count('u0') + a011.count('u0') + a012.count('u0') +
       a013.count('u0') + a014.count('u0') + a016.count('u0')) 

x61 = (a010.count('u1') + a011.count('u1') + a012.count('u1') +
       a013.count('u1') + a014.count('u1') + a016.count('u1'))

x62 = (a010.count('u2') + a011.count('u2') + a012.count('u2') +
       a013.count('u2') + a014.count('u2') + a016.count('u2'))

x63 = (a010.count('d0') + a011.count('d0') + a012.count('d0') +
       a013.count('d0') + a014.count('d0') + a016.count('d0'))

x64 = (a010.count('d1') + a011.count('d1') + a012.count('d1') +
       a013.count('d1') + a014.count('d1') + a016.count('d1'))

x65 = (a010.count('d2') + a011.count('d2') + a012.count('d2') +
       a013.count('d2') + a014.count('d2') + a016.count('d2'))

a123 = [0]
a123.insert(0, x60)
a123.insert(1, x61)
a123.insert(2, x62)
a123.insert(3, x63)
a123.insert(4, x64)
a123.insert(5, x65)
a123.pop(6)

x70 = (a015.count('u0') + a017.count('u0')) 
x71 = (a015.count('u1') + a017.count('u1'))
x72 = (a015.count('u2') + a017.count('u2'))
x73 = (a015.count('d0') + a017.count('d0'))
x74 = (a015.count('d1') + a017.count('d1'))
x75 = (a015.count('d2') + a017.count('d2'))

a124 = [0]
a124.insert(0, x70)
a124.insert(1, x71)
a124.insert(2, x72)
a124.insert(3, x73)
a124.insert(4, x74)
a124.insert(5, x75)
a124.pop(6)


x80 = a018.count('u0') 
x81 = a018.count('u1')
x82 = a018.count('u2')
x83 = a018.count('d0')
x84 = a018.count('d1')
x85 = a018.count('d2')

a125 = [0]
a125.insert(0, x80)
a125.insert(1, x81)
a125.insert(2, x82)
a125.insert(3, x83)
a125.insert(4, x84)
a125.insert(5, x85)
a125.pop(6)

x90 = (a030.count('u0') + a031.count('u0') + a032.count('u0') +
       a033.count('u0') + a034.count('u0') + a036.count('u0')) 

x91 = (a030.count('u1') + a031.count('u1') + a032.count('u1') +
       a033.count('u1') + a034.count('u1') + a036.count('u1'))

x92 = (a030.count('u2') + a031.count('u2') + a032.count('u2') +
       a033.count('u2') + a034.count('u2') + a036.count('u2'))

x93 = (a030.count('d0') + a031.count('d0') + a032.count('d0') +
       a033.count('d0') + a034.count('d0') + a036.count('d0'))

x94 = (a030.count('d1') + a031.count('d1') + a032.count('d1') +
       a033.count('d1') + a034.count('d1') + a036.count('d1'))

x95 = (a030.count('d2') + a031.count('d2') + a032.count('d2') +
       a033.count('d2') + a034.count('d2') + a036.count('d2'))

a126 = [0]
a126.insert(0, x90)
a126.insert(1, x91)
a126.insert(2, x92)
a126.insert(3, x93)
a126.insert(4, x94)
a126.insert(5, x95)
a126.pop(6)

x100 = (a035.count('u0') + a037.count('u0')) 
x101 = (a035.count('u1') + a037.count('u1'))
x102 = (a035.count('u2') + a037.count('u2'))
x103 = (a035.count('d0') + a037.count('d0'))
x104 = (a035.count('d1') + a037.count('d1'))
x105 = (a035.count('d2') + a037.count('d2'))

a127 = [0]
a127.insert(0, x100)
a127.insert(1, x101)
a127.insert(2, x102)
a127.insert(3, x103)
a127.insert(4, x104)
a127.insert(5, x105)
a127.pop(6)

x110 = a038.count('u0') 
x111 = a038.count('u1')
x112 = a038.count('u2')
x113 = a038.count('d0')
x114 = a038.count('d1')
x115 = a038.count('d2')

a128 = [0]
a128.insert(0, x110)
a128.insert(1, x111)
a128.insert(2, x112)
a128.insert(3, x113)
a128.insert(4, x114)
a128.insert(5, x115)
a128.pop(6)

a13 = [a123, a124, a125, a126, a127, a128]
a13 = array(a13)

"""It looks visually.
a13 = array( [[2.0, 1.0, 1.0, 1.0, 1.0, 0.0], 
              [0.0, 1.0, 0.0, 0.0, 0.0, 1.0], 
              [0.0, 0.0, 1.0, 0.0, 0.0, 0.0], 
              [1.0, 1.0, 0.0, 2.0, 1.0, 1.0], 
              [0.0, 0.0, 1.0, 0.0, 1.0, 0.0], 
              [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])"""

# V = 4/3πR**3
# The formula can be changed after updating the geometric shape of a quark core
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
    
bq = array ([Algorithm.SHELLP0, Algorithm.SHELLP1, Algorithm.SHELLP2, 
             Algorithm.SHELLN0, Algorithm.SHELLN1, Algorithm.SHELLN2])
    
bv = array ([vrpc, vrpi, vrpo, vrnc, vrni, vrno])
    
bm = array ([mpc, mpi, mpo, mnc, mni, mno])

# Calculation of the charge in the electric charges of an electron for the
# core and shells of the "u" and "d" quarks, and their twins.
# The numbers from [0] to [2] refer to the "u" quark, and his twin.
# The numbers from [3] to [5] refer to the "d" quark, and his twin.
xq02 = linalg.solve(a02, bq)
xq13 = linalg.solve(a13, bq)

# Calculation of volume for core and shells of the "u" and "d" quarks, 
# and their twins.
# The numbers from [0] to [2] refer to the "u" quark, and his twin.
# The numbers from [3] to [5] refer to the "d" quark, and his twin.
xv02 = linalg.solve(a02, bv)
xv13 = linalg.solve(a13, bv)

# Calculation of mass for core and shells of the "u" and "d" quarks, 
# and their twins.
# The numbers from [0] to [2] refer to the "u" quark, and his twin.
# The numbers from [3] to [5] refer to the "d" quark, and his twin.
xm02 = linalg.solve(a02, bm)
xm13 = linalg.solve(a13, bm)

# Calculation of the charge for the core and shells of the "u" and "d" 
#quarks, and their twins.
# The numbers from [0] to [2] refer to the "u" quark, and his twin.
# The numbers from [3] to [5] refer to the "d" quark, and his twin.
for i, item in enumerate(xq02):
    xq02[i] *= Algorithm.qe
    
for i, item in enumerate(xq13):
    xq13[i] *= Algorithm.qe

unit = Algorithm(xq02, xq13, xv02, xv13, xm02, xm13)

"""Calculation of values for "u" and "d" quarks, their twins"""

class Wave():
    constant = unit.constanth/unit.constantc
    def __init__ (self, wave0, wave1):
        self.wave0 = wave0
        self.wave1 = wave1
        
wave0 = 1/unit.xm02
wave1 = 1/unit.xm13

for i, item in enumerate(wave0): 
    wave0[i] *= Wave.constant
for i, item in enumerate(wave1): 
    wave1[i] *= Wave.constant
    
unit2 = Wave(wave0, wave1)

class Electromagnetism():
    constant1 = 1/(2 * unit.constantε0 * unit.constanth * unit.constantc)
    constant2 = 1/(2 * unit.constantε02 * unit.constanth * unit.constantc)
    def __init__ (self, electro0, electro1, electro02, electro12):
        self.electro0 = electro0
        self.electro1 = electro1
        self.electro02 = electro02
        self.electro12 = electro12
        
electro0 = unit.xq02 ** 2
electro1 = unit.xq13 ** 2

electro02 = unit.xq02 ** 2
electro12 = unit.xq13 ** 2

for i, item in enumerate(electro0): 
    electro0[i] *= Electromagnetism.constant1
for i, item in enumerate(electro1): 
    electro1[i] *= Electromagnetism.constant1
    
for i, item in enumerate(electro02): 
    electro02[i] *= Electromagnetism.constant2
for i, item in enumerate(electro12): 
    electro12[i] *= Electromagnetism.constant2
    
unit3 = Electromagnetism(electro0, electro1, electro02, electro12)

class Gravity():
    constant3 = 2 * unit.π * unit.constantg/(unit.constanth * unit.constantc)
    constant4 = 2 * unit.π * unit.constantg2/(unit.constanth * unit.constantc)
    def __init__ (self, grav0, grav1, grav02, grav12):
        self.grav0 = grav0
        self.grav1 = grav1
        self.grav02 = grav02
        self.grav12 = grav12
        
grav0 = unit.xm02 ** 2
grav1 = unit.xm13 ** 2

grav02 = unit.xm02 ** 2
grav12 = unit.xm13 ** 2

for i, item in enumerate(grav0): 
    grav0[i] *= Gravity.constant3
for i, item in enumerate(grav1): 
    grav1[i] *= Gravity.constant3
    
for i, item in enumerate(grav02): 
    grav02[i] *= Gravity.constant4
for i, item in enumerate(grav12): 
    grav12[i] *= Gravity.constant4
        
unit4 = Gravity(grav0, grav1, grav02, grav12)

class Frequency():
    def __init__ (self, frequency0, frequency1):
        self.frequency0 = frequency0
        self.frequency1 = frequency1
        
frequency0 = 1/unit2.wave0
frequency1 = 1/unit2.wave1

for i, item in enumerate(frequency0):
       frequency0[i] *= unit.constantc
for i, item in enumerate(frequency1):
       frequency1[i] *= unit.constantc
        
unit5 = Frequency(frequency0, frequency1)

"""Formation of a data set for a proton, neutron, their twins,""" 
"""a carrier of an electromagnetic field. """    

class Particles():
    def __init__ (self, proton0, proton1, neutron0, neutron1, melectron_charge,
                 melectron_mass, melectron_volume):
        self.proton0 = proton0
        self.proton1 = proton1
        self.neutron0 = neutron0
        self.neutron1 = neutron1
        self.melectron_charge = melectron_charge
        self.melectron_mass = melectron_mass
        self.melectron_volume = melectron_volume        
        
# Matrices from a0 to a3 from the Algorithm class are used to form a data set
# for a proton, neutron, and their twins. 
# The x...02 values are used for the matrices a0 and a2.
# The x...13 values are used for the matrices a1 and a3.
Newproton = namedtuple('Newproton', 'name1 charge name2 mass name3 volume')

proton0 = [[1, 'pq1', unit.xq02[0], 'pm1', unit.xm02[0], 'pv1', unit.xv02[0]],
           [2, 'pq2', unit.xq02[1], 'pm2', unit.xm02[1], 'pv2', unit.xv02[1]],
           [3, 'pq3', unit.xq02[0], 'pm3', unit.xm02[0], 'pv3', unit.xv02[0]],
           [4, 'pq4', unit.xq02[2], 'pm4', unit.xm02[2], 'pv4', unit.xv02[2]],
           [5, 'pq5', unit.xq02[1], 'pm5', unit.xm02[1], 'pv5', unit.xv02[1]],           
           [6, 'pq6', unit.xq02[3], 'pm6', unit.xm02[3], 'pv6', unit.xv02[3]],           
           [7, 'pq7', unit.xq02[2], 'pm7', unit.xm02[2], 'pv7', unit.xv02[2]],           
           [8, 'pq8', unit.xq02[4], 'pm8', unit.xm02[4], 'pv8', unit.xv02[4]],
           [9, 'pq9', unit.xq02[5], 'pm9', unit.xm02[5], 'pv9', unit.xv02[5]]] 

table3 = PrettyTable(['#', 'Charge sym.', 'Charge in Cl', 'Mass sym.',
                      'Mass in kg.', 'Volume sym.', 'Volume in cbm'])

for rec in proton0:
    table3.add_row(rec) 

Proton = namedtuple('Proton', 'name1 charge name2 mass name3 volume')
proton1 = [[1, 'pq1', unit.xq13[0], 'pm1', unit.xm13[0], 'pv1', unit.xv13[0]], 
           [2, 'pq2', unit.xq13[1], 'pm2', unit.xm13[1], 'pv2', unit.xv13[1]], 
           [3, 'pq3', unit.xq13[3], 'pm3', unit.xm13[3], 'pv3', unit.xv13[3]],
           [4, 'pq4', unit.xq13[2], 'pm4', unit.xm13[2], 'pv4', unit.xv13[2]],
           [5, 'pq6', unit.xq13[4], 'pm6', unit.xm13[4], 'pv6', unit.xv13[4]],           
           [6, 'pq5', unit.xq13[0], 'pm5', unit.xm13[0], 'pv5', unit.xv13[0]],
           [7, 'pq7', unit.xq13[5], 'pm7', unit.xm13[5], 'pv7', unit.xv13[5]],           
           [8, 'pq8', unit.xq13[1], 'pm8', unit.xm13[1], 'pv8', unit.xv13[1]],
           [9, 'pq9', unit.xq13[2], 'pm9', unit.xm13[2], 'pv9', unit.xv13[2]]] 

table4 = PrettyTable(['#', 'Charge sym.', 'Charge in Cl', 'Mass sym.',
                      'Mass in kg.', 'Volume sym.', 'Volume in cbm'])
for rec in proton1:
    table4.add_row(rec)

Newneutron = namedtuple('Newneutron', 'name1 charge name2 mass name3 volume')
neutron0 = [[1, 'nq1', unit.xq02[3], 'nm1', unit.xm02[3], 'nv1', unit.xv02[3]], 
            [2, 'nq2', unit.xq02[4], 'nm2', unit.xm02[4], 'nv2', unit.xv02[4]],
            [3, 'nq5', unit.xq02[3], 'nm5', unit.xm02[3], 'nv5', unit.xv02[3]],
            [4, 'nq4', unit.xq02[5], 'nm4', unit.xm02[5], 'nv4', unit.xv02[5]],
            [5, 'nq8', unit.xq02[4], 'nm8', unit.xm02[4], 'nv8', unit.xv02[4]],            
            [6, 'nq3', unit.xq02[0], 'nm3', unit.xm02[0], 'nv3', unit.xv02[0]],
            [7, 'nq8', unit.xq02[4], 'nm8', unit.xm02[4], 'nv8', unit.xv02[4]],            
            [8, 'nq6', unit.xq02[1], 'nm6', unit.xm02[1], 'nv6', unit.xv02[1]],
            [9, 'nq7', unit.xq02[2], 'nm7', unit.xm02[2], 'nv7', unit.xv02[2]]]

table5 = PrettyTable(['#', 'Charge sym.', 'Charge in Cl', 'Mass sym.',
                      'Mass in kg.', 'Volume sym.', 'Volume in cbm'])
for rec in neutron0:
    table5.add_row(rec)

Neutron = namedtuple('Neutron', 'name1 charge name2 mass name3 volume')    
neutron1 = [[1, 'nq1', unit.xq13[3], 'nm1', unit.xm13[3], 'nv1', unit.xv13[3]], 
            [2, 'nq2', unit.xq13[4], 'nm2', unit.xm13[4], 'nv2', unit.xv13[4]], 
            [3, 'nq3', unit.xq13[0], 'nm3', unit.xm13[0], 'nv3', unit.xv13[0]],
            [4, 'nq4', unit.xq13[5], 'nm4', unit.xm13[5], 'nv4', unit.xv13[5]],
            [5, 'nq6', unit.xq13[1], 'nm6', unit.xm13[1], 'nv6', unit.xv13[1]],
            [6, 'nq5', unit.xq13[3], 'nm5', unit.xm13[3], 'nv5', unit.xv13[3]],            
            [7, 'nq7', unit.xq13[2], 'nm7', unit.xm13[2], 'nv7', unit.xv13[2]],
            [8, 'nq8', unit.xq13[4], 'nm8', unit.xm13[4], 'nv8', unit.xv13[4]],
            [9, 'nq9', unit.xq13[5], 'nm9', unit.xm13[5], 'nv9', unit.xv13[5]]]

table6 = PrettyTable(['#', 'Charge sym.', 'Charge in Cl', 'Mass sym.',
                      'Mass in kg.', 'Volume sym.', 'Volume in cbm'])
for rec in neutron1:
    table6.add_row(rec)

proton0 = list(zip(* proton0))
proton1 = list(zip(* proton1))
neutron0 = list(zip(* neutron0))
neutron1 = list(zip(* neutron1))

proton0_min_charge = min((proton0)[2], key=abs)
proton1_min_charge = min((proton1)[2], key=abs)
neutron0_min_charge = min((neutron0)[2], key=abs)
neutron1_min_charge = min((neutron1)[2], key=abs)

# Let's compare the minimum values of charges in a new proton, new neutron, 
# neutron, proton and find the value of a point charged particle

if (proton0_min_charge == neutron0_min_charge and  
    proton1_min_charge == neutron1_min_charge and 
    proton0_min_charge == proton1_min_charge):
    
# Algorithm for finding GCD by subtraction
   
    a = proton0_min_charge
    b = unit.qe2
    while a != b:
        if a > b:
            a = a - b
        else:
            b = b - a

# A point charged particle(melectron) - a carrier of an electromagnetic field.
melectron_charge = a
        
# Find the mass of a point charged particle.
melectron_mass = unit.me/(unit.qe2/melectron_charge)
 
# Minimum_volume particles. 
melectron_volume = ve/(unit.qe2/melectron_charge)

# Let's define a proton, neutron, new proton, new neutron through the definition
# of electric charge.

if (sum(neutron0[2]) > sum(neutron1[2]) and sum(proton0[2]) == sum(proton1[2])): 
    neutron = neutron1 
    new_neutron = neutron0 
    proton = proton1 
    new_proton = proton0
else:
    print('Algorithm requires verification.')     
                         
unit6 = Particles(proton0, proton1, neutron0, neutron1, melectron_charge,
                 melectron_mass, melectron_volume)

"""Calculation of wave, frequency, gravitational values for"""
"""proton and neutron, their twins."""
# new neutron
nnq = unit6.neutron0[2]
nnq = array(nnq)
# neutron
nq = unit6.neutron1[2]
nq = array(nq)
# new proton
npq = unit6.proton0[2]
npq = array(npq)
# proton
pq = unit6.proton1[2]
pq = array(pq)

# new neutron
nnm = unit6.neutron0[4]
nnm = array(nnm)
# neutron
nm = unit6.neutron1[4]
nm = array(nm)
# new proton
npm = unit6.proton0[4]
npm = array(npm)
# proton
pm = unit6.proton1[4]
pm = array(pm)

class Wavep():
    constant = unit.constanth/unit.constantc
    def __init__ (self, wave01, wave11, wave21, wave31):
        self.wave01 = wave01
        self.wave11 = wave11
        self.wave21 = wave21
        self.wave31 = wave31
# new neutron        
wave01 = 1/nnm
# neutron
wave11 = 1/nm
# new proton
wave21 = 1/npm
# proton
wave31 = 1/pm

for i, item in enumerate(wave01): 
    wave01[i] *= Wavep.constant
for i, item in enumerate(wave11): 
    wave11[i] *= Wavep.constant
for i, item in enumerate(wave21): 
    wave21[i] *= Wavep.constant
for i, item in enumerate(wave31): 
    wave31[i] *= Wavep.constant
    
unit7 = Wavep(wave01, wave11, wave21, wave31)

class Electromagnetismp():
    constant1 = 1/(2 * unit.constantε0 * unit.constanth * unit.constantc)
    constant2 = 1/(2 * unit.constantε02 * unit.constanth * unit.constantc)
    def __init__ (self, electro01, electro11, electro03, electro13,
                  electro02, electro14, electro04, electro15):
        self.electro01 = electro01
        self.electro11 = electro11
        self.electro03 = electro03
        self.electro13 = electro13
        self.electro02 = electro02
        self.electro14 = electro14
        self.electro04 = electro04
        self.electro15 = electro15
# new neutron         
electro01 = nnq ** 2
electro11 = nnq ** 2
# neutron 
electro03 = nq ** 2
electro13 = nq ** 2
# new proton
electro02 = npq ** 2
electro14 = npq ** 2
# proton
electro04 = pq ** 2
electro15 = pq ** 2

for i, item in enumerate(electro01): 
    electro01[i] *= Electromagnetismp.constant1
for i, item in enumerate(electro11): 
    electro11[i] *= Electromagnetismp.constant2
    
for i, item in enumerate(electro03): 
    electro03[i] *= Electromagnetismp.constant1
for i, item in enumerate(electro13): 
    electro13[i] *= Electromagnetismp.constant2
    
for i, item in enumerate(electro02): 
    electro02[i] *= Electromagnetismp.constant1
for i, item in enumerate(electro14): 
    electro14[i] *= Electromagnetismp.constant2
    
for i, item in enumerate(electro04): 
    electro04[i] *= Electromagnetismp.constant1
for i, item in enumerate(electro15): 
    electro15[i] *= Electromagnetismp.constant2
    
unit8 = Electromagnetismp(electro01, electro11, electro03, electro13,
                         electro02, electro14, electro04, electro15)

class Gravityp():
    constant3 = 2 * unit.π * unit.constantg/(unit.constanth * unit.constantc)
    constant4 = 2 * unit.π * unit.constantg2/(unit.constanth * unit.constantc)
    def __init__ (self, grav01, grav11, grav03, grav13, grav02, grav14,
                 grav04, grav15):
        self.grav01 = grav01
        self.grav11 = grav11
        self.grav03 = grav03
        self.grav13 = grav13
        self.grav02 = grav02
        self.grav14 = grav14
        self.grav04 = grav04
        self.grav15 = grav15
# new neutron        
grav01 = nnm ** 2
grav11 = nnm ** 2
# neutron
grav03 = nm ** 2
grav13 = nm ** 2
# new proton
grav02 = npm ** 2
grav14 = npm ** 2
# proton
grav04 = pm ** 2
grav15 = pm ** 2

for i, item in enumerate(grav01): 
    grav01[i] *= Gravityp.constant3
for i, item in enumerate(grav11): 
    grav11[i] *= Gravityp.constant4
    
for i, item in enumerate(grav03): 
    grav03[i] *= Gravityp.constant3
for i, item in enumerate(grav13): 
    grav13[i] *= Gravityp.constant4
    
for i, item in enumerate(grav02): 
    grav02[i] *= Gravityp.constant3
for i, item in enumerate(grav14): 
    grav14[i] *= Gravityp.constant4
    
for i, item in enumerate(grav04): 
    grav04[i] *= Gravityp.constant3
for i, item in enumerate(grav15): 
    grav15[i] *= Gravityp.constant4
        
unit9 = Gravityp(grav01, grav11, grav03, grav13, grav02, grav14,
                 grav04, grav15)

class Frequencyp():
    def __init__ (self, frequency02, frequency03, 
                  frequency04, frequency05):
        self.frequency02 = frequency02
        self.frequency03 = frequency03
        self.frequency04 = frequency04
        self.frequency05 = frequency05
# new neutron        
frequency02 = 1/unit7.wave01
# neutron
frequency03 = 1/unit7.wave11
# new proton
frequency04 = 1/unit7.wave21
# proton
frequency05 = 1/unit7.wave31

for i, item in enumerate(frequency02):
       frequency02[i] *= unit.constantc
for i, item in enumerate(frequency03):
       frequency03[i] *= unit.constantc
        
for i, item in enumerate(frequency04):
       frequency04[i] *= unit.constantc
for i, item in enumerate(frequency05):
       frequency05[i] *= unit.constantc
        
unit10 = Frequencyp(frequency02, frequency03, 
                    frequency04, frequency05)

class Practice():
# The relationship of frequency and wavelength at different speeds.
# Speed - v

    def __init__ (self, frequenn, frequennn, frequenp,
                  frequennp, wavelengthn, wavelengthnn,
                  wavelengthp, wavelengthnp, v):
        self.frequenn = frequenn
        self.frequennn = frequennn
        self.frequenp = frequenp
        self.frequennp = frequennp
        self.wavelengthn = wavelengthn
        self.wavelengthnn = wavelengthnn
        self.wavelengthp = wavelengthp
        self.wavelengthnp = wavelengthnp
        self.v = v
        
frequenn = [1/unit7.wave11[0], 1/unit7.wave11[1], 
            1/unit7.wave11[2], 1/unit7.wave11[3], 
            1/unit7.wave11[4], 1/unit7.wave11[5], 
            1/unit7.wave11[6], 1/unit7.wave11[7], 
            1/unit7.wave11[8]]

frequennn = [1/unit7.wave01[0], 1/unit7.wave01[1], 
             1/unit7.wave01[2], 1/unit7.wave01[3], 
             1/unit7.wave01[4], 1/unit7.wave01[5], 
             1/unit7.wave01[6], 1/unit7.wave01[7], 
             1/unit7.wave01[8]]

frequenp = [1/unit7.wave31[0], 1/unit7.wave31[1], 
            1/unit7.wave31[2], 1/unit7.wave31[3], 
            1/unit7.wave31[4], 1/unit7.wave31[5], 
            1/unit7.wave31[6], 1/unit7.wave31[7], 
            1/unit7.wave31[8]]

frequennp = [1/unit7.wave21[0], 1/unit7.wave21[1], 
             1/unit7.wave21[2], 1/unit7.wave21[3], 
             1/unit7.wave21[4], 1/unit7.wave21[5], 
             1/unit7.wave21[6], 1/unit7.wave21[7], 
             1/unit7.wave21[8]]

# neutron        
for j in range(8):
    wavelengthn = np.random.uniform(unit7.wave11[j], unit7.wave11[j+1], size=(1, 8000))
    frequenn = 1/wavelengthn
    v = np.random.uniform(- Algorithm.constantc, Algorithm.constantc, size=(1, 8000))

for i, item in enumerate(frequenn):
    frequenn[i] *= v[i]  
# new neutron 
for j in range(8):
    wavelengthnn = np.random.uniform(unit7.wave01[j], unit7.wave01[j+1], size=(1, 8000))
    frequennn = 1/wavelengthnn        
for i, item in enumerate(frequennn):
    frequennn *= v[i]
# proton     
for j in range(8):
    wavelengthp = np.random.uniform(unit7.wave31[j], unit7.wave31[j+1], size=(1, 8000))
    frequenp = 1/wavelengthp        
for i, item in enumerate(frequenp):
    frequenp[i] *= v[i]
# new proton    
for j in range(8):
    wavelengthnp = np.random.uniform(unit7.wave21[j], unit7.wave21[j+1], size=(1, 8000))
    frequennp = 1/wavelengthnp        
for i, item in enumerate(frequennp):
    frequennp[i] *= v[i]
        
unit11 = Practice(frequenn, frequennn, frequenp,
                  frequennp, wavelengthn, wavelengthnn,
                  wavelengthp, wavelengthnp, v)

Basic_design_data = [[1, 'Proton core volume \n', vrpc],
                    [2, 'Neutron core volume \n', vrnc],
                    [3, 'The volume of the inner shell of the proton \n', vrpi],
                    [4, 'The volume of the outer shell of the proton \n', vrpo],
                    [5, 'The volume of the inner shell of the neutron \n', vrni],
                    [6, 'The volume of the outer shell of the neutron \n', vrno],
                    [7, 'Proton core mass \n', mpc],
                    [8, 'Neutron core mass \n', mnc],
                    [9, 'The mass of the inner shell of the proton \n', mpi],
                    [10, 'The mass of the inner shell of the neutron \n', mni],
                    [11, 'The mass of the outer shell of the proton \n', mpo],
                    [12, 'The mass of the outer shell of the neutron \n', mno],
                    [13, 'The electric charge of the new quark core "u" \n', unit.xq02[0]],
                    [14, 'The electric charge of the quark nucleus "u" \n', unit.xq13[0]],
                    [15, 'The electric charge of the inner new quark \n'
                     'shell "u" \n', unit.xq02[1]],
                    [16, 'The electric charge of the inner quark \n'
                     'shell "u" \n',unit.xq13[1]],
                    [17, 'The electric charge of the outer new quark \n'
                     'shell "u" \n', unit.xq02[2]],
                    [18, 'The electric charge of the outer quark \n'
                     'shell "u" \n',unit.xq13[2]],
                    [19, 'The mass of the new quark core "u" \n', unit.xm02[0]],
                    [20, 'The mass of the quark nucleus "u" \n', unit.xm13[0]],
                    [21, 'The mass of the inner new quark shell "u" \n', unit.xm02[1]],
                    [22, 'The mass of the inner quark shell "u" \n',unit.xm13[1]],
                    [23, 'The mass of the outer new quark shell "u" \n', unit.xm02[2]],
                    [24, 'The mass of the outer quark shell "u" \n',unit.xm13[2]],
                    [25, 'The volume of the new quark core "u" \n', unit.xv02[0]],
                    [26, 'The volume of the quark nucleus "u" \n', unit.xv13[0]],
                    [27, 'The volume of the inner new quark shell "u" \n', unit.xv02[1]],
                    [28, 'The volume of the inner quark shell "u" \n',unit.xv13[1]],
                    [29, 'The volume of the outer new quark shell "u" \n', unit.xv02[2]],
                    [30, 'The volume of the outer quark shell "u" \n',unit.xv13[2]],                     
                    [31, 'Matrix of the proton + neutron for\n'
                     'calculating the values of the quarks \n'
                     '"u" and "d".\n \n'                     
                     'The first column shows the quantity \n'
                     'for the "u" quark core. \n'
                     'The second column shows the amount for \n'
                     'the inner \n'
                     'shell of the  quark "u". \n'
                     'The third column shows the amount for \n'
                     'the outer shell of the quark "u".\n'
                     'The fourth column shows the quantity \n'
                     'for the "d" quark core. \n'
                     'The fifth column shows the amount \n'
                     'for the inner shell of the quark "d". \n'
                     'The sixth column shows the amount \n'
                     'for the outer \n'
                     'shell of the quark "d". \n',a13],                     
                    [32, 'Matrix of the new proton + \n'
                     'neutron for\n'
                     'calculating the values of the new \n'
                     'quarks\n'
                     '"u" and "d".\n'
                     'The first line is the core of the \n'
                     'new proton.\n'
                     'The second line is the inner \n'
                     'shell of the new proton.\n'
                     'The third line is the outer \n'
                     'shell of the new proton.\n'
                     'The fourth line is the core of the\n'
                     'new neutron.\n'
                     'The fifth line is the inner shell of the\n'
                     'new neutron.\n'
                     'The sixth line is the outer shell of the\n'
                     'new neutron.\n', a02]]

table2 = PrettyTable(['#', 'Description', 'Design data'])
for rec in Basic_design_data:
    table2.add_row(rec)
    
Smallest_electrical_charge = [[1, 'Electric charge \n', unit6.melectron_charge],
                              [2, 'Mass \n', unit6.melectron_mass],
                              [3, 'Volume', unit6.melectron_volume]]
table7 = PrettyTable(['#', 'Description', 'Design data'])

for rec in Smallest_electrical_charge:
                              table7.add_row(rec)
        
# Visualization of the distribution of electric charge
# in a proton, neutron and their twins along the shells.

x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])

# proton free state
yp  = np.array([unit.xq13[0], unit.xq13[1], unit.xq13[3], unit.xq13[2], unit.xq13[4],
                unit.xq13[0], unit.xq13[5], unit.xq13[1], unit.xq13[2]])  

# new proton free state
znp = np.array([unit.xq02[0], unit.xq02[1], unit.xq02[0], unit.xq02[2],
                unit.xq02[1], unit.xq02[3], unit.xq02[2], unit.xq02[4], unit.xq02[5]])  

xx = np.linspace(x.min(),x.max(), 1000)
fig, axs = plt.subplots(1, 1, figsize=(14, 11))

itp1 = PchipInterpolator(x,yp)
itp2 = PchipInterpolator(x,znp)
window_size, poly_order = 57, 2

ypyp_sg = savgol_filter(itp1(xx), window_size, poly_order)
znpznp_sg = savgol_filter(itp2(xx), window_size, poly_order)

axs.plot(x, yp, 'gs', label= 'The proton')
axs.plot(xx, ypyp_sg, 'green', label= "Smoothed curve")

axs.plot(x, znp, 'bs', label= 'The new proton')
axs.plot(xx, znpznp_sg, 'b', label= "Smoothed curve")

# neutron free state
yn  = np.array([unit.xq13[3], unit.xq13[4], unit.xq13[0], unit.xq13[5], unit.xq13[1],
                unit.xq13[3], unit.xq13[2], unit.xq13[4], unit.xq13[5]])  

# new neutron free state
znn = np.array([unit.xq02[3], unit.xq02[4], unit.xq02[3], unit.xq02[5],
                unit.xq02[4], unit.xq02[0], unit.xq02[4], unit.xq02[1], unit.xq02[2]])  

itp3 = PchipInterpolator(x,yn)
itp4 = PchipInterpolator(x,znn)

ynyn_sg = savgol_filter(itp3(xx), window_size, poly_order)
znnznn_sg = savgol_filter(itp4(xx), window_size, poly_order)

axs.plot(x, yn, 'ms', label= 'The neutron')
axs.plot(xx, ynyn_sg, 'm', label= "Smoothed curve")

axs.plot(x, znn, 'ys', label= 'The new neutron')
axs.plot(xx, znnznn_sg, 'y', label= "Smoothed curve")

# or fit to a global function for neutron and new neutron 
def func(x, A, B, x0, sigma):
    return abs(A)+B*np.tanh((x-x0)/sigma)

fit, _ = curve_fit(func, x, yn)
ynyn_fit = func(xx, *fit)

fit, _ = curve_fit(func, x, znn)
znnznn_fit = func(xx, *fit)

axs.plot(xx, ynyn_fit, 'm--', 
         label=r"$f(xn) = |A| + B \tanh\left(\frac{x-x_0}{\sigma}\right)$")

axs.plot(xx, znnznn_fit, 'y--', 
         label=r"$f(xp) = |A| + B \tanh\left(\frac{x-x_0}{\sigma}\right)$")

plt.ylabel('The amount of charge \n \n in Cl х Е-19', fontsize=15)
plt.xlabel('Shell number', fontsize=15)

yticks(fontsize=12)
plt.legend(loc='upper left', fontsize=16)
grid()         
plt.title('THE DISTRIBUTION ELECTRIC CHARGE OVER SHELLS \n Graph#1.\n', fontsize=17)

# Visualization of the distribution of volume
# in a proton, neutron and their twins along the shells.

x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])

# proton free state
ypv  = np.array([unit.xv13[0], unit.xv13[1], unit.xv13[3], unit.xv13[2], unit.xv13[4],
                unit.xv13[0], unit.xv13[5], unit.xv13[1], unit.xv13[2]])  

# new proton free state
znpv = np.array([unit.xv02[0], unit.xv02[1], unit.xv02[0], unit.xv02[2],
                unit.xv02[1], unit.xv02[3], unit.xv02[2], unit.xv02[4], unit.xv02[5]])  

xx = np.linspace(x.min(),x.max(), 1000)
fig, axs = plt.subplots(1, 1, figsize=(14, 11))

itp1 = PchipInterpolator(x,ypv)
itp2 = PchipInterpolator(x,znpv)
window_size, poly_order = 57, 2

ypvypv_sg = savgol_filter(itp1(xx), window_size, poly_order)
znpvznpv_sg = savgol_filter(itp2(xx), window_size, poly_order)

axs.plot(x, ypv, 'gs', label= 'The proton')
axs.plot(xx, ypvypv_sg, 'green', label= "Smoothed curve")

axs.plot(x, znpv, 'bs', label= 'The new proton')
axs.plot(xx, znpvznpv_sg, 'b', label= "Smoothed curve")

# neutron free state
ynv  = np.array([unit.xv13[3], unit.xv13[4], unit.xv13[0], unit.xv13[5], unit.xv13[1],
                unit.xv13[3], unit.xv13[2], unit.xv13[4], unit.xv13[5]])  

# new neutron free state
znnv = np.array([unit.xv02[3], unit.xv02[4], unit.xv02[3], unit.xv02[5],
                unit.xv02[4], unit.xv02[0], unit.xv02[4], unit.xv02[1], unit.xv02[2]])  

itp3 = PchipInterpolator(x,ynv)
itp4 = PchipInterpolator(x,znnv)

ynvynv_sg = savgol_filter(itp3(xx), window_size, poly_order)
znnvznnv_sg = savgol_filter(itp4(xx), window_size, poly_order)

axs.plot(x, ynv, 'ms', label= 'The neutron')
axs.plot(xx, ynvynv_sg, 'm', label= "Smoothed curve")

axs.plot(x, znnv, 'ys', label= 'The new neutron')
axs.plot(xx, znnvznnv_sg, 'y', label= "Smoothed curve")

plt.ylabel('The amount of volume \n \n in cbm х Е-45', fontsize=15)
plt.xlabel('Shell number', fontsize=15)

yticks(fontsize=12)
plt.legend(loc='upper left', fontsize=16)
grid()         
plt.title('THE DISTRIBUTION VOLUME OVER SHELLS \n Graph#2.\n', fontsize=17)

# The interrelation of frequency, electromagnetic, gravitational 
# characteristics of a quark "u" and "d"., 3D graph.

fig = plt.figure()
ax = Axes3D(fig)
X = ([unit5.frequency1[0], unit5.frequency1[1], unit5.frequency1[2]])
Y = ([unit4.grav1[0], unit4.grav1[1], unit4.grav1[2]])
X,Y = np.meshgrid(X,Y)
def f(x,y):
    return (sin(x) + cos(y))
ax.plot_surface(X,Y,f(X,Y), rstride=1, cstride=1, cmap=cm.coolwarm)
ax.view_init(elev=10,azim=155)
ax.set_xlabel('\n \n \n The frequency \n', fontsize = 15)
ax.set_ylabel('\n \n \n Gravitational', fontsize = 15)
ax.text2D(0.2, 0.95, 
          '\nThe interrelation of frequency,\n'
          'gravitational characteristics of \n'
          'a quark "u" at v = c.\n'
          'Graph#3.\n', 
          transform=ax.transAxes, fontsize = 16)
fig = plt.figure()
ax = Axes3D(fig)
X = ([unit5.frequency0[0], unit5.frequency0[1], unit5.frequency0[2]])
Y = ([unit4.grav0[0], unit4.grav0[1], unit4.grav0[2]])
X,Y = np.meshgrid(X,Y)
def f(x,y):
    return (sin(x) + cos(y))
ax.plot_surface(X,Y,f(X,Y), rstride=1, cstride=1, cmap=cm.coolwarm)
ax.view_init(elev=10,azim=155)
ax.set_xlabel('\n \n \n The frequency \n', fontsize = 15)
ax.set_ylabel('\n \n \n Gravitational', fontsize = 15)
ax.text2D(0.2, 0.95, 
          '\nThe interrelation of frequency,\n'
          'gravitational characteristics of \n'
          'a new quark "u" at v = c.\n'
          'Graph#4.\n', 
          transform=ax.transAxes, fontsize = 16)
fig = plt.figure()
ax = Axes3D(fig)
X = ([unit3.electro1[0], unit3.electro1[1], unit3.electro1[2]])
Y = ([unit4.grav1[0], unit4.grav1[1], unit4.grav1[2]])
X,Y = np.meshgrid(X,Y)
def f(x,y):
    return (sin(x) + cos(y))
ax.plot_surface(X,Y,f(X,Y), rstride=1, cstride=1, cmap=cm.coolwarm)
ax.view_init(elev=60,azim=70)
ax.set_xlabel('\n \n \n Electromagnetic \n', fontsize = 15)
ax.set_ylabel('\n \n \n Gravitational', fontsize = 15)

ax.text2D(0.2, 0.95, 
          '\nThe interrelation of electromagnetic,\n'
          'gravitational characteristics of \n'
          'a quark "u" at v = c.\n'
          'Graph#5.\n', 
          transform=ax.transAxes, fontsize = 16)
fig = plt.figure()
ax = Axes3D(fig)
X = ([unit3.electro0[0], unit3.electro0[1], unit3.electro0[2]])
Y = ([unit4.grav0[0], unit4.grav0[1], unit4.grav0[2]])
X,Y = np.meshgrid(X,Y)
def f(x,y):
    return (sin(x) + cos(y))
ax.plot_surface(X,Y,f(X,Y), rstride=1, cstride=1, cmap=cm.coolwarm)
ax.view_init(elev=60,azim=70)
ax.set_xlabel('\n \n \n Electromagnetic \n', fontsize = 15)
ax.set_ylabel('\n \n \n Gravitational', fontsize = 15)

ax.text2D(0.2, 0.95, 
          '\nThe interrelation of electromagnetic,\n'
          'gravitational characteristics of \n'
          'a new quark "u" at v = c.\n'
          'Graph#6.\n', 
          transform=ax.transAxes, fontsize = 16)
fig = plt.figure()
ax = Axes3D(fig)
X = ([unit5.frequency1[3], unit5.frequency1[4], unit5.frequency1[5]])
Y = ([unit4.grav1[3], unit4.grav1[4], unit4.grav1[5]])
X,Y = np.meshgrid(X,Y)
def f(x,y):
    return (sin(x) + cos(y))
ax.plot_surface(X,Y,f(X,Y), rstride=1, cstride=1, cmap=cm.coolwarm)
ax.view_init(elev=10,azim=100)
ax.set_xlabel('\n \n \n The frequency \n', fontsize = 15)
ax.set_ylabel('\n \n \n Gravitational', fontsize = 15)
ax.text2D(0.2, 0.95, 
          '\nThe interrelation of frequency,\n'
          'gravitational characteristics of \n'
          'a quark "d" at v = c.\n'
          'Graph#7.\n', 
          transform=ax.transAxes, fontsize = 16)
fig = plt.figure()
ax = Axes3D(fig)
X = ([unit5.frequency0[3], unit5.frequency0[4], unit5.frequency0[5]])
Y = ([unit4.grav0[3], unit4.grav0[4], unit4.grav0[5]])
X,Y = np.meshgrid(X,Y)
def f(x,y):
    return (sin(x) + cos(y))
ax.plot_surface(X,Y,f(X,Y), rstride=1, cstride=1, cmap=cm.coolwarm)
ax.view_init(elev=10,azim=60)
ax.set_xlabel('\n \n \n The frequency \n', fontsize = 15)
ax.set_ylabel('\n \n \n Gravitational', fontsize = 15)

ax.text2D(0.2, 0.95, 
          '\nThe interrelation of frequency,\n'
          'gravitational characteristics of \n'
          'a new quark "d" at v = c.\n'
          'Graph#8.\n', 
          transform=ax.transAxes, fontsize = 16)
fig = plt.figure()
ax = Axes3D(fig)
X = ([unit3.electro1[3], unit3.electro1[4], unit3.electro1[5]])
Y = ([unit4.grav1[3], unit4.grav1[4], unit4.grav1[5]])
X,Y = np.meshgrid(X,Y)
def f(x,y):
    return (sin(x) + cos(y))
ax.plot_surface(X,Y,f(X,Y), rstride=1, cstride=1, cmap=cm.coolwarm)
ax.view_init(elev=75,azim=65)
ax.set_xlabel('\n \n \n Electromagnetic \n', fontsize = 15)
ax.set_ylabel('\n \n \n Gravitational', fontsize = 15)
ax.text2D(0.2, 0.95, 
          '\nThe interrelation of electromagnetic,\n'
          'gravitational characteristics of \n'
          'a quark "d" at v = c.\n'
          'Graph#9.\n', 
          transform=ax.transAxes, fontsize = 16)
fig = plt.figure()
ax = Axes3D(fig)
X = ([unit3.electro0[3], unit3.electro0[4], unit3.electro0[5]])
Y = ([unit4.grav0[3], unit4.grav0[4], unit4.grav0[5]])
X,Y = np.meshgrid(X,Y)
def f(x,y):
    return (sin(x) + cos(y))
ax.plot_surface(X,Y,f(X,Y), rstride=1, cstride=1, cmap=cm.coolwarm)
ax.view_init(elev=75,azim=65)
ax.set_xlabel('\n \n \n Electromagnetic \n', fontsize = 15)
ax.set_ylabel('\n \n \n Gravitational', fontsize = 15)
ax.text2D(0.2, 0.95, 
          '\nThe interrelation of electromagnetic,\n'
          'gravitational characteristics of \n'
          'a new quark "d" at v = c.\n'
          'Graph#10.\n', 
          transform=ax.transAxes, fontsize = 16)
plt.show()

# The interrelation of frequency, electromagnetic 
# characteristics of a proton, neutron and their twins 3D graph.
fig = plt.figure()
ax = Axes3D(fig)
X = ([unit10.frequency02[0], unit10.frequency02[1], unit10.frequency02[2],
      unit10.frequency02[3], unit10.frequency02[4], unit10.frequency02[5], 
      unit10.frequency02[6], unit10.frequency02[7], unit10.frequency02[8]])
Y = ([unit8.electro01[0], unit8.electro01[1], unit8.electro01[2], 
      unit8.electro01[3], unit8.electro01[4], unit8.electro01[5], 
      unit8.electro01[6], unit8.electro01[7], unit8.electro01[8]])
X,Y = np.meshgrid(X,Y)
def f(x,y):
    return (sin(x) + cos(y))
ax.plot_surface(X,Y,f(X,Y), rstride=1, cstride=1, cmap=plt.cm.hot)
ax.view_init(elev=20,azim=85)
ax.set_xlabel('\n \n \n The frequency \n', fontsize = 15)
ax.set_ylabel('\n \n \n Electromagnetic', fontsize = 15)
ax.text2D(0.2, 0.95, 
          'The interrelation of frequency,\n'
          'electromagnetic characteristics\n'
          'of a new neutron at v = c \n'
          'by shells. Graph#11.\n', 
          transform=ax.transAxes, fontsize = 16)
fig = plt.figure()
ax = Axes3D(fig)
X = ([unit10.frequency03[0], unit10.frequency03[1], unit10.frequency03[2], 
      unit10.frequency03[3], unit10.frequency03[4], unit10.frequency03[5], 
      unit10.frequency03[6], unit10.frequency03[7], unit10.frequency03[8]])
Y = ([unit8.electro03[0], unit8.electro03[1], unit8.electro03[2], 
      unit8.electro03[3], unit8.electro03[4], unit8.electro03[5], 
      unit8.electro03[6], unit8.electro03[7], unit8.electro03[8]])
X,Y = np.meshgrid(X,Y)
def f(x,y):
    return (sin(x) + cos(y))
ax.plot_surface(X,Y,f(X,Y), rstride=1, cstride=1, cmap=plt.cm.hot)
ax.view_init(elev=20,azim=85)
ax.set_xlabel('\n \n \n The frequency \n', fontsize = 15)
ax.set_ylabel('\n \n \n Electromagnetic', fontsize = 15)
ax.text2D(0.2, 0.95, 
          'The interrelation of frequency,\n'
          'electromagnetic characteristics\n'
          'of a neutron at v = c \n'
          'by shells. Graph#12.\n', 
          transform=ax.transAxes, fontsize = 16)
fig = plt.figure()
ax = Axes3D(fig)
X = ([unit10.frequency05[0], unit10.frequency05[1], unit10.frequency05[2], 
      unit10.frequency05[3], unit10.frequency05[4], unit10.frequency05[5], 
      unit10.frequency05[6], unit10.frequency05[7], unit10.frequency05[8]])
Y = ([unit8.electro04[0], unit8.electro04[1], unit8.electro04[2], 
      unit8.electro04[3], unit8.electro04[4], unit8.electro04[5], 
      unit8.electro04[6], unit8.electro04[7], unit8.electro04[8]])
X,Y = np.meshgrid(X,Y)
def f(x,y):
    return (sin(x) + cos(y))
ax.plot_surface(X,Y,f(X,Y), rstride=1, cstride=1, cmap=plt.cm.hot)
ax.view_init(elev=20,azim=85)
ax.set_xlabel('\n \n \n The frequency \n', fontsize = 15)
ax.set_ylabel('\n \n \n Electromagnetic', fontsize = 15)
ax.text2D(0.2, 0.95, 
          'The interrelation of frequency,\n'
          'electromagnetic characteristics\n'
          'of a proton at v = c \n'
          'by shells. Graph#13.\n', 
          transform=ax.transAxes, fontsize = 16)
fig = plt.figure()
ax = Axes3D(fig)
X = ([unit10.frequency04[0], unit10.frequency04[1], unit10.frequency04[2], 
      unit10.frequency04[3], unit10.frequency04[4], unit10.frequency04[5], 
      unit10.frequency04[6], unit10.frequency04[7], unit10.frequency04[8]])
Y = ([unit8.electro02[0], unit8.electro02[1], unit8.electro02[2], 
      unit8.electro02[3], unit8.electro02[4], unit8.electro02[5], 
      unit8.electro02[6], unit8.electro02[7], unit8.electro02[8]])
X,Y = np.meshgrid(X,Y)
def f(x,y):
    return (sin(x) + cos(y))
ax.plot_surface(X,Y,f(X,Y), rstride=1, cstride=1, cmap=plt.cm.hot)
ax.view_init(elev=20,azim=85)
ax.set_xlabel('\n \n \n The frequency \n', fontsize = 15)
ax.set_ylabel('\n \n \n Electromagnetic', fontsize = 15)
ax.text2D(0.2, 0.95, 
          'The interrelation of frequency,\n'
          'electromagnetic characteristics\n'
          'of a new proton at v = c \n'
          'by shells. Graph#14.\n', 
          transform=ax.transAxes, fontsize = 16)
plt.show()

# The interrelation of frequency, electromagnetic, wavelength 
# characteristics of a proton, 3D graph.
fig = plt.figure(figsize=plt.figaspect(0.3))
ax = fig.add_subplot(1, 2, 1, projection='3d')
# proton
Xpu = ([unit10.frequency05[0], unit10.frequency05[1], 
        unit10.frequency05[2], unit10.frequency05[3], unit10.frequency05[4], 
        unit10.frequency05[5], unit10.frequency05[6], unit10.frequency05[7], 
        unit10.frequency05[8]])
Ypu = ([unit8.electro15[0], unit8.electro15[1], 
        unit8.electro15[2], unit8.electro15[3], unit8.electro15[4], 
        unit8.electro15[5], unit8.electro15[6], unit8.electro15[7], 
        unit8.electro15[8]])
Zpu = ([unit7.wave31[0], unit7.wave31[1], unit7.wave31[2], 
        unit7.wave31[3], unit7.wave31[4], unit7.wave31[5], 
        unit7.wave31[6], unit7.wave31[7], unit7.wave31[8]])

ax.plot(Xpu,Ypu,Zpu) 

ax.set_xlabel('\n \n \n The frequency \n', fontsize = 15)
ax.set_zlabel('\n \n \n Electromagnetic', fontsize = 15)
ax.set_ylabel('\n \n \n Wavelength', fontsize = 15)

ax.text2D(0.2, 0.95, 
          'The interrelation of frequency,\n'
          'electromagnetic, wavelength\n'
          'characteristics of a proton\n'
          'at v = c by shells. Graph#15.\n', 
          transform=ax.transAxes, fontsize = 16)

# The interrelation of frequency, electromagnetic, wavelength 
# characteristics of a new proton, 3D graph.

ax = fig.add_subplot(1, 2, 2, projection='3d')

Xku = ([unit10.frequency04[0], unit10.frequency04[1], 
        unit10.frequency04[2], unit10.frequency04[3], unit10.frequency04[4], 
        unit10.frequency04[5], unit10.frequency04[6], unit10.frequency04[7], 
        unit10.frequency04[8]])
Yku = ([unit8.electro14[0], unit8.electro14[1], 
        unit8.electro14[2], unit8.electro14[3], unit8.electro14[4], 
        unit8.electro14[5], unit8.electro14[6], unit8.electro14[7], 
        unit8.electro14[8]])
Zku = ([unit7.wave21[0], unit7.wave21[1], unit7.wave21[2], 
        unit7.wave21[3], unit7.wave21[4], unit7.wave21[5], 
        unit7.wave21[6], unit7.wave21[7], unit7.wave21[8]])

ax.plot(Xku,Yku,Zku) 

ax.set_xlabel('\n \n \n The frequency \n', fontsize = 15)
ax.set_zlabel('\n \n \n Electromagnetic', fontsize = 15)
ax.set_ylabel('\n \n \n Wavelength', fontsize = 15)

ax.text2D(0.2, 0.95, 
          'The interrelation of frequency,\n'
          'electromagnetic, wavelength\n'
          'characteristics of a new proton \n'
          'at v = c by shells. Graph#16.\n', 
          transform=ax.transAxes, fontsize = 16)

# The interrelation of frequency, electromagnetic, wavelength 
# characteristics of a neutron, 3D graph.
fig = plt.figure(figsize=plt.figaspect(0.3))
ax = fig.add_subplot(1, 2, 1, projection='3d')
# neutron
Xpu = ([unit10.frequency03[0], unit10.frequency03[1], 
        unit10.frequency03[2], unit10.frequency03[3], unit10.frequency03[4], 
        unit10.frequency03[5], unit10.frequency03[6], unit10.frequency03[7], 
        unit10.frequency03[8]])
Ypu = ([unit8.electro03[0], unit8.electro03[1], 
        unit8.electro03[2], unit8.electro03[3], unit8.electro03[4], 
        unit8.electro03[5], unit8.electro03[6], unit8.electro03[7], 
        unit8.electro03[8]])
Zpu = ([unit7.wave11[0], unit7.wave11[1], unit7.wave11[2], 
        unit7.wave11[3], unit7.wave11[4], unit7.wave11[5], 
        unit7.wave11[6], unit7.wave11[7], unit7.wave11[8]])

ax.plot(Xpu,Ypu,Zpu) 

ax.set_xlabel('\n \n \n The frequency \n', fontsize = 15)
ax.set_zlabel('\n \n \n Electromagnetic', fontsize = 15)
ax.set_ylabel('\n \n \n Wavelength', fontsize = 15)

ax.text2D(0.2, 0.95, 
          'The interrelation of frequency,\n'
          'electromagnetic, wavelength\n'
          'characteristics of a neutron\n'
          'at v = c by shells. Graph#17.\n', 
          transform=ax.transAxes, fontsize = 16)

# The interrelation of frequency, electromagnetic, wavelength 
# characteristics of a new neutron, 3D graph.

ax = fig.add_subplot(1, 2, 2, projection='3d')

Xku = ([unit10.frequency02[0], unit10.frequency02[1], 
        unit10.frequency02[2], unit10.frequency02[3], unit10.frequency02[4], 
        unit10.frequency02[5], unit10.frequency02[6], unit10.frequency02[7], 
        unit10.frequency02[8]])
Yku = ([unit8.electro01[0], unit8.electro01[1], 
        unit8.electro01[2], unit8.electro01[3], unit8.electro01[4], 
        unit8.electro01[5], unit8.electro01[6], unit8.electro01[7], 
        unit8.electro01[8]])
Zku = ([unit7.wave01[0], unit7.wave01[1], unit7.wave01[2], 
        unit7.wave01[3], unit7.wave01[4], unit7.wave01[5], 
        unit7.wave01[6], unit7.wave01[7], unit7.wave01[8]])

ax.plot(Xku,Yku,Zku) 

ax.set_xlabel('\n \n \n The frequency \n', fontsize = 15)
ax.set_zlabel('\n \n \n Electromagnetic', fontsize = 15)
ax.set_ylabel('\n \n \n Wavelength', fontsize = 15)

ax.text2D(0.2, 0.95, 
          'The interrelation of frequency,\n'
          'electromagnetic, wavelength\n'
          'characteristics of a new neutron \n'
          'at v = c by shells. Graph#18.\n', 
          transform=ax.transAxes, fontsize = 16)

# neutron
# The more N, then to generate a dense grid.
N = 20
# Random data is generated in class 11.
volume = np.random.rand(N, N, N)
# The x, y, and z coordinate arrays taken from class 11.  
x = unit11.frequenn
y = unit11.wavelengthn
z = unit11.v
x, y, z = np.broadcast_arrays(x, y, z)
# The volumetric data turn into an RGB array.
c = np.tile(volume.ravel()[:, None], [1, 3])
fig = plt.figure(figsize=(12,12))
ax = fig.gca(projection='3d')
ax.scatter(x.ravel(),
           y.ravel(),
           z.ravel(),
           c=c)
ax.set_xlabel('\n \n \n Frequency \n', fontsize = 15)
ax.set_ylabel('\n \n \n Wavelength\n', fontsize = 15)
ax.set_zlabel('\n \n \n Velocity\n', fontsize = 15)
ax.text2D(0.2, 0.95, 
          '\nThe interrelation of frequency,wavelength and\n'
          'velocity of a neutron. Graph#19.', 
          transform=ax.transAxes, fontsize = 16)
fig = plt.figure()
# new neutron
# The more N, then to generate a dense grid.
N = 20
# Random data is generated in class 11.
volume = np.random.rand(N, N, N)
# The x, y, and z coordinate arrays taken from class 11.  
y = unit11.frequennn
x = unit11.wavelengthnn
z = unit11.v
x, y, z = np.broadcast_arrays(x, y, z)
# The volumetric data turn into an RGB array.
c = np.tile(volume.ravel()[:, None], [1, 3])
fig = plt.figure(figsize=(12,12))
ax = fig.gca(projection='3d')
ax.scatter(x.ravel(),
           y.ravel(),
           z.ravel(),
           c=c)
ax.set_xlabel('\n \n \n Wavelength \n', fontsize = 15)
ax.set_ylabel('\n \n \n Frequency\n', fontsize = 15)
ax.set_zlabel('\n \n \n Velocity\n', fontsize = 15)
ax.text2D(0.2, 0.95, 
          '\nThe interrelation of frequency,wavelength and\n'
          'velocity of a new neutron. Graph#20.', 
          transform=ax.transAxes, fontsize = 16)
fig = plt.figure()
# proton
# The more N, then to generate a dense grid.
N = 20
# Random data is generated in class 11.
volume = np.random.rand(N, N, N)
# The x, y, and z coordinate arrays taken from class 11.  
z = unit11.frequenp
y = unit11.wavelengthp
x = unit11.v
x, y, z = np.broadcast_arrays(x, y, z)
# The volumetric data turn into an RGB array.
c = np.tile(volume.ravel()[:, None], [1, 3])
fig = plt.figure(figsize=(12,12))
ax = fig.gca(projection='3d')
ax.scatter(x.ravel(),
           y.ravel(),
           z.ravel(),
           c=c)
ax.set_xlabel('\n \n \n Velocity  \n', fontsize = 15)
ax.set_ylabel('\n \n \n Wavelength\n', fontsize = 15)
ax.set_zlabel('\n \n \n Frequency\n', fontsize = 15)
ax.text2D(0.2, 0.95, 
          '\nThe interrelation of frequency,wavelength and\n'
          'velocity of a proton. Graph#21.', 
          transform=ax.transAxes, fontsize = 16)
fig = plt.figure() 
# new proton
# The more N, then to generate a dense grid.
N = 20
# Random data is generated in class 11.
volume = np.random.rand(N, N, N)
# The x, y, and z coordinate arrays taken from class 11.  
x = unit11.frequennp
z = unit11.wavelengthnp
y = unit11.v
x, y, z = np.broadcast_arrays(x, y, z)
# The volumetric data turn into an RGB array.
c = np.tile(volume.ravel()[:, None], [1, 3])
fig = plt.figure(figsize=(12,12))
ax = fig.gca(projection='3d')
ax.scatter(x.ravel(),
           y.ravel(),
           z.ravel(),
           c=c)
ax.set_xlabel('\n \n \n Frequency \n', fontsize = 15)
ax.set_ylabel('\n \n \n Velocity \n', fontsize = 15)
ax.set_zlabel('\n \n \n Wavelength \n', fontsize = 15)
ax.text2D(0.2, 0.95, 
          '\nThe interrelation of frequency,wavelength and\n'
          'velocity of a new proton. Graph#22.', 
          transform=ax.transAxes, fontsize = 16)
fig = plt.figure() 
print('\n Initial conditions. Table 1.\n')
print(table1)
print('\nData obtained during the implementation of the algorithm. Table 2.\n')
print(table2)        
print('\nValues of electric charge, mass, volume by shells for the new proton.\n'
     'Table 3.')
print(table3)
print('\nValues of electric charge, mass, volume by shells for the proton.\n'
     'Table 4.')
print(table4)
print('\nValues of electric charge, mass, volume by shells for the new neutron.\n'
     'Table 5.')
print(table5)
print('\nValues of electric charge, mass, volume by shells for the neutron.\n'
     'Table 6.')
print(table6)
print('\n Mass, electric charge and volume of the \n' 
      'smallest particle carrying electric charge.\n'
      'Table 7.')
print(table7)

