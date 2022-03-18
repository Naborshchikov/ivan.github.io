#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from prettytable import PrettyTable
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
print(table1)

