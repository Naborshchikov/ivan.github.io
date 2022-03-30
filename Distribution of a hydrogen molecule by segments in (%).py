#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""Distribution of a hydrogen molecule by segments in (%)"""
# The app serves for simplified visualization.
# For scientific analysis, use data from program version 6.2.
# For example MolecularH2_1_Present_1e and MolecularH2_1_Present_2e are two different segments.

import matplotlib.pyplot as plt
import matplotlib as mpl

# Distribution of the volume of a hydrogen molecule by segments in (%)
# for the first molecule

fig = plt.figure(figsize=plt.figaspect(0.3))
data_names = ['Present, volume "+" \n', 
              'Past, volume "-"  \n',
              'Future, volume "-"']
data_values = [(unit02.MolecularH2_1_Present_1e[2] + unit02.MolecularH2_1_Present_2e[2]) * 10e46,
               -(unit02.MolecularH2_1_Past_1e[2]) * 10e46, 
               -(unit02.MolecularH2_1_Future_1e[2] + unit02.MolecularH2_1_Future_2e[2]) * 10e46]
ax = fig.add_subplot(2, 3, 1)

mpl.rcParams.update({'font.size': 14})

plt.title('Distribution for hydrogen molecule #1 by segments in %: \n \n'
          'volume, Graph#1')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the electric charge of a hydrogen molecule by segments in (%) 
# for the first molecule

data_names = ['Present, charge "+" \n', 
              'Past, charge "-" \n',
              'Future, charge "-"']
data_values = [(unit02.MolecularH2_1_Present_1e[0] + unit02.MolecularH2_1_Present_2e[0]) * 10e20,
               -(unit02.MolecularH2_1_Past_1e[0]) * 10e20, 
               -(unit02.MolecularH2_1_Future_1e[0] + unit02.MolecularH2_1_Future_2e[0]) * 10e20]
ax = fig.add_subplot(2, 3, 2)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n electric charge, Graph#2')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the mass of a hydrogen molecule by segments in (%) 
# for the first molecule

data_names = ['Present time \n', 
              'Past time \n',
              'Future time']
data_values = [(unit02.MolecularH2_1_Present_1e[1] + unit02.MolecularH2_1_Present_2e[1]) * 10e29,
               (unit02.MolecularH2_1_Past_1e[1]) * 10e29, 
               (unit02.MolecularH2_1_Future_1e[1] + unit02.MolecularH2_1_Future_2e[1]) * 10e29]
ax = fig.add_subplot(2, 3, 3)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n mass, Graph#3')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the volume of a hydrogen molecule by segments in (%)
# for the 2 molecule

data_names = ['Present, volume "+" \n', 
              'Past, volume "-"  \n',
              'Future, volume "-"']
data_values = [(unit02.MolecularH2_2_Present_1e[2] + unit02.MolecularH2_2_Present_2e[2]) * 10e46,
               -(unit02.MolecularH2_2_Past_1e[2]) * 10e46, 
               -(unit02.MolecularH2_2_Future_1e[2] + unit02.MolecularH2_2_Future_2e[2]) * 10e46]
ax = fig.add_subplot(2, 3, 4)
mpl.rcParams.update({'font.size': 14})

plt.title('Distribution for hydrogen molecule #2 by segments in %: \n \n'
          'volume, Graph#4')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the electric charge of a hydrogen molecule by segments in (%) 
# for the 2 molecule

data_names = ['Present, charge "+" \n', 
              'Past, charge "-" \n',
              'Future, charge "-"']
data_values = [(unit02.MolecularH2_2_Present_1e[0] + unit02.MolecularH2_2_Present_2e[0]) * 10e20,
               -(unit02.MolecularH2_2_Past_1e[0]) * 10e20, 
               -(unit02.MolecularH2_2_Future_1e[0] + unit02.MolecularH2_2_Future_2e[0]) * 10e20]
ax = fig.add_subplot(2, 3, 5)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n electric charge, Graph#5')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)] )
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the mass of a hydrogen molecule by segments in (%) 
# for the 2 molecule
data_names = ['Present time\n', 
              'Past time \n',
              'Future time']
data_values = [(unit02.MolecularH2_2_Present_1e[1] + unit02.MolecularH2_2_Present_2e[1]) * 10e29,
               (unit02.MolecularH2_2_Past_1e[1]) * 10e29, 
               (unit02.MolecularH2_2_Future_1e[1] + unit02.MolecularH2_2_Future_2e[1]) * 10e29]

ax = fig.add_subplot(2, 3, 6)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n mass, Graph#6')

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

# Distribution of the volume of a hydrogen molecule by segments in (%)
# for the 3 molecule

fig = plt.figure(figsize=plt.figaspect(0.3))
data_names = ['Present, volume "+" \n', 
              'Past, volume "-"  \n',
              'Future, volume "-"']
data_values = [(unit02.MolecularH2_3_Present_1e[2] + unit02.MolecularH2_3_Present_2e[2]) * 10e46,
               -(unit02.MolecularH2_3_Past_1e[2]) * 10e46, 
               -(unit02.MolecularH2_3_Future_1e[2]) * 10e46]
ax = fig.add_subplot(2, 3, 1)
mpl.rcParams.update({'font.size': 14})

plt.title('Distribution for hydrogen molecule #3 by segments in %: \n \n'
          'volume, Graph#7')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the electric charge of a hydrogen molecule by segments in (%) 
# for the 3 molecule

data_names = ['Present, charge "+" \n', 
              'Past, charge "-" \n',
              'Future, charge "-"']
data_values = [(unit02.MolecularH2_3_Present_1e[0] + unit02.MolecularH2_3_Present_2e[0]) * 10e20,
               -(unit02.MolecularH2_3_Past_1e[0]) * 10e20, 
               -(unit02.MolecularH2_3_Future_1e[0]) * 10e20]
ax = fig.add_subplot(2, 3, 2)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n electric charge, Graph#8')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the mass of a hydrogen molecule by segments in (%) 
# for the 3 molecule

data_names = ['Present time \n', 
              'Past time \n',
              'Future time']
data_values = [(unit02.MolecularH2_3_Present_1e[1] + unit02.MolecularH2_3_Present_2e[1]) * 10e29,
               (unit02.MolecularH2_3_Past_1e[1]) * 10e29, 
               (unit02.MolecularH2_3_Future_1e[1]) * 10e29]
ax = fig.add_subplot(2, 3, 3)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n mass, Graph#9')
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

# Distribution of the volume of a hydrogen molecule by segments in (%)
# for the 4 molecule

fig = plt.figure(figsize=plt.figaspect(0.3))
data_names = ['Present, volume "+" \n', 
              'Past, volume "-"  \n',
              'Future, volume "-"']
data_values = [(unit02.MolecularH2_4_Present_1e[2] + 
                unit02.MolecularH2_4_Present_2e[2]) * 10e46,
               -(unit02.MolecularH2_4_Past_1e[2]) * 10e46, 
               -(unit02.MolecularH2_4_Future_1e[2]) * 10e46]
ax = fig.add_subplot(2, 3, 1)

mpl.rcParams.update({'font.size': 14})

plt.title('Distribution for hydrogen molecule #4 by segments in %: \n \n'
          'volume, Graph#10')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the electric charge of a hydrogen molecule by segments in (%) 
# for the 4 molecule

data_names = ['Present, charge "+" \n', 
              'Past, charge "+" \n',
              'Future, charge "+"']
data_values = [(unit02.MolecularH2_4_Present_1e[0] + 
                unit02.MolecularH2_4_Present_2e[0]) * 10e20,
               (unit02.MolecularH2_4_Past_1e[0]) * 10e20, 
               (unit02.MolecularH2_4_Future_1e[0]) * 10e20]
ax = fig.add_subplot(2, 3, 2)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n electric charge, Graph#11')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the mass of a hydrogen molecule by segments in (%) 
# for the 4 molecule

data_names = ['Present time \n', 
              'Past time \n',
              'Future time']
data_values = [(unit02.MolecularH2_4_Present_1e[1] + 
                unit02.MolecularH2_4_Present_2e[1]) * 10e29,
               (unit02.MolecularH2_4_Past_1e[1]) * 10e29, 
               (unit02.MolecularH2_4_Future_1e[1]) * 10e29]
ax = fig.add_subplot(2, 3, 3)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n mass, Graph#12')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

plt.subplots_adjust(left=0.3,
                    bottom=0.001, 
                    right=1.0, 
                    top=2.5, 
                    wspace=0.1, 
                    hspace=0.1)

# Distribution of the volume of a hydrogen molecule by segments in (%)
# for the 5 molecule

data_names = ['Present, volume "+" \n', 
              'Past, volume "+"  \n',
              'Future, volume "-"']
data_values = [(unit02.MolecularH2_5_Present_1e[2] + 
                unit02.MolecularH2_5_Present_2e[2]) * 10e46,
               (unit02.MolecularH2_5_Past_1e[2] + 
                unit02.MolecularH2_5_Past_2e[2]) * 10e46, 
               -(unit02.MolecularH2_5_Future_1e[2]) * 10e46]
ax = fig.add_subplot(2, 3, 4)
mpl.rcParams.update({'font.size': 14})

plt.title('Distribution for hydrogen molecule #5 by segments in % \n'
          '(common part of protons in a molecule: 3.41797682e-19 C,\n'
          '1.16000696e-27 kg, 4.28213181e-45 m.cube): \n'
          'volume, Graph#13')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the electric charge of a hydrogen molecule by segments in (%) 
# for the 5 molecule

data_names = ['Present, charge "+" \n', 
              'Past, charge "+" \n',
              'Future, charge "+"']
data_values = [(unit02.MolecularH2_5_Present_1e[0] + 
                 unit02.MolecularH2_5_Present_2e[0]) * 10e20,
               (unit02.MolecularH2_5_Past_1e[0] + 
                unit02.MolecularH2_5_Past_2e[0]) * 10e20, 
               (unit02.MolecularH2_5_Future_1e[0]) * 10e20]
ax = fig.add_subplot(2, 3, 5)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n electric charge, Graph#14')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)] )
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the mass of a hydrogen molecule by segments in (%) 
# for the 5 molecule
data_names = ['Present time\n', 
              'Past time \n',
              'Future time']
data_values = [(unit02.MolecularH2_5_Present_1e[1] + 
                unit02.MolecularH2_5_Present_2e[1]) * 10e29,
               (unit02.MolecularH2_5_Past_1e[1] + 
                unit02.MolecularH2_5_Past_2e[1]) * 10e29, 
               (unit02.MolecularH2_5_Future_1e[1]) * 10e29]

ax = fig.add_subplot(2, 3, 6)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n mass, Graph#15')

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

# Distribution of the volume of a hydrogen molecule by segments in (%)
# for the 6 molecule

fig = plt.figure(figsize=plt.figaspect(0.3))
data_names = ['Present, volume "+" \n', 
              'Past, volume "+"  \n',
              'Future, volume "-"']
data_values = [(unit02.MolecularH2_6_Present_1e[2]) * 10e46,
               (unit02.MolecularH2_6_Past_1e[2] + 
                 unit02.MolecularH2_6_Past_2e[2] + 
                 unit02.MolecularH2_6_Past_3e[2]) * 10e46, 
               -(unit02.MolecularH2_6_Future_1e[2]) * 10e46]
ax = fig.add_subplot(2, 3, 1)
mpl.rcParams.update({'font.size': 14})

plt.title('Distribution for hydrogen molecule #6 by segments in % \n'
          '(common part of protons in a molecule:-2.42996789e-19 C,\n'
          '1.18744445e-27 kg, 2.34472000e-45 m.cube): \n \n'
          'volume, Graph#16')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the electric charge of a hydrogen molecule by segments in (%) 
# for the 6 molecule

data_names = ['Present, charge "-" \n', 
              'Past, charge "+" \n',
              'Future, charge "+"']
data_values = [-(unit02.MolecularH2_6_Present_1e[0]) * 10e20,
               (unit02.MolecularH2_6_Past_1e[0] + 
                 unit02.MolecularH2_6_Past_2e[0] + 
                 unit02.MolecularH2_6_Past_3e[0]) * 10e20, 
               (unit02.MolecularH2_6_Future_1e[0]) * 10e20]
ax = fig.add_subplot(2, 3, 2)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n electric charge, Graph#17')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the mass of a hydrogen molecule by segments in (%) 
# for the 6 molecule

data_names = ['Present time \n', 
              'Past time \n',
              'Future time']
data_values = [(unit02.MolecularH2_6_Present_1e[1]) * 10e29,
               (unit02.MolecularH2_6_Past_1e[1] + 
                unit02.MolecularH2_6_Past_2e[1] + 
                unit02.MolecularH2_6_Past_3e[1]) * 10e29, 
               (unit02.MolecularH2_6_Future_1e[1]) * 10e29]
ax = fig.add_subplot(2, 3, 3)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n mass, Graph#18')
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

# Distribution of the volume of a hydrogen molecule by segments in (%)
# for the 7 molecule

fig = plt.figure(figsize=plt.figaspect(0.3))
data_names = ['Present, volume "+" \n', 
              'Past, volume "+"  \n',
              'Future, volume "-"']
data_values = [(unit02.MolecularH2_7_Present_1e[2]) * 10e46,
               (unit02.MolecularH2_7_Past_1e[2] + 
                 unit02.MolecularH2_7_Past_2e[2] + 
                 unit02.MolecularH2_7_Past_3e[2]) * 10e46, 
               -(unit02.MolecularH2_7_Future_1e[2]) * 10e46]
ax = fig.add_subplot(2, 3, 1)

mpl.rcParams.update({'font.size': 14})

plt.title('Distribution for hydrogen molecule #7 by segments in %: \n \n'
          'volume, Graph#19')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the electric charge of a hydrogen molecule by segments in (%) 
# for the 7 molecule

data_names = ['Present, charge "+" \n', 
              'Past, charge "+" \n',
              'Future, charge "-"']
data_values = [(unit02.MolecularH2_7_Present_1e[0]) * 10e20,
               (unit02.MolecularH2_7_Past_1e[0] + 
                unit02.MolecularH2_7_Past_2e[0] + 
                unit02.MolecularH2_7_Past_3e[0]) * 10e20, 
               -(unit02.MolecularH2_7_Future_1e[0]) * 10e20]
ax = fig.add_subplot(2, 3, 2)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n electric charge, Graph#20')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the mass of a hydrogen molecule by segments in (%) 
# for the 7 molecule

data_names = ['Present time \n', 
              'Past time \n',
              'Future time']
data_values = [(unit02.MolecularH2_7_Present_1e[1]) * 10e29,
               (unit02.MolecularH2_7_Past_1e[1] + 
                unit02.MolecularH2_7_Past_2e[1] + 
                unit02.MolecularH2_7_Past_3e[1]) * 10e29, 
               (unit02.MolecularH2_7_Future_1e[1]) * 10e29]
ax = fig.add_subplot(2, 3, 3)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n mass, Graph#21')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the volume of a hydrogen molecule by segments in (%)
# for the 8 molecule

data_names = ['Present, volume "+" \n', 
              'Past, volume "+"  \n',
              'Future, volume "-"']
data_values = [(unit02.MolecularH2_8_Present_1e[2] + 
                unit02.MolecularH2_8_Present_2e[2]) * 10e46,
               (unit02.MolecularH2_8_Past_1e[2] + 
                unit02.MolecularH2_8_Past_2e[2]) * 10e46, 
               -(unit02.MolecularH2_8_Future_1e[2]) * 10e46]
ax = fig.add_subplot(2, 3, 4)
mpl.rcParams.update({'font.size': 14})

plt.title('Distribution for hydrogen molecule #8 by segments in %: \n \n'
          'volume, Graph#22')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the electric charge of a hydrogen molecule by segments in (%) 
# for the 8 molecule

data_names = ['Present, charge "+" \n', 
              'Past, charge "-" \n',
              'Future, charge "+"']
data_values = [(unit02.MolecularH2_8_Present_1e[0] + 
                 unit02.MolecularH2_8_Present_2e[0]) * 10e20,
               -(unit02.MolecularH2_8_Past_1e[0] + 
                unit02.MolecularH2_8_Past_2e[0]) * 10e20, 
               (unit02.MolecularH2_8_Future_1e[0]) * 10e20]
ax = fig.add_subplot(2, 3, 5)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n electric charge, Graph#23')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)] )
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the mass of a hydrogen molecule by segments in (%) 
# for the 8 molecule
data_names = ['Present time\n', 
              'Past time \n',
              'Future time']
data_values = [(unit02.MolecularH2_8_Present_1e[1] + 
                unit02.MolecularH2_8_Present_2e[1]) * 10e29,
               (unit02.MolecularH2_8_Past_1e[1] + 
                unit02.MolecularH2_8_Past_2e[1]) * 10e29, 
               (unit02.MolecularH2_8_Future_1e[1]) * 10e29]

ax = fig.add_subplot(2, 3, 6)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n mass, Graph#24')

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

