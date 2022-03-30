#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""Distribution of a He2_1 by segments in (%)""" 
# The app serves for simplified visualization.
# For scientific analysis, use data from program version 6.2.
# For example He2_2_Past_1_e and He2_2_Past_2_e are two different segments.

import matplotlib.pyplot as plt
import matplotlib as mpl
He2_1_Present_1_1e_sum = np.sum(unit03.He2_1_Present_1_1e, axis=0, dtype=None, out=None)
He2_1_Present_2_1e_sum = np.sum(unit03.He2_1_Present_2_1e, axis=0, dtype=None, out=None)
He2_1_Future_2_e_sum = np.sum(unit03.He2_1_Future_2_e, axis=0, dtype=None, out=None)
He2_1_Future_3_e_sum = np.sum(unit03.He2_1_Future_3_e, axis=0, dtype=None, out=None)

He2_2_Present_1_2e_sum = np.sum(unit03.He2_2_Present_1_2e, axis=0, dtype=None, out=None)
He2_2_Present_3e_sum = np.sum(unit03.He2_2_Present_3e, axis=0, dtype=None, out=None)
He2_2_Present_4e_sum = np.sum(unit03.He2_2_Present_4e, axis=0, dtype=None, out=None)
He2_2_Future_1e_sum = np.sum(unit03.He2_2_Future_1e, axis=0, dtype=None, out=None)

fig = plt.figure(figsize=plt.figaspect(0.3))
data_names = ['Present, volume "+" \n', 
              'Past, volume "-"  \n',
              'Future, volume "-"']
data_values = [(He2_1_Present_1_1e_sum[2] + He2_1_Present_2_1e_sum[2]) * 10e46,
               -(unit03.He2_1_Past_e[0, 2]) * 10e46, 
               -(unit03.He2_1_Future_1_e[0, 2] + He2_1_Future_2_e_sum[2] +
                 He2_1_Future_3_e_sum[2]) * 10e46]
ax = fig.add_subplot(2, 3, 1)

mpl.rcParams.update({'font.size': 14})

plt.title('Distribution for He2_1 by segments in %: \n \n'
          'volume, Graph#1')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the electric charge of a He2_1 by segments in (%) 

data_names = ['Present, charge "+" \n', 
              'Past, charge "-" \n',
              'Future, charge "-"']
data_values = [(He2_1_Present_1_1e_sum[0] + He2_1_Present_2_1e_sum[0]) * 10e20,
               -(unit03.He2_1_Past_e[0, 0]) * 10e20, 
               -(unit03.He2_1_Future_1_e[0, 0] + He2_1_Future_2_e_sum[0] +
                 He2_1_Future_3_e_sum[0]) * 10e20]
ax = fig.add_subplot(2, 3, 2)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n electric charge, Graph#2')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the mass of a He2_1 by segments in (%) 

data_names = ['Present time \n', 
              'Past time \n',
              'Future time']
data_values = [(He2_1_Present_1_1e_sum[1] + 
                He2_1_Present_2_1e_sum[1]) * 10e29,
               (unit03.He2_1_Past_e[0, 1]) * 10e29, 
               (unit03.He2_1_Future_1_e[0, 1] + He2_1_Future_2_e_sum[1] +
                 He2_1_Future_3_e_sum[1]) * 10e29]
ax = fig.add_subplot(2, 3, 3)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n mass, Graph#3')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the volume of a He2_2 by segments in (%)

data_names = ['Present, volume "+" \n', 
              'Past, volume "+"  \n',
              'Future, volume "-"']
data_values = [(He2_2_Present_1_2e_sum[2] + unit03.He2_2_Present_2e[0, 2] + 
                He2_2_Present_3e_sum[2] + He2_2_Present_4e_sum[2]) * 10e46, 
               (np.array(unit03.He2_2_Past_1_e)[0, 2] + 
                unit03.He2_2_Past_2_e[0, 2]) * 10e46, 
               -(He2_2_Future_1e_sum[2]) * 10e46]
ax = fig.add_subplot(2, 3, 4)
mpl.rcParams.update({'font.size': 14})

plt.title('Distribution for He2_2 by segments in %: \n \n'
          'volume, Graph#4')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)])
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the electric charge of a He2_2 by segments in (%) 

data_names = ['Present, charge "+" \n', 
              'Past, charge "-" \n',
              'Future, charge "+"']
data_values = [(He2_2_Present_1_2e_sum[0] + unit03.He2_2_Present_2e[0, 0] + 
                He2_2_Present_3e_sum[0] + He2_2_Present_4e_sum[0]) * 10e20,
               -(np.array(unit03.He2_2_Past_1_e)[0, 0] + 
                unit03.He2_2_Past_2_e[0, 0]) * 10e20, 
               (He2_2_Future_1e_sum[0]) * 10e20]
ax = fig.add_subplot(2, 3, 5)
mpl.rcParams.update({'font.size': 14})

plt.title('\n \n electric charge, Graph#5')

xs = range(len(data_names))

plt.pie(data_values, autopct='%.1f', radius = 1.1,
        explode = [0.15] + [0 for _ in range(len(data_names) - 1)] )
plt.legend(bbox_to_anchor = (0.05, -0.4, 0.25, 0.25),
           loc = 'lower left', labels = data_names)

# Distribution of the mass of a He2_2 by segments in (%) 

data_names = ['Present time\n', 
              'Past time \n',
              'Future time']
data_values = [(He2_2_Present_1_2e_sum[1] + unit03.He2_2_Present_2e[0, 1] + 
                He2_2_Present_3e_sum[1] + He2_2_Present_4e_sum[1]) * 10e29,
               (np.array(unit03.He2_2_Past_1_e)[0, 1] + 
                unit03.He2_2_Past_2_e[0, 1]) * 10e29, 
               (He2_2_Future_1e_sum[1]) * 10e29]

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

