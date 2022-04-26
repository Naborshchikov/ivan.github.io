#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from math import *
class He2_properties():
    """The calculation is made for vacuum conditions"""
    
    def __init__ (self, divE_He2_1_Past, divE_He2_1_Present_1, divE_He2_1_Present_2, 
                  divE_He2_1_Future_1, divE_He2_1_Future_2, divE_He2_1_Future_3,                   
                  divE_He2_1_Past_e, divE_He2_1_Present_1_1e, divE_He2_1_Present_2_1e, 
                  divE_He2_1_Future_1_e, divE_He2_1_Future_2_e, divE_He2_1_Future_3_e,
                  
                  divE_He2_2_Past_1, divE_He2_2_Past_2, divE_He2_2_Present_1, 
                  divE_He2_2_Present_2, divE_He2_2_Present_3, divE_He2_2_Present_4, 
                  divE_He2_2_Future_1, divE_He2_2_Past_1_e, divE_He2_2_Past_2_e, 
                  divE_He2_2_Present_1_2e, divE_He2_2_Present_2e, divE_He2_2_Present_3e, 
                  divE_He2_2_Present_4e, divE_He2_2_Future_1e, T16, T17, T18):
        
        self.divE_He2_1_Past = divE_He2_1_Past        
        self.divE_He2_1_Present_1 = divE_He2_1_Present_1
        self.divE_He2_1_Present_2 = divE_He2_1_Present_2
        self.divE_He2_1_Future_1 = divE_He2_1_Future_1        
        self.divE_He2_1_Future_2 = divE_He2_1_Future_2
        self.divE_He2_1_Future_3 = divE_He2_1_Future_3
        
        self.divE_He2_1_Past_e = divE_He2_1_Past_e
        self.divE_He2_1_Present_1_1e = divE_He2_1_Present_1_1e
        self.divE_He2_1_Present_2_1e = divE_He2_1_Present_2_1e
        self.divE_He2_1_Future_1_e = divE_He2_1_Future_1_e
        self.divE_He2_1_Future_2_e = divE_He2_1_Future_2_e
        self.divE_He2_1_Future_3_e = divE_He2_1_Future_3_e
        
        self.divE_He2_2_Past_1 = divE_He2_2_Past_1
        self.divE_He2_2_Past_2 = divE_He2_2_Past_2
        self.divE_He2_2_Present_1 = divE_He2_2_Present_1
        self.divE_He2_2_Present_2 = divE_He2_2_Present_2
        self.divE_He2_2_Present_3 = divE_He2_2_Present_3
        self.divE_He2_2_Present_4 = divE_He2_2_Present_4
        self.divE_He2_2_Future_1 = divE_He2_2_Future_1
        
        self.divE_He2_2_Past_1_e = divE_He2_2_Past_1_e
        self.divE_He2_2_Past_2_e = divE_He2_2_Past_2_e
        self.divE_He2_2_Present_1_2e = divE_He2_2_Present_1_2e
        self.divE_He2_2_Present_2e = divE_He2_2_Present_2e
        self.divE_He2_2_Present_3e = divE_He2_2_Present_3e
        self.divE_He2_2_Present_4e = divE_He2_2_Present_4e
        self.divE_He2_2_Future_1e = divE_He2_2_Future_1e
        
        self.T16 = T16
        self.T17 = T17
        self.T18 = T18
        
        
# He2_1_Past  

#  ρ at the beginning of the symbol indicates the volume charge density
# divE is determined according to the first Maxwell equation

# Nuclear Communications - He2_1_Past
ρ_He2_1_Past = np.array(He2_1_Past)[:, 0]
D = np.array(He2_1_Past)[:, 2]
divE_He2_1_Past = ρ_He2_1_Past/unit.constantε0
for i, item in enumerate(D):
    ρ_He2_1_Past[i] /= D[i]
    divE_He2_1_Past[i] = ρ_He2_1_Past[i]/unit.constantε0

# He2_1_Present_1 

ρ_He2_1_Present_1 = np.array(He2_1_Present_1)[:, 0]
D1 = np.array(He2_1_Present_1)[:, 2]
divE_He2_1_Present_1 = ρ_He2_1_Present_1/unit.constantε0
for i, item in enumerate(D1):
    ρ_He2_1_Present_1[i] /= D1[i]
    divE_He2_1_Present_1[i] = ρ_He2_1_Present_1[i]/unit.constantε0

# He2_1_Present_2 

ρ_He2_1_Present_2 = np.array(He2_1_Present_2)[:, 0]
D2 = np.array(He2_1_Present_2)[:, 2]
divE_He2_1_Present_2 = ρ_He2_1_Present_2/unit.constantε0
for i, item in enumerate(D2):
    ρ_He2_1_Present_2[i] /= D2[i]
    divE_He2_1_Present_2[i] = ρ_He2_1_Present_2[i]/unit.constantε0
    

k1 = 1/(4 * unit.π * unit.constantε0)

# He2_1_Past & (He2_1_Present_1 and He2_1_Present_2) - 
# Assessing the Impact on Nuclear Communications Using Coulomb's Law

R_He2_1_Past = (3/4 * abs(He2_1_Past[0, 2])/unit.π) ** (1./3.)
R_He2_1_Present_1 = (3/4 * (sum(row[2] for row in He2_1_Present_1))/unit.π) ** (1./3.)
R_He2_1_Present_2 = (3/4 * (sum(row[2] for row in He2_1_Present_2))/unit.π) ** (1./3.)

Q_He2_1_Past = He2_1_Past[0, 0]
Q_He2_1_Present_1 = sum(row[0] for row in He2_1_Present_1)
Q_He2_1_Present_2 = sum(row[0] for row in He2_1_Present_2)

if (He2_1_Present_1[0][2] > 0 and He2_1_Present_2[0][2] > 0 
    and He2_1_Past[0, 2] < 0):
    F_He2_1_Past_He2_1_Present_1 = (Q_He2_1_Past * Q_He2_1_Present_1 *
                                    k1/((R_He2_1_Past + R_He2_1_Present_1) ** 2))
    F_He2_1_Past_He2_1_Present_2 = (Q_He2_1_Past * Q_He2_1_Present_2 *
                                    k1/((R_He2_1_Past + R_He2_1_Present_2) ** 2))
    
    F_He2_1_Present_1_He2_1_Present_2 = (Q_He2_1_Present_1 * Q_He2_1_Present_2 * 
                                              k1/((R_He2_1_Present_1 + 
                                                   R_He2_1_Present_2) ** 2))


# He2_1_Future_1 
# Nuclear Communications - He2_1_Future_1

ρ_He2_1_Future_1 = np.array(He2_1_Future_1)[:, 0]
D3 = np.array(He2_1_Future_1)[:, 2]
divE_He2_1_Future_1 = ρ_He2_1_Future_1/unit.constantε0
for i, item in enumerate(D3):
    ρ_He2_1_Future_1[i] /= D3[i]
    divE_He2_1_Future_1[i] = ρ_He2_1_Future_1[i]/unit.constantε0

# He2_1_Future_2 

ρ_He2_1_Future_2 = np.array(He2_1_Future_2)[:, 0]
D4 = np.array(He2_1_Future_2)[:, 2]
divE_He2_1_Future_2 = ρ_He2_1_Future_2/unit.constantε0
for i, item in enumerate(D4):
    ρ_He2_1_Future_2[i] /= D4[i]
    divE_He2_1_Future_2[i] = ρ_He2_1_Future_2[i]/unit.constantε0

# He2_1_Future_3 

ρ_He2_1_Future_3 = np.array(He2_1_Future_3)[:, 0]
D5 = np.array(He2_1_Future_3)[:, 2]
divE_He2_1_Future_3 = ρ_He2_1_Future_3/unit.constantε0
for i, item in enumerate(D5):
    ρ_He2_1_Future_3[i] /= D5[i]
    divE_He2_1_Future_3[i] = ρ_He2_1_Future_3[i]/unit.constantε0
    
    
# He2_1_Future_1&(He2_1_Present_1 and He2_1_Present_2)
# He2_1_Future_1&(He2_1_Future_2 and He2_1_Future_3)
# Assessing the Impact on Nuclear Communications Using Coulomb's Law

R_He2_1_Future_1 = (3/4 * abs(He2_1_Future_1[0, 2])/unit.π) ** (1./3.)
R_He2_1_Future_2 = (3/4 * abs(sum(row[2] for row in He2_1_Future_2))/unit.π) ** (1./3.)
R_He2_1_Future_3 = (3/4 * abs(sum(row[2] for row in He2_1_Future_3))/unit.π) ** (1./3.)

Q_He2_1_Future_1 = He2_1_Future_1[0, 0]
Q_He2_1_Future_2 = sum(row[0] for row in He2_1_Future_2)
Q_He2_1_Future_3 = sum(row[0] for row in He2_1_Future_3)

if ((He2_1_Present_1[0][2] > 0 and He2_1_Present_2[0][2] > 0) 
    and He2_1_Future_1[0, 2] < 0 and (He2_1_Future_2[0][2] < 0 and 
                                      He2_1_Future_3[0][2] < 0)):
    F_He2_1_Future_1_He2_1_Present_1 = (Q_He2_1_Future_1 * Q_He2_1_Present_1 *
                                        k1/((R_He2_1_Future_1 + R_He2_1_Present_1) ** 2))
    
    F_He2_1_Future_1_He2_1_Present_2 = (Q_He2_1_Future_1 * Q_He2_1_Present_2 *
                                        k1/((R_He2_1_Future_1 + R_He2_1_Present_2) ** 2))    
      
    F_He2_1_Future_1_He2_1_Future_2 = (Q_He2_1_Future_1 * Q_He2_1_Future_2 *
                                       k1/((R_He2_1_Future_1 + R_He2_1_Future_2) ** 2))
    
    F_He2_1_Future_1_He2_1_Future_3 = (Q_He2_1_Future_1 * Q_He2_1_Future_3 *
                                       k1/((R_He2_1_Future_1 + R_He2_1_Future_3) ** 2))
    
    F_He2_1_Future_2_He2_1_Future_3 = (Q_He2_1_Future_2 * Q_He2_1_Future_3 * 
                                              k1/((R_He2_1_Future_2 + 
                                                   R_He2_1_Future_3) ** 2))
    
# find the mass of the first proton without taking into account 
# the masses with the second proton
M_He2_1_Present_1 = sum(row[1] for row in He2_1_Present_1)
M_He2_1_Future_2 = sum(row[1] for row in He2_1_Future_2)
M_He2_1 = M_He2_1_Present_1 + M_He2_1_Future_2

# find the mass of the second proton without taking into account 
# the masses shared with the first proton
M_He2_1_Present_2 = sum(row[1] for row in He2_1_Present_2)
M_He2_1_Future_3 = sum(row[1] for row in He2_1_Future_3)
M_He2_11 = M_He2_1_Present_2 + M_He2_1_Future_3


# find the resulting force vector
F_F1_Present_1_Present_2 = (F_He2_1_Future_1_He2_1_Present_1 ** 2 +
                            F_He2_1_Future_1_He2_1_Present_2 ** 2) ** (1./2.)
F_F1_Future_2_Future_3 = (F_He2_1_Future_1_He2_1_Future_2 ** 2 + 
                          F_He2_1_Future_1_He2_1_Future_3 ** 2) ** (1./2.)
F_F1 = (F_F1_Present_1_Present_2 ** 2 + F_F1_Future_2_Future_3 ** 2) ** (1./2.)

# find the disconnection time 2
T11 = abs(unit.constantc/(F_F1/M_He2_1))
T12 = abs(unit.constantc/(F_F1/M_He2_11))
T13 = min(T11, T12)

# find the resulting force vector
F_P_Present_1_Present_2 = (F_He2_1_Past_He2_1_Present_1 ** 2 + 
                           F_He2_1_Present_1_He2_1_Present_2 ** 2 + 
                           2 * F_He2_1_Past_He2_1_Present_1 * 
                           F_He2_1_Present_1_He2_1_Present_2 * 
                           cos(135 * unit.π/180)) ** (1./2.)

F_P_Present_1_Present22 = (F_He2_1_Past_He2_1_Present_2 ** 2 + 
                           F_He2_1_Present_1_He2_1_Present_2 ** 2 + 
                           2 * F_He2_1_Past_He2_1_Present_2 * 
                           F_He2_1_Present_1_He2_1_Present_2 * 
                           cos(135 * unit.π/180)) ** (1./2.)

F_P = (F_P_Present_1_Present_2 ** 2 + F_P_Present_1_Present22 ** 2 + 
       2 * F_P_Present_1_Present_2 * F_P_Present_1_Present22 * 
       cos(45 * unit.π/180)) ** (1./2.)

# find the disconnection time 1
T14 = unit.constantc/(F_P/M_He2_1)
T15 = unit.constantc/(F_P/M_He2_11)
T16 = min(T14, T15)

# Decay time of the first helium nucleus, which decays into two protons
T17 = max(T13, T16)

# Time spent by the first nucleus of a helium atom in the form of 
# an ion of a hydrogen molecule
T18 = T17 - T16    
        
# He2_1_Past_e 

ρ_He2_1_Past_e = np.array(He2_1_Past_e)[:, 0]
D6 = np.array(He2_1_Past_e)[:, 2]
divE_He2_1_Past_e = ρ_He2_1_Past_e/unit.constantε0
for i, item in enumerate(D6):
    ρ_He2_1_Past_e[i] /= D6[i]
    divE_He2_1_Past_e[i] = ρ_He2_1_Past_e[i]/unit.constantε0

# He2_1_Present_1_1e 

ρ_He2_1_Present_1_1e = np.array(He2_1_Present_1_1e)[:, 0]
D7 = np.array(He2_1_Present_1_1e)[:, 2]
divE_He2_1_Present_1_1e = ρ_He2_1_Present_1_1e/unit.constantε0
for i, item in enumerate(D7):
    ρ_He2_1_Present_1_1e[i] /= D7[i]
    divE_He2_1_Present_1_1e[i] = ρ_He2_1_Present_1_1e[i]/unit.constantε0

# He2_1_Present_2_1e 

ρ_He2_1_Present_2_1e = np.array(He2_1_Present_2_1e)[:, 0]
D8 = np.array(He2_1_Present_2_1e)[:, 2]
divE_He2_1_Present_2_1e = ρ_He2_1_Present_2_1e/unit.constantε0
for i, item in enumerate(D8):
    ρ_He2_1_Present_2_1e[i] /= D8[i]
    divE_He2_1_Present_2_1e[i] = ρ_He2_1_Present_2_1e[i]/unit.constantε0

# He2_1_Future_1_e 

ρ_He2_1_Future_1_e = np.array(He2_1_Future_1_e)[:, 0]
D9 = np.array(He2_1_Future_1_e)[:, 2]
divE_He2_1_Future_1_e = ρ_He2_1_Future_1_e/unit.constantε0
for i, item in enumerate(D9):
    ρ_He2_1_Future_1_e[i] /= D9[i]
    divE_He2_1_Future_1_e[i] = ρ_He2_1_Future_1_e[i]/unit.constantε0

field_He2_1_Future_2_e = []

ρ_He2_1_Future_2_e = np.array(He2_1_Future_2_e)[:, 0]
D10 = np.array(He2_1_Future_2_e)[:, 2]
divE_He2_1_Future_2_e = ρ_He2_1_Future_2_e/unit.constantε0
for i, item in enumerate(D10):
    ρ_He2_1_Future_2_e[i] /= D10[i]
    divE_He2_1_Future_2_e[i] = ρ_He2_1_Future_2_e[i]/unit.constantε0

# He2_1_Future_3_e 

ρ_He2_1_Future_3_e = np.array(He2_1_Future_3_e)[:, 0]
D11 = np.array(He2_1_Future_3_e)[:, 2]
divE_He2_1_Future_3_e = ρ_He2_1_Future_3_e/unit.constantε0
for i, item in enumerate(D11):
    ρ_He2_1_Future_3_e[i] /= D11[i]
    divE_He2_1_Future_3_e[i] = ρ_He2_1_Future_3_e[i]/unit.constantε0

        
# He2_2_Past_1 

ρ_He2_2_Past_1 = np.array(He2_2_Past_1)[:, 0]
D12 = np.array(He2_2_Past_1)[:, 2]
divE_He2_2_Past_1 = ρ_He2_2_Past_1/unit.constantε0
for i, item in enumerate(D12):
    ρ_He2_2_Past_1[i] /= D12[i]
    divE_He2_2_Past_1[i] = ρ_He2_2_Past_1[i]/unit.constantε0

# He2_2_Past_2 
# Nuclear Communications - He2_2_Past_2

ρ_He2_2_Past_2 = np.array(He2_2_Past_2)[:, 0]
D13 = np.array(He2_2_Past_2)[:, 2]
divE_He2_2_Past_2 = ρ_He2_2_Past_2/unit.constantε0
for i, item in enumerate(D13):
    ρ_He2_2_Past_2[i] /= D13[i]
    divE_He2_2_Past_2[i] = ρ_He2_2_Past_2[i]/unit.constantε0


# He2_2_Present_1 

ρ_He2_2_Present_1 = np.array(He2_2_Present_1)[:, 0]
D14 = np.array(He2_2_Present_1)[:, 2]
divE_He2_2_Present_1 = ρ_He2_2_Present_1/unit.constantε0
for i, item in enumerate(D14):
    ρ_He2_2_Present_1[i] /= D14[i]
    divE_He2_2_Present_1[i] = ρ_He2_2_Present_1[i]/unit.constantε0

# He2_2_Present_2 
# Nuclear Communications - He2_2_Present_2

ρ_He2_2_Present_2 = np.array(He2_2_Present_2)[:, 0]
D15 = np.array(He2_2_Present_2)[:, 2]
divE_He2_2_Present_2 = ρ_He2_2_Present_2/unit.constantε0
for i, item in enumerate(D15):
    ρ_He2_2_Present_2[i] /= D15[i]
    divE_He2_2_Present_2[i] = ρ_He2_2_Present_2[i]/unit.constantε0

# He2_2_Present_3 

ρ_He2_2_Present_3 = np.array(He2_2_Present_3)[:, 0]
D16 = np.array(He2_2_Present_3)[:, 2]
divE_He2_2_Present_3 = ρ_He2_2_Present_3/unit.constantε0
for i, item in enumerate(D16):
    ρ_He2_2_Present_3[i] /= D16[i]
    divE_He2_2_Present_3[i] = ρ_He2_2_Present_3[i]/unit.constantε0

# He2_2_Present_4 

ρ_He2_2_Present_4 = np.array(He2_2_Present_4)[:, 0]
D17 = np.array(He2_2_Present_4)[:, 2]
divE_He2_2_Present_4 = ρ_He2_2_Present_4/unit.constantε0
for i, item in enumerate(D17):
    ρ_He2_2_Present_4[i] /= D17[i]
    divE_He2_2_Present_4[i] = ρ_He2_2_Present_4[i]/unit.constantε0

# He2_2_Future_1 

ρ_He2_2_Future_1 = np.array(He2_2_Future_1)[:, 0]
D18 = np.array(He2_2_Future_1)[:, 2]
divE_He2_2_Future_1 = ρ_He2_2_Future_1/unit.constantε0
for i, item in enumerate(D18):
    ρ_He2_2_Future_1[i] /= D18[i]
    divE_He2_2_Future_1[i] = ρ_He2_2_Future_1[i]/unit.constantε0
        
# He2_2_Past_1_e 

ρ_He2_2_Past_1_e = np.array(He2_2_Past_1_e)[:, 0]
D19 = np.array(He2_2_Past_1_e)[:, 2]
divE_He2_2_Past_1_e = ρ_He2_2_Past_1_e/unit.constantε0
for i, item in enumerate(D19):
    ρ_He2_2_Past_1_e[i] /= D19[i]
    divE_He2_2_Past_1_e[i] = ρ_He2_2_Past_1_e[i]/unit.constantε0

# He2_2_Past_2_e 

ρ_He2_2_Past_2_e = np.array(He2_2_Past_2_e)[:, 0]
D20 = np.array(He2_2_Past_2_e)[:, 2]
divE_He2_2_Past_2_e = ρ_He2_2_Past_2_e/unit.constantε0
for i, item in enumerate(D20):
    ρ_He2_2_Past_2_e[i] /= D20[i]
    divE_He2_2_Past_2_e[i] = ρ_He2_2_Past_2_e[i]/unit.constantε0

# He2_2_Present_1_2e 

ρ_He2_2_Present_1_2e = np.array(He2_2_Present_1_2e)[:, 0]
D21 = np.array(He2_2_Present_1_2e)[:, 2]
divE_He2_2_Present_1_2e = ρ_He2_2_Present_1_2e/unit.constantε0
for i, item in enumerate(D21):
    ρ_He2_2_Present_1_2e[i] /= D21[i]
    divE_He2_2_Present_1_2e[i] = ρ_He2_2_Present_1_2e[i]/unit.constantε0

# He2_2_Present_2e 

ρ_He2_2_Present_2e = np.array(He2_2_Present_2e)[:, 0]
D22 = np.array(He2_2_Present_2e)[:, 2]
divE_He2_2_Present_2e = ρ_He2_2_Present_2e/unit.constantε0
for i, item in enumerate(D22):
    ρ_He2_2_Present_2e[i] /= D22[i]
    divE_He2_2_Present_2e[i] = ρ_He2_2_Present_2e[i]/unit.constantε0

# He2_2_Present_3e 

ρ_He2_2_Present_3e = np.array(He2_2_Present_3e)[:, 0]
D23 = np.array(He2_2_Present_3e)[:, 2]
divE_He2_2_Present_3e = ρ_He2_2_Present_3e/unit.constantε0
for i, item in enumerate(D23):
    ρ_He2_2_Present_3e[i] /= D23[i]
    divE_He2_2_Present_3e[i] = ρ_He2_2_Present_3e[i]/unit.constantε0

# He2_2_Present_4e 

ρ_He2_2_Present_4e = np.array(He2_2_Present_4e)[:, 0]
D24 = np.array(He2_2_Present_4e)[:, 2]
divE_He2_2_Present_4e = ρ_He2_2_Present_4e/unit.constantε0
for i, item in enumerate(D24):
    ρ_He2_2_Present_4e[i] /= D24[i]
    divE_He2_2_Present_4e[i] = ρ_He2_2_Present_4e[i]/unit.constantε0

# He2_2_Future_1e 

ρ_He2_2_Future_1e = np.array(He2_2_Future_1e)[:, 0]
D25 = np.array(He2_2_Future_1e)[:, 2]
divE_He2_2_Future_1e = ρ_He2_2_Future_1e/unit.constantε0
for i, item in enumerate(D25):
    ρ_He2_2_Future_1e[i] /= D25[i]
    divE_He2_2_Future_1e[i] = ρ_He2_2_Future_1e[i]/unit.constantε0


unit06 = He2_properties(divE_He2_1_Past, divE_He2_1_Present_1, divE_He2_1_Present_2, 
                  divE_He2_1_Future_1, divE_He2_1_Future_2, divE_He2_1_Future_3,                   
                  divE_He2_1_Past_e, divE_He2_1_Present_1_1e, divE_He2_1_Present_2_1e, 
                  divE_He2_1_Future_1_e, divE_He2_1_Future_2_e, divE_He2_1_Future_3_e,
                  
                  divE_He2_2_Past_1, divE_He2_2_Past_2, divE_He2_2_Present_1, 
                  divE_He2_2_Present_2, divE_He2_2_Present_3, divE_He2_2_Present_4, 
                  divE_He2_2_Future_1, divE_He2_2_Past_1_e, divE_He2_2_Past_2_e, 
                  divE_He2_2_Present_1_2e, divE_He2_2_Present_2e, divE_He2_2_Present_3e, 
                  divE_He2_2_Present_4e, divE_He2_2_Future_1e, T16, T17, T18)        

