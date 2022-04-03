#!/usr/bin/env python
# coding: utf-8

# In[ ]:


class Deltam():
    def __init__ (self, NEUTRON_Delta, NEUTRON_Delta2):
        
            self.NEUTRON_Delta = NEUTRON_Delta
            self.NEUTRON_Delta2 = NEUTRON_Delta2        
        
PROTON = np.sum(unit01.PROTON_Past_M + unit01.PROTON_Present_M + unit01.PROTON_Future_M)
PROTON2 = np.sum(unit01.PROTON2_Past_M + unit01.PROTON2_Present_M + unit01.PROTON2_Future_M)

NEUTRON = np.sum(unit01.NEUTRON_Past_M + unit01.NEUTRON_Present_M + unit01.NEUTRON_Future_M)
NEUTRON2 = np.sum(unit01.NEUTRON2_Past_M + unit01.NEUTRON2_Present_M + unit01.NEUTRON2_Future_M)

if (PROTON - PROTON2 == 0):
    NEUTRON_Delta = NEUTRON - NEUTRON2 
if (NEUTRON_Delta != 0):
    NEUTRON_Delta2 = - NEUTRON_Delta
    
if (NEUTRON_Delta > 0):
    print('NEUTRON_Delta is used for exothermic reactions\n')
if (NEUTRON_Delta2 < 0):
    print('NEUTRON_Delta2 is used in endothermic reactions')
        
unit04 = Deltam(NEUTRON_Delta, NEUTRON_Delta2)

