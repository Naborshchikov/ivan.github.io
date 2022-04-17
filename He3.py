#!/usr/bin/env python
# coding: utf-8

# In[ ]:


class He3():
    
# He3 Atom

    def __init__ (self, He3_0_Past_0_e, He3_0_Past_1_e, He3_0_Past_2_e, He3_0_Past_3_e, 
                  He3_0_Present_0_e, He3_0_Present_1_e, He3_0_Present_2_e, 
                  He3_0_Present_3_e, He3_0_Future_0_e, He3_0_Future_1_e, 
                  
                  He3_1_Past_0_e, He3_1_Past_1_e, He3_1_Past_2_e, He3_1_Present_0_e,
                  He3_1_Present_1_e, He3_1_Present_2_e, He3_1_Present_3_e, 
                  He3_1_Future_0_e,  
                  
                  He3_2_Past_0_e, He3_2_Past_1_e, He3_2_Present_0_e, He3_2_Present_1_e, 
                  He3_2_Present_2_e, He3_2_Present_3_e, He3_2_Future_0_e, He3_2_Future_1_e, 
                  He3_2_Future_2_e, He3_2_Future_3_e, He3_2_Future_4_e, 
                  
                  He3_3_Past_0_e, He3_3_Past_1_e, He3_3_Past_2_e, He3_3_Past_3_e, 
                  He3_3_Present_0_e, He3_3_Present_1_e, He3_3_Present_2_e, 
                  He3_3_Present_3_e, He3_3_Future_0_e, He3_3_Future_1_e, He3_3_Future_2_e, 
                  
                  He3_4_Past_0_e, He3_4_Past_1_e, He3_4_Past_2_e, He3_4_Present_0_e, 
                  He3_4_Present_1_e, He3_4_Present_2_e, He3_4_Present_3_e, 
                  He3_4_Future_0_e, He3_4_Future_1_e):
        
        self.He3_0_Past_0_e = He3_0_Past_0_e
        self.He3_0_Past_1_e = He3_0_Past_1_e
        self.He3_0_Past_2_e = He3_0_Past_2_e
        self.He3_0_Past_3_e = He3_0_Past_3_e
        
        self.He3_0_Present_0_e = He3_0_Present_0_e
        self.He3_0_Present_1_e = He3_0_Present_1_e
        self.He3_0_Present_2_e = He3_0_Present_2_e
        self.He3_0_Present_3_e = He3_0_Present_3_e
        
        self.He3_0_Future_0_e = He3_0_Future_0_e
        self.He3_0_Future_1_e = He3_0_Future_1_e
                
        self.He3_1_Past_0_e = He3_1_Past_0_e
        self.He3_1_Past_1_e = He3_1_Past_1_e
        self.He3_1_Past_2_e = He3_1_Past_2_e
                
        self.He3_1_Present_0_e = He3_1_Present_0_e
        self.He3_1_Present_1_e = He3_1_Present_1_e
        self.He3_1_Present_2_e = He3_1_Present_2_e
        self.He3_1_Present_3_e = He3_1_Present_3_e
        
        self.He3_1_Future_0_e = He3_1_Future_0_e
        
        self.He3_2_Past_0_e = He3_2_Past_0_e
        self.He3_2_Past_1_e = He3_2_Past_1_e
                
        self.He3_2_Present_0_e = He3_2_Present_0_e
        self.He3_2_Present_1_e = He3_2_Present_1_e
        self.He3_2_Present_2_e = He3_2_Present_2_e
        self.He3_2_Present_3_e = He3_2_Present_3_e
        
        self.He3_2_Future_0_e = He3_2_Future_0_e
        self.He3_2_Future_1_e = He3_2_Future_1_e
        self.He3_2_Future_2_e = He3_2_Future_2_e
        self.He3_2_Future_3_e = He3_2_Future_3_e
        
        self.He3_2_Future_4_e = He3_2_Future_4_e
                        
        self.He3_3_Past_0_e = He3_3_Past_0_e
        self.He3_3_Past_1_e = He3_3_Past_1_e
        self.He3_3_Past_2_e = He3_3_Past_2_e
        self.He3_3_Past_3_e = He3_3_Past_3_e
        
        self.He3_3_Present_0_e = He3_3_Present_0_e
        self.He3_3_Present_1_e = He3_3_Present_1_e
        self.He3_3_Present_2_e = He3_3_Present_2_e
        self.He3_3_Present_3_e = He3_3_Present_3_e
        
        self.He3_3_Future_0_e = He3_3_Future_0_e
        self.He3_3_Future_1_e = He3_3_Future_1_e
        self.He3_3_Future_2_e = He3_3_Future_2_e
        
        self.He3_4_Past_0_e = He3_4_Past_0_e
        self.He3_4_Past_1_e = He3_4_Past_1_e
        self.He3_4_Past_2_e = He3_4_Past_2_e
                
        self.He3_4_Present_0_e = He3_4_Present_0_e
        self.He3_4_Present_1_e = He3_4_Present_1_e
        self.He3_4_Present_2_e = He3_4_Present_2_e
        self.He3_4_Present_3_e = He3_4_Present_3_e
        
        self.He3_4_Future_0_e = He3_4_Future_0_e
        self.He3_4_Future_1_e = He3_4_Future_1_e              

# He3_0
# Core
He3_0_Past_0 = unit01.PROTON2_Future_matrix
He3_0_Past_1 = unit03.He3_intersection
He3_0_Past_2 = unit01.PROTON2_Present_matrix[:4]
He3_0_Past_3 = unit01.NEUTRON_Past_matrix[:1]

He3_0_Present_0 = unit03.He3_intersection1
He3_0_Present_1 = unit01.PROTON_Present_matrix[:3]
He3_0_Present_2 = unit01.PROTON2_Past_matrix[:2]
He3_0_Present_3 = unit01.NEUTRON_Present_matrix[:4]

He3_0_Future_0 = unit01.PROTON_Future_matrix
He3_0_Future_1 = unit01.NEUTRON_Future_matrix

# Atom
He3_0_Past_0_e = unit01.PROTON2_Future_matrix
He3_0_Past_1_e = unit03.He3_intersection
He3_0_Past_2_e = unit01.PROTON2_Present_matrix[:4]
He3_0_Past_3_e = unit01.NEUTRON_Past_matrix[:1]

He3_0_Present_0_e = unit03.He3_intersection1

He3_0_Present_1_e = He3_0_Present_1[2] + unit02.e
He3_0_Present_1_e = np.vstack((np.array(He3_0_Present_1[:2]), np.array(He3_0_Present_1_e)))

He3_0_Present_2_e = unit01.PROTON2_Past_matrix[:2]

He3_0_Present_3_e = He3_0_Present_3[3] + unit02.e
He3_0_Present_3_e = np.vstack((np.array(He3_0_Present_3[:3]), np.array(He3_0_Present_3_e)))

He3_0_Future_0_e = unit01.PROTON_Future_matrix
He3_0_Future_1_e = unit01.NEUTRON_Future_matrix

# He3_1
# Core
He3_1_Past_0 = unit01.PROTON_Past_matrix
He3_1_Past_1 = unit03.He3_intersection2
He3_1_Past_2 = unit01.NEUTRON_Past_matrix[:1]

He3_1_Present_0 = unit03.He3_intersection3
He3_1_Present_1 = unit01.NEUTRON_Future_matrix[:1]
He3_1_Present_2 = unit01.PROTON_Present_matrix[:3]
He3_1_Present_3 = unit01.PROTON2_Past_matrix[:2]

He3_1_Future_0 = unit01.PROTON_Future_matrix

# Atom
He3_1_Past_0_e = unit01.PROTON_Past_matrix
He3_1_Past_1_e = unit03.He3_intersection2
He3_1_Past_2_e = unit01.NEUTRON_Past_matrix[:1]

He3_1_Present_0_e = unit03.He3_intersection3

He3_1_Present_1_e = He3_1_Present_1 + unit02.e

He3_1_Present_2_e = He3_1_Present_2[2] + unit02.e
He3_1_Present_2_e = np.vstack((np.array(He3_1_Present_2[:2]), np.array(He3_1_Present_2_e))) 

He3_1_Present_3_e = unit01.PROTON2_Past_matrix[:2]

He3_1_Future_0_e = unit01.PROTON_Future_matrix

# He3_2
# Core
He3_2_Past_0 = unit01.PROTON2_Future_matrix
He3_2_Past_1 = unit01.PROTON_Past_matrix

He3_2_Present_0 = unit03.He3_intersection4
He3_2_Present_1 = unit01.PROTON_Present_matrix[:3]
He3_2_Present_2 = unit01.PROTON2_Present_matrix[:4]
He3_2_Present_3 = unit01.NEUTRON_Past_matrix[:1]

He3_2_Future_0 = unit03.He3_intersection5
He3_2_Future_1 = unit01.PROTON_Future_matrix[:3]
He3_2_Future_2 = unit01.PROTON2_Past_matrix[:2]
He3_2_Future_3 = unit01.NEUTRON_Present_matrix[:4]

He3_2_Future_4 = unit01.NEUTRON_Future_matrix

# Atom
He3_2_Past_0_e = unit01.PROTON2_Future_matrix
He3_2_Past_1_e = unit01.PROTON_Past_matrix

He3_2_Present_0_e = unit03.He3_intersection4

He3_2_Present_1_e = He3_2_Present_1[2] + unit02.e
He3_2_Present_1_e = np.vstack((np.array(He3_2_Present_1[:2]), np.array(He3_2_Present_1_e)))

He3_2_Present_2_e = He3_2_Present_2[3] + unit02.e
He3_2_Present_2_e = np.vstack((np.array(He3_2_Present_2[:3]), np.array(He3_2_Present_2_e)))

He3_2_Present_3_e = unit01.NEUTRON_Past_matrix[:1]

He3_2_Future_0_e = unit03.He3_intersection5
He3_2_Future_1_e = unit01.PROTON_Future_matrix[:3]
He3_2_Future_2_e = unit01.PROTON2_Past_matrix[:2]
He3_2_Future_3_e = unit01.NEUTRON_Present_matrix[:4]

He3_2_Future_4_e = unit01.NEUTRON_Future_matrix

# He3_3
# Core
He3_3_Past_0 = unit01.PROTON2_Future_matrix

He3_3_Past_1 = unit01.PROTON2_Present_matrix
He3_3_Past_2 = unit01.PROTON_Past_matrix
He3_3_Past_3 = unit01.NEUTRON2_Future_matrix

He3_3_Present_0 = unit03.He3_intersection6
He3_3_Present_1 = unit01.PROTON_Present_matrix[:3]
He3_3_Present_2 = unit01.PROTON2_Past_matrix[:2]
He3_3_Present_3 = unit01.NEUTRON2_Present_matrix[:2]

He3_3_Future_0 = unit03.He3_intersection7
He3_3_Future_1 = unit01.PROTON_Future_matrix[:3]
He3_3_Future_2 = unit01.NEUTRON2_Past_matrix[:2]

# Atom
He3_3_Past_0_e = unit01.PROTON2_Future_matrix

He3_3_Past_1_e = unit01.PROTON2_Present_matrix
He3_3_Past_2_e = unit01.PROTON_Past_matrix
He3_3_Past_3_e = unit01.NEUTRON2_Future_matrix

He3_3_Present_0_e = unit03.He3_intersection6

He3_3_Present_1_e = He3_3_Present_1[2] + unit02.e
He3_3_Present_1_e = np.vstack((np.array(He3_3_Present_1[:2]), np.array(He3_3_Present_1_e)))

He3_3_Present_2_e = unit01.PROTON2_Past_matrix[:2]

He3_3_Present_3_e = He3_3_Present_3[1] + unit02.e
He3_3_Present_3_e = np.vstack((np.array(He3_3_Present_3[:1]), np.array(He3_3_Present_3_e)))

He3_3_Future_0_e = unit03.He3_intersection7
He3_3_Future_1_e = unit01.PROTON_Future_matrix[:3]
He3_3_Future_2_e = unit01.NEUTRON2_Past_matrix[:2]

# He3_4
# Core
He3_4_Past_0 = unit01.NEUTRON2_Future_matrix

He3_4_Past_1 = unit03.He3_intersection8
He3_4_Past_2 = unit01.NEUTRON2_Present_matrix[:2]

He3_4_Present_0 = unit03.He3_intersection9
He3_4_Present_1 = unit01.PROTON_Present_matrix[:3]
He3_4_Present_2 = unit01.PROTON2_Present_matrix[:4]
He3_4_Present_3 = unit01.NEUTRON2_Past_matrix[:2]

He3_4_Future_0 = unit01.PROTON_Future_matrix
He3_4_Future_1 = unit01.PROTON2_Past_matrix

# Atom
He3_4_Past_0_e = unit01.NEUTRON2_Future_matrix

He3_4_Past_1_e = unit03.He3_intersection8
He3_4_Past_2_e = unit01.NEUTRON2_Present_matrix[:2]

He3_4_Present_0_e = unit03.He3_intersection9

He3_4_Present_1_e = He3_4_Present_1[2] + unit02.e
He3_4_Present_1_e = np.vstack((np.array(He3_4_Present_1[:2]), np.array(He3_4_Present_1_e)))

He3_4_Present_2_e = He3_4_Present_2[3] + unit02.e
He3_4_Present_2_e = np.vstack((np.array(He3_4_Present_2[:3]), np.array(He3_4_Present_2_e)))

He3_4_Present_3_e = unit01.NEUTRON2_Past_matrix[:2]

He3_4_Future_0_e = unit01.PROTON_Future_matrix
He3_4_Future_1_e = unit01.PROTON2_Past_matrix
                  
unit05 = He3(He3_0_Past_0_e, He3_0_Past_1_e, He3_0_Past_2_e, He3_0_Past_3_e, 
                  He3_0_Present_0_e, He3_0_Present_1_e, He3_0_Present_2_e, 
                  He3_0_Present_3_e, He3_0_Future_0_e, He3_0_Future_1_e, 
                  
                  He3_1_Past_0_e, He3_1_Past_1_e, He3_1_Past_2_e, He3_1_Present_0_e,
                  He3_1_Present_1_e, He3_1_Present_2_e, He3_1_Present_3_e, 
                  He3_1_Future_0_e,  
                  
                  He3_2_Past_0_e, He3_2_Past_1_e, He3_2_Present_0_e, He3_2_Present_1_e, 
                  He3_2_Present_2_e, He3_2_Present_3_e, He3_2_Future_0_e, He3_2_Future_1_e, 
                  He3_2_Future_2_e, He3_2_Future_3_e, He3_2_Future_4_e, 
                  
                  He3_3_Past_0_e, He3_3_Past_1_e, He3_3_Past_2_e, He3_3_Past_3_e, 
                  He3_3_Present_0_e, He3_3_Present_1_e, He3_3_Present_2_e, 
                  He3_3_Present_3_e, He3_3_Future_0_e, He3_3_Future_1_e, He3_3_Future_2_e, 
                  
                  He3_4_Past_0_e, He3_4_Past_1_e, He3_4_Past_2_e, He3_4_Present_0_e, 
                  He3_4_Present_1_e, He3_4_Present_2_e, He3_4_Present_3_e, 
                  He3_4_Future_0_e, He3_4_Future_1_e)

