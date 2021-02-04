# ivan.github.io
<html>
<head>
<title>
Page README <br>
Project name: <br> 
<h1> Keywords: python, environmentally friendly energy, new materials and equipment, exotic materials, quarks, proton, neutron, nuclear, chemical, pharmaceutical, thermoelectric generator. </h1>
</title>
</head>
<body>
<h2> <p>Physics and programming as a source of optimization of existing technologies, obtaining environmentally friendly energy, new materials and equipment. </p>
Description: </h2>
<p> The Project has implemented an algorithm for a part of the physical-mathematical model. The physical-mathematical model allows: </p>

<p> To calculate the optimal characteristics for work in the nuclear, chemical, pharmaceutical industries. </p>

<p> Receiving exotic materials such as metallic hydrogen as a consequence of strong environmentally clean power sources. </p>

<p> To begin the development of new medical equipment based on non-contact point destruction of harmful bacteria, viruses, and mutated cells in the human body. </p>

<p> The physical-mathematical model allows you to calculate the optimal characteristics for obtaining new ultra-light, super-strong materials. </p>

<p> <h3> Project implemented in python ver. 3.7.6, with using anaconda 2020 02. </h3> <br>
To implement the project used: <br>
1. The libraries: pandas; numpy; matplotlib; scipy; prettytable; collections; plotly. <br>
2. The data: Nobel laureates; constants recommended by NIST USA (Table#1). <br>
3. The author's algorithm. </p>

<p><h3> The program includes 23 classes: </h3> </p>
<p> <h4> 1. class Preliminary() </h4>  - Calculation of additional data obtained based on data from Table#1. The table contains links to the data obtained by the Nobel laureates, the data of constants recommended by the US NIST. </p>

<p> <h4> 2. class Newnewproton()</h4>   - Calculation of three shells for new quarks "u" and "d", nine shells for a new proton, a new neutron in terms of mass, charge, volume. The calculation is implemented according to the author's algorithm. The applied algorithm allows to calculate any known particles and find new ones. This algorithm is only one of three parts of the physical-mathematical model. Code for class Newnewproton(): </p>
<p> class Newnewproton():
# The magnitude of the charge of the core, shells in the newproton, newneutron, respectively<br>
#Robert Hofstadter the Nobel laureate<br>
    SHELLSP1 = 0.35<br>
    SHELLSP2 = 0.5<br>
    SHELLSP3 = 0.15<br>
    SHELLSN1 = 0.35<br>
    SHELLSN2 = -0.5<br>
    SHELLSN3 = 0.15<br>    
# The mass of the core, shells in the new proton, new neutron, respectively<br>
    shellsmp1 = 1.67262192369E-27 * 0.91<br>   
    shellsmp2 = Preliminary.mps1<br>
    shellsmp3 = Preliminary.mps2<br>    
    shellsmn1 = 1.67492749804E-27 * 0.91<br>   
    shellsmn2 = Preliminary.mns1<br>
    shellsmn3 = Preliminary.mns2<br>     
    def __init__ (self, array): <br>
        self.array = array<br>
# Array input according to the matrix proposed by the author<br>
a1 = array ([[2.0 , 1.0, 1.0, 1.0, 1.0, 0.0], <br> 
             [0.0, 1.0, 0.0, 0.0, 0.0, 0.0], <br> 
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0], <br>
             [1.0, 1.0, 0.0, 2.0, 1.0, 1.0], <br>
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0], <br>
             [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]) <br>
unit = Newnewproton(a1) <br>
unit.array</p>

<p> <h4> 3. class Pseudonewneutron() </h4>   - Calculation of three shells for quarks "u" and "d", nine shells for a proton, a neutron in terms of mass, charge, volume. The calculation is implemented according to the author's algorithm. According to the calculation, the proton and neutron have three wormholes each. Table 10 contains data on three negative values of the proton shells. Table 11 contains data on three negative values of the neutron shells. A negative volume is that part of the volume that is in another dimension, that is, a part of a wormhole located in another dimension. </p>
<p> Code for class Pseudonewneutron():  </p>
<p> class Pseudonewneutron():<br>
# The difference from the class new proton in the matrix<br>    
    def __init__ (self, arr): <br>
        self.arr = arr<br>        
# Array input according to the matrix proposed by the author<br>
a2 = array ([[2.0 , 2.0, 1.0, 1.0, 0.0, 0.0], <br>
             [0.0, 0.0, 1.0, 0.0, 1.0, 0.0], <br>
             [0.0, 0.0, 0.0, 0.0, 0.0, 1.0], <br> 
             [1.0, 0.0, 0.0, 2.0, 2.0, 1.0], <br>
             [0.0, 1.0, 0.0, 0.0, 0.0, 1.0], <br>
             [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]]) <br>
uni = Pseudonewneutron(a2) <br>
uni.arr</p>
<p> Graphs 11, 17, 23 demonstrate the presence of shells with negative volume, that is, with volume placed in another dimension. </p>
<p>  <img src=" https://github.com/Naborshchikov/ivan.github.io/blob/main/Graph%2011.png " align="left">
 
</p>
<p> <img src=" https://github.com/Naborshchikov/ivan.github.io/blob/main/Graph%2017.png " align="left">

</p>
<p>
<p> <img src=" https://github.com/Naborshchikov/ivan.github.io/blob/main/Graph%2023.png" align="left">

 </p>
<p> According to the article of famous physicists A. Yu. Khlestkov & Yu. A. Khlestkov from Russian National Research Nuclear University "MEPhI" (https://link.springer.com/article/10.1007/s11182-019-01709-9), the neutron has a wormhole, which can be considered as proof of the author's algorithm. </p>
<p> <h4> 4. class Wavep(), class Waven(), class Wavepsn(), class Wavepsp() </h4>   - Classes for calculating the Compton wavelength from shells for a new proton, new neutron, proton, neutron. Code for class Wavep(): </p>
<p> class Wavep():<br>
# Planck's constant<br>
    CONSTANTH = 6.62607015E-34<br>    
# The speed of light in a vacuum<br>
    CONSTANTC = 299792458<br> 
        
# The ratio of Planck's constant to the speed of light in a vacuum, D<br>
    D = CONSTANTH/CONSTANTC  <br>    
    def __init__ (self, newcomptonlp): <br>
        self.newcomptonlp = newcomptonlp<br>
numbers = [1/newnewprotons[0].mass, 1/newnewprotons[1].mass, <br>
           1/newnewprotons[2].mass, 1/newnewprotons[3].mass, <br> 
           1/newnewprotons[4].mass, 1/newnewprotons[5].mass, <br> 
           1/newnewprotons[6].mass, 1/newnewprotons[7].mass, <br> 
           1/newnewprotons[8].mass] <br>
for i, item in enumerate(numbers): <br>
       numbers[i] *= Wavep.D<br>        
unit2 = Wavep(numbers) </p>
<p> <h4> 5. class ElectricWavep(), class ElectricWaven(), class ElectricWavepsn(), class ElectricWavepsp() </h4>  - Classes for calculating the electrical ε for each shell at the speed of light.  Code for class ElectricWavep(): </p>
<p> class ElectricWavep():<br>
# Planck's constant<br>
    CONSTANTH = 6.62607015E-34<br>    
# The speed of light in a vacuum<br>
    CONSTANTC = 299792458<br>
# The electrical constant ε  <br> 
    CONSTANTE0 = 8.85418781762039E-12<br>        
# The ratio to the unit of the doubled product of the electrical constant, <br> 
# Planck's constant, the speed of light in vacuum<br>
    constantd = 1/(2 * CONSTANTE0 * CONSTANTH * CONSTANTC) <br>    
    def __init__ (self, newelektromagnetikp): <br>
        self.newelektromagnetikp = newelektromagnetikp<br>
numbers = [newnewprotons[0].charge **2, newnewprotons[1].charge **2, <br>
           newnewprotons[2].charge **2, <br>
           newnewprotons[3].charge **2, newnewprotons[4].charge **2, <br>
           newnewprotons[5].charge **2, <br> 
           newnewprotons[6].charge **2, newnewprotons[7].charge **2, <br> 
           newnewprotons[8].charge **2] <br>

for i, item in enumerate(numbers): <br>
       numbers[i] *= ElectricWavep.constantd<br>        
unit6 = ElectricWavep(numbers) </p>
<p> <h4> 6. class GravityWavep(), class GravityWaven(), class GravityWavepsn(), class GravityWavepsp() </h4>   - Calculation of gravity characteristics for shells. Code for class GravityWavep():</p>
<p> class GravityWavep():<br>
# Planck's constant<br>
    CONSTANTH = 6.62607015E-34<br>    
# The speed of light in a vacuum<br>
    CONSTANTC = 299792458<br>
# The Gravitational constant<br>  
    CONSTANTEG = 6.67448478E-11<br>    
    π = 3.14159265358979<br>        
# The ratio of the doubled product of pi and the gravitational constant to<br> 
# Planck's constant, the speed of light in vacuum<br>
    constantg = 2 * π * CONSTANTEG/(CONSTANTH * CONSTANTC) <br>    
    def __init__ (self, newgravp): <br>
        self.newgravp = newgravp<br>
numbers = [newnewprotons[0].mass **2, newnewprotons[1].mass **2, <br>
           newnewprotons[2].mass **2, newnewprotons[3].mass **2, <br> 
           newnewprotons[4].mass **2, newnewprotons[5].mass **2, <br>
           newnewprotons[6].mass **2, newnewprotons[7].mass **2, <br> 
           newnewprotons[8].mass **2] <br>
for i, item in enumerate(numbers): <br>
       numbers[i] *= GravityWavep.constantg<br>        
unit10 = GravityWavep(numbers) </p>
<p> <h4> 7. class Frequencyp(), class Frequencyn(), class Frequencypsn(), class Frequencypsp() </h4>   - Calculation of frequency characteristics for shells taking into account the speed of light in vacuum. Code for class Frequencyp(): </p>
<p> class Frequencyp():<br>
# The speed of light in a vacuum<br>
    CONSTANTC = 299792458<br>  


    def __init__ (self, newfrequencep): <br>
        self.newfrequencep = newfrequencep<br>
numbers = [1/unit2.newcomptonlp[0], <br>  
           1/unit2.newcomptonlp[1], 1/unit2.newcomptonlp[2], <br>
           1/unit2.newcomptonlp[3], <br>
           1/unit2.newcomptonlp[4], 1/unit2.newcomptonlp[5], <br>  
           1/unit2.newcomptonlp[6], <br>
           1/unit2.newcomptonlp[7], 1/unit2.newcomptonlp[8]] <br>
for i, item in enumerate(numbers): <br>
       numbers[i] *= Frequencyp.CONSTANTC<br>        
unit14 = Frequencyp(numbers) </p>
<p> <h4> 8. class Melectron() </h4>   - Calculation of the existence of a point charged particle with  a charge 60 times less than the charge of an electron. Particle electric charge microminus (author's name) - -2.6702944e-21C. Microplus particle electric charge (author's name) -  2.6702944e-21C. Microplus and microminus mass (Kg): 1.518230622602297e-32. Volume microplus and microminus (cbm): 8.726646292652101e-66. Code for class Melectron(): </p>
<p> class Melectron():<br>          
     def __init__ (self, melectron_charge_amount, melectron_mass_amount, <br>
                  melectron_volume_amount): <br>  
        self.melectron_charge_amount = melectron_charge_amount<br>
        self.melectron_mass_amount = melectron_mass_amount<br>
        self.melectron_volume_amount = melectron_volume_amount<br>        
# Let's compare the minimum values of charges in a newproton, newneutron, <br> 
# pseudo new neutron, pseudo new proton and<br>
# find the value of a point charged particle<br>
# Let's name the new particle - melectron<br>
if (newnewprotons_min_charge == newnewneutrons_min_charge and newnewprotons_min_charge < <br> 
    protopsns_min_charge == neutropsns_min_charge): <br>
    for i in range(9): <br>
        melectron_charge_amount = (newnewprotons[i].charge/newnewprotons_min_charge - <br>
                                   newnewprotons[i].charge//newnewprotons_min_charge) <br>
        melectron_charge_amount = melectron_charge_amount * newnewprotons_min_charge<br>        
# find the mass of a point charged particle - melectron<br>
# Qe - electron charge modulo<br>
melectron_mass_amount = Preliminary.me/(Qe/melectron_charge_amount) <br>
# minimum_volume_amount particles - melectron<br>
melectron_volume_amount = Preliminary.Ve/(Qe/melectron_charge_amount) <br>     
unit18 = Melectron(melectron_charge_amount, melectron_mass_amount, melectron_volume_amount) </p>
<p> <h4> 9. class NewnewprotonCycles(), class NewneutronCycles(), class PsnewneutronCycles(), class PsnewnewprotonsCycles() </h4>   - Calculation of the characteristics of a new proton, new neutron, neutron, proton, by phases within the cycle. The new proton has five phases, the rest have four phases in a cycle. It is suggested that you familiarize yourself with these classes in the program itself. <br> </p>
<p>In the Project, the algorithm was tested for a proton, a neutron in terms of electric charge after the calculations. The result of testing the algorithm showed no error. </p>

<p>The result of the assessment using the algorithm of a new proton, a new neutron showed the presence of a total for shells electric charge in 2.403264951e-20 C for a new neutron and revealed differences between a proton and a new proton. The presence of a small charge in the new neutron and the absence of a total of shells electric charge in the neutron speaks of the validity of the Dirac equation without exception. </p>
<p> Graphs 26, 32 show the distribution of the electric charge over the shells for the neutron and new neutron over the cycle phases, respectively. </p>
<p>
<img src=" https://github.com/Naborshchikov/ivan.github.io/blob/main/Graph%2026.png" align="left">
 </p>
<p>
<img src=" https://github.com/Naborshchikov/ivan.github.io/blob/main/Graph%2032.png" align="left">
 </p>

<p> To visualize the data obtained, 45 different graphs are presented. Key data including conclusions are presented in 12 tables.  Graphs 34 through 45 show the relationship of fundamental forces, mass, charge, volume. The result shows that the five fundamental forces, mass, charge, volume, represent a pyramidal dependence of the 8th order. To visualize imagine a square pyramid, each edge of which represents one of the fundamental forces, mass, charge, volume. </p>
<p>The data obtained by the Project allows the user to use it to calculate the General field.  Using the data obtained by the Project and the capabilities of the function “from scipy import integrate” the user can calculate the characteristics required to optimize the technological process. Data arrays for each indicator, taking into account the identified linear patterns, can be both integrated and differentiated depending on the task at hand. This will allow you to calculate the characteristics required for a particular case, for example, field strength, current strength, and so on. </p>
<p>These results simplify the control of reactions, enable users to receive cumulative benefits, reduce energy and time costs; receive new materials. </p>

<h2> Installation: </h2>

<p> The project does not require any special skills to install. It will work on Jupiter Notebook. If you are missing some packages, then they are easily loaded using pip. </p>

<h2> Usage:  </h2>

<p> To start the Project, it is enough to use Run. Everything will be done automatically. The program contains all the necessary comments.  </p>
<p> Therefore, if you need to use the received data in another program, you can easily find a place to insert the script.  </p>
<p> The project allows developing new industrial technologies, including obtaining exotic materials, assessing the impact of fields on the human condition, and obtaining new environmentally friendly energy sources. </p>
<p> It is possible to obtain exotic materials with the following characteristics: </p>

<p> Radiation resistant from neutrons pcs/m3, respectively ≥ 1.5978313081543083*e48. Note, that the decay of one ton of uranium 92U235 will produce 6.12*e27 (approximately) neutrons. Thus, one cubic meter of material will absorb neutrons from the decay 6.6*e20 tons of uranium 92U235. </p>

<p> Temperature resistant from photons, 0C, respectively  ≤ 81*e18. Note, that the temperature at the center of a nuclear explosion at 0C is e7 (approximately), that is, the material of the EMH absorbs heat from e11 nuclear explosions, and the protective field will withstand temperatures exceeding the temperature at the center of a nuclear explosion 81*e11 times.  </p>
<p> Knowledge of the structure of the proton, the cyclicity of its phases, the relationship of fundamental forces, mass, charge, volume allows the use of their mutual reinforcement in the production of metallic hydrogen. Imagine a "positive resonance" that, with little external stimulus, produces a significant positive effect. </p>
<p> In terms of use for the development of new equipment for human treatment. Each person hears a lot of what he considers to be noise. A similar thing happens in the human body, each particle has its own characteristics and causes noise characteristic only for it. But just as each person recognizes another by their voice, so the equipment recognizes a particle, bacterium, molecule by its characteristics. To do this, the device needs to know its characteristics. The computer into which the proposed physical and mathematical model is loaded will calculate such characteristics. </p>
<p> Today, the location of any smartphone user can be found. The device will be able to find a particle, bacterium, virus in the human body in the same way. The knowledge obtained with the help of a physical-mathematical model, which, on the one hand, allows creating atoms, particles with minimal energy costs, will, on the other hand, also destroy viruses and bacteria harmful to the body. </p>
<p> Thus, any person will be able to recover from a large number of diseases in a hospital in a short time and return to the family full of strength and health. </p>
<p> The project allows the development of a new generation thermoelectric generator. It will work for tens of years, providing electricity to a house, car, plane, and spacecraft. It will be eco-friendly, safe. </p>
<p> The project allows the production of an ultra-light, ultra-strong material from hydrogen. Density: 287.5 kg/m. cube. The number of atoms/ions in the lattice 3/1. Heat resistance 2.3xe25 J/kg.  </p>
<p> Only part of the Project is posted on the site in the form of program code due to the specifics of its capabilities. </p>
<p> For example, the points of small extrema on the graph demonstrate the possibility of creating a measuring device for studying the structure of matter, and the point of maximum extremum will allow the user to create a neutron cutter. Used 512 dots. </p>
<p>
<img src="https://github.com/Naborshchikov/ivan.github.io/blob/main/Graph%20n.png" align="left">

</p>

<p> The graph below allows the user to better understand the effect of gravity, frequencies at different speeds on neutron shells at different speeds. Used 64000 dots. </p>
<p>
<img src="https://github.com/Naborshchikov/ivan.github.io/blob/main/Graph%20n3.png" align="left">
 
</p>

<h2> Credits: </h2>

<p> The author of the project is <h3> Ivan Naborshchikov</h3>.  </p>
<p> If you have any questions, you can ask them via: </p>
<p> Information about me is located at this address: <br>
<h4>  http://www.linkedin.com/in/ivan-naborshchikov/ </h4> </p>
 <p> The project is located at this address: <br>
 <h4>   https://github.com/Naborshchikov/ivan.github.io.git/  </h4> </p>

</body>
</html>
