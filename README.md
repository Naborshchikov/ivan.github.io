# ivan.github.io
<!DOCTYPE html>
<html lang="en">
<head>
<title> Physics and programming as a source of optimization of existing technologies, obtaining environmentally friendly energy, new materials and equipment.</title>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
body {
  font-family: Arial, Helvetica, sans-serif;
}
</style>
</head>
<body>

<div class="header">
<h1> Physics and programming as a source of optimization of existing technologies, obtaining environmentally friendly energy, new materials and equipment.</h1>
</div>
.header {
  padding: 80px; /* some padding */
  text-align: center; /* center the text */
  background: #1abc9c; /* green background */
  color: white; /* white text color */
}

/* Increase the font size of the <h1> element */
.header h1 {
  font-size: 40px;
}
<div class="navbar">
  <a href=" The project: ">https://github.com/Naborshchikov/ivan.github.io.git/</a>
  <a href=" Information about me: ">http://www.linkedin.com/in/ivan-naborshchikov/</a>
  <a href="#" class="right">Link</a>
</div>
/* Style the top navigation bar */
.navbar {
  overflow: hidden; /* Hide overflow */
  background-color: #333; /* Dark background color */
}

/* Style the navigation bar links */
.navbar a {
  float: left; /* Make sure that the links stay side-by-side */
  display: block; /* Change the display to block, for responsive reasons (see below) */
  color: white; /* White text color */
  text-align: center; /* Center the text */
  padding: 14px 20px; /* Add some padding */
  text-decoration: none; /* Remove underline */
}

/* Right-aligned link */
.navbar a.right {
  float: right; /* Float a link to the right */
}

/* Change color on hover/mouse-over */
.navbar a:hover {
  background-color: #ddd; /* Grey background color */
  color: black; /* Black text color */
}
<div class="row">
  <div class="side">
<h2> Credits: </h2>
<p> The author of the project is <h3> Ivan Naborshchikov</h3>.  </p>
<p> If you have any questions, you can ask them via my pages on Linkedin, Github: </p> </div>
  <div class="main">
<h2> Description: </h2>
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

<p> <h4> 2. class Newnewproton()</h4>   - Calculation of three shells for new quarks "u" and "d", nine shells for a new proton, a new neutron in terms of mass, charge, volume. The calculation is implemented according to the author's algorithm. The applied algorithm allows to calculate any known particles and find new ones. This algorithm is only one of three parts of the physical-mathematical model. </p>

<p> <h4> 3. class Pseudonewneutron() </h4>   - Calculation of three shells for quarks "u" and "d", nine shells for a proton, a neutron in terms of mass, charge, volume. The calculation is implemented according to the author's algorithm. According to the calculation, the proton and neutron have three wormholes each. Table 10 contains data on three negative values of the proton shells. Table 11 contains data on three negative values of the neutron shells. A negative volume is that part of the volume that is in another dimension, that is, a part of a wormhole located in another dimension. </p>

<p> According to the article of famous physicists A. Yu. Khlestkov & Yu. A. Khlestkov from Russian National Research Nuclear University "MEPhI" (https://link.springer.com/article/10.1007/s11182-019-01709-9), the neutron has a wormhole, which can be considered as proof of the author's algorithm. </p>
<p> <h4> 4. class Wavep(), class Waven(), class Wavepsn(), class Wavepsp() </h4>   - Classes for calculating the Compton wavelength from shells for a new proton, new neutron, proton, neutron. </p>
<p> <h4> 5. class ElectricWavep(), class ElectricWaven(), class ElectricWavepsn(), class ElectricWavepsp() </h4>  - Classes for calculating the electrical ε for each shell at the speed of light. </p>
<p> <h4> 6. class GravityWavep(), class GravityWaven(), class GravityWavepsn(), class GravityWavepsp() </h4>   - Calculation of gravity characteristics for shells. </p>
<p> <h4> 7. class Frequencyp(), class Frequencyn(), class Frequencypsn(), class Frequencypsp() </h4>   - Calculation of frequency characteristics for shells taking into account the speed of light in vacuum. </p>
<p> <h4> 8. class Melectron() </h4>   - Calculation of the existence of a point charged particle with  a charge 60 times less than the charge of an electron. Particle electric charge microminus (author's name) - -2.6702944e-21C. Microplus particle electric charge (author's name) -  2.6702944e-21C. Microplus and microminus mass (Kg): 1.518230622602297e-32. Volume microplus and microminus (cbm): 8.726646292652101e-66. </p>
<p> <h4> 9. class NewnewprotonCycles(), class NewneutronCycles(), class PsnewneutronCycles(), class PsnewnewprotonsCycles() </h4>   - Calculation of the characteristics of a new proton, new neutron, neutron, proton, by phases within the cycle. The new proton has five phases, the rest have four phases in a cycle. <br> </p>
<p>In the Project, the algorithm was tested for a proton, a neutron in terms of electric charge after the calculations. The result of testing the algorithm showed no error. </p>

<p>The result of the assessment using the algorithm of a new proton, a new neutron showed the presence of an electric charge in 2.403264951e-20 C for a new neutron and revealed differences between a proton and a new proton. The presence of a small charge in the new neutron and the absence of an electric charge in the neutron speaks of the validity of the Dirac equation without exception. </p>

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
</div>
</div>
/* Ensure proper sizing */
* {
  box-sizing: border-box;
}

/* Column container */
.row {
  display: flex;
  flex-wrap: wrap;
}

/* Create two unequal columns that sits next to each other */
/* Sidebar/left column */
.side {
  flex: 20%; /* Set the width of the sidebar */
  background-color: #f1f1f1; /* Grey background color */
  padding: 20px; /* Some padding */
}

/* Main column */
.main {
  flex: 80%; /* Set the width of the main content */
  background-color: white; /* White background color */
  padding: 20px; /* Some padding */
}
/* Responsive layout - when the screen is less than 700px wide, make the two columns stack on top of each other instead of next to each other */
@media screen and (max-width: 700px) {
  .row {
    flex-direction: column;
  }
}

/* Responsive layout - when the screen is less than 400px wide, make the navigation links stack on top of each other instead of next to each other */
@media screen and (max-width: 400px) {
  .navbar a {
    float: none;
    width: 100%;
  }
}
<div class="footer">
  <h2><p> Only part of the Project is posted on the site in the form of program code due to the specifics of its capabilities. </p>
</h2>
</div>
.footer {
  padding: 20px; /* Some padding */
  text-align: center; /* Center text*/
  background: #ddd; /* Grey background */
}
</body>
</html>


