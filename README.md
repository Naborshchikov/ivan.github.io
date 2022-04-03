# ivan.github.io
<html>
<head>
<title>
Page README <br>
Project name: <br> 
<h1> Keywords: python3, tachyon, proton, neutron, molecular hydrogen, hydrogen bond, helium 2, helium 3, time, space-temporal, neutron radiation, thermonuclear fusion. </h1>
</title>
</head>
<body>
<h2> <p> Python 3 program for segmentation of protons and neutrons by physical properties to build models of atoms and molecules using hydrogen and helium as an example. Ver. 6.2 </p>
Description: </h2>
<p> The python 3 program made it possible to identify two types of protons, neutrons and to single out three segments in each of them for electric charge, mass, volume. </p>
<p> The program uses publicly available data that is contained on the NIST (US National Institute of Standards and Technology) website. Sometimes validated scientific data from open sources were used in addition to NIST data. </p>
<p> The program contains the mathematical conversion of a neutron to a proton. This can be used to create a program to reduce energy costs for modifying a substance, eliminating or reducing neutron radiation. </p>

<p> The computed structures protons and neutrons indicate that: a proton and a neutron move along the time axis; a proton 2 and a neutron 2 move in the opposite direction. On the one hand, it contributes to the formation of atoms, molecules, on the other hand, it contributes to their destruction. </p>
<p> Protons, neutrons consist of three segments located in the future, present, and past time. As proof of this, one can take the fact that scientists have confirmed that the proton and neutron have a passage to another dimension: https://link.springer.com/article/10.1007/s11182-019-01709-9 </p>
<p> The presented segmentation by the time of the proton and the neutron explains the reason for the decay of the neutron, the stability of the proton. This approach can be extrapolated to modeling materials by calculating their service life in advance, which will improve the durability and quality of the developed materials. </p>
<p> The project has produced two types of neutrons with different characteristics. This explains why the "bottle" and "beam" methods for estimating the neutron lifetime give scientists different results. This can be viewed as proof of the obtained results. </p>
<p> Protons, neutrons have the form of a dipole. This has been confirmed by a large group of scientists experimentally: https://doi.org/10.1063/1.4967465 The project contains a theoretical basis for this phenomenon. The calculation showed the presence of “complex” dipoles. </p>
<p> The use of the wile cycle made it possible to reveal an ultra-small particle - presumably a tachyon. </p>
<p> Characteristics: </p>
<p> Electric charge: 9.027796614315168e-36 C </p>
<p> Mass: 5.1328712618966e-47 kg </p>
<p> Volume: 2.9503259453108798e-80 m.cube </p>
<p> This information can be used to develop computer models that analyze the processes occurring between quarks, aging and rejuvenation, and the possibility of temporal transitions. </p>
<p> The implementation of the program showed that there are 8 hydrogen molecules with different characteristics, two different helium atoms 2, five helium atoms 3. </p>
<p> It was also found that one common overlap in proton shells leads to the formation of a molecular bond, while two common overlaps in proton shells lead to the formation of an atom nucleus. </p>
<p> The results obtained allow us to better understand: </p>
<p> In the field of chemistry: </p>
<p> 1. The diversity of hydrogen bonds, intra- and intermolecular interactions, the physical properties of water and many organic liquids (alcohols, carboxylic acids, amides of carboxylic acids, esters). </p>
<p> 2. The presented knowledge makes it possible to solve one of the important problems of modern chemistry - the creation of an efficient and safe hydrogen storage system, to develop a system that combines a high hydrogen storage density with low costs for its release. </p>
<p> 3. We will be able to better understand the properties of biologically important substances such as proteins and nucleic acids. In particular, elements of the secondary structure (for example, α-helices, β-folds) and tertiary structure in protein, RNA and DNA molecules are stabilized by hydrogen bonds. In these macromolecules, hydrogen bonds hold parts of the same macromolecule together, causing it to fold into a particular shape. </p>
<p> 4. Many polymers are reinforced with hydrogen bonds in their backbones. Among synthetic polymers, the most famous example is nylon, where hydrogen bonds play a major role in the crystallization of the material. Hydrogen bonds are also important in the structure of artificially produced polymers (such as cellulose) and in many different forms in nature, such as wood, cotton, and linen. </p>
<p> For each of these directions in the field of chemistry, new, more efficient computer models can be developed. </p>
<p> The results obtained allow us to assert that the calculated He2_1 will decay into two protons, and He2_2 will give deuterium upon decay. This is clearly seen from the obtained structures and their characteristics. Detailed script can be executed additionally. </p>
<p> The presence of five different types of helium 3 atoms allows us to state that there are at least five different technological regimes for producing helium 3, for the process of thermonuclear fusion with obtaining helium 3. Computer programs for these technologies can be developed. </p>
<p> The obtained results will reduce the temperature required for thermonuclear fusion hundreds of times, achieve a neutron-free reaction. This opens up the possibility of developing an individual environmentally-friendly electric fusion generator for vehicles, buildings. </p>

<p> The error of the data obtained using the program version 6.2 tends to zero in comparison with the publicly available experimental data. For example, if we compare the data on quarks "u", we get an error equal to zero.  For example, the value of the electric charge for quarks “u” and “u2” according to the program is: </p>
<p> 1.0681177559999997e-19 C (Tuxq02 from class Algorithm()) </p>
<p> 1.068117756e-19 C (Tuxq13 from class Algorithm()) </p>

<p> Combinatorics, matrix calculus, cycles, logical calculus are important parts of the program code. </p>
<p> Two scripts are attached for visualization:  </p>
<p> Distribution of a hydrogen molecule by segments in (%). </p>
<p> Distribution of a He2_1 by segments in (%). </p>
<p> For scientific analysis, use data from program version 6.2 please. For example, He2_2_Past_1_e and He2_2_Past_2_e are two different segments.  </p>
<p> The laws of nature obey mathematical laws. Therefore, they can be digitized and represented in the form of computer code. The more correct the algorithm is chosen for the computer code, the smaller the difference between the calculated data and the well-known experimental data used will be, which is typical for the proposed software code.  </p>
<p> <h3> Project implemented in python ver. 3, with using anaconda 2020 02. </h3> <br>
To implement the project used: <br>
The libraries: numpy; scipy; prettytable; collections. For the attached separately scripts you will need: matplotlib.pyplot, matplotlib <br>  </p>
<p><h3> The program includes 5 classes, table with the original data placed separately: </h3> </p>
<p> The program code contains detailed comments and explanations. </p>
<p> <h4> Table with the original data </h4> Table#1. The table contains the data recommended by the US NIST and data of various groups of scientists. The table is made in the form of a separate script. </p>
<p> <h4> 1. class Algorithm () </h4>   - The 'class Algorithm()' class is used to calculate nuclei, shells of quarks 'u', 'd', and contains the data of constants used in the program. </p>
<p> <h4> 2. class Particles () </h4>   - This class forms the data set for protons, neutrons, calculates the characteristics of the tachyon. </p>
<p> <h4> 3. class Segments() </h4> - This class defines segments of protons, neutrons conditionally located in the past, present and future time. The values that can be obtained from this class are very simple and do not require additional comments. </p>
<p> For example: PROTON_Present_Q means it's a proton, the segment is present, and Q is the electric charge for this segment. The value of M at the end will talk about the mass. The V value at the end will indicate volume. The “matrix” value at the end will indicate the content of the matrix of values for the shells for the specified segment.  </p>
<p> <h4> 4. class Molecularhydrogen() </h4> - This class calculates the electric charge, mass, volume of eight different hydrogen molecules segment by segment. </p>
<p> Consider the notation using the MolecularH2_1_Past_1e as an example: </p>
<p> MolecularH2 - hydrogen molecule H2 </p>
<p> 1 - hydrogen molecule of the first type </p>
<p> Past - segment Past </p>
<p> 1 - first layer </p>
<p> e - the presence of an electron is taken into account </p>
<p> The absence of “e” at the end of the symbol means the absence of an electron. </p>
<p> <h4> 5. class He() </h4> - This class calculates the electric charge, mass, volume of two helium atoms 2 by segments; the presence of five helium atoms 3 and characteristics of common areas for protons, neutrons. The calculation was made taking into account each shell. The “intersection” extension denotes that the value refers to a common area of particles. </p>
<p>  He2_2_Present_4e means that this is helium 2, the second kind of atom, the fourth part of the segment is “present”, the data on the electron are taken into account. </p>
<p> <h4> 6. class Deltam() </h4> - This class explains "bond energy". It is taken from two types of neutrons having different weights. This class also introduces: </p>
<p> NEUTRON_Delta is used for exothermic reactions </p>
<p> NEUTRON_Delta2 is used in endothermic reactions </p>
<p> The nucleus of an atom, a molecule has binding energies, the formation of which or their destruction leads to the release of energy or its absorption. Chemical and nuclear reactions are like a handshake. For molecules, this is a "handshake with one hand", for particles in the nucleus, this is a "handshake with two hands". The destruction or formation of these bonds and causes heat or cold. The transition from "handshake with two hands" to "handshake with one hand" is the transition from nuclear to chemical bonding. </p>
<h2> Installation: </h2>
<p> The project does not require any special skills to install. It will work on Jupiter Notebook. If you are missing some packages, then they are easily loaded using pip. </p>
<h2> Usage:  </h2>
<p> To start the Project, it is enough to use Run. Everything will be done automatically. </p>
<h2> Credits: </h2>
<p> The author of the project is <h3> Ivan Naborshchikov</h3>.  </p>
<p> If you have any questions, you can ask them via: </p>
<p> e-mail: ivannaborshchikov@yahoo.com </p>
<p> Information about me is located at this address: <br>
<h4>  http://www.linkedin.com/in/ivan-naborshchikov/ </h4> </p>
 </body>
</html>
