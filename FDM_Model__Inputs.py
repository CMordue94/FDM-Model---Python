"""
--------------------------------------INTRODUCTION--------------------------------------

NAME:
FDM MODEL - INPUTS

AUTHOR:
Christopher Mordue

DESCRIPTION:
This model utilises the FDM to find the vertical displacement, rotational displacement and spring constant of nonprismatic
beams. The model can be split into two parts: Mechanical & Thermal. The former calculates the displacements and spring constant
for external force application to the beam. The latter calculates the displacements from external temperature change (this occurs
due to different materials present in the beam with different thermal expansion coefficients). Each of these two parts draw upon 
the 'FDM_Model__Inputs' file, but they each have particular features enabling complex, non-prismatic beams to be analysed in various
conditions. These are hence listed below:

Mechanical Features:
- Substrate Linear Width Taper
- Substrate Gap/Width Change (x2)
- Composite Beams (x2 Materials)
- Cantilever Angled Orientation - Constant Angle 
- Cantilever Angled Orientation - Changinging Angle/Pre-Bend
- Cantilever Part Angled Orientation
- Offset Force 
- Variable Mesh (x2 Meshes)
- Sensitivity Analysis

Thermal Features: 
- Temperature Change Distribution Array
- Coating Linear Width Taper
- Coating Gap/Width Change x2
- Cantilever Angled Orientation - Constant Angle
- Cantilever Angled Orientation - Changing Angle/Pre-Bend
- Cantilever Part Angled Orientation
- Variable Temperature Distribution 
- Variable Mesh (x2 Meshes)
- Sensitivity Analysis
"""



"""
--------------------------------------IMPORTS--------------------------------------
"""
import numpy as np



"""
--------------------------------------OPERATIONS--------------------------------------
"""


"""Spring Constant Output""" #("Y" or "N")
Spring="Y"
 
"""Displacement""" #("Y" or "N")
Displacement="Y"

"""Thermal (Bimetallic) Displacement""" #("Y" or "N")
Th_Displacement="N" 

"""Sensitivity Analysis""" #("Y" or "N")
SensitivityAnalysis="N" 

"""IO Save""" #("Y" or "N")
I_O="N"


"""
--------------------------------------SENSITIVITY ANALYSIS INPUTS--------------------------------------
"""

#Analysed Input (Equals the Algebraic Term with Quotes)
SI="_"                 
#Start Value
SS=0 
#End Value
SE=0
#Increments
Si=100



"""
--------------------------------------OVERALL ANALYSIS INPUTS--------------------------------------
"""

"""General Inputs"""
#Number of elements
Ne=200
#Number of Nodes/Length of Arrays (=len(array[Ne]))
Nn=Ne+1

"""Advanced Inputs"""
"This creates two different meshes along the cantilever where one is stated, and the other calculated according to the remaining elements and length"
###Variable Mesh ("Y" or "N")
VM="N"
#Stated Mesh Start X-Point
SM_S=0
#Stated Mesh End X-Point
SM_E=90e-6
#Staed Mesh Element Number
SM_Ne=20



"""
--------------------------------------SUBSTRATE INPUTS--------------------------------------
"""

"""Substrate General Inputs"""
#Young's Modulus (Pa)
E=2.315e+11 
#Poisson's Ratio
vs=0.255
#Substrate X Length (m)
Ls_x=150e-6 
#Substrate thickness (m)
ts=5.0e-7 
#Substrate Width (m)
ws=120.0e-6 

"""Substrate Tapered Width Inputs"""
#Tapered Substrate ("Y" or "N")
T="Y"
#Start of Taper on X-Axis (m)
TS=90e-6
#End of Taper on X-Axis (m)
TE=Ls_x
#Width at End of Taper (m)
Tw=0
#Taper Gradient
if T != "Y":
    TM=0
else:
    TM=(ws-Tw)/(TS-TE)

"""Substrate Gap/Width Change Inputs 1"""
#Gap/Width Change ("Y" or "N")
G1="N"
#Start of Gap/Width Change on X-Axis (m)
G1S=0e-6
#End of Gap/Width Change on X-Axis (m)
G1E=80e-6
#Gap/Width Change at Start (m)
G1Sw=10e-6
#Gap/Width Change at End (m)
G1Ew=10E-6
#Gap Gradient
if G1 != "Y":
    G1M=0
else:
    G1M=(G1Ew-G1Sw)/(G1E-G1S)

"""Substrate Gap/Width Change Inputs 2"""
#Gap/Width Change ("Y" or "N")
G2="N"
#Start of Gap/Width Change on X-Axis (m)
G2S=0
#End of Gap/Width Change on X-Axis (m)
G2E=0
#Gap/Width Change at Start (m)
G2Sw=0
#Gap/Width Change at End (m)
G2Ew=0
#Gap Gradient
if G2 != "Y":
    G2M=0
else:
    G2M=(G2Ew-G2Sw)/(G2E-G2S)



"""
--------------------------------------COATING INPUTS--------------------------------------
"""

"""Coating General Inputs"""
#Coating ("Y" or "N")
Coating="Y"
#Young's Modulus (Pa)
Ec=8e+10
#Poisson's Ratio
vc=0.42
#Coating Start X-Point (m)
cS=0
#Coating End X-Point (m)
cE=90e-6
#Coating Thickness (m)
tc=150e-9
#Coating Width (m)
wc=50e-6

"""Coating Tapered Width Inputs"""
#Tapered Coating ("Y" or "N")
Tc="N"
#Start of Taper on X-Axis (m)
TcS=90E-6
#End of Taper on X-Axis (m)
TcE=Ls_x
#Width at End of Taper (m)
Tcw=0
#Taper Gradient
if Tc != "Y":
    TcM=0
else:
    TcM=(wc-Tcw)/(TcS-TcE)

"""Coating Gap/Width Change Inputs 1"""
#Gap/Width Change ("Y" or "N")
Gc1="N"
#Start of Gap/Width Change on X-Axis (m)
Gc1S=0
#End of Gap/Width Change on X-Axis (m)
Gc1E=90e-6
#Gap/Width Change at Start (m)
Gc1Sw=40e-6
#Gap/Width Change at End (m)
Gc1Ew=40e-6
#Gap Gradient
if Gc1 != "Y":
    Gc1M=0
else:
    Gc1M=(Gc1Ew-Gc1Sw)/(Gc1E-Gc1S)

"""Coating Gap/Width Change Inputs 2"""
#Gap/Width Change ("Y" or "N")
Gc2="N"
#Start of Gap/Width Change on X-Axis (m)
Gc2S=110E-6
#End of Gap/Width Change on X-Axis (m)
Gc2E=148E-6
#Gap/Width Change at Start (m)
Gc2Sw=30E-6
#Gap/Width Change at End (m)
Gc2Ew=0
#Gap Gradient
if Gc2 != "Y":
    Gc2M=0
else:
    Gc2M=(Gc2Ew-Gc2Sw)/(Gc2E-Gc2S)



"""
--------------------------------------CANTILEVER INPUTS--------------------------------------
"""

"""Angled Inputs - Constant Angle"""
#Angled Cantilever ("Y" or "N")
AC="N"
#Angle (Degrees)
Ad=13
#Angle (Radians)
Ar=(np.deg2rad(Ad) if AC == "Y" else 0)

"""Angled Inputs - Changing Angle/Pre-Bend"""
#Angled Cantilever ("Y" or "N")
ACCh="N"
#Changing Angle Function (state function as a string)
AChr="(1.4*n_x**0.4)"                                       #Made a string as n_x has not been defined yet

"""Part Angled Inputs"""
#Part Angled Cantilever ("Y" or "N")
PAC="Y"
#Start of Angled Part on X-Axis
PS=144E-6
#End of Angled Part on X-Axis
PE=Ls_x
#Angle (Degree)
PAd=56

PAr=np.deg2rad(PAd)


"""
--------------------------------------MECHANICAL INPUTS--------------------------------------
"""

#Force (N)
F=100.0e-9
#Offset Force ("Y" or "N")
OF="N"
#Force Application X-Point (m)
OF_Px=145E-6



"""
--------------------------------------THERMAL INPUTS--------------------------------------
"""

"""General"""
#Substarte Thermal Expansion Coefficient (K-1)
ThEs=2.6e-6
#Coating Thermal Expansion Coefficient (K-1)
ThEc=1.43e-5
#Temperature Change (K)
dT=10

"""Advanced"""
#Variable Temperature Distribution ("Y" or "N")
VTD="N"
#Distribution Function (state function as a string)
VTD_F="0"










