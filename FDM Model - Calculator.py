"""
--------------------------------------INTRODUCTION--------------------------------------

NAME:
FDM MODEL - CALCULATOR

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

print 'FDM Model - Calculator\n'

"""
--------------------------------------IMPORTS--------------------------------------
"""

from FDM_Model__Inputs import *
import shutil as sh
import os
import time as t
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from IPython import get_ipython as g_ipy
g_ipy().run_line_magic('matplotlib', 'qt5')
#from matplotlib import rc
#rc('font',**{'family':'serif','serif':['Computer Modern Roman']}) #Alters matplotlib's default font to one used in Latex
#rc('text', usetex=True) #Executes the alteration



"""
--------------------------------------INPUT CALCULATIONS--------------------------------------
"""

"""OVERALL ANALYSIS INPUT CALCULATIONS"""
"Basic"
#Base Element X-Length (m)
Le_x=Ls_x/Ne
#Node X-Coordinates (m)
n_x=(np.arange(0,Nn,1)*Le_x if VM != "Y" else np.zeros(Nn))
#Node Elements Centroid X-Coordinate (m)
n_e_xc=((n_x)-(Le_x/2) if VM != "Y" else np.zeros(Nn))
n_e_xc[0]=0                                                                                 #Removes coordinate for non-existent element
#Node Along-Cantilever Coordinate (m)
n_A=(n_x if PAC != "Y" else np.zeros(Nn))
#Base Element X-Length Array
Le_xa=np.ones(Nn)*Le_x
Le_xa[0]=0
#Base Cantilever Length (m)
Ls_A=Ls_x
"Advanced"
if VM == "Y":
    #Remaining Mesh Elements
    Ner=Ne-SM_Ne
    #Remaining Mesh X-Length (m)
    Ls_xr=Ls_x-(SM_E-SM_S)
    #Stated Mesh Element X-Length (m)
    SM_Le_x=(SM_E-SM_S)/SM_Ne
    #Remaining Element X-Length (m)
    RM_Le_x=Ls_xr/Ner
    #Node X-Coordinates (m) & Node Elements Centroid X-Coordinate (m)
    for i in range(1,Nn):
        if (i-1)*SM_Le_x >= SM_S and (i-1)*SM_Le_x < SM_E:
            n_x[i]=n_x[i-1]+SM_Le_x
            n_e_xc[i]=n_x[i]-(SM_Le_x/2)
        else:
            n_x[i]=n_x[i-1]+RM_Le_x
            n_e_xc[i]=n_x[i]-(RM_Le_x/2)
    #Variable Element X-Length Array (m)
    Le_xa=np.asarray([SM_Le_x if n_e_xc[i]>SM_S and n_e_xc[i]<=SM_E else RM_Le_x for i in range(0,Nn)])
    Le_xa[0]=0
    #Element Change Indicator
    VM_N=([i for i in range(1,Nn) if n_e_xc[i-1]<=SM_E and n_e_xc[i]>SM_E])[0]


"""SUBSTRATE INPUT CALCULATIONS"""                                                            
"Indicators"
#Taper Indicator
TS_N=((([i for i in range(0,Ne) if n_e_xc[i]<=TS and n_e_xc[i+1]>TS])[0]) if T == "Y" else 0)
#Gap 1 Indicator
G1S_N=((([i for i in range(0,Ne) if n_e_xc[i]<=G1S and n_e_xc[i+1]>G1S])[0]) if G1 == "Y" else 0)
G1E_N=((([i if n_e_xc[i-1]<=G1E and n_e_xc[i]>G1E else Nn for i in range(G1S_N,Nn)])[0]) if G1 == "Y" else 0)                                                    
#Gap 2 Indicator
G2S_N=((([i for i in range(0,Ne) if n_e_xc[i]<=G2S and n_e_xc[i+1]>G2S])[0]) if G2 == "Y" else 0)
G2E_N=((([i if n_e_xc[i-1]<=G2E and n_e_xc[i]>G2E else Nn for i in range(G2S_N,Nn)])[0]) if G2 == "Y" else 0)                                                    
"Substrate Element Widths"
#Pristmatic Widths (m)
n_e_ws=np.ones(Nn)*ws                                                     #Removes width for non-existent element
n_e_ws[0]=0      
#+Tapered Widths (m)             
if T == "Y":
    for i in range(TS_N,Nn):
        if n_e_xc[i]>TS and n_e_xc[i]<=TE:
            n_e_ws[i]=(TM*(n_e_xc[i]-TS))+ws
        elif n_e_xc[i]>TE:
            n_e_ws[i]=Tw 
#+Gap 1 Widths (m)
if G1 == "Y":
    for i in range (G1S_N,G1E_N):
        if n_e_xc[i]>G1S and n_e_xc[i]<=G1E:
           n_e_ws[i]=n_e_ws[i]-G1Sw-(G1M*(n_e_xc[i]-G1S))
#+Gap 2 Widths (m)
if G2 == "Y":
    for i in range (G2S_N,G2E_N):
        if n_e_xc[i]>G2S and n_e_xc[i]<=G2E:
           n_e_ws[i]=n_e_ws[i]-G2Sw-(G2M*(n_e_xc[i]-G2S))
"Substrate Element Second Moment of Area (m^4)"
n_e_I=n_e_ws*ts**3/12. 
n_e_I[0]=np.Infinity                                                                        #Removes I for non-existent element


"""COATING INPUT CALCULATIONS"""
if Coating == "Y":
    "Indicators"
    #Taper Approx. Indicator
    TcS_N=((([i for i in range(0,Ne) if n_e_xc[i]<=TcS and n_e_xc[i+1]>TcS])[0]) if Tc == "Y" else 0)
    #Gap 1 Approx. Indicator
    Gc1S_N=((([i for i in range(0,Ne) if n_e_xc[i]<=Gc1S and n_e_xc[i+1]>Gc1S])[0]) if Gc1 == "Y" else 0)
    Gc1E_N=((([i if n_e_xc[i-1]<=Gc1E and n_e_xc[i]>Gc1E else Nn for i in range(Gc1S_N,Nn)])[0]) if Gc1 == "Y" else 0) 
    #Gap 2 Approx. Indicator
    Gc2S_N=((([i for i in range(1,Nn) if n_e_xc[i-1]<=Gc2S and n_e_xc[i]>=Gc2S])[0]) if Gc2 == "Y" else 0)
    Gc2E_N=((([i if n_e_xc[i-1]<=Gc2E and n_e_xc[i]>Gc2E else Nn for i in range(Gc2S_N,Nn)])[0]) if Gc2 == "Y" else 0) 
    "Coating Node Widths"
    #Pristmatic Widths (m)
    n_e_wc=np.zeros(Nn)
    n_e_wc=np.asarray([wc if round(n_e_xc[i],8)>cS and round(n_e_xc[i],8)<cE else 0 for i in range(0,Nn)])                                                                                                                                
    #+Tapered Widths (m)              
    if Tc == "Y":
        for i in range(TcS_N,Nn):
            if n_e_xc[i]>TcS and n_e_xc[i]<=TcE:
                n_e_wc[i]=(TcM*(n_e_xc[i]-TcS))+wc
            elif n_e_xc[i]>TcE:
                n_e_wc[i]=Tcw 
    #+Gap 1 Widths (m)
    if Gc1 == "Y":
        for i in range (Gc1S_N,Gc1E_N):
            if n_e_xc[i]>Gc1S and n_e_xc[i]<=Gc1E:
               n_e_wc[i]=n_e_wc[i]-Gc1Sw-(Gc1M*(n_e_xc[i]-Gc1S))
    #+Gap 2 Widths (m)
    if Gc2 == "Y":
        for i in range (Gc2S_N,Gc2E_N):
            if n_e_xc[i]>Gc2S and n_e_xc[i]<=Gc2E:
               n_e_wc[i]=n_e_wc[i]-Gc2Sw-(Gc2M*(n_e_xc[i]-Gc2S))
    "Composite Mechanical Calculations" 
    if Spring == "Y" or Displacement == "Y": 
        #Young's Modulus Ratio
        n=Ec/E
        #Bending Axis Location (m)
        n_e_xb=np.zeros(Nn)
        n_e_xb[1:Nn+1]=((n_e_ws[1:Nn+1]*ts**2/2)+(n*n_e_wc[1:Nn+1]*tc*ts)+(n*n_e_wc[1:Nn+1]*tc**2/2))/((n_e_ws[1:Nn+1]*ts)+(n*n_e_wc[1:Nn+1]*tc))
        #Second Moment of Area
        n_e_I=(n_e_ws*n_e_xb**3/3)+(n_e_ws*(ts-n_e_xb)**3/3)+((n*n_e_wc*tc**3/12)+(n*n_e_wc*tc*((tc/2)+ts-n_e_xb)**2))
        n_e_I[0]=np.Infinity                                                                    #Removes I for non-existent element
 
    
"""CANTILEVER INPUT CALCULATIONS"""
"Angled Inputs - Changing Angle/Pre-Bend"
AChr=(eval(AChr) if ACCh == "Y" else np.zeros(Nn))
"Part Angled Inputs"
PAr=np.zeros(Nn)                                                                        #Generate base Angle Array which is always needed
if PAC == "Y":                                           
    #Total Cantilever Along-Length (m)
    Ls_A=(Ls_x-(PE-PS))+((PE-PS)/np.cos(np.deg2rad(PAd)))
    #Part Angled End Along-Length (m)
    Ls_EA=PS+((PE-PS)/np.cos(np.deg2rad(PAd)))                                          #This may rectify the limitation of this code being for the part andgle change dynamic element length only workign for the part angle occuring to the end of the strucutre
    #Part Angled Substrate Element X-Length (m)
    PS_N=[i for i in range(1,Nn) if n_e_xc[i-1]<=PS and n_e_xc[i]>PS][0]                
    PE_N=[i if n_e_xc[i-1]<=PE and n_e_xc[i]>PE else Nn for i in range(PS_N,Nn)][0]
    Pe_x=(Ls_EA-n_x[PS_N-1])/(PE_N-PS_N)
    #Part Angled Element X-Length Array
    Le_xa=np.asarray([Pe_x if n_e_xc[i]>PS and n_e_xc[i]<=PE else Le_xa[i] for i in range(0,Nn)])
    #Part Angled Angle Array
    PAr=np.asarray([np.deg2rad(PAd) if n_e_xc[i]>=PS and n_e_xc[i]<=PE else 0 for i in range(0,Nn)])                                                              
    #Part Angled Node Along-Cantilever Coordinates (m)
    for i in range(1,Nn):                                                               #Currently unsure if fully working, may be better to scrap the above and make the element number equal to the point in change on probes (i.e. make n_x generate a node point on the point of change)!
        n_A[i]=n_A[i-1]+Le_xa[i]                                                        #This element length improves end displacement accuracy when angled to the end (which is always the case for our probes) as it alters the element length to the known geometry end. However, it does not improve the end displacement calculation if not angled to the end and should be borne in mind if not the case! 


"""MECHANICAL INPUT CALCULATIONS"""                                                                                 
if Spring == "Y" or Displacement == "Y":
    "Offset Force Calculations"
    if OF == "Y":
        #Offset Force Application Along-Cantilever point (m)
        OF_PA=((PS+((OF_Px-PS)/np.cos(PAr))) if PAC == "Y" and OF_Px>PS and OF_Px<=PE else OF_Px)
        #Offset Force Zero-ing Array
        OF_0=[0 for i in range(0,Nn) if n_x[i]>OF_Px]                                   #CHECK THIS CHANGE!!!
    "Force Application Point"
    #Force Application X-Point
    F_Px=(OF_Px if OF == "Y" else Ls_x)
    #Force Application Along-Point
    F_PA=(OF_PA if OF == "Y" else Ls_A)
    #Force Application Point Array
    F_PAa=np.asarray([F_PA if PAC == "Y" and n_e_xc[i]>PS and n_e_xc[i]<=PE else F_Px for i in range(0,Nn)])


"""THERMAL INPUT CALCULATIONS"""
if Th_Displacement == "Y":
    #Total Cantielver Thickness
    tT=ts+tc
    #Substrate Biaxial Young's Modulus (Pa)
    EsB=E/(1-vs)
    #Coating Biaxial Young's Modulus (Pa)
    EcB=Ec/(1-vc)
    #Biaxial Modulus Transformation Ratio (Coating to Substrate)
    nB=EcB/EsB
    #Bending Axis Location from Underside (m)
    Bax=((ws*ts**2/2)+(nB*wc*tc*ts)+(nB*wc*tc**2/2))/((ws*ts)+(nB*wc*tc))
    #Flexure Rigidity (Nm^2)
    D=EsB*((Bax**3/3)+((ts-Bax)**3/3)+((nB*tc**3/12)+(nB*tc*((tc/2)+ts-Bax)**2)))
    #Effective Young's Modulus for Different Widths (Pa)
    Eeff=12*D/tT**3
    #Effective Biaxial Modulus Transformation Ratio (Coating to Substrate)
    neff=EcB/Eeff
    #Temperature Distribution Array (K)
    Ta=(np.ones(Nn)*dT if VTD != "Y" else eval(VTD_F))
    Ta[0]=0                                                                                     #Removes temperature in non-existent element                                                 



"""
--------------------------------------MECHANICAL OUTPUT CALCULATION FUNCTIONS--------------------------------------
"""

"""SPRING CONSTANTS"""
#Spring Constant Component of Y Displacement
def k_UY():
    k_UY=(Le_xa**2/(E*n_e_I)*((Le_xa/3)+((F_PAa-n_A)/2)))*np.cos(PAr+Ar+AChr)**2
    if OF == "Y":
        k_UY[Nn-len(OF_0):Nn]=OF_0
    return k_UY
#Spring Constant Component for Rotational Y Displacement
def k_R_UY():
    k_RXY=(Le_xa**2/(E*n_e_I)*((Le_xa/2)+F_PAa-n_A))*np.cos(PAr+Ar+AChr)**2
    if OF == "Y":
        k_RXY[Nn-len(OF_0):Nn]=OF_0
    k_R_UY=np.zeros(Nn)
    for i in range(1,Nn):
        if VM == "Y" and Le_xa[i] != SM_Le_x:
            k_R_UY[i]=(np.sum(k_RXY[0:VM_N])/SM_Le_x*RM_Le_x)+np.sum(k_RXY[VM_N:i])
        else:
            k_R_UY[i]=np.sum(k_RXY[0:i])                           #This code appears to produce ever so slightly different values to excel, however it is close to negligible difference! This difference could be due to the lsight difference in the carrying of the numbers mathematically which has been relaised in the past.
    return k_R_UY
#Total Spring Constant
def ks_e():    
    kT=k_UY()+k_R_UY()
    kT[0]=np.Infinity                                                       #Removes division by zero error
    k_e=1/kT
    k_e[0]=np.Infinity                                                      #Remove division by zero error
    return k_e
#Along Substrate/Cantilever Spring Constant
def ks_Along():
    ks_a=np.zeros(Nn)
    for i in range(1,Nn+1):
        ke=1/ks_e()
        ke[i:Nn]=0
        ks_a[0]=0                                                           #Remove division by zero error
        if np.sum(ke) == 0:                                                 
            ks_a[i-1]=0                                                     #Remove division by zero error
        else:
            ks_a[i-1]=1/np.sum(ke)                                          
    return ks_a
#Substrate/Cantilever Spring Constant
def ks_Tot(): 
    ks_Tot=1/np.sum(1/ks_e())
    return ks_Tot


"""DISPLACEMENTS"""
#Y Displacement
def UY():
    UY=((F*Le_xa**3/(3*E*n_e_I))+(F*(F_PAa-n_A)*Le_xa**2/(2*E*n_e_I)))*np.cos(PAr+Ar+AChr)**2
    if OF == "Y":
        UY[Nn-len(OF_0):Nn]=OF_0
    return UY
#Rotational Displacement 
def RXY():
    R=((F*Le_xa**2/(2*E*n_e_I))+(F*(F_PAa-n_A)*Le_xa/(E*n_e_I)))*np.cos(PAr+Ar+AChr)**2
    if OF == "Y":
        R[Nn-len(OF_0):Nn]=OF_0
    return R
#Rotational Y Displacement
def R_UY(): 
    R=RXY()
    R_UY=np.zeros(Nn)
    UY=np.tan(R)*Le_xa
    for i in range(1,Nn):
        if VM == "Y" and Le_xa[i] != SM_Le_x:
            R_UY[i]=(np.sum(UY[0:VM_N])/SM_Le_x*RM_Le_x)+np.sum(UY[VM_N:i])
        else:
            R_UY[i]=np.sum(UY[0:i])        
    return R_UY
#Total Y Displacement
def UY_Tot(): 
    UY_T=UY()+R_UY()
    return UY_T
#Along Substrate/Cantilever Y Displacement
def UY_Along():
    UY_A=np.zeros(Nn)  
    for i in range(1,Nn):
        UY_A[i]=np.sum(UY_Tot()[0:i+1])
    return UY_A
#Along Substrate/Cantilever Rotational Displacement
def RXY_Along():
    R_A=np.zeros(Nn)
    for i in range(1,Nn):
        R_A[i]=np.sum(RXY()[0:i+1])
    return R_A
#End Y Displacement
def UY_End(): 
    UY_End=np.sum(UY_Tot())
    return UY_End



"""
--------------------------------------THERMAL OUTPUT CALCULATION FUNCTIONS--------------------------------------
"""

#Y Displacement
def Th_UY():
    UY=np.zeros(Nn)                                                         #This method is to get around the divide by zero error and run-time issue
    UY[1:Nn+1]=((3*neff*n_e_ws[1:Nn+1]*n_e_wc[1:Nn+1]*ts*tc*tT*(ThEs-ThEc)*Ta[1:Nn+1])/((n_e_ws[1:Nn+1]**2*ts**4)+(neff**2*n_e_wc[1:Nn+1]**2*tc**4)+(neff*n_e_ws[1:Nn+1]*n_e_wc[1:Nn+1]*ts*tc*((6*ts*tc)+(4*ts**2)+(4*tc**2)))))*Le_xa[1:Nn+1]**2*np.cos(PAr[1:Nn+1]+Ar+AChr[1:Nn+1])
    return UY
#Rotational Displacement
def Th_RXY():
    RXY=np.zeros(Nn)                                                        #This method is to get around the divide by zero error and run-time issue
    RXY[1:Nn+1]=Th_UY()[1:Nn+1]*2./Le_xa[1:Nn+1]
    return RXY
#Rotational Y Displacement
def Th_R_UY():
    R=Th_RXY()
    R_UY=np.zeros(Nn)
    UY=np.tan(R)*Le_xa
    for i in range(1,Nn):  
        if VM == "Y" and Le_xa[i] != SM_Le_x:
            R_UY[i]=(np.sum(UY[0:VM_N])/SM_Le_x*RM_Le_x)+np.sum(UY[VM_N:i])
        else:
            R_UY[i]=np.sum(UY[0:i])
    return R_UY
#Total Y Displacement
def Th_UY_Tot():
    UY_T=Th_UY()+Th_R_UY()
    return UY_T
#Along Cantilever Y Displacement
def Th_UY_Along():
    UY_A=np.zeros(Nn)
    for i in range(1,Nn):
        UY_A[i]=np.sum(Th_UY_Tot()[0:i+1])
    return UY_A
#Along Cantilever Rotational Displacement
def Th_RXY_Along():
    R_A=np.zeros(Nn)
    for i in range(1,Nn):
        R_A[i]=np.sum(Th_RXY()[0:i+1])
    return R_A
#End Y Displacement
def Th_UY_End():
    UY_End=np.sum(Th_UY_Tot())
    return UY_End



"""
--------------------------------------ADDITIONAL FUNCTIONS--------------------------------------
"""



#Ytitles={"ks_e()":"Elemental Spring Constant (N/m)","ks_Along()":"Along Cantilever Spring Constant (N/m)","Y1k":"Spring Constant (N/m)","Y2u":"Y Displacement (um)"}
#Xtitles={SI:SI,"n_x":"Cantilever X-Length (m)"}
#
#"""2 Plot Window"""
#def WP2(X,Y1,Y2,c):
#    Grid=plt.figure(figsize=(11.7, 4.15))
#    Spec=gridspec.GridSpec(ncols=2, nrows=1, wspace=0.4, hspace=0.4)
#    L=Grid.add_subplot(Spec[0,0]).plot(X,Y1,color=c)
#    L=plt.xlim(0)
#    L=plt.title(Ytitles[Y1]+" vs"+Xtitles[X])
#    L=plt.xlabel(Xtitles[X])
#    L=plt.ylabel(Ytitles[Y1])
#    R=Grid.add_subplot(Spec[0,0]).plot(X,Y2,color=c)
#    R=plt.xlim(0)
#    R=plt.title(Ytitles[Y2]+" vs"+Xtitles[X])
#    R=plt.xlabel(Xtitles[X])
#    R=ply.ylabel(Ytitles[Y1])
#    #Grid.text(0.008,0.9,"Spring Constant = {0:.6f} N/m".format(ks_Tot()),fontsize=13,bbox=dict(facecolor="orange",alpha=0.5))
#    L=plt.subplots_adjust(top=0.77)
#    return
#
#"""4 Plot Window"""
#def WP4():
#    return

"""
--------------------------------------OUTPUTS--------------------------------------
"""


"""INPUTS SAVE"""
if I_O == "Y":
    DateTime=t.strftime('%Y%m%d-%Hhr%Mm')
    InputStamp=("k"if Spring=="Y"else"")+("U"if Displacement=="Y"else"")+("Th"if Th_Displacement=="Y"else"")+("SA"if SensitivityAnalysis=="Y"else"")+"_L"+str(int(Ls_x*1e+6))+("-C" if Coating=="Y"else"")+("-OF" if OF=="Y"else"")+("-AC" if AC=="Y"else"")+("-ACCh" if ACCh=="Y"else"")+("-PAC" if PAC=="Y"else"")
    Foldername=DateTime+"__"+InputStamp
    os.makedirs('Results/'+Foldername)
    sh.copy2('FDM_Model__Inputs.py','Results/'+Foldername+'/Inputs.txt')


"""SENSITIVITY ANALYSIS OUTPUTS"""
if SensitivityAnalysis == "Y":
    X=np.zeros(Si+1)
    Y1k=np.zeros(Si+1)
    Y2u=np.zeros(Si+1)
    for i in range(0,Si+1):
        Change=(SE-SS)/Si
        globals()[SI]=SS+(Change*i)
        #Outputs
        X[i]=globals()[SI]
        Y1k[i]=ks_Tot()
        Y2u[i]=UY_End()     
    SAGrid1=plt.figure(figsize=(11.7, 4.15))
    Spec=gridspec.GridSpec(ncols=2, nrows=1, wspace=0.4, hspace=0.4)
    TL=SAGrid1.add_subplot(Spec[0,0]).plot(X,Y1k,color="orange")
    TL=plt.xlim(0)
    TL=plt.title("Cantilever Spring Constant vs "+SI)
    TL=plt.xlabel(SI)
    TL=plt.ylabel("Spring Constant (N/m)")
    TR=SAGrid1.add_subplot(Spec[0, 1]).plot(X,Y2u*1e+6,color="orange")
    TR=plt.xlim(0)
    TR=plt.title("End Y Displacement vs "+SI)
    TR=plt.xlabel(SI)
    TR=plt.ylabel("End Y Displacement (um)")
    TL=plt.subplots_adjust(top=0.77)
    if I_O == "Y":
        SAGrid1.savefig('Results/'+Foldername+'/Outputs-Sensitivity Analysis Plots.svg')
    
    
"""MECHANICAL OUTPUTS"""
"Spring Constant"
if Spring == "Y" and SensitivityAnalysis != "Y":
    print "\nSpring Constant = {0:.6f} N/m\n".format(ks_Tot())
    kGrid1=plt.figure(figsize=(11.7, 4.15))
    Spec=gridspec.GridSpec(ncols=2, nrows=1, wspace=0.4, hspace=0.4)
    TL=kGrid1.add_subplot(Spec[0,0]).plot(n_x*1e+6,ks_e(),color="orange")
    TL=plt.xlim(0)
    TL=plt.yscale('log')
    TL=plt.title("Element Spring Constant vs X-Length")
    TL=plt.xlabel("Cantilever X-Length (um)")
    TL=plt.ylabel("Spring Constant (N/m)")
    TR=kGrid1.add_subplot(Spec[0, 1]).plot(n_x*1e+6,ks_Along(),color="orange")
    TR=plt.xlim(0)
    TR=plt.yscale('log')
    TR=plt.title("Along Cantilever Spring Constant vs X-Length")
    TR=plt.xlabel("Cantilever X-Length (um)")
    TR=plt.ylabel("Spring Constant (N/m)")
    kGrid1.text(0.008,0.9,"Spring Constant = {0:.6f} N/m".format(ks_Tot()),fontsize=13,bbox=dict(facecolor="orange",alpha=0.5))
    TL=plt.subplots_adjust(top=0.77)
    if I_O == "Y":
        kGrid1.savefig('Results/'+Foldername+'/Outputs-Spring Constant Plots.svg')
    

"Displacements"
if Displacement == "Y" and SensitivityAnalysis != "Y":
    #End Displacement
    print "\nEnd Y Displacement = {0:.6f} um\n".format(UY_End()*1e+6)
    #Plots
    DGrid1=plt.figure(figsize=(11.7, 8.3))
    Spec=gridspec.GridSpec(ncols=2, nrows=2, wspace=0.4, hspace=0.4)
    TL=DGrid1.add_subplot(Spec[0,0]).plot(n_x*1e+6,UY_Tot()*1e+6)
    TL=plt.xlim(0)
    TL=plt.ylim(0)
    TL=plt.title("Element Total Y Displacement vs X-Length")
    TL=plt.xlabel("Cantilever X-Length (um)")
    TL=plt.ylabel("Y Displacement (um)")
    TR=DGrid1.add_subplot(Spec[0, 1]).plot(n_x*1e+6,RXY())
    TR=plt.xlim(0)
    TR=plt.ylim(0)
    TR=plt.title("Element Rotational Displacement vs X-Length")
    TR=plt.xlabel("Cantilever X-Length (um)")
    TR=plt.ylabel("Rotational Displacement (rad)")
    BL=DGrid1.add_subplot(Spec[1, 0]).plot(n_x*1e+6,UY_Along()*1e+6)
    BL=plt.xlim(0)
    BL=plt.ylim(0)
    BL=plt.title("Along Cantilever Y Displacement vs X-Length")
    BL=plt.xlabel("Cantilever X-Length (um)")
    BL=plt.ylabel("Y Displacement (um)")
    BR=DGrid1.add_subplot(Spec[1, 1]).plot(n_x*1e+6,RXY_Along())
    BR=plt.xlim(0)
    BR=plt.xlim(0)
    BR=plt.title("Along Cantilever Rotational Displacement vs X-Length")
    BR=plt.xlabel("Cantilever X-Length (um)")
    BR=plt.ylabel("Rotational Displacement (rad)")
    DGrid1.text(0.008,0.95,"End Y Displacement = {0:.6f} um".format(UY_End()*1e+6),fontsize=13,bbox=dict(facecolor="blue",alpha=0.5))
    if I_O == "Y":
        DGrid1.savefig('Results/'+Foldername+'/Outputs-Displacement Plots.svg')
  
"""THERMAL OUTPUTS"""
if Th_Displacement == "Y" and SensitivityAnalysis != "Y":   
    print "\nThermal End Y Displacement = {0:.6f} um\n".format(Th_UY_End()*1e+6)
    ThGrid1=plt.figure(figsize=(11.7, 8.3))
    Spec=gridspec.GridSpec(ncols=2, nrows=2, wspace=0.4, hspace=0.4)
    TL=ThGrid1.add_subplot(Spec[0,0]).plot(n_x*1e+6,Th_UY_Tot()*1e+6,color='red')
    TL=plt.xlim(0)
    TL=plt.ylim(None,0)
    TL=plt.title("Element Total Thermal Y Displacement vs X-Length")
    TL=plt.xlabel("Cantilever X-Length (um)")
    TL=plt.ylabel("Y Displacement (um)")
    TR=ThGrid1.add_subplot(Spec[0, 1]).plot(n_x*1e+6,Th_RXY(),color='red')
    TR=plt.xlim(0)
    TR=plt.ylim(None,0)
    TR=plt.title("Element Thermal Rotational Displacement vs X-Length")
    TR=plt.xlabel("Cantilever X-Length (um)")
    TR=plt.ylabel("Rotational Displacement (rad)")
    BL=ThGrid1.add_subplot(Spec[1, 0]).plot(n_x*1e+6,Th_UY_Along()*1e+6,color='red')
    BL=plt.xlim(0)
    BL=plt.ylim(None,0)
    BL=plt.title("Along Cantilever Thermal Y Displacement vs X-Length")
    BL=plt.xlabel("Cantilever X-Length (um)")
    BL=plt.ylabel("Y Displacement (um)")
    BR=ThGrid1.add_subplot(Spec[1, 1]).plot(n_x*1e+6,Th_RXY_Along(),color='red')
    BR=plt.xlim(0)
    BR=plt.ylim(None,0)
    BR=plt.title("Along Cantilever Thermal Rotational Displacement vs X-Length")
    BR=plt.xlabel("Cantilever X-Length (um)")
    BR=plt.ylabel("Rotational Displacement (rad)")
    ThGrid1.text(0.008,0.95,"End Y Displacement = {0:.6f} um".format(Th_UY_End()*1e+6),fontsize=13,bbox=dict(facecolor="yellow",alpha=0.5))
    if I_O == "Y":
        ThGrid1.savefig('Results/'+Foldername+'/Outputs-Thermal Plots.svg')


"""
--------------------------------------MISCELLANEOUS CODE--------------------------------------

##########
Node Element Tapered Widths (m)
if T == "Y":
    a=np.asarray([(TM*(n_e_xc[i-1]-TS))+ws for i in range(int(np.floor(TS/Le_x))+1,int(np.ceil(TE/Le_x))+2) if n_e_xc[i-1]>TS and n_e_xc[i-1]<TE])
    b=np.insert(a,0,n_e_ws[0:int(np.floor(TS/Le_x))+1])
    if Ne-int(np.ceil(TE/Le_x)) >0:
        c=np.ones(Ne-int(np.ceil(TE/Le_x)))*Tw
        n_e_ws=np.append(b,c)
    else:
        n_e_ws=b
print n_e_ws
##########


"""