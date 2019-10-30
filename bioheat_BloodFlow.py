#> \author Elias Ghadam Soltani
#> \brief This is an example program to solve a diffusion equation using OpenCMISS calls.
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s):
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. If you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. If you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>





#================================================================================================================================
#  Initialise OpenCMISS
#================================================================================================================================

# Import the libraries (OpenCMISS,python,numpy,scipy)
import numpy,sys,os,time,math,cmath,csv
from shutil import copyfile
from opencmiss.iron import iron

#coupled1D0D
from scipy.fftpack import fft,ifft
from scipy.sparse  import linalg
from scipy.linalg  import inv,eig
from scipy.special import jn

t=time.time()
# Diagnostics
#iron.DiagnosticsSetOn(iron.DiagnosticTypes.IN,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])
#iron.DiagnosticsSetOn(iron.DiagnosticTypes.ALL,[1,2,3,4,5],"Diagnostics",[""])
#iron.ErrorHandlingModeSet(iron.ErrorHandlingModes.TRAP_ERROR)
#iron.OutputSetOn("Testing")

# Set the OpenCMISS random seed so that we can test this example by using the
# same parallel decomposition.
numberOfRandomSeeds = iron.RandomSeedsSizeGet()
randomSeeds = [0]*numberOfRandomSeeds
randomSeeds[0] = 100
iron.RandomSeedsSet(randomSeeds)

# Get the computational nodes info
#computationEnvironment = iron.ComputationEnvironment()
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber    = iron.ComputationalNodeNumberGet()

print(' ')

#================================================================================================================================
# Modules
#================================================================================================================================
# Copying current example file into output directory

if computationalNodeNumber==0:
  if os.path.exists("./output"):
    input('\033[1;31m'+'Output directory already exists. RENAME IT first, if required.'+'\033[0m')
  else:
    os.makedirs("./output/Details")
  copyfile(os.path.basename(__file__), './output/Details/file.py')

def module_exists(module_name):
    try:
        __import__(module_name)
    except ImportError:
        return False
    else:
        return True

if computationalNodeNumber==0 and module_exists('DescriptionGenerator'):
  import DescriptionGenerator

#================================================================================================================================
#  Problem Control Panel
#================================================================================================================================
numberOfDimensionsEnergy     = 1  #(One-dimensional)
#
# derivIdx = 1
#
ProgressDiagnostics = False   # Set to diagnostics
#
# energySolve       = False # Set to solve energy equation
# tissueSolve       = False # Set to solve energy equation
# tissueEnergySolve = True  # Set to solve energy equation
meshOrigin = [0.0,0.0,0.0]


# =================================
# F L O W
#if (CoupledBioheatFlow or TestFlow):
# =================================
# Only one of these could be true.
TestFlow = True
Bioheat = False
CoupledBioheatFlow = False
if (CoupledBioheatFlow):
  Bioheat = False
  TestFlow = False
# ENERGY=1
# TISSUE=2
# FLOW=3

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  numberOfDimensionsFlow     = 1  #(One-dimensional)
  numberOfComponentsFlow     = 2  #(Flow & Area)
  numberOfInflowNodes    = 0
  numberOfInputNodes     = 0
  numberOfBifurcations   = 0
  numberOfTrifurcations  = 0
  numberOfTerminalNodes  = 0
  inputNodeNumber        = []
  bifurcationNodeNumber  = []
  trifurcationNodeNumber = []
  coupledNodeNumber      = []

  # Set the user number
  compIdx    = 1
  derivIdx   = 1
  versionIdx = 1

  # Set the flags
  Heart               = False   # Set to use coupled 0D Windkessel models (from CellML) at model inlet boundaries
  HeartInput          = True
  RCRBoundaries       = False   # Set to use coupled 0D Windkessel models (from CellML) at model outlet boundaries
  nonReflecting       = True   # Set to use non-reflecting outlet boundaries
  streeBoundaries     = False   # Set to use structured tree outlet boundaries
  coupledAdvection    = False   # Set to solve a coupled advection problem
  timestepStability   = False   # Set to do a basic check of the stability of the hyperbolic problem based on the timestep size
  initialiseFromFile  = False   # Set to initialise values
  ProgressDiagnostics = False   # Set to diagnostics
# =================================

#================================================================================================================================
#  Start Program
#================================================================================================================================

controlLoopNode                        = 0
CoordinateSystemUserNumber             = 1
regionUserNumberEnergy                 = 2
regionUserNumberTissue                 = 3
basisUserNumberEnergy                  = 4
BasisUserNumberTissue                  = 5
meshUserNumberEnergy                   = 6
meshUserNumberTissue                   = 7
decompositionUserNumberEnergy          = 8
decompositionUserNumberTissue          = 9
geometricFieldUserNumberEnergy         = 10
geometricFieldUserNumberTissue         = 11
equationsSetUserNumberEnergy           = 12
equationsSetUserNumberTissue           = 13
equationsSetFieldUserNumberEnergy      = 14
equationsSetFieldUserNumberTissue      = 15
dependentFieldUserNumberEnergy         = 16
dependentFieldUserNumberTissue         = 17
materialsFieldUserNumberEnergy         = 18
materialsFieldUserNumberTissue         = 19
sourceFieldUserNumberEnergy            = 20
sourceFieldUserNumberTissue            = 21
independentFieldUserNumberEnergy       = 22
independentFieldUserNumberTissue       = 23
problemUserNumber                      = 24
cellMLUserNumberTissue                 = 25
cellMLModelsFieldUserNumberTissue      = 26  
cellMLParametersFieldUserNumberTissue  = 27
CellMLIntermediateFieldUserNumberTissue = 28

if (CoupledBioheatFlow or TestFlow):
  # Set program variables
  # CoordinateSystemUserNumber                = 1
  BasisUserNumberSpace                      = 29
  BasisUserNumberTime                       = 30
  RegionUserNumber                          = 31
  RegionUserNumber2                         = 32
  MeshUserNumber                            = 33
  MeshUserNumber2                           = 34
  DecompositionUserNumber                   = 35
  DecompositionUserNumber2                  = 36
  GeometricFieldUserNumber                  = 37
  GeometricFieldUserNumber2                 = 38
  EquationsSetFieldUserNumberStree          = 39
  EquationsSetFieldUserNumberCharacteristic = 40
  EquationsSetFieldUserNumberNavierStokes   = 41
  EquationsSetFieldUserNumberAdvection      = 42
  DependentFieldUserNumber                  = 43
  DependentFieldUserNumber2                 = 44
  DependentFieldUserNumber3                 = 45
  MaterialsFieldUserNumber                  = 46
  MaterialsFieldUserNumber2                 = 47
  IndependentFieldUserNumber                = 48
  EquationsSetUserNumberStree               = 49
  EquationsSetUserNumberCharacteristic      = 50
  EquationsSetUserNumberNavierStokes        = 51
  EquationsSetUserNumberAdvection           = 52
  # ProblemUserNumber                         = 26
  CellMLUserNumber                          = 53
  CellMLModelsFieldUserNumber               = 54
  CellMLStateFieldUserNumber                = 55
  CellMLIntermediateFieldUserNumber         = 56
  CellMLParametersFieldUserNumber           = 57
  MaterialsFieldUserNumberCellML            = 58
  AnalyticFieldUserNumber                   = 59
  # Solver user numbers
  SolverDAEUserNumber                       = 1
  SolverStreeUserNumber                     = 1
  SolverCharacteristicUserNumber            = 2
  SolverNavierStokesUserNumber              = 3
  SolverAdvectionUserNumber                 = 4
  # Materials constants
  MaterialsFieldUserNumberMu     = 1
  MaterialsFieldUserNumberRho    = 2
  MaterialsFieldUserNumberAlpha  = 3
  MaterialsFieldUserNumberPext   = 4
  MaterialsFieldUserNumberLs     = 5
  MaterialsFieldUserNumberTs     = 6
  MaterialsFieldUserNumberMs     = 7
  MaterialsFieldUserNumberG0     = 8
  MaterialsFieldUserNumberFr     = 9
  MaterialsFieldUserNumberD      = 1
  # Materials variables
  MaterialsFieldUserNumberA0     = 1
  MaterialsFieldUserNumberE      = 2
  MaterialsFieldUserNumberH      = 3
  MaterialsFieldUserNumberkp     = 4
  MaterialsFieldUserNumberk1     = 5
  MaterialsFieldUserNumberk2     = 6
  MaterialsFieldUserNumberk3     = 7
  MaterialsFieldUserNumberb1     = 8

#================================================================================================================================
#  Mesh Reading
#================================================================================================================================

# LENGTH = 0.5 # m
# numberGlobalXElements = 192
numberOfNodesEnergy    = 26
numberOfElementsEnergy = 25

#ELEMENTS
elementNodesEnergy= (numberOfElementsEnergy)*[2*[0]]
#Brachial artery
elementNodesEnergy[0]=[1,2]
elementNodesEnergy[1]=[2,3]
elementNodesEnergy[2]=[3,4]
elementNodesEnergy[3]=[4,5]
elementNodesEnergy[4]=[5,6]

#Radial
elementNodesEnergy[5] =[6,7]
elementNodesEnergy[6] =[7,8]
elementNodesEnergy[7] =[8,9]
elementNodesEnergy[8] =[9,10]
elementNodesEnergy[9]=[10,11]
elementNodesEnergy[10]=[11,12]
elementNodesEnergy[11]=[12,13]
elementNodesEnergy[12]=[13,14]
elementNodesEnergy[13]=[14,15]
elementNodesEnergy[14]=[15,16]

#Ulnar
elementNodesEnergy[15]=[6,17]
elementNodesEnergy[16]=[17,18]
elementNodesEnergy[17]=[18,19]
elementNodesEnergy[18]=[19,20]
elementNodesEnergy[19]=[20,21]
elementNodesEnergy[20]=[21,22]
elementNodesEnergy[21]=[22,23]
elementNodesEnergy[22]=[23,24]
elementNodesEnergy[23]=[24,25]
elementNodesEnergy[24]=[25,26]

#NODES
xValues = numpy.zeros((numberOfNodesEnergy,1),dtype = numpy.float)
yValues = numpy.zeros((numberOfNodesEnergy,1),dtype = numpy.float)
zValues = numpy.zeros((numberOfNodesEnergy,1),dtype = numpy.float)
radius= (numberOfNodesEnergy)*[1*[0]]
ra=(numberOfElementsEnergy)*[1*[0]]
Tt= (numberOfNodesEnergy)*[1*[0]]
#brachial artery
xValues[0]=152.12 #mm
xValues[1]=164.04
xValues[2]=179
xValues[3]=187.37
xValues[4]=195.89
xValues[5]=215.5

yValues[0]=-80.904
yValues[1]=-76.816
yValues[2]=-76.099
yValues[3]=-82.713
yValues[4]=-86.659
yValues[5]=-85.043

zValues[0]=1243.4
zValues[1]=1214.5
zValues[2]=1171
zValues[3]=1123
zValues[4]=1073.7
zValues[5]=1022.5

radius[0]=3.3062 #mm
radius[1]=3.1344
radius[2]=2.8703
radius[3]=2.5923
radius[4]=2.2908
radius[5]=1.9026

Tt[0]=37.0 #C
#Tt[1]=36.77
Tt[1]=37.0
# Tt[2]=36.67
Tt[2]=37.0
#Tt[3]=36.57
Tt[3]=37.0
# Tt[4]=36.47
Tt[4]=37.0
# Tt[5]=36.37
Tt[5]=37.0

#Radial
# le   = xValues[1][0]
xValues[6] =227.3
xValues[7] =239.54
xValues[8] =250.22
xValues[9] =256.49
xValues[10]=256.77
xValues[11]=258.42
xValues[12]=262.93
xValues[13]=267.01
xValues[14]=270.44
xValues[15]=278.41

yValues[6] =-90.085
yValues[7] =-93.529
yValues[8] =-93.818
yValues[9] =-90.243
yValues[10]=-85.607
yValues[11]=-88.336
yValues[12]=-100.68
yValues[13]=-115.58
yValues[14]=-121.21
yValues[15]=-119.14

zValues[6] =1008.6
zValues[7] =992.83
zValues[8] =975.58
zValues[9] =954.23
zValues[10]=930.66
zValues[11]=905.41
zValues[12]=877.61
zValues[13]=852.92
zValues[14]=823.34
zValues[15]=792.14

radius[6] =1.8158 #mm
radius[7] =1.7652
radius[8] =1.7176
radius[9] =1.6593
radius[10]=1.6072
radius[11]=1.5359
radius[12]=1.4673
radius[13]=1.3942
radius[14]=1.305
radius[15]=0.93851

#Tt[6] =36.26
# Tt[6]=37.0
# Tt[7] =36.2
# Tt[8] =36.14
# Tt[9] =36.08
# Tt[10]=36.02
# Tt[11]=35.96
# Tt[12]=35.90
# Tt[13]=35.84
# Tt[14]=35.76
# Tt[15]=35.70

#Ulnar
xValues[16]=213.85
xValues[17]=212.98
xValues[18]=213.25
xValues[19]=214.23
xValues[20]=215.15
xValues[21]=221.02
xValues[22]=228.77
xValues[23]=235.44
xValues[24]=239.83
xValues[25]=242.49

yValues[16]=-80.244
yValues[17]=-78.828
yValues[18]=-85.493
yValues[19]=-96.217
yValues[20]=-101.56
yValues[21]=-105.22
yValues[22]=-111.98
yValues[23]=-121.22
yValues[24]=-128.52
yValues[25]=-128.92

zValues[16]=1009.5
zValues[17]=999.14
zValues[18]=982.26
zValues[19]=962.45
zValues[20]=928.98
zValues[21]=893.59
zValues[22]=868.46
zValues[23]=844.1
zValues[24]=820.37
zValues[25]=795.86

radius[16]=1.8299
radius[17]=1.797
radius[18]=1.7584
radius[19]=1.6931
radius[20]=1.6111
radius[21]=1.5121
radius[22]=1.4446
radius[23]=1.3709
radius[24]=1.3085
radius[25]=1.1698

# Tt[16]=36.26
# Tt[17]=36.20
# Tt[18]=36.14
# Tt[19]=36.08
# Tt[20]=36.02
# Tt[21]=35.96
# Tt[22]=35.90
# Tt[23]=35.84
# Tt[24]=35.76
# Tt[25]=35.70

for nodeIdx in range(numberOfNodesEnergy):
    Tt[nodeIdx]=37 #C

arteriesElements=range(1,numberOfElementsEnergy+1)

for elemIdx in range(numberOfElementsEnergy):
    ra[elemIdx]=(radius[elemIdx]+radius[elemIdx+1])/2

ra[15]=(radius[5]+radius[16])/2 #ulnar nodes in bifurcation is different

# print(ra)
def elementLength(p1,p2):
    length=math.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2)
    return length

le= (numberOfElementsEnergy)*[1*[0]]
for i in range(numberOfElementsEnergy):
    p1=[xValues[i],yValues[i],zValues[i]]
    p2=[xValues[i+1],yValues[i+1],zValues[i+1]]
    le[i]=elementLength(p1,p2) #mm

brachialElements=[1,2,3,4,5]
radialElements  =[6,7,8,9,10,11,12,13,14,15]
ulnarElements   =[16,17,18,19,20,21,22,23,24,25]
xValuesEnergy=xValues
yValuesEnergy=yValues
zValuesEnergy=zValues
# # ================ VTK to csv file =================
# def vtkTocsv():
#   vtkFile='body.1.vtk'
#   nodesFile='nodes.csv'
#   elementsFile='elements.csv'
#   nodesFound=False
#   with open(vtkFile,'r') as vtk, open(nodesFile,'w') as nodes, open(elementsFile,'w') as elements:
#     for line in vtk:
#       line=line.strip().split()
#       if line[0].lower()=='points':
#         numberOfNodesTissue=int(line[1])
#         nodesFound=True
#       if nodesFound:
        

# ==================================================

# Tissue mesh
# vtkTocsv()
# FileName_ele = "MaxVol1000/WholeLeftArm.1.ele"
numberOfNodesTissue    = 58557
numberOfElementsTissue = 253489
numberOfLocalNodes = 4
# offset = 1 # offset is 1 if nodes and elements number begin from 0, offset for default = 0

# Printing the head and tail to see if we are reading the information correctly.
printHead = True
printTail = True
head = range(5)
tail = range(numberOfElementsTissue-4,numberOfElementsTissue+1)

# Label for parts
muscleRegionLabel      = 2
leftRadiusRegionLabel  = 3
leftHumerusRegionLabel = 4
leftUlnaRegionLabel    = 5

kidneyLeft             = -20
kidneyRight            = -30

muscleElements      = []
leftRadiusElements  = []
leftHumerusElements = []
leftUlnaElements    = []

kidneyElementsLeft  = []
kidneyElementsRight = []

localNodes=[0]*numberOfLocalNodes

elementNumber=0
totalNumberOfElements=1



# elementNodes = False
#with open("upperlimb_refined2.exelem", "r") as f:

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  if (ProgressDiagnostics):
      print( " == >> Reading geometry from files... << == ")

  # Read the node file
  with open('input/Node.csv','r') as csvfile:
      reader = csv.reader(csvfile, delimiter=',')
  #    workSpace(vars())
      rownum = 0
      for row in reader:
          if (rownum == 0):
              # Read the header row
              header = row
          else:
              # Read the number of nodes
              if (rownum == 1):
                  numberOfNodesSpace = int(row[10])
                  totalNumberOfNodes = numberOfNodesSpace*3
                  xValues = numpy.zeros((numberOfNodesSpace+1,4),dtype = numpy.float)
                  yValues = numpy.zeros((numberOfNodesSpace+1,4),dtype = numpy.float)
                  zValues = numpy.zeros((numberOfNodesSpace+1,4),dtype = numpy.float)
                  A0      = numpy.zeros((numberOfNodesSpace+1,4),dtype = numpy.float)  # Area        (m2)
                  H       = numpy.zeros((numberOfNodesSpace+1,4),dtype = numpy.float)  # Thickness   (m)
                  E       = numpy.zeros((numberOfNodesSpace+1,4),dtype = numpy.float)  # Elasticity  (Pa)
                  kp      = numpy.zeros((numberOfNodesSpace+1,4),dtype = numpy.float)
                  k1      = numpy.zeros((numberOfNodesSpace+1,4),dtype = numpy.float)
                  k2      = numpy.zeros((numberOfNodesSpace+1,4),dtype = numpy.float)
                  k3      = numpy.zeros((numberOfNodesSpace+1,4),dtype = numpy.float)
                  b1      = numpy.zeros((numberOfNodesSpace+1,4),dtype = numpy.float)
                  Tp      = [0]*(numberOfNodesSpace+1)                                 # Type (Artery or Vein)
                  segment = [0]*(numberOfNodesSpace+1)
                  Sigma = 0.5                                                          # Poisson Ratio
              # Initialise the coordinates
              segment[rownum] = row[0]
              xValues[rownum][0] = float(row[1])
              yValues[rownum][0] = float(row[2])
              zValues[rownum][0] = float(row[3])
              A0[rownum][0] = float(row[4])
              E [rownum][0] = float(row[5])
              H [rownum][0] = float(row[6])
              Tp[rownum] = row[8]
              # Read the input nodes
              if (row[9] == 'input'):
                  numberOfInputNodes = numberOfInputNodes+1
                  inputNodeNumber.append(rownum)
              # Read the bifurcation nodes
              elif (row[9] == 'bifurcation'):
                  numberOfBifurcations+=1
                  bifurcationNodeNumber.append(rownum)
                  xValues[rownum][1] = float(row[1])
                  yValues[rownum][1] = float(row[2])
                  zValues[rownum][1] = float(row[3])
                  xValues[rownum][2] = float(row[1])
                  yValues[rownum][2] = float(row[2])
                  zValues[rownum][2] = float(row[3])
              # Read the trifurcation nodes
              elif (row[9] == 'trifurcation'):
                  numberOfTrifurcations+=1
                  trifurcationNodeNumber.append(rownum)
                  xValues[rownum][1] = float(row[1])
                  yValues[rownum][1] = float(row[2])
                  zValues[rownum][1] = float(row[3])
                  xValues[rownum][2] = float(row[1])
                  yValues[rownum][2] = float(row[2])
                  zValues[rownum][2] = float(row[3])
                  xValues[rownum][3] = float(row[1])
                  yValues[rownum][3] = float(row[2])
                  zValues[rownum][3] = float(row[3])
              # Read the terminal nodes
              elif (row[9] == 'terminal'):
                  numberOfTerminalNodes+=1
                  coupledNodeNumber.append(rownum)
              if (Tp[rownum] == 'Artery' or Tp[rownum] == 'Vein'):
                  #Arterial system (Sherwin 2003)
                  kp[rownum][0] = (math.pi**0.5)/(1.0-Sigma**2)
                  k1[rownum][0] = -0.5
                  k2[rownum][0] = 1.0
                  k3[rownum][0] = 1.0
                  b1[rownum][0] = 0.5
              #elif (Tp[rownum] == 'Vein'):
                  #Venous system (Keijsers 2015, Pedley 1996)
              #    kp[rownum][0] = (math.pi**1.5)/(12.0*(1.0-Sigma**2)**0.5)
              #    k1[rownum][0] = -1.5
              #    k2[rownum][0] = 1.0
              #    k3[rownum][0] = 3.0
              #    b1[rownum][0] = 10.0
          # Next line
          rownum+=1

  #------------------

  # Read the element file
  with open('input/Element.csv','r') as csvfile:
      reader = csv.reader(csvfile, delimiter=',')
      rownum = 0
      i = 0
      k = 0
      for row in reader:
          if (rownum == 0):
              # Read the header row
              header = row
          else:
              # Read the number of elements
              if (rownum == 1):
                  totalNumberOfElements = int(row[11])
                  elementNodes          = (totalNumberOfElements+1)*[3*[0]]
                  bifurcationElements   = (numberOfBifurcations+1)*[3*[0]]
                  trifurcationElements  = (numberOfTrifurcations+1)*[4*[0]]
              # Read the element nodes
              elementNodes[rownum] = [int(row[1]),int(row[2]),int(row[3])]
              # Read the bifurcation elements
              if (row[4]):
                  i+=1
                  bifurcationElements[i] = [int(row[4]),int(row[5]),int(row[6])]
              # Read the trifurcation elements
              elif (row[7]):
                  k+=1
                  trifurcationElements[k] = [int(row[7]),int(row[8]),int(row[9]),int(row[10])]
          # Next line
          rownum+=1

  if (ProgressDiagnostics):
      print( " Input at nodes: " + str(inputNodeNumber))
      print( " Bifurcations at nodes: " + str(bifurcationNodeNumber))
      print( " Trifurcations at nodes: " + str(trifurcationNodeNumber))
      print( " Terminal at nodes: " + str(coupledNodeNumber))
      print( " == >> Finished reading geometry... << == ")
# =================================

#================================================================================================================================
#  Initial Data & Default values
#================================================================================================================================

# ===============================
# List of what I have changed in the code that might change in the future.
# I do not want to look for these parameters. The parameters that is possible to change later on in the code.
# For example I might need to change the blood density. So I do not want to look deep in the source code. So it should not be
# HARD CODED. 
# Therefore, the structure should be like this. we create the code layer by layer. the deepest layer is the layer that we most likely
# do not change them. and the first layer is the user layer that most likey will change. like the stop time. Other things would probably
# not change by the mesh or geometry.
# Wrtie down the list of the parameters that should be apart this thermoregulation code. So we read a file to control our problem.
# So as the time goes this source code will work for a general problem and later on I will only have a file to change the parameters.
# I will have a main file to read the setup of the thermoregulation parameters and I will not change this source code anymore. Only the
# control panel/setting file would change that makes it very easy for us. Also we can write some codes to show us the parameters we would 
# most probably check to test if the code is working or not.
# The best way is to create a software dedicated for thermoregulation. For example a software that uses OpenCMISS and example file as a
# library and then you can be focused on the results. For example you easily say compare the results for effect of kidney. so you can 
# easily toggle kidney and see the results.
# ===============================

# Artery =========
k_bl               = 0.5*1e-3      # W/mm.K blood conductivity.
rho_bl             = 1069.0*1e-9   # kg/mm3 blood density
c_bl               = 3650.0        # J/Kg.K blood specific heat
Alpha              = k_bl/(rho_bl*c_bl)       # mm2/s Diffusivity
# r                  = 1.5           # mm, inner radius of the artery
# CArtery            = 4*Alpha/(r*r) # 0.27911 1/s
# Tt                 = 37.0          # C

# Tissue ==========

rho_ms             = 1085.0*1e-9   # kg/mm3   muscle density
c_ms               = 3768.0        # J/Kg.K   muscle specific heat
rho_bn             = 1357.0*1e-9   # kg/mm3    bone density
c_bn               = 1700.0        # J/Kg.K   bone specific heat
rho_sk             = 1085.0*1e-9   # kg/mm3    skin density
c_sk               = 3680.0        # J/Kg.K   skin specific heat

k_ms               = 0.42*1e-3     # W/mm.K muscle conductivity.
k_bn               = 0.75*1e-3     # W/mm.K bone conductivity.
k_sk               = 0.47*1e-3     # W/mm.K skin conductivity.

h_conv             = 2.0*1e-6      # W/mm2.K
#h_conv            = 200.0*1e-6    # W/mm2.K for water
hr_rad             = 5.9*1e-6      # W/mm2.K See example 3.12 Incropera

# R_arm              = 0.03          # m

Tb                 = 37.0          # C blood temeprature
Tair               = 17.7          # C

w                  = 5e-4          # 1/s terminal blood flow per volume of tissue.

cMuscle            = rho_bl*c_bl/(rho_ms*c_ms) *w   # 4.51128e-4 1/s

cBone              = 0.0           #\see equations section

# Set the time parameters
#timeIncrement   = 0.5
timeIncrementBioheat   = 10
startTimeBioheat       = 0.0
stopTimeBioheat        = 8200

# Set the output parameters
DYNAMIC_SOLVER_DIFFUSION_OUTPUT_FREQUENCY = 100

# Set the solver parameters
#relativeTolerance   = 1.0E-05  # default: 1.0E-05
#absoluteTolerance   = 1.0E-08  # default: 1.0E-10
#DIVERGENCE_TOLERANCE = 1.0e+10  # default: 1.0e+05
MAXIMUM_ITERATIONS   = 1000   # default: 100000
#RESTART_VALUE        = 3000     # default: 30

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  # Set the material parameters
  Rho  = 1050.0                                 # Rho         (kg/m3)
  Mu   = 0.004                                  # Mu          (Pa.s)
  Fr   = 8.0*math.pi*Mu/Rho                     # Viscous resistance per unit length (m2/s)
  D    = 2.5                                    # Diffusivity (m2/s)
  G0   = 0.0                                    # Gravitational acceleration (m/s2)
  Pext = 0.0                                    # External pressure (Pa)
  Pv   = 2667.0                                 # Venous pressure = 20.0 mmHg (Pa)
  dt   = [0]*(numberOfNodesSpace+1)             # TimeStep    (s)
  eig  = [0]*(numberOfNodesSpace+1)             # Eigenvalues

  # Material parameter scaling factors
  Ls = 1000.0              # Length   (m -> mm)
  Ts = 1000.0              # Time     (s -> ms)
  Ms = 1000.0              # Mass     (kg -> g)

  Alpha = 1.3              # Flow profile
  Qs    = (Ls**3.0)/Ts     # Flow             (m3/s)
  As    = Ls**2.0          # Area             (m2)
  Hs    = Ls               # Vessel thickness (m)
  Es    = Ms/(Ls*Ts**2.0)  # Elasticity Pa    (kg/(ms2) --> g/(mm.ms^2)
  Rhos  = Ms/(Ls**3.0)     # Density          (kg/m3)
  Mus   = Ms/(Ls*Ts)       # Viscosity        (kg/(ms))
  Ps    = Ms/(Ls*Ts**2.0)  # Pressure         (kg/(ms2))
  Fs    = Mus/Rhos         # Viscous R        (m2/s)
  Ds    = (Ls**2.0)/Ts     # Diffusivity      (m2/s)
  Zs    = Ps/Qs            # Impedance        (pa/(m3/s))
  Gs    = Ls/(Ts**2.0)     # Acceleration     (m/s2)

  Rho = Rho*Rhos
  Mu  = Mu*Mus
  P   = Pext*Ps
  A0  = A0*As
  E   = E*Es
  H   = H*Hs
  Fr  = Fr*Fs
  D   = D*Ds
  G0  = G0*Gs

  Q  = numpy.zeros((numberOfNodesSpace+1,4))
  A  = numpy.zeros((numberOfNodesSpace+1,4))
  dQ = numpy.zeros((numberOfNodesSpace+1,4))
  dA = numpy.zeros((numberOfNodesSpace+1,4))

  # Set A0 for branch nodes
  for bifIdx in range(1,numberOfBifurcations+1):
      nodeIdx = bifurcationNodeNumber[bifIdx-1]
      for versionIdx in range(1,3):
          if (Tp[nodeIdx] == 'Artery'):
              A0[nodeIdx][versionIdx] = 2*A0[elementNodes[bifurcationElements[bifIdx][versionIdx]][1]][0]-A0[elementNodes[bifurcationElements[bifIdx][versionIdx]][2]][0]
              if (A0[nodeIdx][versionIdx] <= 0):
                  A0[nodeIdx][versionIdx] = A0[elementNodes[bifurcationElements[bifIdx][versionIdx]][1]][0]
          elif (Tp[nodeIdx] == 'Vein'):
              A0[nodeIdx][versionIdx] = A0[elementNodes[bifurcationElements[bifIdx][versionIdx]][1]][0]
          E [nodeIdx][versionIdx] = E [elementNodes[bifurcationElements[bifIdx][versionIdx]][1]][0]
          H [nodeIdx][versionIdx] = H [elementNodes[bifurcationElements[bifIdx][versionIdx]][1]][0]
          kp[nodeIdx][versionIdx] = kp[elementNodes[bifurcationElements[bifIdx][versionIdx]][1]][0]
          k1[nodeIdx][versionIdx] = k1[elementNodes[bifurcationElements[bifIdx][versionIdx]][1]][0]
          k2[nodeIdx][versionIdx] = k2[elementNodes[bifurcationElements[bifIdx][versionIdx]][1]][0]
          k3[nodeIdx][versionIdx] = k3[elementNodes[bifurcationElements[bifIdx][versionIdx]][1]][0]
          b1[nodeIdx][versionIdx] = b1[elementNodes[bifurcationElements[bifIdx][versionIdx]][1]][0]

  # Area for circle of willis
  #A0[422][1]=2
  #A0[422][2]=2
  #A0[425][1]=2
  #A0[428][1]=2
  #A0[425][2]=0.8
  #A0[428][2]=0.8

  for trifIdx in range(1,numberOfTrifurcations+1):
      nodeIdx = trifurcationNodeNumber[trifIdx-1]
      for versionIdx in range(1,4):
          if (Tp[nodeIdx] == 'Artery'):
              A0[nodeIdx][versionIdx] = 2*A0[elementNodes[trifurcationElements[trifIdx][versionIdx]][1]][0]-A0[elementNodes[trifurcationElements[trifIdx][versionIdx]][2]][0]
              if (A0[nodeIdx][versionIdx] <= 0):
                  A0[nodeIdx][versionIdx] = A0[elementNodes[trifurcationElements[trifIdx][versionIdx]][1]][0]
          elif (Tp[nodeIdx] == 'Vein'):
              A0[nodeIdx][versionIdx] = A0[elementNodes[trifurcationElements[trifIdx][versionIdx]][1]][0]
          E [nodeIdx][versionIdx] = E [elementNodes[trifurcationElements[trifIdx][versionIdx]][1]][0]
          H [nodeIdx][versionIdx] = H [elementNodes[trifurcationElements[trifIdx][versionIdx]][1]][0]
          kp[nodeIdx][versionIdx] = kp[elementNodes[trifurcationElements[trifIdx][versionIdx]][1]][0]
          k1[nodeIdx][versionIdx] = k1[elementNodes[trifurcationElements[trifIdx][versionIdx]][1]][0]
          k2[nodeIdx][versionIdx] = k2[elementNodes[trifurcationElements[trifIdx][versionIdx]][1]][0]
          k3[nodeIdx][versionIdx] = k3[elementNodes[trifurcationElements[trifIdx][versionIdx]][1]][0]
          b1[nodeIdx][versionIdx] = b1[elementNodes[trifurcationElements[trifIdx][versionIdx]][1]][0]

  # Start with Q=0, A=A0 state
  A = A0

  # Or initialise from init file
  if (initialiseFromFile):
      init = numpy.zeros([numberOfNodesSpace+1,4,4])
      init = numpy.load('./input/init.npy')
      Q [1:numberOfNodesSpace+1,:] = init[:,0,:]
      A [1:numberOfNodesSpace+1,:] = init[:,1,:]
      dQ[1:numberOfNodesSpace+1,:] = init[:,2,:]
      dA[1:numberOfNodesSpace+1,:] = init[:,3,:]

  # Set the boundary conditions flag
  InletBoundaryConditionType = iron.BoundaryConditionsTypes.FIXED_INLET
  if (nonReflecting):
      OutletBoundaryConditionType = iron.BoundaryConditionsTypes.FIXED_NONREFLECTING
  elif (RCRBoundaries):
      OutletBoundaryConditionType = iron.BoundaryConditionsTypes.FIXED_CELLML
  elif (streeBoundaries):
      OutletBoundaryConditionType = iron.BoundaryConditionsTypes.FIXED_STREE
  else:
      OutletBoundaryConditionType = iron.BoundaryConditionsTypes.FIXED_OUTLET

  # Set the output parameters
  # (NONE/PROGRESS/TIMING/SOLVER/MATRIX)
  DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE    = iron.SolverOutputTypes.NONE
  NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE  = iron.SolverOutputTypes.NONE
  NONLINEAR_SOLVER_CHARACTERISTIC_OUTPUT_TYPE = iron.SolverOutputTypes.NONE
  LINEAR_SOLVER_CHARACTERISTIC_OUTPUT_TYPE    = iron.SolverOutputTypes.NONE
  LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE     = iron.SolverOutputTypes.NONE
  # (NONE/TIMING/SOLVER/MATRIX)
  CMISS_SOLVER_OUTPUT_TYPE = iron.SolverOutputTypes.NONE
  DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY = 10

  # Set the time parameters
  numberOfPeriods = 5.0
  timePeriod      = 800.0
  timeIncrement   = 0.1
  startTime       = 0.0
  stopTime  = numberOfPeriods*timePeriod
  dynamicSolverNavierStokesTheta = [1.0]
  dynamicSolverAdvectionTheta    = [0.5]

  # Set the solver parameters
  relativeToleranceNonlinearNavierStokes   = 1.0E-05  # default: 1.0E-05
  absoluteToleranceNonlinearNavierStokes   = 1.0E-08  # default: 1.0E-10
  solutionToleranceNonlinearNavierStokes   = 1.0E-05  # default: 1.0E-05
  relativeToleranceLinearNavierStokes      = 1.0E-05  # default: 1.0E-05
  absoluteToleranceLinearNavierStokes      = 1.0E-08  # default: 1.0E-10
  relativeToleranceNonlinearCharacteristic = 1.0E-05  # default: 1.0E-05
  absoluteToleranceNonlinearCharacteristic = 1.0E-08  # default: 1.0E-10
  solutionToleranceNonlinearCharacteristic = 1.0E-05  # default: 1.0E-05
  relativeToleranceLinearCharacteristic    = 1.0E-05  # default: 1.0E-05
  absoluteToleranceLinearCharacteristic    = 1.0E-08  # default: 1.0E-10

  DIVERGENCE_TOLERANCE = 1.0e+10  # default: 1.0e+05
  MAXIMUM_ITERATIONS   = 100000   # default: 100000
  RESTART_VALUE        = 3000     # default: 30

  # N-S/C coupling tolerance
  couplingTolerance1D = 1.0E+6
  # 1D-0D coupling tolerance
  couplingTolerance1D0D = 0.001
  
  # Check the CellML flag
  if (RCRBoundaries or Heart):
      if (coupledAdvection):
          # Navier-Stokes solver
          EquationsSetSubtype = iron.EquationsSetSubtypes.COUPLED1D0D_ADV_NAVIER_STOKES
          # Characteristic solver
          EquationsSetCharacteristicSubtype = iron.EquationsSetSubtypes.CHARACTERISTIC
          # Advection solver
          EquationsSetAdvectionSubtype = iron.EquationsSetSubtypes.ADVECTION
          if (TestFlow):
            ProblemSubtype = iron.ProblemSubtypes.COUPLED1D0D_ADV_NAVIER_STOKES
      else:
          # Navier-Stokes solver
          EquationsSetSubtype = iron.EquationsSetSubtypes.COUPLED1D0D_NAVIER_STOKES
          # Characteristic solver
          EquationsSetCharacteristicSubtype = iron.EquationsSetSubtypes.CHARACTERISTIC
          if (TestFlow):
            ProblemSubtype = iron.ProblemSubtypes.COUPLED1D0D_NAVIER_STOKES
  elif (streeBoundaries):
      if (coupledAdvection):
          # Navier-Stokes solver
          EquationsSetSubtype = iron.EquationsSetSubtypes.COUPLED1D0D_ADV_NAVIER_STOKES
          # Characteristic solver
          EquationsSetCharacteristicSubtype = iron.EquationsSetSubtypes.CHARACTERISTIC
          # Stree solver
          EquationsSetStreeSubtype = iron.EquationsSetSubtypes.STREE1D0D_ADV
          # Advection solver
          EquationsSetAdvectionSubtype = iron.EquationsSetSubtypes.ADVECTION
          if (TestFlow):
            ProblemSubtype = iron.ProblemSubtypes.STREE1D0DAdv_NAVIER_STOKES
      else:
          # Navier-Stokes solver
          EquationsSetSubtype = iron.EquationsSetSubtypes.COUPLED1D0D_NAVIER_STOKES
          # Characteristic solver
          EquationsSetCharacteristicSubtype = iron.EquationsSetSubtypes.CHARACTERISTIC
          # Stree solver
          EquationsSetStreeSubtype = iron.EquationsSetSubtypes.STREE1D0D
          if (TestFlow):
            ProblemSubtype = iron.ProblemSubtypes.STREE1D0D
  else:
      if (coupledAdvection):
          # Navier-Stokes solver
          EquationsSetSubtype = iron.EquationsSetSubtypes.OnedTransientAdv_NAVIER_STOKES
          # Characteristic solver
          EquationsSetCharacteristicSubtype = iron.EquationsSetSubtypes.CHARACTERISTIC
          # Advection solver
          EquationsSetAdvectionSubtype = iron.EquationsSetSubtypes.ADVECTION
          if (TestFlow):
            ProblemSubtype = iron.ProblemSubtypes.TRANSIENT1D_ADV_NAVIER_STOKES
      else:
          # Navier-Stokes solver
          EquationsSetSubtype = iron.EquationsSetSubtypes.TRANSIENT1D_NAVIER_STOKES
          # Characteristic solver
          EquationsSetCharacteristicSubtype = iron.EquationsSetSubtypes.CHARACTERISTIC
          if (TestFlow):
            ProblemSubtype = iron.ProblemSubtypes.TRANSIENT1D_NAVIER_STOKES
# =================================
if (Bioheat):
  # Navier-Stokes solver
  equationsSetEnergySubtype = iron.EquationsSetSubtypes.ADVECTION_DIFFUSION
  equationsSetTissueSubtype = iron.EquationsSetSubtypes.LINEAR_SOURCE_DIFFUSION
  ProblemSubtype      = iron.ProblemSubtypes.THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION

if (CoupledBioheatFlow):
  ProblemType    = iron.ProblemTypes.PROBLEM_NAVIER_STOKES_DIFFUSION_ADVECTION_DIFFUSION_TYPE
  ProblemSubtype = iron.ProblemSubtypes.PROBLEM_COUPLED_BIOHEAT_NAVIERSTOKES_DIFF_ADV_DIFF_SUBTYPE
#================================================================================================================================
#  Coordinate System
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> COORDINATE SYSTEM << == ")

# Start the creation of RC coordinate system
CoordinateSystem = iron.CoordinateSystem()
CoordinateSystem.CreateStart(CoordinateSystemUserNumber)
CoordinateSystem.DimensionSet(3)
CoordinateSystem.CreateFinish()

print('\033[1;32m'+'Coordinate System COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Region
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> COORDINATE SYSTEM << == ")

# Start the creation of Energy region
regionEnergy = iron.Region()
regionEnergy.CreateStart(regionUserNumberEnergy,iron.WorldRegion)
regionEnergy.LabelSet("Energy")
regionEnergy.CoordinateSystemSet(CoordinateSystem)
regionEnergy.CreateFinish()

# Start the creation of Tissue region
regionTissue = iron.Region()
regionTissue.CreateStart(regionUserNumberTissue,iron.WorldRegion)
regionTissue.LabelSet("Tissue")
regionTissue.CoordinateSystemSet(CoordinateSystem)
regionTissue.CreateFinish()

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  # Start the creation of SPACE region
  Region = iron.Region()
  Region.CreateStart(RegionUserNumber,iron.WorldRegion)
  Region.label = "ArterialSystem"
  Region.coordinateSystem = CoordinateSystem
  Region.CreateFinish()

  if (streeBoundaries):
      # Start the creation of TIME region
      RegionStree = iron.Region()
      RegionStree.CreateStart(RegionUserNumber2,iron.WorldRegion)
      RegionStree.label = "StructuredTree"
      RegionStree.coordinateSystem = CoordinateSystem
      RegionStree.CreateFinish()
# =================================

print('\033[1;32m'+'Region            COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Bases
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> BASIS << == ")

# Start the creation of bases for blood energy equaiton
basisXiGaussSpace = 3
basisEnergy = iron.Basis()
basisEnergy.CreateStart(basisUserNumberEnergy)
basisEnergy.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
basisEnergy.NumberOfXiSet(numberOfDimensionsEnergy)
basisEnergy.InterpolationXiSet([iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE])
basisEnergy.QuadratureNumberOfGaussXiSet([basisXiGaussSpace])
basisEnergy.CreateFinish()


# Create a tri-linear Simplex basis
basisTissue = iron.Basis()
basisTissue.CreateStart(BasisUserNumberTissue)
basisTissue.TypeSet(iron.BasisTypes.SIMPLEX)
basisTissue.NumberOfXiSet(3)
basisTissue.InterpolationXiSet([iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX]*3)
basisTissue.QuadratureOrderSet(2)
basisTissue.CreateFinish()

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  # Start the creation of SPACE bases
  basisXiGaussSpace = 3
  BasisSpace = iron.Basis()
  BasisSpace.CreateStart(BasisUserNumberSpace)
  BasisSpace.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
  BasisSpace.numberOfXi = numberOfDimensionsFlow
  BasisSpace.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]
  BasisSpace.quadratureNumberOfGaussXi = [basisXiGaussSpace]
  BasisSpace.CreateFinish()

  if (streeBoundaries):
      # Start the creation of TIME bases
      basisXiGaussSpace = 3
      BasisTime = iron.Basis()
      BasisTime.CreateStart(BasisUserNumberTime)
      BasisTime.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
      BasisTime.numberOfXi = numberOfDimensionsFlow
      BasisTime.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]
      BasisTime.quadratureNumberOfGaussXi = [basisXiGaussSpace]
      BasisTime.CreateFinish()
# =================================

print('\033[1;32m'+'Bases             COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Nodes
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> NODES << == ")
#
# # Start the creation of mesh nodes
NodesEnergy = iron.Nodes()
NodesEnergy.CreateStart(regionEnergy,numberOfNodesEnergy)
NodesEnergy.CreateFinish()

# Define nodes for the tissue mesh
nodesTissue = iron.Nodes()
nodesTissue.CreateStart(regionTissue,numberOfNodesTissue)
nodesTissue.CreateFinish()

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  # Start the creation of mesh nodes
  Nodes = iron.Nodes()
  Nodes.CreateStart(Region,totalNumberOfNodes)
  Nodes.CreateFinish()

  if (streeBoundaries):
      NodesStree = iron.Nodes()
      NodesStree.CreateStart(RegionStree,timePeriod+1)
      NodesStree.CreateFinish()
# =================================

print('\033[1;32m'+'Nodes             COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Mesh
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> MESH << == ")


# Start the creation of artery mesh for energy =======CellMLShiv
# Create a generated mesh
# generatedMesh = iron.GeneratedMesh()
# generatedMesh.CreateStart(meshUserNumberEnergy,regionEnergy)
# generatedMesh.TypeSet(iron.GeneratedMeshTypes.REGULAR)
# generatedMesh.BasisSet([basisEnergy])
# generatedMesh.ExtentSet([0,LENGTH,0])
# generatedMesh.NumberOfElementsSet([numberGlobalXElements])
#
# meshEnergy = iron.Mesh()
# generatedMesh.CreateFinish(meshUserNumberEnergy,meshEnergy)
meshEnergy = iron.Mesh()
meshEnergy.CreateStart(meshUserNumberEnergy,regionEnergy,3)
meshEnergy.NumberOfElementsSet(numberOfElementsEnergy)
meshEnergy.NumberOfComponentsSet(1)

# Start the creation of mesh elements for energy
meshElementsEnergy  = iron.MeshElements()
meshComponentNumber = 1
meshElementsEnergy.CreateStart(meshEnergy,meshComponentNumber,basisEnergy)

for elemIdx in range(0,numberOfElementsEnergy):
    meshElementsEnergy.NodesSet(elemIdx+1,elementNodesEnergy[elemIdx])
meshElementsEnergy.CreateFinish()

meshEnergy.CreateFinish()


# Start the creation of tissue mesh
meshTissue = iron.Mesh()
meshTissue.CreateStart(meshUserNumberTissue,regionTissue,3)
meshTissue.origin=meshOrigin
meshTissue.NumberOfComponentsSet(1)
meshTissue.NumberOfElementsSet(numberOfElementsTissue)

# Start the creation of mesh elements for tissue
meshElementsTissue  = iron.MeshElements()
meshComponentNumber = 1
meshElementsTissue.CreateStart(meshTissue, meshComponentNumber, basisTissue)


print( "Elapsed time before reading ele file is: ", time.time()-t)
# reading elements and its local Nodes and setting elements.Nodes
with open('elements2.csv') as elementscsv:
  reader = csv.reader(elementscsv)
  # next(elementscsv)
  elementNumber=0
  for row in reader:
      elementNumber=elementNumber+1
      
      for elemIdx in range(numberOfLocalNodes):
        localNodes[elemIdx]=int(row[elemIdx])+1
      materialType = int(row[numberOfLocalNodes])
      if materialType == kidneyLeft:
        kidneyElementsLeft.append(elementNumber)
      elif materialType == kidneyRight:
        kidneyElementsRight.append(elementNumber)
      meshElementsTissue.NodesSet(elementNumber,localNodes)
print(kidneyElementsLeft[1],kidneyElementsRight[1],'Hi',len(kidneyElementsLeft),len(kidneyElementsRight))
# elementscsv.close
# with open(FileName_ele, "r") as f:
#     target=f.readlines()
#     for lineNum,line in enumerate(target):
#         target[lineNum] = target[lineNum].rstrip('\n\r').replace("\t"," ").replace(","," ").split()
#     for lineNum,line in enumerate(target):
#         if (lineNum !=0 and lineNum<=numberOfElementsTissue):
# # reading elements and its localNodes
#             elementNumber = int(target[lineNum][0])+offset
#             localNodes[0] = int(target[lineNum][1])+offset
#             localNodes[1] = int(target[lineNum][2])+offset
#             localNodes[2] = int(target[lineNum][3])+offset
#             localNodes[3] = int(target[lineNum][4])+offset
#             materialType  = int(target[lineNum][5])
# # Giving elements material label
#             if materialType == muscleRegionLabel:
#                 muscleElements.append(elementNumber)
#             elif materialType == leftRadiusRegionLabel:
#                 leftRadiusElements.append(elementNumber)
#             elif materialType == leftHumerusRegionLabel:
#                 leftHumerusElements.append(elementNumber)
#             elif materialType == leftUlnaRegionLabel:
#                 leftUlnaElements.append(elementNumber)
# # printing the head and the tail of the elements
#             if elementNumber in head and printHead:
#               print(elementNumber,localNodes,materialType)
#             elif elementNumber in tail and printTail:
#               print(elementNumber,localNodes,materialType)
#
#             meshElementsTissue.NodesSet(elementNumber,localNodes)

# print( "number of muscle Elements = %d\nnumber of left Ulna elements = %d\ntotal number of elements = %d\n"%(len(muscleElements)
# ,len(leftUlnaElements),len(muscleElements)+len(leftRadiusElements)+len(leftHumerusElements)+len(leftUlnaElements)))

print( "Elapsed time after reading ele file is: ", time.time()-t)
#input("Press Enter to continue...")


meshElementsTissue.CreateFinish()
meshTissue.CreateFinish()

# Obtain boundary Nodes
# meshNodesTissue=iron.MeshNodes()
# meshTissue.NodesGet(1,meshNodesTissue)

# with open('boundary_nodes.csv','w') as csvFile:
#   writer=csv.writer(csvFile)
#   for nodeNumber in range(1,numberOfNodesTissue+1):
#     if meshNodesTissue.NodeOnBoundaryGet(nodeNumber) :
#       writer.writerow([nodeNumber])

# input("Stop Here")

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  if (ProgressDiagnostics):
      print( " == >> MESH << == ")

  # Start the creation of SPACE mesh
  Mesh = iron.Mesh()
  Mesh.CreateStart(MeshUserNumber,Region,numberOfDimensionsFlow)
  Mesh.NumberOfElementsSet(totalNumberOfElements)
  if (coupledAdvection):
      meshNumberOfComponents = 2
      Mesh.NumberOfComponentsSet(meshNumberOfComponents)
      # Specify the mesh components
      MeshElementsSpace = iron.MeshElements()
      MeshElementsConc  = iron.MeshElements()
      meshComponentNumberSpace = 1
      meshComponentNumberConc  = 2
  else:
      meshNumberOfComponents = 1
      # Specify the mesh components
      Mesh.NumberOfComponentsSet(meshNumberOfComponents)
      # Specify the mesh components
      MeshElementsSpace = iron.MeshElements()
      meshComponentNumberSpace = 1

  #------------------

  # Specify the SPACE mesh component
  MeshElementsSpace.CreateStart(Mesh,meshComponentNumberSpace,BasisSpace)
  for elemIdx in range(1,totalNumberOfElements+1):
      MeshElementsSpace.NodesSet(elemIdx,elementNodes[elemIdx])
  for bifIdx in range(1,numberOfBifurcations+1):
      nodeIdx = bifurcationNodeNumber[bifIdx-1]
      if (Tp[nodeIdx] == 'Artery'):
          MeshElementsSpace.LocalElementNodeVersionSet(int(bifurcationElements[bifIdx][0]),1,1,3)
          MeshElementsSpace.LocalElementNodeVersionSet(int(bifurcationElements[bifIdx][1]),2,1,1)
          MeshElementsSpace.LocalElementNodeVersionSet(int(bifurcationElements[bifIdx][2]),3,1,1)
      elif (Tp[nodeIdx] == 'Vein'):
          MeshElementsSpace.LocalElementNodeVersionSet(int(bifurcationElements[bifIdx][0]),1,1,1)
          MeshElementsSpace.LocalElementNodeVersionSet(int(bifurcationElements[bifIdx][1]),2,1,3)
          MeshElementsSpace.LocalElementNodeVersionSet(int(bifurcationElements[bifIdx][2]),3,1,3)
  for trifIdx in range(1,numberOfTrifurcations+1):
      nodeIdx = trifurcationNodeNumber[trifIdx-1]
      if (Tp[nodeIdx] == 'Artery'):
          MeshElementsSpace.LocalElementNodeVersionSet(int(trifurcationElements[trifIdx][0]),1,1,3)
          MeshElementsSpace.LocalElementNodeVersionSet(int(trifurcationElements[trifIdx][1]),2,1,1)
          MeshElementsSpace.LocalElementNodeVersionSet(int(trifurcationElements[trifIdx][2]),3,1,1)
          MeshElementsSpace.LocalElementNodeVersionSet(int(trifurcationElements[trifIdx][3]),4,1,1)
      elif (Tp[nodeIdx] == 'Vein'):
          MeshElementsSpace.LocalElementNodeVersionSet(int(trifurcationElements[trifIdx][0]),1,1,1)
          MeshElementsSpace.LocalElementNodeVersionSet(int(trifurcationElements[trifIdx][1]),2,1,3)
          MeshElementsSpace.LocalElementNodeVersionSet(int(trifurcationElements[trifIdx][2]),3,1,3)
          MeshElementsSpace.LocalElementNodeVersionSet(int(trifurcationElements[trifIdx][3]),4,1,3)
  MeshElementsSpace.CreateFinish()

  #------------------

  # Specify the CONCENTRATION mesh component
  if (coupledAdvection):
      MeshElementsConc.CreateStart(Mesh,meshComponentNumberConc,BasisSpace)
      for elemIdx in range(1,totalNumberOfElements+1):
          MeshElementsConc.NodesSet(elemIdx,elementNodes[elemIdx])
      MeshElementsConc.CreateFinish()

  # Finish the creation of the mesh
  Mesh.CreateFinish()

  #------------------

  if (streeBoundaries):
      # Start the creation of TIME mesh
      MeshTime = iron.Mesh()
      MeshTime.CreateStart(MeshUserNumber2,RegionStree,numberOfDimensionsFlow)
      MeshTime.NumberOfElementsSet(timePeriod)
      MeshTime.NumberOfComponentsSet(1)
      # Specify the mesh components
      MeshElementsTime = iron.MeshElements()
      meshComponentNumberTime = 1
      MeshElementsTime.CreateStart(MeshTime,meshComponentNumberTime,BasisTime)
      for elemIdx in range(1,timePeriod+1):
          MeshElementsTime.NodesSet(elemIdx,[elemIdx,elemIdx+1])
      MeshElementsTime.CreateFinish()
      MeshTime.CreateFinish()
# =================================

print('\033[1;32m'+'Mesh              COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Decomposition
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> MESH DECOMPOSITION << == ")

# Create a decomposition for artery energy mesh
decompositionEnergy = iron.Decomposition()
decompositionEnergy.CreateStart(decompositionUserNumberEnergy,meshEnergy)
decompositionEnergy.TypeSet(iron.DecompositionTypes.CALCULATED)
decompositionEnergy.NumberOfDomainsSet(numberOfComputationalNodes)
# decompositionEnergy.CalculateFacesSet(True)
decompositionEnergy.CreateFinish()

# Start the creation of Tissue mesh decomposition =======
decompositionTissue = iron.Decomposition()
decompositionTissue.CreateStart(decompositionUserNumberTissue,meshTissue)
decompositionTissue.TypeSet(iron.DecompositionTypes.CALCULATED)
decompositionTissue.NumberOfDomainsSet(numberOfComputationalNodes)
decompositionTissue.CalculateFacesSet(True)
decompositionTissue.CreateFinish()

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  if (ProgressDiagnostics):
    print( " == >> MESH DECOMPOSITION << == ")

  # Start the creation of SPACE mesh decomposition
  Decomposition = iron.Decomposition()
  Decomposition.CreateStart(DecompositionUserNumber,Mesh)
  Decomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
  Decomposition.NumberOfDomainsSet(numberOfComputationalNodes)
  Decomposition.CreateFinish()

  #------------------

  if (streeBoundaries):
      # Start the creation of TIME mesh decomposition
      DecompositionTime = iron.Decomposition()
      DecompositionTime.CreateStart(DecompositionUserNumber2,MeshTime)
      DecompositionTime.TypeSet(iron.DecompositionTypes.CALCULATED)
      DecompositionTime.NumberOfDomainsSet(numberOfComputationalNodes)
      DecompositionTime.CreateFinish()
# =================================

print('\033[1;32m'+'Decomposition     COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Geometric Field
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> GEOMETRIC FIELD << == ")

# Start the creation of geometric field for artery energy
geometricFieldEnergy = iron.Field()
geometricFieldEnergy.CreateStart(geometricFieldUserNumberEnergy,regionEnergy)
geometricFieldEnergy.LabelSet('Geometric Field')
geometricFieldEnergy.MeshDecompositionSet(decompositionEnergy)
geometricFieldEnergy.NumberOfVariablesSet(1)
geometricFieldEnergy.VariableLabelSet(iron.FieldVariableTypes.U,'Artery Coordinates')
geometricFieldEnergy.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricFieldEnergy.ScalingTypeSet(iron.FieldScalingTypes.NONE)

for componentNumber in range(1,CoordinateSystem.dimension+1):
    geometricFieldEnergy.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentNumber,1)
# geometricFieldEnergy.CreateFinish()
# generatedMesh.GeometricParametersCalculate(geometricFieldEnergy)

geometricFieldEnergy.CreateFinish()

# Set the geometric field values for version 1
versionIdx = 1
derivIdx = 1
for nodeIdx in range(0,numberOfNodesEnergy):
    nodeDomain = decompositionEnergy.NodeDomainGet(nodeIdx+1,1)
    if (nodeDomain == computationalNodeNumber):
        geometricFieldEnergy.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
         versionIdx,derivIdx,nodeIdx+1,1,xValuesEnergy[nodeIdx][0])
        geometricFieldEnergy.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
         versionIdx,derivIdx,nodeIdx+1,2,yValuesEnergy[nodeIdx][0])
        geometricFieldEnergy.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
         versionIdx,derivIdx,nodeIdx+1,3,zValuesEnergy[nodeIdx][0])

# Finish the parameter update
geometricFieldEnergy.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
geometricFieldEnergy.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)


# Start the creation of geometric field for Tissue
geometricFieldTissue = iron.Field()
geometricFieldTissue.CreateStart(geometricFieldUserNumberTissue,regionTissue)
geometricFieldTissue.LabelSet('Geometric Field')
geometricFieldTissue.MeshDecompositionSet(decompositionTissue)
geometricFieldTissue.NumberOfVariablesSet(1)
geometricFieldTissue.VariableLabelSet(iron.FieldVariableTypes.U,'Tissue Coordinates')
geometricFieldTissue.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricFieldTissue.ScalingTypeSet(iron.FieldScalingTypes.NONE)
for componentNumber in range(1,CoordinateSystem.dimension+1):
    geometricFieldTissue.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentNumber,1)
geometricFieldTissue.CreateFinish()



# Get nodes
nodes = iron.Nodes()
regionTissue.NodesGet(nodes)
numberOfNodes = nodes.numberOfNodes
print( numberOfNodes)
# Get or calculate geometric Parameters

print( "Elapsed time before reading node file is: ", time.time()-t)

FileName_node = "nodes.csv"

# X input file is mm. if you want to keep it that way you need to multiply k and rho*c by factors of 10^-3 and 10^-9 respectively.
Units = 1e0
kFactor= 1#1e-3
rhoCFactor = 1#1e-9

printHead = True
printTail = True
head = range(5)
tail = range(numberOfNodes-4,numberOfNodes+1)

#nodePositions = False
boundaryMarker=0
nodeNumber = 0

boundaryTissue = []
skinMarker = 2

with open(FileName_node) as nodescsv:
  reader = csv.reader(nodescsv)
  # next(elementscsv)
  nodeNumber=0
  for row in reader:
      # for nodeIdx in range(numberOfLocalNodes):
      nodeNumber=nodeNumber+1
      # if(nodeNumber not in [23646,25193]):
      x=float(row[0])
      y=float(row[1])
      z=float(row[2])
      nodeDomain = decompositionTissue.NodeDomainGet(nodeNumber,1)
      if nodeDomain == computationalNodeNumber:
        geometricFieldTissue.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,
        nodeNumber,1,x)
        geometricFieldTissue.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,
        nodeNumber,2,y)
        geometricFieldTissue.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,
        nodeNumber,3,z)



with open('boundary_nodes.csv') as csvFile:
  reader = csv.reader(csvFile)
  # next(elementscsv)
  for row in reader:
      # for nodeIdx in range(numberOfLocalNodes):
    boundaryTissue.append(int(row[0]))



# with open(FileName_node, "r") as f:
#     target=f.readlines()
#     for lineNum,line in enumerate(target):
#         target[lineNum] = target[lineNum].rstrip('\n\r').replace("\t"," ").replace(","," ").split()
#     for lineNum,line in enumerate(target):
#         if lineNum !=0:
#           nodeNumber = int(target[lineNum][0])+offset
#           x = float(target[lineNum][1]) * Units
#           y = float(target[lineNum][2]) * Units
#           z = float(target[lineNum][3]) * Units
#           boundaryMarker = int(target[lineNum][4])
#           if boundaryMarker == skinMarker:
#             boundaryTissue.append(nodeNumber)
#
#           if nodeNumber in head and printHead:
#             print( nodeNumber,x,y,z,boundaryMarker)
#           elif nodeNumber in tail and printTail:
#             print( nodeNumber,x,y,z,boundaryMarker)
#
#           nodeDomain = decompositionTissue.NodeDomainGet(nodeNumber,1)
#           if nodeDomain == computationalNodeNumber:
#             geometricFieldTissue.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,nodeNumber,1,x)
#             geometricFieldTissue.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,nodeNumber,2,y)
#             geometricFieldTissue.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,nodeNumber,3,z)

print( "number of boundary nodes = %d"%len(boundaryTissue))

print( "Elapsed time after reading node file is: ", time.time()-t)

# Update the geometric field
geometricFieldTissue.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
geometricFieldTissue.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  if (ProgressDiagnostics):
    print( " == >> GEOMETRIC FIELD << == ")

  # Start the creation of SPACE geometric field
  GeometricField = iron.Field()
  GeometricField.CreateStart(GeometricFieldUserNumber,Region)
  GeometricField.NumberOfVariablesSet(1)
  GeometricField.VariableLabelSet(iron.FieldVariableTypes.U,'Coordinates')
  GeometricField.TypeSet = iron.FieldTypes.GEOMETRIC
  GeometricField.meshDecomposition = Decomposition
  GeometricField.ScalingTypeSet = iron.FieldScalingTypes.NONE
  # Set the mesh component to be used by the geometric field components
  for componentNumber in range(1,CoordinateSystem.dimension+1):
      GeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentNumber,
      meshComponentNumberSpace)
  GeometricField.CreateFinish()

  # Set the geometric field values for version 1
  versionIdx = 1
  for nodeIdx in range(1,numberOfNodesSpace+1):
      nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberSpace)
      if (nodeDomain == computationalNodeNumber):
          GeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
          versionIdx,derivIdx,nodeIdx,1,xValues[nodeIdx][0])
          GeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
          versionIdx,derivIdx,nodeIdx,2,yValues[nodeIdx][0])
          GeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
          versionIdx,derivIdx,nodeIdx,3,zValues[nodeIdx][0])
  # Set the geometric field for bifurcation
  for bifIdx in range (1,numberOfBifurcations+1):
      nodeIdx = bifurcationNodeNumber[bifIdx-1]
      nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberSpace)
      if (nodeDomain == computationalNodeNumber):
          for versionNumber in range(2,4):
              GeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              versionNumber,derivIdx,nodeIdx,1,xValues[nodeIdx][versionNumber-1])
              GeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              versionNumber,derivIdx,nodeIdx,2,yValues[nodeIdx][versionNumber-1])
              GeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              versionNumber,derivIdx,nodeIdx,3,zValues[nodeIdx][versionNumber-1])
  # Set the geometric field for trifurcation
  for trifIdx in range (1,numberOfTrifurcations+1):
      nodeIdx = trifurcationNodeNumber[trifIdx-1]
      nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberSpace)
      if nodeDomain == computationalNodeNumber:
          for versionNumber in range(2,5):
              GeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              versionNumber,derivIdx,nodeIdx,1,xValues[nodeIdx][versionNumber-1])
              GeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              versionNumber,derivIdx,nodeIdx,2,yValues[nodeIdx][versionNumber-1])
              GeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              versionNumber,derivIdx,nodeIdx,3,zValues[nodeIdx][versionNumber-1])

  # Finish the parameter update
  GeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
  GeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

  #------------------

  if (streeBoundaries):
      # Start the creation of TIME geometric field
      GeometricFieldTime = iron.Field()
      GeometricFieldTime.CreateStart(GeometricFieldUserNumber2,RegionStree)
      GeometricFieldTime.NumberOfVariablesSet(1)
      GeometricFieldTime.VariableLabelSet(iron.FieldVariableTypes.U,'Time')
      GeometricFieldTime.TypeSet = iron.FieldTypes.GEOMETRIC
      GeometricFieldTime.meshDecomposition = DecompositionTime
      GeometricFieldTime.ScalingTypeSet = iron.FieldScalingTypes.NONE
      # Set the mesh component to be used by the geometric field components
      for componentNumber in range(1,CoordinateSystem.dimension+1):
          GeometricFieldTime.ComponentMeshComponentSet(iron.FieldVariableTypes.U,
          componentNumber,meshComponentNumberTime)
      GeometricFieldTime.CreateFinish()
# =================================

print('\033[1;32m'+'Geometric Field   COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Equations Sets
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> EQUATIONS SET << == ")

if (not TestFlow):
  # Create the equations set for advection diffusion in arteries
  # dT/dt+u dT/dx-alpha d2T/dx2-(b-cT)=0
  equationsSetEnergy = iron.EquationsSet()
  equationsSetFieldEnergy = iron.Field()
  # Set the equations set to be a dynamic linear problem
  equationsSetEnergySpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
    iron.EquationsSetTypes.ADVECTION_EQUATION,
    equationsSetEnergySubtype]
  equationsSetEnergy.CreateStart(equationsSetUserNumberEnergy,regionEnergy,geometricFieldEnergy,
  equationsSetEnergySpecification,equationsSetFieldUserNumberEnergy,equationsSetFieldEnergy)
  equationsSetEnergy.LabelSet('Advec Diffusion Equation')
  equationsSetEnergy.CreateFinish()


  # Create standard Diffusion equations set
  # dT/dt-div(alpha grad(T))-(b-cT)=0
  equationsSetTissue = iron.EquationsSet()
  equationsSetFieldTissue = iron.Field()
  equationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
    iron.EquationsSetTypes.DIFFUSION_EQUATION,
    equationsSetTissueSubtype]
  equationsSetTissue.CreateStart(equationsSetUserNumberTissue,regionTissue,geometricFieldTissue,
          equationsSetSpecification,equationsSetFieldUserNumberTissue,equationsSetFieldTissue)
  equationsSetTissue.LabelSet('Diffusion Equation')
  equationsSetTissue.CreateFinish()

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  if (ProgressDiagnostics):
    print( " == >> EQUATIONS SET << == ")

  # Create the equations set for STREE
  if (streeBoundaries):
      EquationsSetStree = iron.EquationsSet()
      EquationsSetFieldStree = iron.Field()
      StreeEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                iron.EquationsSetTypes.STREE_EQUATION,
                EquationsSetStreeSubtype]
      # Set the equations set to be a dynamic linear problem
      EquationsSetStree.CreateStart(EquationsSetUserNumberStree,RegionStree,GeometricFieldTime,
          StreeEquationsSetSpecification,EquationsSetFieldUserNumberStree,EquationsSetFieldStree)
      EquationsSetStree.CreateFinish()

  #------------------

  # Create the equations set for CHARACTERISTIC
  EquationsSetCharacteristic = iron.EquationsSet()
  EquationsSetFieldCharacteristic = iron.Field()
  # Set the equations set to be a static nonlinear problem
  CharacteristicEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
            iron.EquationsSetTypes.CHARACTERISTIC_EQUATION,
            EquationsSetCharacteristicSubtype]
  EquationsSetCharacteristic.CreateStart(EquationsSetUserNumberCharacteristic,Region,GeometricField,
      CharacteristicEquationsSetSpecification,EquationsSetFieldUserNumberCharacteristic,EquationsSetFieldCharacteristic)
  EquationsSetCharacteristic.CreateFinish()

  # Create the equations set for NAVIER-STOKES
  EquationsSetNavierStokes = iron.EquationsSet()
  EquationsSetFieldNavierStokes = iron.Field()
  # Set the equations set to be a dynamic nonlinear problem
  NavierStokesEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                EquationsSetSubtype]
  EquationsSetNavierStokes.CreateStart(EquationsSetUserNumberNavierStokes,Region,GeometricField,
      NavierStokesEquationsSetSpecification,EquationsSetFieldUserNumberNavierStokes,EquationsSetFieldNavierStokes)
  EquationsSetNavierStokes.CreateFinish()

  #------------------

  # Create the equations set for ADVECTION
  if (coupledAdvection):
      EquationsSetAdvection = iron.EquationsSet()
      EquationsSetFieldAdvection = iron.Field()
      # Set the equations set to be a dynamic linear problem
      AdvectionEquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                    iron.EquationsSetTypes.ADVECTION_EQUATION,
                EquationsSetAdvectionSubtype]
      EquationsSetAdvection.CreateStart(EquationsSetUserNumberAdvection,Region,GeometricField,
          AdvectionEquationsSetSpecification,EquationsSetFieldUserNumberAdvection,EquationsSetFieldAdvection)
      EquationsSetAdvection.CreateFinish()
# =================================

print('\033[1;32m'+'Equations Set     COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Dependent Field
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> DEPENDENT FIELD << == ")

if (not TestFlow):
  # Create the equations set dependent field variables
  dependentFieldEnergy = iron.Field()
  equationsSetEnergy.DependentCreateStart(dependentFieldUserNumberEnergy,dependentFieldEnergy)
  dependentFieldEnergy.LabelSet('Dependent Field')
  dependentFieldEnergy.VariableLabelSet(iron.FieldVariableTypes.U,'Blood Temperature')
  dependentFieldEnergy.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'Blood Temperature Gradient')
  dependentFieldEnergy.DOFOrderTypeSet(iron.FieldVariableTypes.U,iron.FieldDOFOrderTypes.SEPARATED)
  dependentFieldEnergy.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN,iron.FieldDOFOrderTypes.SEPARATED)
  equationsSetEnergy.DependentCreateFinish()

  # Initialise dependent field
  dependentFieldEnergy.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,37.0)

  dependentFieldEnergy.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
  dependentFieldEnergy.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)


  # Create dependent field
  dependentFieldTissue = iron.Field()
  equationsSetTissue.DependentCreateStart(dependentFieldUserNumberTissue,dependentFieldTissue)
  dependentFieldTissue.VariableLabelSet(iron.FieldVariableTypes.U,'Tissue Temperature')
  dependentFieldTissue.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'Tissue Temperature Gradient')
  dependentFieldTissue.DOFOrderTypeSet(iron.FieldVariableTypes.U,iron.FieldDOFOrderTypes.SEPARATED)
  dependentFieldTissue.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN,iron.FieldDOFOrderTypes.SEPARATED)
  equationsSetTissue.DependentCreateFinish()

  # Initialise dependent field
  dependentFieldTissue.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,36.3)

  dependentFieldTissue.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
  dependentFieldTissue.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  if (ProgressDiagnostics):
    print( " == >> DEPENDENT FIELD << == ")

  # STREE
  if (streeBoundaries):
      # Create the equations set dependent field variables
      DependentFieldStree = iron.Field()
      EquationsSetStree.DependentCreateStart(DependentFieldUserNumber3,DependentFieldStree)
      DependentFieldStree.VariableLabelSet(iron.FieldVariableTypes.U,'Stree 1st Variable')
      DependentFieldStree.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'Stree 2nd Variable')
      EquationsSetStree.DependentCreateFinish()

  #------------------

  # CHARACTERISTIC
  # Create the equations set dependent field variables
  DependentFieldNavierStokes = iron.Field()
  EquationsSetCharacteristic.DependentCreateStart(DependentFieldUserNumber,DependentFieldNavierStokes)
  DependentFieldNavierStokes.VariableLabelSet(iron.FieldVariableTypes.U,'General')
  DependentFieldNavierStokes.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'Derivatives')
  DependentFieldNavierStokes.VariableLabelSet(iron.FieldVariableTypes.V,'Characteristics')
  if (RCRBoundaries or Heart):
      DependentFieldNavierStokes.VariableLabelSet(iron.FieldVariableTypes.U1,'CellML Q and P')
  DependentFieldNavierStokes.VariableLabelSet(iron.FieldVariableTypes.U2,'Pressure')
  # Set the mesh component to be used by the field components.
  # Flow & Area
  DependentFieldNavierStokes.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,meshComponentNumberSpace)
  DependentFieldNavierStokes.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,meshComponentNumberSpace)
  # Derivatives
  DependentFieldNavierStokes.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,1,meshComponentNumberSpace)
  DependentFieldNavierStokes.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,2,meshComponentNumberSpace)
  # Riemann
  DependentFieldNavierStokes.ComponentMeshComponentSet(iron.FieldVariableTypes.V,1,meshComponentNumberSpace)
  DependentFieldNavierStokes.ComponentMeshComponentSet(iron.FieldVariableTypes.V,2,meshComponentNumberSpace)
  # qCellML & pCellml
  if (RCRBoundaries or Heart):
      DependentFieldNavierStokes.ComponentMeshComponentSet(iron.FieldVariableTypes.U1,1,meshComponentNumberSpace)
      DependentFieldNavierStokes.ComponentMeshComponentSet(iron.FieldVariableTypes.U1,2,meshComponentNumberSpace)
  # Pressure
  DependentFieldNavierStokes.ComponentMeshComponentSet(iron.FieldVariableTypes.U2,1,meshComponentNumberSpace)
  DependentFieldNavierStokes.ComponentMeshComponentSet(iron.FieldVariableTypes.U2,2,meshComponentNumberSpace)

  EquationsSetCharacteristic.DependentCreateFinish()

  #------------------

  # NAVIER-STOKES
  EquationsSetNavierStokes.DependentCreateStart(DependentFieldUserNumber,DependentFieldNavierStokes)
  EquationsSetNavierStokes.DependentCreateFinish()

  DependentFieldNavierStokes.ParameterSetCreate(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.PREVIOUS_VALUES)

  # Initialise the dependent field variables
  for nodeIdx in range (1,numberOfNodesSpace+1):
      nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberSpace)
      if (nodeDomain == computationalNodeNumber):
          if (nodeIdx in trifurcationNodeNumber):
              versions = [1,2,3,4]
          elif (nodeIdx in bifurcationNodeNumber):
              versions = [1,2,3]
          else:
              versions = [1]
          for versionIdx in versions:
              # U variables
              DependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              versionIdx,derivIdx,nodeIdx,1,Q[nodeIdx][versionIdx-1])
              DependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              versionIdx,derivIdx,nodeIdx,2,A[nodeIdx][versionIdx-1])
              DependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.PREVIOUS_VALUES,
              versionIdx,derivIdx,nodeIdx,1,Q[nodeIdx][versionIdx-1])
              DependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.PREVIOUS_VALUES,
              versionIdx,derivIdx,nodeIdx,2,A[nodeIdx][versionIdx-1])
              # delUdelN variables
              DependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.DELUDELN,iron.FieldParameterSetTypes.VALUES,
              versionIdx,derivIdx,nodeIdx,1,dQ[nodeIdx][versionIdx-1])
              DependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.DELUDELN,iron.FieldParameterSetTypes.VALUES,
              versionIdx,derivIdx,nodeIdx,2,dA[nodeIdx][versionIdx-1])

  # revert default version to 1
  versionIdx = 1

  # Finish the parameter update
  DependentFieldNavierStokes.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
  DependentFieldNavierStokes.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

  #------------------

  # ADVECTION
  if (coupledAdvection):
      # Create the equations set dependent field variables
      DependentFieldAdvection = iron.Field()
      EquationsSetAdvection.DependentCreateStart(DependentFieldUserNumber2,DependentFieldAdvection)
      DependentFieldAdvection.VariableLabelSet(iron.FieldVariableTypes.U,'Concentration')
      DependentFieldAdvection.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'Deriv')
      # Set the mesh component to be used by the field components.
      DependentFieldAdvection.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,meshComponentNumberConc)
      DependentFieldAdvection.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,1,meshComponentNumberConc)
      EquationsSetAdvection.DependentCreateFinish()

      # Initialise the dependent field variables
      for inputIdx in range (1,numberOfInputNodes+1):
          nodeIdx = inputNodeNumber[inputIdx-1]
          nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberConc)
          if (nodeDomain == computationalNodeNumber):
              DependentFieldAdvection.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              versionIdx,derivIdx,nodeIdx,1,0.0)

      # Finish the parameter update
      DependentFieldAdvection.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
      DependentFieldAdvection.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
# =================================

print('\033[1;32m'+'Dependent Field   COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Materials Field
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> MATERIALS FIELD << == ")

if (not TestFlow):
  # Create the equations set material field variables for artery energy
  materialsFieldEnergy = iron.Field()
  equationsSetEnergy.MaterialsCreateStart(materialsFieldUserNumberEnergy,materialsFieldEnergy)
  materialsFieldEnergy.VariableLabelSet(iron.FieldVariableTypes.U,'Blood Properties')
  materialsFieldEnergy.ComponentLabelSet(iron.FieldVariableTypes.U,1,'Blood Diffusivity')
  materialsFieldEnergy.ComponentLabelSet(iron.FieldVariableTypes.U,2,'Source Tb coeff.')
  equationsSetEnergy.MaterialsCreateFinish()

  # Initialise the properties and source values
  diffusivity=Alpha #+U*beta*le/2 #U*beta*le/2=0.000416667 almost 3000 times of the real diffusivity Pe=Ule/2a=0.2*0.05/12/2/0.0004=1
  materialsFieldEnergy.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,diffusivity)

  for elemIdx in arteriesElements:
      elemDomain = decompositionEnergy.ElementDomainGet(elemIdx)
      if elemDomain == computationalNodeNumber:
          # ra=(radius[elemIdx-1]+radius[elemIdx])/2
          cArtery=4*Alpha/ra[elemIdx-1]**2
          materialsFieldEnergy.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
    elemIdx,2, cArtery)

  materialsFieldEnergy.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
  materialsFieldEnergy.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)



  # Create the equations set material field variables for tissue
  materialsFieldTissue = iron.Field()
  equationsSetTissue.MaterialsCreateStart(materialsFieldUserNumberTissue,materialsFieldTissue)
  materialsFieldTissue.VariableLabelSet(iron.FieldVariableTypes.U,'Materials')
  materialsFieldTissue.ComponentLabelSet(iron.FieldVariableTypes.U,1,'Diffusivity 1')
  materialsFieldTissue.ComponentLabelSet(iron.FieldVariableTypes.U,2,'Diffusivity 2')
  materialsFieldTissue.ComponentLabelSet(iron.FieldVariableTypes.U,3,'Diffusivity 3')
  materialsFieldTissue.ComponentLabelSet(iron.FieldVariableTypes.U,4,'Source T coeff.')
  equationsSetTissue.MaterialsCreateFinish()


  # for elementNumber in muscleElements:
  #     elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
  #     if elementDomain == computationalNodeNumber:
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,1, k_ms*kFactor/(rho_ms*c_ms*rhoCFactor)) #0.42*1e-3
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,2, k_ms*kFactor/(rho_ms*c_ms*rhoCFactor))
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,3, k_ms*kFactor/(rho_ms*c_ms*rhoCFactor))
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,4, cMuscle)  #4.0e6*1e-9 0.42/4.0e6

  # for elementNumber in leftRadiusElements:
  #     elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
  #     if elementDomain == computationalNodeNumber:
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,1, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,2, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,3, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,4, cBone)

  # for elementNumber in leftHumerusElements:
  #     elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
  #     if elementDomain == computationalNodeNumber:
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,1, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,2, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,3, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,4, cBone)

  # for elementNumber in leftUlnaElements:
  #     elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
  #     if elementDomain == computationalNodeNumber:
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,1, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,2, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,3, k_bn*kFactor/(rho_bn*c_bn*rhoCFactor))
  #       materialsFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #   elementNumber,4, cBone)

  materialsFieldTissue.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
  k_ms*kFactor/(rho_ms*c_ms*rhoCFactor))
  materialsFieldTissue.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
  k_ms*kFactor/(rho_ms*c_ms*rhoCFactor))
  materialsFieldTissue.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
  k_ms*kFactor/(rho_ms*c_ms*rhoCFactor))
  materialsFieldTissue.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,cMuscle)

  materialsFieldTissue.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
  materialsFieldTissue.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  if (ProgressDiagnostics):
    print( " == >> MATERIALS FIELD << == ")

  # STREE
  if (streeBoundaries):
      # Create the equations set materials field variables
      MaterialsFieldStree = iron.Field()
      EquationsSetStree.MaterialsCreateStart(MaterialsFieldUserNumber2,MaterialsFieldStree)
      MaterialsFieldStree.VariableLabelSet(iron.FieldVariableTypes.U,'Stree Impedance')
      MaterialsFieldStree.VariableLabelSet(iron.FieldVariableTypes.V,'Stree Flow')
      EquationsSetStree.MaterialsCreateFinish()

  #------------------

  # CHARACTERISTIC
  # Create the equations set materials field variables
  MaterialsFieldNavierStokes = iron.Field()
  EquationsSetCharacteristic.MaterialsCreateStart(MaterialsFieldUserNumber,MaterialsFieldNavierStokes)
  MaterialsFieldNavierStokes.VariableLabelSet(iron.FieldVariableTypes.U,'MaterialsConstants')
  MaterialsFieldNavierStokes.VariableLabelSet(iron.FieldVariableTypes.V,'MaterialsVariables')
  # Set the mesh component to be used by the field components.
  for componentNumber in range(1,4):
      MaterialsFieldNavierStokes.ComponentMeshComponentSet(iron.FieldVariableTypes.V,componentNumber,meshComponentNumberSpace)
  EquationsSetCharacteristic.MaterialsCreateFinish()

  # NAVIER-STOKES
  EquationsSetNavierStokes.MaterialsCreateStart(MaterialsFieldUserNumber,MaterialsFieldNavierStokes)
  EquationsSetNavierStokes.MaterialsCreateFinish()

  # Set the materials field constants
  MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
  iron.FieldParameterSetTypes.VALUES,MaterialsFieldUserNumberMu,Mu)
  MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
  iron.FieldParameterSetTypes.VALUES,MaterialsFieldUserNumberRho,Rho)
  MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
  iron.FieldParameterSetTypes.VALUES,MaterialsFieldUserNumberAlpha,Alpha)
  MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
  iron.FieldParameterSetTypes.VALUES,MaterialsFieldUserNumberPext,Pext)
  MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
  iron.FieldParameterSetTypes.VALUES,MaterialsFieldUserNumberLs,Ls)
  MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
  iron.FieldParameterSetTypes.VALUES,MaterialsFieldUserNumberTs,Ts)
  MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
  iron.FieldParameterSetTypes.VALUES,MaterialsFieldUserNumberMs,Ms)
  MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
  iron.FieldParameterSetTypes.VALUES,MaterialsFieldUserNumberG0,G0)
  #MaterialsFieldNavierStokes.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
  # iron.FieldParameterSetTypes.VALUES,MaterialsFieldUserNumberFr,Fr)

  # Initialise the materials field variables (A0,E,H,kp,k1,k2,k3,b1)
  bifIdx = 0
  trifIdx = 0
  for nodeIdx in range(1,numberOfNodesSpace+1,1):
      nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberSpace)
      if (nodeDomain == computationalNodeNumber):
          if (nodeIdx in trifurcationNodeNumber):
              versions = [1,2,3,4]
          elif (nodeIdx in bifurcationNodeNumber):
              versions = [1,2,3]
          else:
              versions = [1]
          for versionIdx in versions:
              MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,
              versionIdx,derivIdx,nodeIdx,MaterialsFieldUserNumberA0,A0[nodeIdx][versionIdx-1])
              MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,
              versionIdx,derivIdx,nodeIdx,MaterialsFieldUserNumberE,E[nodeIdx][versionIdx-1])
              MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,
              versionIdx,derivIdx,nodeIdx,MaterialsFieldUserNumberH,H[nodeIdx][versionIdx-1])
  #            MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,
  #             versionIdx,derivIdx,nodeIdx,MaterialsFieldUserNumberkp,kp[nodeIdx][versionIdx-1])
  #            MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,
  #             versionIdx,derivIdx,nodeIdx,MaterialsFieldUserNumberk1,k1[nodeIdx][versionIdx-1])
  #            MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,
  #             versionIdx,derivIdx,nodeIdx,MaterialsFieldUserNumberk2,k2[nodeIdx][versionIdx-1])
  #            MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,
  #             versionIdx,derivIdx,nodeIdx,MaterialsFieldUserNumberk3,k3[nodeIdx][versionIdx-1])
  #            MaterialsFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,
  #             versionIdx,derivIdx,nodeIdx,MaterialsFieldUserNumberb1,b1[nodeIdx][versionIdx-1])

  # Finish the parameter update
  MaterialsFieldNavierStokes.ParameterSetUpdateStart(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES)
  MaterialsFieldNavierStokes.ParameterSetUpdateFinish(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES)

  #------------------

  # ADVECTION
  if (coupledAdvection):
      # Create the equations set materials field variables
      MaterialsFieldAdvection = iron.Field()
      EquationsSetAdvection.MaterialsCreateStart(MaterialsFieldUserNumber2,MaterialsFieldAdvection)
      MaterialsFieldAdvection.VariableLabelSet(iron.FieldVariableTypes.U,'Diffusivity')
      EquationsSetAdvection.MaterialsCreateFinish()
      # Set the materials field constant
      MaterialsFieldAdvection.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
      MaterialsFieldUserNumberD,D)
# =================================

print('\033[1;32m'+'Materials Field   COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Source Field
#================================================================================================================================

if (not TestFlow):
  #------------------------------
  #creating ulnar artery
  #-----------------------------

  #FileName_cell = "MaxVol1000/ulnarArteryFromWholeUpLi"

  #ulnarElementList=[]

  #with open(FileName_cell, "r") as f:
  #    target=f.readlines()
  #    for lineNum,line in enumerate(target):
  #        target[lineNum] = target[lineNum].rstrip('\n\r').replace("\t"," ").replace(","," ").split()
  #    for lineNum,line in enumerate(target):
  #        if lineNum !=0:
  #          ulnarElementList.append(int(target[lineNum][0])+1) # plus one because we are taking vtk numbers and we need to add 1.

  if (ProgressDiagnostics):
      print( " == >> SOURCE FIELD << == ")

  # Create source field for artery energy
  sourceFieldEnergy = iron.Field()
  equationsSetEnergy.SourceCreateStart(sourceFieldUserNumberEnergy,sourceFieldEnergy)
  equationsSetEnergy.SourceCreateFinish()

  # Source is b-cT; e.g. for my case C(Tt-T), b=CTt, c=C. Because source field is scalar type I cannot define 2 components.
  # So it is in 2nd component of the materials field.

  #sourceFieldEnergy.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,cArtery*Tt)
  for elemIdx in arteriesElements:
      elemDomain = decompositionEnergy.ElementDomainGet(elemIdx)
      if elemDomain == computationalNodeNumber:
          # ra=(radius[elemIdx-1]+radius[elemIdx])/2
          cArtery=4*Alpha/ra[elemIdx-1]**2
          sourceFieldEnergy.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elemIdx,1,cArtery*Tt[elemIdx])

  sourceFieldEnergy.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
  sourceFieldEnergy.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

  # Create source field for tissue
  sourceFieldTissue = iron.Field()
  equationsSetTissue.SourceCreateStart(sourceFieldUserNumberTissue,sourceFieldTissue)
  equationsSetTissue.SourceCreateFinish()

  # for elementNumber in muscleElements:
  #     elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
  #     if elementDomain == computationalNodeNumber:
  #       sourceFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementNumber,1,700.0e-9/(rho_ms*c_ms)+cMuscle*Tb) # source=qm/rhoC+CTb-CT
  # for elementNumber in leftRadiusElements:
  #     elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
  #     if elementDomain == computationalNodeNumber:
  #       sourceFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementNumber,1,0.0)
  # for elementNumber in leftHumerusElements:
  #     elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
  #     if elementDomain == computationalNodeNumber:
  #       sourceFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementNumber,1,0.0)
  # for elementNumber in leftUlnaElements:
  #     elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
  #     if elementDomain == computationalNodeNumber:
  #       sourceFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementNumber,1,0.0)

  sourceFieldTissue.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
  700.0e-9/(rho_ms*c_ms)+cMuscle*Tb)  # source=qm/rhoC+CTb-CT

  #for elementNumber in ulnarElementList:
  #    elementDomain = decomposition.ElementDomainGet(elementNumber)
  #    if elementDomain == computationalNodeNumber:
  ##      sourceField.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementNumber,1,2.5e6/2.0/4.0e6)
  #      sourceField.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementNumber,1,0.0)
  # Modifying source term for the element with the artery passing through.
  # for elementNumber in range(1,81):
  #   sourceFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementNumber,1,
  #     700/(rho_t*c_t)+CWall*Tb*volumeCorrection) # Because element volume here is much smaller than the artery passing it. We could consider more elements as well.

  sourceFieldTissue.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
  sourceFieldTissue.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# =================================
# F L O W
#if (CoupledBioheatFlow or TestFlow):
# =================================

print('\033[1;32m'+'Source Field      COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
# Independent Field
#================================================================================================================================

if (ProgressDiagnostics):
    print (" == >> INDEPENDENT FIELD << == ")

if (not TestFlow):
  # Create the equations set independent field variables
  IndependentFieldEnergy = iron.Field()
  #IndependentFieldEnergy.VariableLabelSet(iron.FieldVariableTypes.U,'flow velocity')
  # Set the mesh component to be used by the field components.
  #IndependentFieldEnergy.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)

  # NAVIER-STOKES
  equationsSetEnergy.IndependentCreateStart(independentFieldUserNumberEnergy,IndependentFieldEnergy)
  equationsSetEnergy.IndependentCreateFinish()

  # Set the velocity
  IndependentFieldEnergy.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
    1,1.0)

  # Finish the parameter update
  IndependentFieldEnergy.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
  IndependentFieldEnergy.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

  # Tissue independent field =================
  IndependentFieldTissue = iron.Field()
  IndependentFieldTissue.CreateStart(independentFieldUserNumberTissue,regionTissue)
  IndependentFieldTissue.VariableLabelSet(iron.FieldVariableTypes.U,'dynamic parameters')
  # Set type to general
  IndependentFieldTissue.TypeSet(iron.FieldTypes.GENERAL)
  # label the field
  IndependentFieldTissue.LabelSet('Control Parameters')
  # define new created field to be independent
  IndependentFieldTissue.DependentTypeSet(iron.FieldDependentTypes.INDEPENDENT)
  # Define decomposition
  IndependentFieldTissue.MeshDecompositionSet(decompositionTissue)
  # point new field to geometric field
  IndependentFieldTissue.GeometricFieldSet(geometricFieldTissue)
  # Set number of variables to 1
  IndependentFieldTissue.NumberOfVariablesSet(1)
  # Set the variable type to U
  IndependentFieldTissue.VariableTypesSet([iron.FieldVariableTypes.U])
  # Set the dimension to scalar
  IndependentFieldTissue.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.VECTOR)
  # Set number of components 2 for T_skin and T_core
  IndependentFieldTissue.NumberOfComponentsSet(iron.FieldVariableTypes.U,4)
  # Set the mesh component to be used by the field components.
  IndependentFieldTissue.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
  IndependentFieldTissue.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
  IndependentFieldTissue.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
  IndependentFieldTissue.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,1)

  #Set the label for them
  IndependentFieldTissue.ComponentLabelSet(iron.FieldVariableTypes.U,1,'T_skin')
  IndependentFieldTissue.ComponentLabelSet(iron.FieldVariableTypes.U,2,'T_core')
  IndependentFieldTissue.ComponentLabelSet(iron.FieldVariableTypes.U,3,'organ')
  IndependentFieldTissue.ComponentLabelSet(iron.FieldVariableTypes.U,4,'volume')

  # Set interpolation type
  IndependentFieldTissue.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.ELEMENT_BASED)
  IndependentFieldTissue.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.ELEMENT_BASED)
  IndependentFieldTissue.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.ELEMENT_BASED)
  IndependentFieldTissue.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.ELEMENT_BASED)
  # Set scaling type to none
  IndependentFieldTissue.ScalingTypeSet(iron.FieldScalingTypes.NONE)

  IndependentFieldTissue.CreateFinish()




  equationsSetTissue.IndependentCreateStart(independentFieldUserNumberTissue,IndependentFieldTissue)
  equationsSetTissue.IndependentCreateFinish()

  # Set the parameters
  IndependentFieldTissue.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
    1,37.0)
  IndependentFieldTissue.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
    2,37.0)
  IndependentFieldTissue.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
    3,1)

  # for elementNumber in kidneyElementsLeft:
  #   elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
  #   if elementDomain == computationalNodeNumber:
  #     IndependentFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #       elementNumber,4,kidneyLeft) 

  # for elementNumber in kidneyElementsRight:
  #   elementDomain = decompositionTissue.ElementDomainGet(elementNumber)
  #   if elementDomain == computationalNodeNumber:
  #     IndependentFieldTissue.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
  #       elementNumber,4,kidneyRight) 

  IndependentFieldTissue.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
    4,1)

  # Finish the parameter update
  IndependentFieldTissue.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
  IndependentFieldTissue.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  if (ProgressDiagnostics):
    print( " == >> INDEPENDENT FIELD << == ")

  # CHARACTERISTIC
  # Create the equations set independent field variables
  IndependentFieldNavierStokes = iron.Field()
  EquationsSetCharacteristic.IndependentCreateStart(IndependentFieldUserNumber,IndependentFieldNavierStokes)
  IndependentFieldNavierStokes.VariableLabelSet(iron.FieldVariableTypes.U,'Normal Wave Direction')
  # Set the mesh component to be used by the field components.
  IndependentFieldNavierStokes.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,meshComponentNumberSpace)
  EquationsSetCharacteristic.IndependentCreateFinish()

  #------------------

  # NAVIER-STOKES
  EquationsSetNavierStokes.IndependentCreateStart(IndependentFieldUserNumber,IndependentFieldNavierStokes)
  EquationsSetNavierStokes.IndependentCreateFinish()

  # Set the normal wave direction for arteries
  for nodeIdx in range(1,numberOfNodesSpace+1,1):
      nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberSpace)
      if (nodeDomain == computationalNodeNumber):
          if (Tp[nodeIdx] == 'Artery'):
              # Incoming
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              versionIdx,derivIdx,nodeIdx,compIdx,1.0)
          elif(Tp[nodeIdx] == 'Vein'):
              # Outgoing
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              versionIdx,derivIdx,nodeIdx,compIdx,-1.0)

  # Set the normal wave direction for bifurcation
  for bifIdx in range (1,numberOfBifurcations+1):
      nodeIdx = bifurcationNodeNumber[bifIdx-1]
      nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberSpace)
      if (nodeDomain == computationalNodeNumber):
          if (Tp[nodeIdx] == 'Artery'):
              # Incoming(parent)
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              1,derivIdx,nodeIdx,compIdx,1.0)
              # Outgoing(branches)
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              2,derivIdx,nodeIdx,compIdx,-1.0)
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              3,derivIdx,nodeIdx,compIdx,-1.0)
          elif(Tp[nodeIdx] == 'Vein'):
              # Outgoing(parent)
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              1,derivIdx,nodeIdx,compIdx,-1.0)
              # Incoming(branches)
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              2,derivIdx,nodeIdx,compIdx,1.0)
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              3,derivIdx,nodeIdx,compIdx,1.0)

  # Set the normal wave direction for trifurcation
  for trifIdx in range (1,numberOfTrifurcations+1):
      nodeIdx = trifurcationNodeNumber[trifIdx-1]
      nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberSpace)
      if (nodeDomain == computationalNodeNumber):
          if (Tp[nodeIdx] == 'Artery'):
              # Incoming(parent)
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              1,derivIdx,nodeIdx,compIdx,1.0)
              # Outgoing(branches)
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              2,derivIdx,nodeIdx,compIdx,-1.0)
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              3,derivIdx,nodeIdx,compIdx,-1.0)
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              4,derivIdx,nodeIdx,compIdx,-1.0)
          elif (Tp[nodeIdx] == 'Vein'):
              # Outgoing(parent)
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              1,derivIdx,nodeIdx,compIdx,-1.0)
              # Incoming(branches)
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              2,derivIdx,nodeIdx,compIdx,1.0)
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              3,derivIdx,nodeIdx,compIdx,1.0)
              IndependentFieldNavierStokes.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              4,derivIdx,nodeIdx,compIdx,1.0)

  # Finish the parameter update
  IndependentFieldNavierStokes.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
  IndependentFieldNavierStokes.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

  #------------------

  # ADVECTION
  if (coupledAdvection):
      # Create the equations set independent field variables
      IndependentFieldAdvection = iron.Field()
      EquationsSetAdvection.IndependentCreateStart(DependentFieldUserNumber,DependentFieldNavierStokes)
      EquationsSetAdvection.IndependentCreateFinish()
# =================================

print('\033[1;32m'+'Independent Field COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))

#================================================================================================================================
#  CellML fields: Shivering Model
#================================================================================================================================

if (not TestFlow):
  # Create the CellML environment
  CellMLShiv = iron.CellML()
  CellMLShiv.CreateStart(cellMLUserNumberTissue,regionTissue)

  # Shivering model
  ShiveringIdx = CellMLShiv.ModelImport("Shivering3.cellml")

  # Flag for known and wanted variables
  # Variables imported from OpenCMISS
  CellMLShiv.VariableSetAsKnown(ShiveringIdx, "Shivering/Tskin")
  CellMLShiv.VariableSetAsKnown(ShiveringIdx, "Shivering/Tcore")
  CellMLShiv.VariableSetAsKnown(ShiveringIdx, "thermophysicalProperties/OrganType")
  CellMLShiv.VariableSetAsKnown(ShiveringIdx, "Source/vol")

  # Variables to get from the CellML
  CellMLShiv.VariableSetAsWanted(ShiveringIdx, "Source/source")

  # Finish the CellML environment creation
  CellMLShiv.CreateFinish()

  # Create the CellML models field
  CellMLShiv.FieldMapsCreateStart()

  # Map the components
  CellMLShiv.CreateFieldToCellMLMap(IndependentFieldTissue, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                                ShiveringIdx, "Shivering/Tskin", iron.FieldParameterSetTypes.VALUES)
  CellMLShiv.CreateFieldToCellMLMap(IndependentFieldTissue, iron.FieldVariableTypes.U, 2, iron.FieldParameterSetTypes.VALUES,
                                ShiveringIdx, "Shivering/Tcore", iron.FieldParameterSetTypes.VALUES)
  CellMLShiv.CreateFieldToCellMLMap(IndependentFieldTissue, iron.FieldVariableTypes.U, 3, iron.FieldParameterSetTypes.VALUES,
                                ShiveringIdx, "thermophysicalProperties/OrganType", iron.FieldParameterSetTypes.VALUES)                              
  CellMLShiv.CreateFieldToCellMLMap(IndependentFieldTissue, iron.FieldVariableTypes.U, 4, iron.FieldParameterSetTypes.VALUES,
                                ShiveringIdx, "Source/vol", iron.FieldParameterSetTypes.VALUES)                              

  CellMLShiv.CreateCellMLToFieldMap(ShiveringIdx, "Source/source", iron.FieldParameterSetTypes.VALUES, sourceFieldTissue,
                                iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

  # Finish the creation of CellML <--> OpenCMISS field maps 
  CellMLShiv.FieldMapsCreateFinish()


  # Map fields

  # Create the CellML models field
  CellMLModelsField = iron.Field()
  CellMLShiv.ModelsFieldCreateStart(cellMLModelsFieldUserNumberTissue, CellMLModelsField)
  CellMLShiv.ModelsFieldCreateFinish()

  # Update cellmlModelsfield parameter set with the T_core?

  # Create the CellML state field
  # CellMLStateField = iron.Field()
  # CellML.StateFieldCreateStart(cellMLStateFieldUserNumber, CellMLStateField)
  # CellML.StateFieldCreateFinish()

  # Create the CellML parameters field
  CellMLParametersField = iron.Field()
  CellMLShiv.ParametersFieldCreateStart(cellMLParametersFieldUserNumberTissue, CellMLParametersField)
  CellMLShiv.ParametersFieldCreateFinish()
  # Create the CellML intermediate field
  CellMLIntermediateField = iron.Field()
  CellMLShiv.IntermediateFieldCreateStart(CellMLIntermediateFieldUserNumberTissue, CellMLIntermediateField)
  CellMLShiv.IntermediateFieldCreateFinish()

  # Finish the parameter update
  IndependentFieldTissue.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
  IndependentFieldTissue.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)  


#================================================================================================================================
# Analytic Field
#================================================================================================================================
# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):

  if (HeartInput):
      if (ProgressDiagnostics):
          print( " == >> ANALYTIC FIELD << == ")

      AnalyticFieldNavierStokes = iron.Field()
      EquationsSetNavierStokes.AnalyticCreateStart(iron.NavierStokesAnalyticFunctionTypes.FLOWRATE_AORTA,AnalyticFieldUserNumber,
      AnalyticFieldNavierStokes) # SplintFromFile,FlowrateAorta,FlowrateOlufsen
      AnalyticFieldNavierStokes.VariableLabelSet(iron.FieldVariableTypes.U,'Input Flow')
      EquationsSetNavierStokes.AnalyticCreateFinish()
  else:
      Q[1][0]=0.0

  #DOC-START cellml define field maps
  #================================================================================================================================
  #  RCR CellML Model Maps
  #================================================================================================================================

  if (RCRBoundaries):

      #----------------------------------------------------------------------------------------------------------------------------
      # Description
      #----------------------------------------------------------------------------------------------------------------------------
      # A CellML OD model is used to provide the impedance from the downstream vascular bed beyond the termination
      # point of the 1D model. This is iteratively coupled with the the 1D solver. In the case of a simple resistance
      # model, P=RQ, which is analogous to Ohm's law: V=IR. A variable map copies the guess for the FlowRate, Q at
      # the boundary from the OpenCMISS Dependent Field to the CellML equation, which then returns presssure, P.
      # The initial guess value for Q is taken from the previous time step or is 0 for t=0. In OpenCMISS this P value is
      # then used to compute a new Area value based on the P-A relationship and the Riemann variable W_2, which gives a
      # new value for Q until the values for Q and P converge within tolerance of the previous value.
      #----------------------------------------------------------------------------------------------------------------------------

      if (ProgressDiagnostics):
          print( " == >> RCR CELLML MODEL << == ")

      qCellMLComponent = 1
      pCellMLComponent = 2

      # Create the CellML environment
      CellML = iron.CellML()
      CellML.CreateStart(CellMLUserNumber,Region)
      # Number of CellML models
      CellMLModelIndex = [0]*(numberOfTerminalNodes+1)

      # Windkessel Model
      for terminalIdx in range (1,numberOfTerminalNodes+1):
          nodeIdx = coupledNodeNumber[terminalIdx-1]
          nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberSpace)
          print(('reading model: ' + "./input/CellMLModels/outlet/"+str(segment[nodeIdx])+"/ModelRCR.cellml"))
          if (nodeDomain == computationalNodeNumber):
              CellMLModelIndex[terminalIdx] = CellML.ModelImport("./input/CellMLModels/outlet/"+str(segment[nodeIdx])+"/ModelRCR.cellml")
              # known (to OpenCMISS) variables
              CellML.VariableSetAsKnown(CellMLModelIndex[terminalIdx],"Circuit/Qin")
              # to get from the CellML side
              CellML.VariableSetAsWanted(CellMLModelIndex[terminalIdx],"Circuit/Pout")
      CellML.CreateFinish()

      # Start the creation of CellML <--> OpenCMISS field maps
      CellML.FieldMapsCreateStart()

      # ModelIndex
      for terminalIdx in range (1,numberOfTerminalNodes+1):
          nodeIdx = coupledNodeNumber[terminalIdx-1]
          nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberSpace)
          if (nodeDomain == computationalNodeNumber):
              # Now we can set up the field variable component <--> CellML model variable mappings.
              # Map the OpenCMISS boundary flow rate values --> CellML
              # Q is component 1 of the DependentField
              CellML.CreateFieldToCellMLMap(DependentFieldNavierStokes,iron.FieldVariableTypes.U,1,
              iron.FieldParameterSetTypes.VALUES,CellMLModelIndex[terminalIdx],"Circuit/Qin",iron.FieldParameterSetTypes.VALUES)
              # Map the returned pressure values from CellML --> CMISS
              # pCellML is component 1 of the Dependent field U1 variable
              CellML.CreateCellMLToFieldMap(CellMLModelIndex[terminalIdx],"Circuit/Pout",iron.FieldParameterSetTypes.VALUES,
              DependentFieldNavierStokes,iron.FieldVariableTypes.U1,pCellMLComponent,iron.FieldParameterSetTypes.VALUES)

      # Finish the creation of CellML <--> OpenCMISS field maps
      CellML.FieldMapsCreateFinish()

      CellMLModelsField = iron.Field()
      CellML.ModelsFieldCreateStart(CellMLModelsFieldUserNumber,CellMLModelsField)
      CellML.ModelsFieldCreateFinish()

      # Set the models field at boundary nodes
      for terminalIdx in range (1,numberOfTerminalNodes+1):
          nodeIdx = coupledNodeNumber[terminalIdx-1]
          nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberSpace)
          if (nodeDomain == computationalNodeNumber):
              print(("Terminal node: " + str(nodeIdx) + " - " + str(segment[nodeIdx])))
              CellMLModelsField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              versionIdx,derivIdx,nodeIdx,1,CellMLModelIndex[terminalIdx])

      CellMLStateField = iron.Field()
      CellML.StateFieldCreateStart(CellMLStateFieldUserNumber,CellMLStateField)
      CellML.StateFieldCreateFinish()

      CellMLParametersField = iron.Field()
      CellML.ParametersFieldCreateStart(CellMLParametersFieldUserNumber,CellMLParametersField)
      CellML.ParametersFieldCreateFinish()

      CellMLIntermediateField = iron.Field()
      CellML.IntermediateFieldCreateStart(CellMLIntermediateFieldUserNumber,CellMLIntermediateField)
      CellML.IntermediateFieldCreateFinish()

      # Finish the parameter update
      DependentFieldNavierStokes.ParameterSetUpdateStart(iron.FieldVariableTypes.U1,iron.FieldParameterSetTypes.VALUES)
      DependentFieldNavierStokes.ParameterSetUpdateFinish(iron.FieldVariableTypes.U1,iron.FieldParameterSetTypes.VALUES)
  # DOC-END cellml define field maps

  #================================================================================================================================
  #  Heart CellML Model Maps
  #================================================================================================================================

  if (Heart):

      #----------------------------------------------------------------------------------------------------------------------------
      # Description
      #----------------------------------------------------------------------------------------------------------------------------
      # A CellML OD model is used to provide Heart model for the 1D model. A variable map copies the guess for the Pressure, P at
      # the inlet node from the OpenCMISS Dependent Field to the CellML equation, which then returns flow, Q. The initial guess
      # value for P is taken from the previous time step or is 0 for t=0. In OpenCMISS this Q value is then imposed to the inlet.
      #----------------------------------------------------------------------------------------------------------------------------

      if (ProgressDiagnostics):
          print( " == >> HEART CELLML MODEL << == ")

      qCellMLComponent = 1
      pCellMLComponent = 2

      # Create the CellML environment
      CellML = iron.CellML()
      CellML.CreateStart(CellMLUserNumber,Region)
      # Number of CellML models
      CellMLModelIndex = [0]*(numberOfInputNodes+1)

      # Heart Model
      for inputIdx in range (1,numberOfInputNodes+1):
          nodeIdx = inputNodeNumber[inputIdx-1]
          nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberSpace)
          print(('reading model: ' + "./input/CellMLModels/inlet/"+str(segment[nodeIdx])+"/Heart.cellml"))
          if (nodeDomain == computationalNodeNumber):
              CellMLModelIndex[inputIdx] = CellML.ModelImport("./input/CellMLModels/inlet/"+str(segment[nodeIdx])+"/Heart.cellml")
              # known (to OpenCMISS) variables
              CellML.VariableSetAsKnown(CellMLModelIndex[inputIdx],"Heart/P_art")
              # to get from the CellML side
              CellML.VariableSetAsWanted(CellMLModelIndex[inputIdx],"Heart/Q_art")
      CellML.CreateFinish()

      # Start the creation of CellML <--> OpenCMISS field maps
      CellML.FieldMapsCreateStart()

      # ModelIndex
      for inputIdx in range (1,numberOfInputNodes+1):
          nodeIdx = inputNodeNumber[inputIdx-1]
          nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberSpace)
          if (nodeDomain == computationalNodeNumber):
              # Now we can set up the field variable component <--> CellML model variable mappings.
              # Map the OpenCMISS boundary flow rate values --> CellML
              # P is component 1 of the Dependent field U2 variable
              CellML.CreateFieldToCellMLMap(DependentFieldNavierStokes,iron.FieldVariableTypes.U2,1,
              iron.FieldParameterSetTypes.VALUES,CellMLModelIndex[inputIdx],"Heart/P_art",iron.FieldParameterSetTypes.VALUES)
              # Map the returned pressure values from CellML --> CMISS
              # qCellML is component 1 of the Dependent field U1 variable
              CellML.CreateCellMLToFieldMap(CellMLModelIndex[inputIdx],"Heart/Q_art",iron.FieldParameterSetTypes.VALUES,
              DependentFieldNavierStokes,iron.FieldVariableTypes.U1,qCellMLComponent,iron.FieldParameterSetTypes.VALUES)

      # Finish the creation of CellML <--> OpenCMISS field maps
      CellML.FieldMapsCreateFinish()

      CellMLModelsField = iron.Field()
      CellML.ModelsFieldCreateStart(CellMLModelsFieldUserNumber,CellMLModelsField)
      CellML.ModelsFieldCreateFinish()

      # Set the models field at inlet boundary nodes
      for inputIdx in range (1,numberOfInputNodes+1):
          nodeIdx = inputNodeNumber[inputIdx-1]
          nodeDomain = Decomposition.NodeDomainGet(nodeIdx,meshComponentNumberSpace)
          if (nodeDomain == computationalNodeNumber):
              print(("Input node: " + str(nodeIdx) + " - " + str(segment[nodeIdx])))
              CellMLModelsField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
              versionIdx,derivIdx,nodeIdx,1,CellMLModelIndex[inputIdx])

      CellMLStateField = iron.Field()
      CellML.StateFieldCreateStart(CellMLStateFieldUserNumber,CellMLStateField)
      CellML.StateFieldCreateFinish()

      CellMLParametersField = iron.Field()
      CellML.ParametersFieldCreateStart(CellMLParametersFieldUserNumber,CellMLParametersField)
      CellML.ParametersFieldCreateFinish()

      CellMLIntermediateField = iron.Field()
      CellML.IntermediateFieldCreateStart(CellMLIntermediateFieldUserNumber,CellMLIntermediateField)
      CellML.IntermediateFieldCreateFinish()

      # Finish the parameter update
      DependentFieldNavierStokes.ParameterSetUpdateStart(iron.FieldVariableTypes.U1,iron.FieldParameterSetTypes.VALUES)
      DependentFieldNavierStokes.ParameterSetUpdateFinish(iron.FieldVariableTypes.U1,iron.FieldParameterSetTypes.VALUES)
# =================================

#================================================================================================================================
#  Equations
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> EQUATIONS << == ")

if (not TestFlow):
  # Create equations for artery energy
  equationsEnergy = iron.Equations()
  equationsSetEnergy.EquationsCreateStart(equationsEnergy)
  equationsEnergy.sparsityType = iron.EquationsSparsityTypes.SPARSE
  equationsEnergy.outputType = iron.EquationsOutputTypes.NONE
  equationsSetEnergy.EquationsCreateFinish()

    # I want to solve this type of equation, dT/dt-div(Sigma grad(T))-(b-cT)=0.
    # Sigma in my case is Sigma=k/rhoC.
    # q=(Sigma grad(T)).n which in my 1D case is q=k/rhoC dT/dx which in the boundary is q=-h/rhoC (T-Tinf)
    # in my case b=qm/rhoC+CTb and c=C

  # Create equations for tissue
  equationsTissue = iron.Equations()
  equationsSetTissue.EquationsCreateStart(equationsTissue)
  equationsTissue.sparsityType = iron.EquationsSparsityTypes.SPARSE
  equationsTissue.outputType = iron.EquationsOutputTypes.NONE
  equationsSetTissue.EquationsCreateFinish()

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  if (ProgressDiagnostics):
    print( " == >> EQUATIONS << == ")

  # 1th Equations Set - STREE
  if (streeBoundaries):
      EquationsStree = iron.Equations()
      EquationsSetStree.EquationsCreateStart(EquationsStree)
      EquationsStree.sparsityType = iron.EquationsSparsityTypes.SPARSE
      # (NONE/TIMING/MATRIX/ELEMENT_MATRIX/NODAL_MATRIX)
      EquationsStree.outputType = iron.EquationsOutputTypes.NONE
      EquationsSetStree.EquationsCreateFinish()

  #------------------

  # 2nd Equations Set - CHARACTERISTIC
  EquationsCharacteristic = iron.Equations()
  EquationsSetCharacteristic.EquationsCreateStart(EquationsCharacteristic)
  EquationsCharacteristic.sparsityType = iron.EquationsSparsityTypes.SPARSE
  # (NONE/TIMING/MATRIX/ELEMENT_MATRIX/NODAL_MATRIX)
  EquationsCharacteristic.outputType = iron.EquationsOutputTypes.NONE
  EquationsSetCharacteristic.EquationsCreateFinish()

  #------------------

  # 3rd Equations Set - NAVIER-STOKES
  EquationsNavierStokes = iron.Equations()
  EquationsSetNavierStokes.EquationsCreateStart(EquationsNavierStokes)
  EquationsNavierStokes.sparsityType = iron.EquationsSparsityTypes.FULL
  EquationsNavierStokes.lumpingType = iron.EquationsLumpingTypes.UNLUMPED
  # (NONE/TIMING/MATRIX/ELEMENT_MATRIX/NODAL_MATRIX)
  EquationsNavierStokes.outputType = iron.EquationsOutputTypes.NONE
  EquationsSetNavierStokes.EquationsCreateFinish()

  #------------------

  # 4th Equations Set - ADVECTION
  if (coupledAdvection):
      EquationsAdvection = iron.Equations()
      EquationsSetAdvection.EquationsCreateStart(EquationsAdvection)
      EquationsAdvection.sparsityType = iron.EquationsSparsityTypes.SPARSE
      # (NONE/TIMING/MATRIX/ELEMENT_MATRIX/NODAL_MATRIX)
      EquationsAdvection.outputType = iron.EquationsOutputTypes.NONE
      EquationsSetAdvection.EquationsCreateFinish()
# =================================

print('\033[1;32m'+'equations         COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Problems
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> PROBLEM << == ")
if(not TestFlow):
  # Start the creation of a problem.
  problem = iron.Problem()
  problemSpecification = [iron.ProblemClasses.MULTI_PHYSICS,
                          iron.ProblemTypes.DIFFUSION_ADVECTION_DIFFUSION,
                          ProblemSubtype]
  problem.CreateStart(problemUserNumber,problemSpecification)
  problem.CreateFinish()

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  if (TestFlow):
    if (ProgressDiagnostics):
        print( " == >> PROBLEM << == ")

    # Start the creation of a problem.
    problem = iron.Problem()
    problemSpecification = [iron.ProblemClasses.FLUID_MECHANICS,
                            iron.ProblemTypes.NAVIER_STOKES_EQUATION,
                            ProblemSubtype]
    problem.CreateStart(problemUserNumber,problemSpecification)
    problem.CreateFinish()
# =================================

print('\033[1;32m'+'Problems          COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Control Loops
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> PROBLEM CONTROL LOOP << == ")

# Create control loops
problem.ControlLoopCreateStart()
TimeLoop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],TimeLoop)
TimeLoop.LabelSet('Time Loop')
if (TestFlow):
  TimeLoop.TimesSet(startTime,stopTime,timeIncrement)
  TimeLoop.TimeOutputSet(DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY)
  TimeLoop.OutputTypeSet(iron.ControlLoopOutputTypes.TIMING)
else:
  TimeLoop.TimesSet(startTimeBioheat,stopTimeBioheat,timeIncrementBioheat)
  TimeLoop.TimeOutputSet(DYNAMIC_SOLVER_DIFFUSION_OUTPUT_FREQUENCY)


# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  if (ProgressDiagnostics):
      print( " == >> PROBLEM CONTROL LOOP << == ")

  '''
    Solver Control Loops

                    L1                                 L2                        L3

  1D-Stree
  ------


                                                        | 1) 0D Simple subloop   | 1) 0D/Structured tree Solver
                                                        |
      Time Loop, L0  | 1) 1D-0D Iterative Coupling, L1  | 2) 1D NS/C coupling:   | 1) Characteristic Nonlinear Solver
                    |    Convergence Loop (while loop) |    (while loop)        | 2) 1DNavierStokes Transient Solver
                    |
                    | 2) (optional) Simple subloop     | 1) Advection Linear Solver


  1D0D
  ------


                                                        | 1) 0D Simple subloop   | 1) 0D/CellML DAE Solver
                                                        |
      Time Loop, L0  | 1) 1D-0D Iterative Coupling, L1  | 2) 1D NS/C coupling:   | 1) Characteristic Nonlinear Solver
                    |    Convergence Loop (while loop) |    (while loop)        | 2) 1DNavierStokes Transient Solver
                    |
                    | 2) (optional) Simple subloop     | 1) Advection Linear Solver


  1D
  ------


      Time Loop, L0  | 1) 1D NS/C coupling subloop      | 1) Characteristic Nonlinear Solver
                    |    (while loop)                  | 2) 1DNavierStokes Transient Solver
                    |
                    | 2) (optional) Simple subloop     | 1) Advection Linear Solver


  '''

  # Order of solvers within their respective subloops
  SolverCharacteristicUserNumber = 1
  SolverNavierStokesUserNumber   = 2
  SolverAdvectionUserNumber      = 1
  SolverCellmlUserNumber         = 1
  if (RCRBoundaries or streeBoundaries or Heart):
    Iterative1d0dControlLoopNumber   = 1
    SimpleAdvectionControlLoopNumber = 2
    Simple0DControlLoopNumber        = 1
    Iterative1dControlLoopNumber     = 2
  else:
    Iterative1dControlLoopNumber     = 1
    SimpleAdvectionControlLoopNumber = 2

  # Start the creation of the problem control loop
  # TimeLoop = iron.ControlLoop()
  # Problem.ControlLoopCreateStart()
  # Problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],TimeLoop)
  # TimeLoop.LabelSet('Time Loop')
  # TimeLoop.TimesSet(startTime,stopTime,timeIncrement)
  # TimeLoop.TimeOutputSet(DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY)
  # TimeLoop.OutputTypeSet(iron.ControlLoopOutputTypes.TIMING)

  # Set tolerances for iterative convergence loops
  if (RCRBoundaries or streeBoundaries or Heart):
      Iterative1DCouplingLoop = iron.ControlLoop()
      problem.ControlLoopGet([Iterative1d0dControlLoopNumber,Iterative1dControlLoopNumber,
      iron.ControlLoopIdentifiers.NODE],Iterative1DCouplingLoop)
      Iterative1DCouplingLoop.AbsoluteToleranceSet(couplingTolerance1D)
      Iterative1D0DCouplingLoop = iron.ControlLoop()
      problem.ControlLoopGet([Iterative1d0dControlLoopNumber,iron.ControlLoopIdentifiers.NODE],
      Iterative1D0DCouplingLoop)
      Iterative1D0DCouplingLoop.AbsoluteToleranceSet(couplingTolerance1D0D)
  else:
      Iterative1DCouplingLoop = iron.ControlLoop()
      if (TestFlow):
        problem.ControlLoopGet([Iterative1dControlLoopNumber,iron.ControlLoopIdentifiers.NODE],Iterative1DCouplingLoop)
        Iterative1DCouplingLoop.AbsoluteToleranceSet(couplingTolerance1D)

  # problem.ControlLoopCreateFinish()
# =================================
problem.ControlLoopCreateFinish()

print('\033[1;32m'+'Control Loops     COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Solvers
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> SOLVERS << == ")

if (TestFlow):
  if (ProgressDiagnostics):
    print( " == >> SOLVERS << == ")

  # Start the creation of the problem solvers
  DynamicSolverNavierStokes     = iron.Solver()
  NonlinearSolverNavierStokes   = iron.Solver()
  LinearSolverNavierStokes      = iron.Solver()
  NonlinearSolverCharacteristic = iron.Solver()
  LinearSolverCharacteristic    = iron.Solver()
  if (streeBoundaries):
      LinearSolverStree         = iron.Solver()
  if (coupledAdvection):
      DynamicSolverAdvection    = iron.Solver()
      LinearSolverAdvection     = iron.Solver()

  problem.SolversCreateStart()

  #------------------

  # 1st Solver, Simple 0D subloop - STREE
  if (streeBoundaries):
      problem.SolverGet([Iterative1d0dControlLoopNumber,Simple0DControlLoopNumber,
      iron.ControlLoopIdentifiers.NODE],SolverStreeUserNumber,LinearSolverStree)
      # Set the nonlinear Jacobian type
      LinearSolverStree.OutputTypeSet(CMISS_SOLVER_OUTPUT_TYPE)

  #------------------

  # 1st Solver, Simple 0D subloop - CellML
  if (RCRBoundaries or Heart):
      CellMLSolver = iron.Solver()
      problem.SolverGet([Iterative1d0dControlLoopNumber,Simple0DControlLoopNumber,
      iron.ControlLoopIdentifiers.NODE],SolverDAEUserNumber,CellMLSolver)
      CellMLSolver.OutputTypeSet(CMISS_SOLVER_OUTPUT_TYPE)

  #------------------

  # 1st Solver, Iterative 1D subloop - CHARACTERISTIC
  if (RCRBoundaries or streeBoundaries or Heart):
      problem.SolverGet([Iterative1d0dControlLoopNumber,Iterative1dControlLoopNumber,
      iron.ControlLoopIdentifiers.NODE],SolverCharacteristicUserNumber,NonlinearSolverCharacteristic)
  else:
      problem.SolverGet([Iterative1dControlLoopNumber,iron.ControlLoopIdentifiers.NODE],
      SolverCharacteristicUserNumber,NonlinearSolverCharacteristic)
  # Set the nonlinear Jacobian type
  NonlinearSolverCharacteristic.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS) #(FD/EQUATIONS)
  NonlinearSolverCharacteristic.OutputTypeSet(NONLINEAR_SOLVER_CHARACTERISTIC_OUTPUT_TYPE)
  # Set the solver settings
  NonlinearSolverCharacteristic.NewtonAbsoluteToleranceSet(absoluteToleranceNonlinearCharacteristic)
  NonlinearSolverCharacteristic.NewtonSolutionToleranceSet(solutionToleranceNonlinearCharacteristic)
  NonlinearSolverCharacteristic.NewtonRelativeToleranceSet(relativeToleranceNonlinearCharacteristic)
  # Get the nonlinear linear solver
  NonlinearSolverCharacteristic.NewtonLinearSolverGet(LinearSolverCharacteristic)
  LinearSolverCharacteristic.OutputTypeSet(LINEAR_SOLVER_CHARACTERISTIC_OUTPUT_TYPE)
  # Set the solver settings
  LinearSolverCharacteristic.LinearTypeSet(iron.LinearSolverTypes.ITERATIVE)
  LinearSolverCharacteristic.LinearIterativeMaximumIterationsSet(MAXIMUM_ITERATIONS)
  LinearSolverCharacteristic.LinearIterativeDivergenceToleranceSet(DIVERGENCE_TOLERANCE)
  LinearSolverCharacteristic.LinearIterativeRelativeToleranceSet(relativeToleranceLinearCharacteristic)
  LinearSolverCharacteristic.LinearIterativeAbsoluteToleranceSet(absoluteToleranceLinearCharacteristic)
  LinearSolverCharacteristic.LinearIterativeGMRESRestartSet(RESTART_VALUE)

  #------------------

  # 2nd Solver, Iterative 1D subloop - NAVIER-STOKES
  if (RCRBoundaries or streeBoundaries or Heart):
      problem.SolverGet([Iterative1d0dControlLoopNumber,Iterative1dControlLoopNumber,
      iron.ControlLoopIdentifiers.NODE],SolverNavierStokesUserNumber,DynamicSolverNavierStokes)
  else:
      problem.SolverGet([Iterative1dControlLoopNumber,iron.ControlLoopIdentifiers.NODE],
      SolverNavierStokesUserNumber,DynamicSolverNavierStokes)
  DynamicSolverNavierStokes.OutputTypeSet(DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE)
  DynamicSolverNavierStokes.DynamicThetaSet(dynamicSolverNavierStokesTheta)
  # Get the dynamic nonlinear solver
  DynamicSolverNavierStokes.DynamicNonlinearSolverGet(NonlinearSolverNavierStokes)
  # Set the nonlinear Jacobian type
  NonlinearSolverNavierStokes.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS) #(FD/EQUATIONS)
  NonlinearSolverNavierStokes.OutputTypeSet(NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE)

  # Set the solver settings
  NonlinearSolverNavierStokes.NewtonAbsoluteToleranceSet(absoluteToleranceNonlinearNavierStokes)
  NonlinearSolverNavierStokes.NewtonSolutionToleranceSet(solutionToleranceNonlinearNavierStokes)
  NonlinearSolverNavierStokes.NewtonRelativeToleranceSet(relativeToleranceNonlinearNavierStokes)
  # Get the dynamic nonlinear linear solver
  NonlinearSolverNavierStokes.NewtonLinearSolverGet(LinearSolverNavierStokes)
  LinearSolverNavierStokes.OutputTypeSet(LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE)
  # Set the solver settings
  LinearSolverNavierStokes.LinearTypeSet(iron.LinearSolverTypes.ITERATIVE)
  LinearSolverNavierStokes.LinearIterativeMaximumIterationsSet(MAXIMUM_ITERATIONS)
  LinearSolverNavierStokes.LinearIterativeDivergenceToleranceSet(DIVERGENCE_TOLERANCE)
  LinearSolverNavierStokes.LinearIterativeRelativeToleranceSet(relativeToleranceLinearNavierStokes)
  LinearSolverNavierStokes.LinearIterativeAbsoluteToleranceSet(absoluteToleranceLinearNavierStokes)
  LinearSolverNavierStokes.LinearIterativeGMRESRestartSet(RESTART_VALUE)

  #------------------

  # 1st Solver, Simple advection subloop - ADVECTION
  if (coupledAdvection):
      problem.SolverGet([SimpleAdvectionControlLoopNumber,iron.ControlLoopIdentifiers.NODE],
      SolverAdvectionUserNumber,DynamicSolverAdvection)
      DynamicSolverAdvection.OutputTypeSet(DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE)
      DynamicSolverAdvection.DynamicThetaSet(dynamicSolverAdvectionTheta)
      # Get the dynamic linear solver
      DynamicSolverAdvection.DynamicLinearSolverGet(LinearSolverAdvection)

  # Finish the creation of the problem solver
  problem.SolversCreateFinish()
else:
  # Create problem solver
  solverEnergy = iron.Solver()
  LinearSolverEnergy = iron.Solver()

  solverTissue = iron.Solver()
  LinearSolverTissue = iron.Solver()

  problem.SolversCreateStart()
  problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,solverEnergy)
  problem.SolverGet([iron.ControlLoopIdentifiers.NODE],3,solverTissue)
  solverEnergy.LabelSet('Arterial Energy Solver')
  solverTissue.LabelSet('Tissue Energy Solver')
  #solver.outputType = iron.SolverOutputTypes.SOLVER
  solverEnergy.DynamicLinearSolverGet(LinearSolverEnergy)
  solverTissue.DynamicLinearSolverGet(LinearSolverTissue)
  #solver.linearType = iron.LinearSolverTypes.ITERATIVE
  #solver.linearIterativeAbsoluteTolerance = 1.0E-12
  #solver.linearIterativeRelativeTolerance = 1.0E-12
  problem.SolversCreateFinish()

# =================================
# F L O W
#if (CoupledBioheatFlow or TestFlow):
# =================================

print('\033[1;32m'+'Solvers           COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Solver Equations
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> SOLVER EQUATIONS << == ")

if (TestFlow):
  if (ProgressDiagnostics):
    print( " == >> SOLVER EQUATIONS << == ")

  # Start the creation of the problem solver equations
  NonlinearSolverCharacteristic = iron.Solver()
  SolverEquationsCharacteristic = iron.SolverEquations()
  DynamicSolverNavierStokes     = iron.Solver()
  SolverEquationsNavierStokes   = iron.SolverEquations()
  if (streeBoundaries):
      LinearSolverStree         = iron.Solver()
      SolverEquationsStree      = iron.SolverEquations()
  if (coupledAdvection):
      DynamicSolverAdvection    = iron.Solver()
      SolverEquationsAdvection  = iron.SolverEquations()

  problem.SolverEquationsCreateStart()

  #------------------

  # STREE Solver
  if (streeBoundaries):
      problem.SolverGet([Iterative1d0dControlLoopNumber,Simple0DControlLoopNumber,
      iron.ControlLoopIdentifiers.NODE],SolverStreeUserNumber,LinearSolverStree)
      LinearSolverStree.SolverEquationsGet(SolverEquationsStree)
      SolverEquationsStree.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
      # Add in the equations set
      EquationsSetStree = SolverEquationsStree.EquationsSetAdd(EquationsSetStree)

  #------------------

  # CellML Solver
  if (RCRBoundaries or Heart):
      CellMLSolver = iron.Solver()
      CellMLEquations = iron.CellMLEquations()
      problem.CellMLEquationsCreateStart()
      problem.SolverGet([Iterative1d0dControlLoopNumber,Simple0DControlLoopNumber,
      iron.ControlLoopIdentifiers.NODE],SolverDAEUserNumber,CellMLSolver)
      CellMLSolver.CellMLEquationsGet(CellMLEquations)
      # Add in the equations set
      CellMLEquations.CellMLAdd(CellML)
      problem.CellMLEquationsCreateFinish()

  #------------------

  # CHARACTERISTIC solver
  if (RCRBoundaries or streeBoundaries or Heart):
      problem.SolverGet([Iterative1d0dControlLoopNumber,Iterative1dControlLoopNumber,
      iron.ControlLoopIdentifiers.NODE],SolverCharacteristicUserNumber,NonlinearSolverCharacteristic)
  else:
      problem.SolverGet([Iterative1dControlLoopNumber,iron.ControlLoopIdentifiers.NODE],
      SolverCharacteristicUserNumber,NonlinearSolverCharacteristic)
  NonlinearSolverCharacteristic.SolverEquationsGet(SolverEquationsCharacteristic)
  SolverEquationsCharacteristic.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
  # Add in the equations set
  EquationsSetCharacteristic = SolverEquationsCharacteristic.EquationsSetAdd(EquationsSetCharacteristic)

  #------------------

  #  NAVIER-STOKES solver
  if (RCRBoundaries or streeBoundaries or Heart):
      problem.SolverGet([Iterative1d0dControlLoopNumber,Iterative1dControlLoopNumber,
      iron.ControlLoopIdentifiers.NODE],SolverNavierStokesUserNumber,DynamicSolverNavierStokes)
  else:
      problem.SolverGet([Iterative1dControlLoopNumber,iron.ControlLoopIdentifiers.NODE],
      SolverNavierStokesUserNumber,DynamicSolverNavierStokes)
  DynamicSolverNavierStokes.SolverEquationsGet(SolverEquationsNavierStokes)
  SolverEquationsNavierStokes.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
  # Add in the equations set
  EquationsSetNavierStokes = SolverEquationsNavierStokes.EquationsSetAdd(EquationsSetNavierStokes)

  #------------------

  # ADVECTION Solver
  if (coupledAdvection):
      problem.SolverGet([SimpleAdvectionControlLoopNumber,iron.ControlLoopIdentifiers.NODE],
      SolverAdvectionUserNumber,DynamicSolverAdvection)
      DynamicSolverAdvection.SolverEquationsGet(SolverEquationsAdvection)
      SolverEquationsAdvection.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
      # Add in the equations set
      EquationsSetAdvection = SolverEquationsAdvection.EquationsSetAdd(EquationsSetAdvection)

  # Finish the creation of the problem solver equations
  problem.SolverEquationsCreateFinish()
else:
  # Create solver equations and add equations set to solver equations
  solverEnergy = iron.Solver()
  solverTissue = iron.Solver()
  cellMLSolver = iron.Solver()

  solverEquationsEnergy = iron.SolverEquations()
  solverEquationsTissue = iron.SolverEquations()
  cellMLEquations = iron.CellMLEquations()

  problem.SolverEquationsCreateStart()
  problem.CellMLEquationsCreateStart()

  problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,cellMLSolver)
  problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,solverEnergy)
  problem.SolverGet([iron.ControlLoopIdentifiers.NODE],3,solverTissue)

  cellMLSolver.CellMLEquationsGet(cellMLEquations)
  solverEnergy.SolverEquationsGet(solverEquationsEnergy)
  solverTissue.SolverEquationsGet(solverEquationsTissue)

  solverEquationsEnergy.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
  equationsSetIndex1 = cellMLEquations.CellMLAdd(CellMLShiv)
  equationsSetIndex2 = solverEquationsEnergy.EquationsSetAdd(equationsSetEnergy)
  equationsSetIndex3 = solverEquationsTissue.EquationsSetAdd(equationsSetTissue)

  problem.SolverEquationsCreateFinish()
  problem.CellMLEquationsCreateFinish()

# =================================
# F L O W
#if (CoupledBioheatFlow or TestFlow):
# =================================

print('\033[1;32m'+'Solver Equations  COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))
#================================================================================================================================
#  Boundary Conditions
#================================================================================================================================

if (ProgressDiagnostics):
    print( " == >> BOUNDARY CONDITIONS << == ")

if (not TestFlow):
  boundaryConditionsEnergy = iron.BoundaryConditions()
  boundaryConditionsTissue = iron.BoundaryConditions()

  solverEquationsEnergy.BoundaryConditionsCreateStart(boundaryConditionsEnergy)

  nodeDomain = decompositionEnergy.NodeDomainGet(1,1)
  if nodeDomain == computationalNodeNumber:
      boundaryConditionsEnergy.SetNode(dependentFieldEnergy,iron.FieldVariableTypes.U,1,1,1,1,
          iron.BoundaryConditionsTypes.FIXED,[37.0])

  dependentFieldEnergy.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
  dependentFieldEnergy.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

  solverEquationsEnergy.BoundaryConditionsCreateFinish()

  #Tissue boundary conditions =============
  solverEquationsTissue.BoundaryConditionsCreateStart(boundaryConditionsTissue)

  nodes = iron.Nodes()
  regionTissue.NodesGet(nodes)


  # Actually in OpenCMISS the equation is divided by rhoC. So q=alpha*gradT.n which alpha=sigma/rhoC.
  # So for Robin BCs you need to pass h/rhoC and q_h/rhoC. Also DON'T FORGET ABOUT THE UNITS
  q_hUnit = 1#1e-6
  hUnit   = 1#1e-6


  Rtot=1/((h_conv+hr_rad)*hUnit)+3/k_sk # units terms has mm2 so no 1e6.

  for nodeNumber in boundaryTissue:
    nodeDomain = decompositionTissue.NodeDomainGet(nodeNumber,1)
    if nodeDomain == computationalNodeNumber:
      dependentFieldTissue.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
        1,1,nodeNumber,1,20.0)
      boundaryConditionsTissue.SetNode(dependentFieldTissue,iron.FieldVariableTypes.DELUDELN,1,1,nodeNumber,1,
        iron.BoundaryConditionsTypes.ROBIN,[1.0/(rho_sk*c_sk*rhoCFactor*Rtot),1.0/(rho_sk*c_sk*rhoCFactor*Rtot)* Tair])


  dependentFieldTissue.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
  dependentFieldTissue.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
  #!  !Finish the creation of the equations set boundary conditions
  #!  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  solverEquationsTissue.BoundaryConditionsCreateFinish()

# =================================
# F L O W
if (CoupledBioheatFlow or TestFlow):
  if (ProgressDiagnostics):
    print( " == >> BOUNDARY CONDITIONS << == ")

  if (streeBoundaries):
      # STREE
      BoundaryConditionsStree = iron.BoundaryConditions()
      SolverEquationsStree.BoundaryConditionsCreateStart(BoundaryConditionsStree)
      SolverEquationsStree.BoundaryConditionsCreateFinish()

  #------------------

  # CHARACTERISTIC
  BoundaryConditionsCharacteristic = iron.BoundaryConditions()
  SolverEquationsCharacteristic.BoundaryConditionsCreateStart(BoundaryConditionsCharacteristic)
  # Area-outlet
  for terminalIdx in range (1,numberOfTerminalNodes+1):
      nodeNumber = coupledNodeNumber[terminalIdx-1]
      nodeDomain = Decomposition.NodeDomainGet(nodeNumber,meshComponentNumberSpace)
      if (nodeDomain == computationalNodeNumber):
          BoundaryConditionsCharacteristic.SetNode(DependentFieldNavierStokes,iron.FieldVariableTypes.U,
          versionIdx,derivIdx,nodeNumber,2,OutletBoundaryConditionType,[A[nodeNumber][0]])
  # Finish the creation of boundary conditions
  SolverEquationsCharacteristic.BoundaryConditionsCreateFinish()

  #------------------

  # NAVIER-STOKES
  BoundaryConditionsNavierStokes = iron.BoundaryConditions()
  SolverEquationsNavierStokes.BoundaryConditionsCreateStart(BoundaryConditionsNavierStokes)
  # Flow-inlet
  for inputIdx in range (1,numberOfInputNodes+1):
      nodeNumber = inputNodeNumber[inputIdx-1]
      nodeDomain = Decomposition.NodeDomainGet(nodeNumber,meshComponentNumberSpace)
      if (nodeDomain == computationalNodeNumber):
          BoundaryConditionsNavierStokes.SetNode(DependentFieldNavierStokes,iron.FieldVariableTypes.U,
          versionIdx,derivIdx,nodeNumber,1,InletBoundaryConditionType,[Q[nodeNumber][0]])
  # Area-outlet
  for terminalIdx in range (1,numberOfTerminalNodes+1):
      nodeNumber = coupledNodeNumber[terminalIdx-1]
      nodeDomain = Decomposition.NodeDomainGet(nodeNumber,meshComponentNumberSpace)
      if (nodeDomain == computationalNodeNumber):
          BoundaryConditionsNavierStokes.SetNode(DependentFieldNavierStokes,iron.FieldVariableTypes.U,
          versionIdx,derivIdx,nodeNumber,2,OutletBoundaryConditionType,[A[nodeNumber][0]])
  # Finish the creation of boundary conditions
  SolverEquationsNavierStokes.BoundaryConditionsCreateFinish()

  #------------------

  # ADVECTION
  if (coupledAdvection):
      BoundaryConditionsAdvection = iron.BoundaryConditions()
      SolverEquationsAdvection.BoundaryConditionsCreateStart(BoundaryConditionsAdvection)
      for inputIdx in range (1,numberOfInputNodes+1):
          nodeNumber = inputNodeNumber[inputIdx-1]
          nodeDomain = Decomposition.NodeDomainGet(nodeNumber,meshComponentNumberConc)
          if (nodeDomain == computationalNodeNumber):
              BoundaryConditionsAdvection.SetNode(DependentFieldAdvection,iron.FieldVariableTypes.U,
              versionIdx,derivIdx,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,[1.0])
      SolverEquationsAdvection.BoundaryConditionsCreateFinish()
# =================================

print('\033[1;32m'+'Boundary Conditions COMPLETED'+'\033[0m',"{0:4.2f}".format(time.time()-t))

#================================================================================================================================
#  Element Length
#================================================================================================================================
if (CoupledBioheatFlow or TestFlow):
  if (timestepStability):
      QMax = 600.0
      # Check the element length
      elementNumber = [0]*(totalNumberOfElements+1)
      elementLength = [0]*(totalNumberOfElements+1)
      for i in range(1,totalNumberOfElements+1):
          Node1 = elementNodes[i][0]
          Node2 = elementNodes[i][1]
          Node3 = elementNodes[i][2]
          Length1 = (((xValues[Node1][0]-xValues[Node2][0])**2)
                    +((yValues[Node1][0]-yValues[Node2][0])**2)
                    +((zValues[Node1][0]-zValues[Node2][0])**2))**0.5
          Length2 = (((xValues[Node2][0]-xValues[Node3][0])**2)
                    +((yValues[Node2][0]-yValues[Node3][0])**2)
                    +((zValues[Node2][0]-zValues[Node3][0])**2))**0.5
          elementNumber[i] = i
          elementLength[i] = Length1 + Length2
          elementLength[0] = elementLength[i]
          print( "Element %1.0f" %elementNumber[i],)
          print( "Length: %1.1f" %elementLength[i],)
          print( "Length1: %1.1f" %Length1,)
          print( "Length2: %1.1f" %Length2)
      maxElementLength = max(elementLength)
      minElementLength = min(elementLength)
      print(("Max Element Length: %1.3f" % maxElementLength))
      print(("Min Element Length: %1.3f" % minElementLength))

      # Check the timestep
      for i in range(1,numberOfNodesSpace+1):
          beta   = (3.0*math.sqrt(math.pi)*H[i,0]*E[i,0])/(4.0*A0[i,0])
          eig[i] = QMax/A0[i,0] + (A0[i,0]**0.25)*(math.sqrt(beta/(2.0*Rho)))
          dt[i]  = ((3.0**(0.5))/3.0)*minElementLength/eig[i]
          dt[0]  = dt[i]
      minTimeStep = min(dt)
      print(("Max allowable timestep:      %3.5f" % minTimeStep ))

#================================================================================================================================
#  Transmission Line Theory
#================================================================================================================================
if (CoupledBioheatFlow or TestFlow):
  '''
  terminal 1 ----- right leg ------- right internal iliac artery
  terminal 2 ----- right leg ------- right deep femoral artery
  terminal 3 ----- right leg ------- right posterior tibial artery
  terminal 4 ----- right leg ------- right anterior tibial artery
  terminal 5 ----- left leg -------- left internal iliac artery
  terminal 6 ----- left leg -------- left deep femoral artery
  terminal 7 ----- left leg -------- left posterior tibial artery
  terminal 8 ----- left leg -------- left anterior tibial artery
  terminal 9 ----- right arm ------- right radial artery
  terminal 10 ---- right arm ------- right ulnar artery
  terminal 11 ---- right arm ------- right interosseous artery
  terminal 12 ---- left arm -------- left radial artery
  terminal 13 ---- left arm -------- left ulnar artery
  terminal 14 ---- left arm -------- left interosseous artery
  terminal 15 ---- brain ----------- right external carotid artery
  terminal 16 ---- brain ----------- left external carotid artery
  terminal 17 ---- right kidney ---- right renal artery 1
  terminal 18 ---- right kidney ---- right renal artery 3
  terminal 19 ---- right kidney ---- right renal artery 4
  terminal 20 ---- left kidney ----- left renal artery 1
  terminal 21 ---- left kidney ----- left renal artery 3
  terminal 22 ---- left kidney ----- left renal artery 4
  terminal 23 ---- right lung ------ right intercostal artery 1
  terminal 24 ---- left lung ------- left intercostal artery 1
  terminal 25 ---- right lung ------ right intercostal artery 2
  terminal 26 ---- left lung ------- left intercostal artery 2
  terminal 27 ---- intestines ------ inferior mesenteric artery
  terminal 28 ---- intestines ------ superior mesenteric artery
  terminal 29 ---- intestines ------ middle colic artery
  terminal 30 ---- intestines ------ jejunal 3 artery
  terminal 31 ---- intestines ------ jejunal 5 artery
  terminal 32 ---- intestines ------ ileocolic artery
  terminal 33 ---- intestines ------ ilieal 3 artery
  terminal 34 ---- intestines ------ ilieal 5 artery
  terminal 35 ---- liver ----------- right hepatic artery
  terminal 36 ---- liver ----------- left hepatic artery
  terminal 37 ---- spline ---------- splenic artery
  terminal 38 ---- stomach --------- left gastric artery
  terminal 39 ---- pancreas -------- pancreatic artery
  terminal 40 ---- brain ----------- right posterior communicating artery 2
  terminal 41 ---- brain ----------- left posterior communicating artery 2
  terminal 42 ---- brain ----------- right middle cerebral artery
  terminal 43 ---- brain ----------- left middle cerebral artery
  terminal 44 ---- brain ----------- right anterior cerebral artery 2
  terminal 45 ---- brain ----------- left anterior cerebral artery 2
  '''

  if (streeBoundaries):
      if (ProgressDiagnostics):
          print( " == >> STREE << == ")

      numberOfTerminalNodes = 27
      # Loop through the terminal nodes
      for terminalIdx in range (1,numberOfTerminalNodes+1):
          # Read the organ node file
          with open('input/stree/'+str(terminalIdx)+'.csv','rb') as csvfile:
              reader = csv.reader(csvfile, delimiter=',')
              rownum = 0
              for row in reader:
                  if (rownum == 0):
                      # Read the header row
                      header = row
                  else:
                      # Read number of nodes
                      if (rownum == 1):
                          numberOfSegments = int(row[6])
                          stree = numpy.zeros((numberOfSegments+1,7,timePeriod+1),dtype = numpy.float)
                      stree[rownum][0] = float(row[0])/Ls  # Length of segment
                      stree[rownum][1] = float(row[1])     # Radius of segment
                      stree[rownum][2] = float(row[2])     # Terminal segment
                      stree[rownum][3] = float(row[3])     # Number of parent segment
                      if (row[4]):
                          stree[rownum][4] = float(row[4]) # Number of daughter segments
                          stree[rownum][5] = float(row[5])
                  # Next line
                  rownum+=1

          # Loop through the segments to calculate each segment impedance
          for idx in range(1,numberOfSegments+1):
              n = numberOfSegments+1-idx                    # Start from last segment
              L = stree[n][0][0]                            # Length of segment
              r = stree[n][1][0]                            # Radius of segment
              term = stree[n][2][0]                         # Terminal segment
              if (term == 0):
                  # Calculate daughter segment impedance
                  zin1=stree[stree[n][4][0]][6][0]
                  zin2=stree[stree[n][5][0]][6][0]

              Ng  = 8                                       # Number of generations
              h   = 0.35*r                                  # Thickness
              E   = 0.4E+6                                  # Elasticity
              A0  = math.pi*(r**2.0)                        # Area at rest
              Qin = 6.5E-6                                  # Input flow
              T   = timePeriod/Ts                           # Time period
              zin = [0]*(timePeriod+1)                      # Impedance
              Cp  = (3.0*A0*(A0/math.pi)**0.5)/(2.0*E*h)    # Vessel wall compliance

              # Non-zero frequency condition
              for k in range(0,timePeriod+1):
                  if (k == 0):
                      # Zero frequency condition
                      # Terminal load
                      if (term == 0):
                          zL = zin1*zin2/(zin1+zin2)
                      else:
                          zL = (Pv/Qin)*(2.0**Ng)
                      # Transfer function
                      zin[k] = 8.0*(Mu/Mus)*L/((A0**2.0)/math.pi)+zL
                  else:
                      # Frequency
                      freq = 2.0*math.pi*k/T
                      # Womersley number
                      w = (A0*freq*(Rho/Rhos)/((Mu/Mus)*math.pi))**0.5
                      w0 = ((1j)**1.5)*w
                      # Bessel function zeroth-order
                      J0 = jn(0,w0)
                      # Bessel function first-order
                      J1 = jn(1,w0)
                      # Bessel function
                      Fj = (2.0*J1)/(w0*J0)
                      # Wave propagation velocity
                      c = cmath.sqrt(A0*(1.0-Fj)/((Rho/Rhos)*Cp))
                      g = c*Cp
                      # Terminal load
                      if (term == 1):
                          zL = 1.0/(c*Cp)
                      else:
                          zL = zin1*zin2/(zin1+zin2)
                      # Transfer function
                      zin[k] = ((1j)*cmath.sin(freq*L/c)/g+zL*cmath.cos(freq*L/c))/(cmath.cos(freq*L/c)+(1j)*g*zL*cmath.sin(freq*L/c))
                  #Saving the line's characteristics
                  stree[n][6][k] = zin[k]
          # Invrese fourier transform
          zt = ifft(stree[1][6])*Zs
          # Set the impedance
          for k in range(0,timePeriod+1):
              MaterialsFieldStree.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
              iron.FieldParameterSetTypes.VALUES,1,1,k+1,terminalIdx,zt[k].real)

#================================================================================================================================
#  Run Solvers
#================================================================================================================================

#!  !Solve the problem
#!  CALL cmfe_Problem_Solve(Problem,Err)
print( "Solving problem...")
start = time.time()
# Solve the problem
problem.Solve()
end = time.time()
elapsed = end - start
print( "Total Number of Elements = %d " %totalNumberOfElements)
print( "Total Number of nodes = %d " %numberOfNodes)
print( "Calculation Time = %3.4f" %elapsed)
print( "Problem solved!")
print( "#")
print( "number of muscle Elements = %d\nnumber of left radius elements = %d\nHumerus=%d\nUlna=%d\ntotal number of elements = %d\n"%(
len(muscleElements),len(leftRadiusElements),len(leftHumerusElements),len(leftUlnaElements),len(muscleElements)+
len(leftRadiusElements)+len(leftHumerusElements)+len(leftUlnaElements)))
# print( "number of boundary nodes = %d"%len(boundary))
print( '\033[1;31m'+"Elapsed time: "+'\033[0m', time.time()-t)
#!# Export results
#!baseName = "laplace"
#!dataFormat = "PLAIN_TEXT"
#!fml = iron.FieldMLIO()
#!fml.OutputCreate(mesh, "", baseName, dataFormat)
#!fml.OutputAddFieldNoType(baseName+".geometric", dataFormat, geometricField,
#!    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
#!fml.OutputAddFieldNoType(baseName+".phi", dataFormat, dependentField,
#!    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
#!fml.OutputWrite("LaplaceExample.xml")
#!fml.Finalise()


#!  !Output Analytic analysis
#!!  CALL cmfe_AnalyticAnalysis_Output(DependentField,"DiffusionAnalytics_x10_y10_z10_L_T1",Err)
#exportField = True
#if (exportField):
#    fields = iron.Fields()
#    fields.Create(Region)
#    fields.Finalise()
#!  EXPORT_FIELD=.TRUE.
#!  IF(EXPORT_FIELD) THEN
#!    CALL cmfe_Fields_Initialise(Fields,Err)
#!    CALL cmfe_Fields_Create(Region,Fields,Err)
#!!    CALL cmfe_Fields_NodesExport(Fields,"DiffusionConstantSourceAnalytic_x10_y10_z10_L_T1","FORTRAN",Err)
#!!    CALL cmfe_Fields_ElementsExport(Fields,"DiffusionConstantSourceAnalytic_x10_y10_z10_L_T1","FORTRAN",Err)
#!    CALL cmfe_Fields_Finalise(Fields,Err)

#!  ENDIF

#!  !CALL cmfe_Finalise(Err)
#!  WRITE(*,'(A)') "Program successfully completed."

#!  STOP
iron.Finalise()
#!END PROGRAM DIFFUSIONCONSTANTSOURCEEXAMPLE


