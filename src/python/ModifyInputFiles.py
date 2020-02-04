# Created by Elias Ghadam Soltani
# This code is written to convert exfile to vtk file and ply to smesh. I use this approach to get the surface mesh from ParaView
# In this version of the code it just produces conectivities and positions on the screen so you need to copy and paste them
# into a vtk template file.
# Don't forget that the order of node numbering for an element in vtk and exfile is different. See Google slides for that

import numpy,time,glob,os,sys,subprocess, pprint
from shutil import copyfile

# ===========================================================================
# making colored prints in console
# ===========================================================================
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# examples
print bcolors.HEADER + "Hi Elias, What color do you see ;)?" + bcolors.ENDC
print bcolors.OKBLUE + "Warning: No active frommets remain. Continue?" + bcolors.ENDC
print bcolors.OKGREEN + "Hi Elias, What color do you see ;)?" + bcolors.ENDC
print bcolors.WARNING + "Warning: No active frommets remain. Continue?" + bcolors.ENDC
print bcolors.FAIL + "Warning: No active frommets remain. Continue?" + bcolors.ENDC
print bcolors.BOLD + "Warning: No active frommets remain. Continue?" + bcolors.ENDC
print bcolors.UNDERLINE + "Warning: No active frommets remain. Continue?" + bcolors.ENDC

def print_format_table():
    """
    prints table of formatted text format options
    """
    for style in range(8):
        for fg in range(30,38):
            s1 = ''
            for bg in range(40,48):
                format = ';'.join([str(style), str(fg), str(bg)])
                s1 += '\x1b[%sm %s \x1b[0m' % (format, format)
            print(s1)
        print('\n')

print_format_table()

print("\x1b[6;36;40m salam \x1b[7;31;47m Elias \x1b[7;32;43m joon \x1b[0m")


# ===========================================================================
# Convert some ply files to smesh
# ===========================================================================
def ply2smesh(FileName,RegionLabel,numberOfLocalNodes,offset=0):
  localNodes = [0]*numberOfLocalNodes
  NodeBegin = False
  nodeNumber = 0
  faceNumber = 0
  FaceBegin = False
  numberOfNodes = -1
  numberOfFaces = -1


  with open(FileName,'r') as f:
    target=f.readlines()
    for lineNum,line in enumerate(target):
      target[lineNum] = target[lineNum].rstrip("\n\r").replace("\t"," ").replace(","," ").split()
    for lineNum,line in enumerate(target):
      if target[lineNum][0] == "element":
        if (target[lineNum][1]=="vertex"):
          numberOfNodes = int(target[lineNum][2])
        elif ((target[lineNum][1]=="face")):
          numberOfFaces = int(target[lineNum][2])
          print numberOfFaces,numberOfNodes
      if NodeBegin:
        x = float(target[lineNum][0])
        y = float(target[lineNum][1])
        z = float(target[lineNum][2])
        print nodeNumber+offset, x, y,z, RegionLabel
        nodeNumber = nodeNumber+1
      if target[lineNum][0] == "end_header":
        NodeBegin = True

      if FaceBegin:
        numberOfLocalNodes = int(target[lineNum][0])
        string=""
        for i in range(numberOfLocalNodes):
          localNodes[i] = int(target[lineNum][i+1])+offset
          string=string+str("localNodes[%d]," %i)
        string="print numberOfLocalNodes,"+string+"RegionLabel"
        exec(string)
#        print numberOfLocalNodes, localNodes[0],localNodes[1],localNodes[2], RegionLabel
        faceNumber = faceNumber+1
      if nodeNumber == numberOfNodes:
        NodeBegin = False
        FaceBegin = True
      if faceNumber == numberOfFaces:
        FaceBegin = False 

#numberOfRegions = 4      
#totalNumberOfNodes = 256
#totalNumberOfFaces = 510
#FileName = "deleteIt.smesh"

# ===========================================================================
# Creating template smesh textfile
# ===========================================================================
def TemplateSmesh(FileName,totalNumberOfNodes,totalNumberOfFaces,numberOfRegions):
  with open(FileName, 'w') as f:
    f.write('# Part 1 - the node list.\n')
    f.write('# The model has %d nodes in 3d, no attributes, with boundary marker.\n' % totalNumberOfNodes)
    f.write('%d  %d  %d  %d\n' %(totalNumberOfNodes,3,0,1))
    f.write('# Skin point clouds\n')
    f.write('%d %f %f %f %d\n' % (29,12.0,2.6,3.5,2))
    f.write('# Left Radius points cloud\n')
    f.write('%d %f %f %f %d\n' % (29,12.0,2.6,3.5,3))
    f.write('# Left Humerus points cloud\n')
    f.write('%d %f %f %f %d\n' % (29,12.0,2.6,3.5,4))
    f.write('# Left Ulna points cloud\n')
    f.write('%d %f %f %f %d\n' % (29,12.0,2.6,3.5,5))
    f.write('# Part 2 - the facet list.\n')
    f.write('# %d facets with boundary markers.\n' % totalNumberOfFaces)
    f.write('%d  %d\n' %(totalNumberOfFaces,1))
    f.write('# Skin facets\n')
    f.write('%d	%d	%d	%d	%d	%d\n' %(4,0,1,2,3,2))
    f.write('# Left Radius facets\n')
    f.write('%d	%d	%d	%d	%d\n' %(3,0,1,2,3))
    f.write('# Left Humerus facets\n')
    f.write('%d	%d	%d	%d	%d\n' %(3,0,1,2,4))
    f.write('# Left Ulna facets\n')
    f.write('%d	%d	%d	%d	%d\n' %(3,0,1,2,5))
    f.write('# Part 3 - the hole list.\n')
    f.write('# There is no hole in regions.\n')
    f.write('%d\n' %0)
    f.write('# Part 4 - the region list.\n')
    f.write('# There are %d regions defined.\n' %numberOfRegions)
    f.write('%d\n' %numberOfRegions)
    f.write('  1 -206.152 -86.7753 978.752 -10 # muscle\n')
    f.write('  2 -251.338 -87.4435 922.511 -20 # bone\n')
    f.write('  3 -251.338 -87.4435 922.511 -20 # bone\n')
    f.write('  4 -251.338 -87.4435 922.511 -20 # bone\n')


# ===========================================================================
# Finds nodes with boundary marker 1 in the final mesh.node for the whole arm
# It uses mesh.face to correct the boundary marker to one of the 2 to 5 numbers according to its region.
# ===========================================================================
def correctBoundaryMarkers(FileName_node,FileName_face):
  t=time.time()
  nodeList = []
  with open(FileName_node, "r") as f:
      target=f.readlines()
      for lineNum,line in enumerate(target):
          target[lineNum] = target[lineNum].rstrip('\n\r').replace("\t"," ").replace(","," ").split()
      for lineNum,line in enumerate(target):
          if (lineNum!=0) and (lineNum!=len(target)-1):
#            nodeNumber=0
            if (int(target[lineNum][4])==1):
              nodeNumber=int(target[lineNum][0])
              nodeList.append(nodeNumber)
#  print nodeList
#  raw_input("stop")
  counter = 0
  ListBoundary = []
  markerList = []
  with open(FileName_face, "r") as f:
      target2=f.readlines()
      for lineNum,line in enumerate(target2):
          target2[lineNum] = target2[lineNum].rstrip('\n\r').replace("\t"," ").replace(","," ").split()
      lastline = len(target2)
      for lineNum,line in enumerate(target2):
          if (lineNum!=0) and (lineNum!=lastline-1):
              node1 = int(target2[lineNum][1])
              node2 = int(target2[lineNum][2])
              node3 = int(target2[lineNum][3])            
              if ((node1 in nodeList) and (node1 not in ListBoundary)):
                  marker = int(target2[lineNum][4])
                  markerList.append(marker)
                  ListBoundary.append(node1)
                  counter=counter+1
#                  print counter,node1,marker
              if ((node2 in nodeList) and (node2 not in ListBoundary)):
                  marker = int(target2[lineNum][4])
                  markerList.append(marker)
                  ListBoundary.append(node2)
                  counter=counter+1
#                  print counter,node2,marker
              if ((node3 in nodeList) and (node3 not in ListBoundary)):
                  marker = int(target2[lineNum][4])
                  markerList.append(marker)
                  ListBoundary.append(node3)
                  counter=counter+1
#                  print counter,node3,marker
  print counter, ListBoundary[-1],markerList[-1],ListBoundary[-2],markerList[-2]
#  with open(FileName_node)
  with open("correctedBoundaryMarker.node",'w') as f:
    for nodeMarkerIncorrect in range(len(ListBoundary)):
      node = ListBoundary[nodeMarkerIncorrect]
      target[node+1][4]=markerList[nodeMarkerIncorrect]
      print node,markerList[nodeMarkerIncorrect],target[node+1][4],len(target),len(markerList)
    for lineNum,line in enumerate(target):
      if lineNum ==0:
        f.write('%d %d %d %d\n' %(int(target[0][0]),int(target[0][1]),int(target[0][2]),int(target[0][3])))
      elif (lineNum !=0) and (lineNum != len(target)-1):
#        print lineNum,target[lineNum][4]
        f.write('%d %f %f %f %d\n' %(int(target[lineNum][0]),float(target[lineNum][1]),float(target[lineNum][2]),float(target[lineNum][3]),int(target[lineNum][4])))





# ===========================================================================        
# converts exfiles to vtk files,
# inputs = exnodeFile, VTKfile before information of exnodefile. 
# ===========================================================================
def exfilesToVTK(exnodeFile,VTKfile,VTKfileCopy,numberOfNodes):
  copyfile(VTKfile, VTKfileCopy)
  with open(exnodeFile, 'r') as f, open(VTKfileCopy, 'a') as g:
      g.write('POINT_DATA %d\n' %numberOfNodes)
      g.write('SCALARS Temperature float 1\n')
      g.write('LOOKUP_TABLE default\n')
      target=f.readlines()
      for lineNum,line in enumerate(target):
          target[lineNum] = target[lineNum].rstrip('\n\r').replace("\t"," ").replace(","," ").split()
      for lineNum,line in enumerate(target):
        if target[lineNum][0]=="Node:":
          temperature = float(target[lineNum+4][0])
          g.write('%f\n'%temperature)




# ===========================================================================
# It gets exfiles as input and creates VTK file from scratch.
# ===========================================================================
# This function creates a template vtk file from the exnode and exelem file.
# Inputs:
#     filePath. This is the folder that the exfiles exist. Does not need to be an absolute address.'
#     numberOfProcessor. 
# Output:
#     A vtk file
# Author: Elias Ghadam Soltani. Aug 2018.
def createTemplate(exnodeFilePath,exelemFilePath,VTKfilePath, numberOfExnodePart=1,numberOfExelemPart=1): 

  #
  #====================== CONSTANTS ========================
  #

  TETRAHEDRA_TYPE = 10 #< Tetrahedra type of the cell/element. \see VTK guide file.
  HEXAHEDRA_OPENCMISS_TYPE = 12 #< Hexahedra type of the cell/element. \see VTK guide file.

  #
  #====================== HELP ======================
  #

  print('\033[96m'+'\nThis code only works for:\n')
  print(' * Tets mesh not mixed mesh')
  print(' * 3 dimension')
  print(' * You need to change INPUT STRUCTURE OF DATA manually')

  #
  #====================== INPUTING STRUCTURE OF DATA ======================
  #

  # Enter only the fields that you are interested to get the values
  #  cellScalars=[(dataName,dataType,numComp,exfileLine)] #the line number for the first component
  dimension = 3 # e.g. 3
  cellType= TETRAHEDRA_TYPE 
  points=[('points','float',3,1),('U','float',1,4)]
  cells=[]
  cellScalars=[('Diffusivity','float',3,0),('Source','float',1,3)]
  pointScalars=[]
  #
  # ===================== GETTING CELL DATA ===============
  #

# creating elem map={elemNumber:(k1,k2,k3,source)} and cellsMap{elemNumber: nodesList}
  elemMap={}
  cellsMap={}

  # Looping over part files.
  for partNumber in range(numberOfExelemPart):
    exelemFile = exelemFilePath+'/MainTime_0.part%d.exelem' %partNumber 
    # Storing data in elemData as a list
    with open(exelemFile, 'r') as fe:
      elemData=fe.readlines()
      for lineNum,line in enumerate(elemData):
          elemData[lineNum] = elemData[lineNum].rstrip('\n\r').replace("\t"," ").replace(","," ").split()
   
  # ===========================    
      # going over the file line by line and extracting data.
      for lineNum, line in enumerate(elemData):
        # Extracting cell scalar data. e.g. diffusivity and source. elemNumber, 
        if elemData[lineNum][0]=="Element:":
          elemNumber = int(elemData[lineNum][1])
          elemList=[]
          for scalar in cellScalars:
            for compNum in range(scalar[2]):
              if scalar[1]=='float':
                elemList.append(float(elemData[lineNum+2+(scalar[3]+compNum)][0]))  # line+2 because we want to skip Values: line
              elif scalar[1]=='int':
                elemList.append(int(elemData[lineNum+2+(scalar[3]+compNum)][0]))  # line+2 because we want to skip Values: line                  

              else:
                raise Exception('This variable type is not implemented')
          elemMap[elemNumber]=elemList
  #=============================
        # Getting nodes on the elements and store it in cellsMap.
        if elemData[lineNum][0]=="Nodes:":
          nodeList=[]
          if cellType==TETRAHEDRA_TYPE:          # celltype 10 means linear tets. See the vtk guide file. 
            for i in range(4):
              nodeList.append(int(elemData[lineNum+1][i])-1)     #lineNum+1 because it is in the next line. 
          elif cellType==HEXAHEDRA_OPENCMISS_TYPE:
            for i in range(8):
              nodeList.append(int(elemData[lineNum+1][i])-1)     #lineNum+1 because it is in the next line.
            SwapList(nodeList)
          else:    
            raise ValueError(' cell Type is not implemented')
 
          cellsMap[elemNumber] = nodeList

      del elemData[:]
  numberOfElements = len(elemMap)
#  print(numberOfElements,elemMap[424])

  #
  # ===================== GETTING POINT DATA ===============
  #

# creating node map={nodeNumber:(temperature,)} and pointsMap={nodeNumber:(x,y,z)}
  nodeMap={}
  pointsMap={}

  # Looping over part files.
  for partNumber in range(numberOfExnodePart):
    exnodeFile = exnodeFilePath+'/MainTime_0.part%d.exnode' %partNumber 
    # Storing data in nodeData as a list
    with open(exnodeFile, 'r') as fn:
      nodeData=fn.readlines()
      for lineNum,line in enumerate(nodeData):
          nodeData[lineNum] = nodeData[lineNum].rstrip('\n\r').replace("\t"," ").replace(","," ").split()
   
  # =============================    
      # going over the file line by line and extracting data.
      for lineNum, line in enumerate(nodeData):
        # Extracting point scalar data. e.g. temperature or source. nodeNumber, and position 
        if nodeData[lineNum][0]=="Node:":

          nodeNumber = int(nodeData[lineNum][1])

          nodeList=[]
          for scalar in pointScalars:
            for compNum in range(scalar[2]):
              if scalar[1]=='float':
                nodeList.append(float(nodeData[lineNum+1+(scalar[3]+compNum)][0]))  # line+1 because we want to skip Node: line
              elif scalar[1]=='int':
                nodeList.append(int(nodeData[lineNum+2+(scalar[3]+compNum)][0]))  # line+2 because we want to skip Values: line                  
              else:
                raise Exception('This variable type is not implemented')
          nodeMap[nodeNumber]=nodeList
  #===============================
        # Getting position of the nodes and store it in pointsMap.
        if nodeData[lineNum][0]=="Node:":
          positionList=[]
          if dimension==3:           
            for i in range(3):
              positionList.append(float(nodeData[lineNum+1+i][0]))     #lineNum+1 because it is in the next line. 
          else:    
            raise ValueError(' only 3 dimension has been implemented')
 
          pointsMap[nodeNumber] = positionList

      del nodeData[:]
  numberOfNodes = len(nodeMap)
#  print(numberOfNodes,pointsMap[2132])



  #
  # ===================== OUTPU VTK TEMPLATE FILE ===============
  #
  #workSpace(vars())
  # writing vtk template file
  VTKfile = VTKfilePath+'/vtk/template.vtk'

  #===============================
  # Writing the positions as x y z in the VTK file.
  with open(VTKfile, 'w') as fv:
    fv.write("# vtk DataFile Version 2.0\nUnstructured Grid\nASCII\nDATASET UNSTRUCTURED_GRID\n")
    fv.write("POINTS %d double\n" %numberOfNodes)  

    if dimension == 3:
      for nodeNumber in range(1,numberOfNodes+1):
        x = pointsMap[nodeNumber][0]
        y = pointsMap[nodeNumber][1]
        z = pointsMap[nodeNumber][2]
        fv.write(" %f %f %f\n" %(x,y,z))
#        if nodeNumber ==30113:
#          print(pointsMap[nodeNumber],x,y,z)
    else:
      raise ValueError(' only 3 dimension has been implemented')
    pointsMap.clear()

  # ==============================
    # Writing the elements and their nodes as 4 10 120 23 50.
    if cellType==TETRAHEDRA_TYPE:
      fv.write("\nCELLS %d %d\n"%(numberOfElements,numberOfElements*5))
    elif cellType==HEXAHEDRA_OPENCMISS_TYPE:
      fv.write("\nCELLS %d %d\n"%(numberOfElements,numberOfElements*9))
    else:
      raise ValueError('Only cell Type of Tets is implemented')

    for i in range(1,numberOfElements+1): #element = [10 20 150 120,...]
      element=cellsMap[i]
      elemLen=len(element)
      fv.write("%d "%elemLen)
      for node in element:
        fv.write("%d "%node)
      fv.write("\n")        
    cellsMap.clear()
  # ==============================
    # Writing the cell types as 10. 
    fv.write("\nCELL_TYPES %d\n" %numberOfElements)
    for i in range(numberOfElements):
      if cellType == TETRAHEDRA_TYPE:
        fv.write("10\n")
      elif cellType == HEXAHEDRA_OPENCMISS_TYPE:
        fv.write("12\n")      
      else:
        raise ValueError('Only cell Type of Tets is implemented')
  # ==============================
    # Writing Cell data information like scalars or vectors,...
    fv.write("\nCELL_DATA %d\n" %numberOfElements)
    # Writing the scalars one by one.
    counterInitial=0
    for scalar in cellScalars:
      fv.write("\nSCALARS %s %s %d\nLOOKUP_TABLE default\n" %(scalar[0],scalar[1],scalar[2]))
      # For the current scalar we go element by element and write the values
      for elemNumber in range(1,numberOfElements+1): 
        elemList=elemMap[elemNumber]
        counter=counterInitial
        for compNum in range(scalar[2]):
          if scalar[1]=='float':
            fv.write("%f\n"%elemList[counter])
          elif scalar[1]=='int':
            fv.write("%d\n"%elemList[counter])
          else:
            raise ValueError('Invalid input') 
          counter+=1
      counterInitial=counter
    elemMap.clear()
  # ==============================
    # Writing Point data information like scalars or vectors,...
    fv.write("\nPOINT_DATA %d\n" %numberOfNodes)
    # Writing the scalars one by one.
    counterInitial=0
    for scalar in pointScalars:
      fv.write("\nSCALARS %s %s %d\nLOOKUP_TABLE default\n" %(scalar[0],scalar[1],scalar[2]))
      # For the current scalar we go point by point and write the values
      for nodeNumber in range(1,numberOfNodes+1): 
        nodeList=nodeMap[nodeNumber]
        counter=counterInitial
        for compNum in range(scalar[2]):
          if scalar[1]=='float':
            fv.write("%f\n"%nodeList[counter])
          elif scalar[1]=='int':
            fv.write("%d\n"%nodeList[counter])
          else:
            raise ValueError('Invalid input') 
          counter+=1
      counterInitial=counter
    nodeMap.clear()




#
# ===========================================================================
# It gets exnodes and vtk template file as input and add temperature data to vtk template file
# ===========================================================================
#
def AddPointData(timeStepNumber,vtkFileTemplate,rootDirectory, exnodeFilePath='./output/output',numberOfExnodePart=1):
  #
  #====================== INPUTING STRUCTURE OF DATA ======================
  #

  # Enter only the fields that you are interested to get the values
  #  cellScalars=[(dataName,dataType,numComp,exfileLine)] #the line number for the first component
  TETRAHEDRA_TYPE = 10
  dimension = 3 # e.g. 3
  cellType= TETRAHEDRA_TYPE # only this type is supported for now
  pointScalars=[('Temperature','float',1,3)]

  #
  # ===================== GETTING POINT DATA ===============
  #

# creating node map={nodeNumber:(temperature,)}
  nodeMap={}

  # Looping over part files.
  for partNumber in range(numberOfExnodePart):
    exnodeFile = exnodeFilePath+'/MainTime_%d.part%d.exnode' %(timeStepNumber,partNumber)  
    # Storing data in nodeData as a list
    with open(exnodeFile, 'r') as fn:
      nodeData=fn.readlines()
      for lineNum,line in enumerate(nodeData):
          nodeData[lineNum] = nodeData[lineNum].rstrip('\n\r').replace("\t"," ").replace(","," ").split()
  # =============================    
      # going over the file line by line and extracting data.
      for lineNum, line in enumerate(nodeData):
        # Extracting point scalar data. e.g. temperature or source. nodeNumber, and position 
        if nodeData[lineNum][0]=="Node:":

          nodeNumber = int(nodeData[lineNum][1])

          nodeList=[]
          for scalar in pointScalars:
            for compNum in range(scalar[2]):
              if scalar[1]=='float':
                nodeList.append(float(nodeData[lineNum+1+(scalar[3]+compNum)][0]))  # line+1 because we want to skip Node: line
              elif scalar[1]=='int':
                nodeList.append(int(nodeData[lineNum+2+(scalar[3]+compNum)][0]))  # line+2 because we want to skip Values: line                  
              else:
                raise Exception('This variable type is not implemented')
          nodeMap[nodeNumber]=nodeList
  #===============================

      del nodeData[:]
  numberOfNodes = len(nodeMap)

  #===============================
  # writing into files one by one
  filePath=rootDirectory
  VTKfile = filePath+'/vtk'+'/MainTime_%d.vtk' %(timeStepNumber)
  copyfile(vtkFileTemplate,VTKfile)
  print(VTKfile)

  with open(VTKfile, 'a') as fv:
    # Writing the scalars one by one.
    counterInitial=0
    for scalar in pointScalars:
      fv.write("\nSCALARS %s %s %d\nLOOKUP_TABLE default\n" %(scalar[0],scalar[1],scalar[2]))
      # For the current scalar we go point by point and write the values
      for nodeNumber in range(1,numberOfNodes+1): 
        nodeList=nodeMap[nodeNumber]
        counter=counterInitial
        for compNum in range(scalar[2]):
          if scalar[1]=='float':
            fv.write("%f\n"%nodeList[counter])
          elif scalar[1]=='int':
            fv.write("%d\n"%nodeList[counter])
          else:
            raise ValueError('Invalid input') 
          counter+=1
      counterInitial=counter
    nodeMap.clear()



#    fv.write("SCALARS Temperature float 1\nLOOKUP_TABLE default\n")
#    for nodeNumber in range(1,numberOfNodes+1):
#        temperature = nodeMap[nodeNumber]
#        fv.write("%f\n"%(temperature))
#  workSpace(vars())
#  input('stop')
#=============================
#  with open("./MainTime_14300.part15.exnode", 'r') as fn, open(storeTemp, 'w') as fs: 
#    nodeData=fn.readlines()
#    for lineNum,line in enumerate(nodeData):
#        nodeData[lineNum] = nodeData[lineNum].rstrip('\n\r').replace("\t"," ").replace(","," ").split()
#    for lineNum, line in enumerate(nodeData):
#      if nodeData[lineNum][0]=="Node:":
#        numberOfValues =int(nodeData[lineNum-1][3])
#    fv.write("# vtk DataFile Version 2.0\nUnstructured Grid\nASCII\nDATASET UNSTRUCTURED_GRID\n")
#    fv.write("POINTS %d double\n" %numberOfNodes)            
#    for lineNum, line in enumerate(nodeData):
#      if nodeData[lineNum][0]=="Node:":
#        fv.write("%f %f %f\n" %(float(nodeData[lineNum+1][0]),float(nodeData[lineNum+2][0]),float(nodeData[lineNum+3][0])))
#    
#del lst[:]
#  subprocess.call("sed -n -e $((17+7*i))p ./MainTime_14300.part15.exnode >> tem",shell=True)



#  with open(vtkFile, 'a') as fv:
#    fv.write("\nSCALARS Source float 1\nLOOKUP_TABLE default\n")
#    temperature=10
#    fv.write("%f\n"%temperature)
#    for i in range(1,5):
#      fv.write("S%d\n"%i)






   
# ===========================================================================
# GENERATING MATLAB-LIKE WORKSPACE
# ===========================================================================
#
# This function gets vars() dictionary as input and can show either the list and dictinaries or all variables.
# You can ask for any variable value by calling its string and either print one of the items in the list or all of them.
# You can also write the list or dictionary in a file with the same name.var.
def workSpace(vari,variableType='All'):
  #==========================
  # printing the lists and dictionaries information. Values are not displayed
  print('\033[1;32m'+'Welcome to Work Space\n')
  for i in range(len(vari)):
    value=vari.values()[i]
    if type(value)==list or type(value)==dict:
      print(vari.keys()[i],type(value),len(value),str(sys.getsizeof(value))+' Bytes')
  #==========================
  # print other types of variables names and values.
  if variableType=='All':
    for i in range(len(vari)):
      value=vari.values()[i]
      if type(value) not in [list,dict]:
        print(vari.keys()[i],vari.values()[i])
  #==========================
  # Printing the value of the user input variable name interactively
  Var='Yes'
  while Var is not 'None':
    try:
      Var=input('\033[35m'+'Enter the variable name. e.g. "dimension". Type "None" to continue : ')
      if Var not in vari.keys() and not 'None':
        raise ValueError 
      if Var=='None':
        raise Exception('Thank you for using the work space ;)')
      VarType=type(vari[Var])
      if VarType==list or VarType==dict:
        fileWrite = input(' Would you like to write the variable in a text file? yes? : ')
        if fileWrite == 'yes':
          with open(Var+'.var','w') as fva:
            fva.write(str(vari[Var]))  
        else:
          index=input('Enter the index you want. If you want all indecis type "all" : ')
          if index=='all':
            print(vari[Var])
          else:
            print(vari[Var][index])
      else:
        print(vari[Var])
    except ValueError:
      print('INVALID INPUT. TRY AGAIN.')  
    except NameError:
      print('The variable name only can be a string type')
    except Exception:
      print('Thank you for using the Work Space ;)'+'\033[96m')

# upgrade:
#   1. If the user asks write the variable in a file





   
# ===========================================================================
# Swap the items of a list  to create a new list
# ===========================================================================
#
# This functions swap the items inside the list and change their order. This has been done originaly to change order 
# of the items from exfile format to vtk format. the vtk reader use the order for generating the connectivity between nodes.
# Inputs: List, 
# Output: list with different order 
def SwapList(List):
  tempList = List[:]
  List[2] = tempList[5] 
  List[3] = tempList[4] 
  List[4] = tempList[2] 
  List[5] = tempList[3] 
  List[6] = tempList[7] 
  List[7] = tempList[6]   

      
# ===========================================================================


#________________________ End of Functions _________________________
#=====================================================================


#RegionLabel = 2
#offset = 0
#FileName = 'fifth layer/DenseSkinSurf.ply'
#numberOfLocalNodes = 3
#ply2smesh(FileName,RegionLabel,numberOfLocalNodes,offset)


#TemplateSmesh("deleteIt.smesh",510,256,4)

#FileName_node = './Cylinder/finerMesh/layer2/cylinderModified1WithBoundaryMarker.1.node'
#FileName_face = './Cylinder/finerMesh/layer2/cylinderModified1WithBoundaryMarker.1.face'
#FileNameCorrected_node = '/Cylinder/CorrectedBoundary.node'
#correctBoundaryMarkers(FileName_node,FileName_face)

#VTKfile = 'WholeLeftArm.1.vtk'
#for i in range(0,3200,100):
#  VTKfileCopy    = 'MainTime_%d.vtk' %i
#  exnodeFile = 'MainTime_%d.part0.exnode' %i
#  exfilesToVTK(exnodeFile,VTKfile,VTKfileCopy,40813)




rootDirectory = 'outputDiffusion'
# get the number of node parts (if we used 4 processors we have 4 exfile parts).
nExnodeParts = int(sys.argv[1])
nExelemParts = int(sys.argv[2])
# if the exfile is passed use the user paths.
if len(sys.argv)>3:
  exnodePath = str(sys.argv[3])
  exelemPath = str(sys.argv[4])
else:
  exnodePath = rootDirectory+'/output'
  exelemPath = rootDirectory+'/exelem'

# create a new folder for vtk files and template file that is used for creating vtk files
if not os.path.exists(rootDirectory+'/vtk'):
  os.makedirs(rootDirectory+'/vtk')
if not os.path.exists(rootDirectory+'/vtk/template.vtk'):
  createTemplate(exnodePath,exelemPath,rootDirectory, numberOfExnodePart=nExnodeParts,numberOfExelemPart=nExelemParts)
else:
  print(rootDirectory+"/vtk/template.vtk already exists")

vtkFileTemplate = rootDirectory+'/vtk/template.vtk'
# Find all the exnode files
exnodeFileList= glob.glob(exnodePath+'/*.part0.exnode')
timeStep=[]

# loop over all exnode files and extract and store the timesteps in a list.
#print(exnodeFileList)
for exnodeFile in exnodeFileList:
  timeStep.append(int(filter(str.isdigit, exnodeFile.split("/")[-1].split(".")[0])))
timeStep.sort()

# Create the vtk files from template file and for all timesteps
#print(timeStep)
for timeStepNumber in timeStep:
  AddPointData(timeStepNumber,vtkFileTemplate,rootDirectory,exnodePath,numberOfExnodePart=nExnodeParts)

# Remove the temporary created template file
os.remove(vtkFileTemplate)
