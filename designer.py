''' 
TODO 

- debug comments
- automate Bonds
- delete "mm", instead, select the units in HFSS
- define in classes
class CircuitElement(object):
	def __init__(self, gates, params):
	
	def get_gates(self):
		returns list of gates (inputs and outpts)

class Capacitor(CircuitElement):
	def draw(self):
		draws the element
		
- LitExp: check if variable exists, updates value if exists, create a new one if not.
- Create drawing script which can start from blank file
- Premesh 

'''

import ScriptEnv

ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
oProject = oDesktop.SetActiveProject("catch_tomo_Q7")
oDesign = oProject.SetActiveDesign("HFSSDesign1")
oEditor = oDesign.SetActiveEditor("3D Modeler")
oModule = oDesign.GetModule("BoundarySetup")
oMesh = oDesign.GetModule("MeshSetup")

'''
Throughout the code, a single frame is used: the default one
+----------+
|	       |
|		   |
^          |
y		   |
+x> -------+
'''

TOP = [0, 1]
DOWN = [0, -1]
RIGHT = [1, 0]
LEFT = [-1, 0]

POS = 0
ORI = 1
TRACK = 2
GAP = 3

trackObjects, gapObjects, bondwireObjects = [], [], []
'''
trackObjects (list of strings) e.g. ['readout_track', 'transmon_pad']
gapObjects (list of strings) e.g. ['readout_box', 'transmon_box']
One starts from a ground plane
1. Substract the gapObject
2. Add the trackObject

throught the code, one populates these two lists with objects, and at the very end
apply the above points 1 and 2
'''

    # def set_variable(self, name, value, postprocessing=False):
        #TODO: check if variable does not exist and quit if it doesn't?
        # if name not in self._design.GetVariables()+self._design.GetPostProcessingVariables():
            # self.create_variable(name, value, postprocessing=postprocessing)
        # else:
            # self._design.SetVariableValue(name, value)
        # return VariableString(name)

    # def get_variable_value(self, name):
        # return self._design.GetVariableValue(name)
        
    # def get_variable_names(self):
        # return [VariableString(s) for s in self._design.GetVariables()+self._design.GetPostProcessingVariables()]  
		
		

		
		
		
def assigneInduc(iName, iIn, iOut, iVal):
	"""
	Inputs:
	-------
	iName (string): name of the object e.g. 'transmon_junction'
	iIn (tuple): (x0, y0) first point of the current line
	iOut (tuple): (x1, y1) last point of the current line
	iVal (float): inductance value in nH e.g. 7

	Returns:
	--------
	takes existing iName object and assigns boundary condition lumped RLC, L=iVal nH, R=0, C=0, along current_line 
	"""
	oModule.AssignLumpedRLC(
		[
			"NAME:LumpRLC_"+iName,
			"Objects:="		, [iName],
			[
				"NAME:CurrentLine",
				"Start:="		, [str(float(iIn.x))+"mm",str(float(iIn.y))+"mm","0mm"],
				"End:="			, [str(float(iOut.x))+"mm",str(float(iOut.y))+"mm","0mm"]
			],
			"UseResist:="		, False,
			"UseInduct:="		, True,
			"Inductance:="		, "("+str(iVal)+")nH",
			"UseCap:="		, False
		])
	

def assigne50Ohm(iName, iIn, iOut):
	"""
	Inputs:
	-------
	iName (string): name of the object e.g. 'transmon_junction'
	iIn (tuple): (x0, y0) first point of the current line
	iOut (tuple): (x1, y1) last point of the current line
	iVal (float): inductance value in nH e.g. 7

	Returns:
	--------
	takes existing iName object and assigns boundary condition lumped RLC, L=iVal nH, R=0, C=0, along current_line 
	"""
	oModule.AssignLumpedRLC(
		[
			"NAME:50Ohm_"+iName,
			"Objects:="		, [iName],
			[
				"NAME:CurrentLine",
				"Start:="		, [str(float(iIn.x))+"mm",str(float(iIn.y))+"mm","0mm"],
				"End:="			, [str(float(iOut.x))+"mm",str(float(iOut.y))+"mm","0mm"]
			],
			"UseResist:="		, True,
			"UseInduct:="		, False,
			"Resistance:="		, "50Ohm",
			"UseCap:="		, False
		])	
		

def assignMeshLength(iName,iVal):
	"""
	Inputs:
	-------
	iName (string): name of the object e.g. 'transmon_junction'
	iVal (float): max length value in mm e.g. 7

	Returns:
	--------
	takes existing iName object and assigns mesh minimum length l=iVal mm
	"""
	oMesh.AssignLengthOp(
		[
			"NAME:Mesh_"+iName,
			"RefineInside:=", False, 
			"Objects:=", [iName],
			"RestrictElem:=", False, 
			"RestrictLength:=",  True,
			"MaxLength:=", "("+str(iVal)+")mm"
		])	

def filet(iObject, iFiletRay, iIndex):
	'''
	Inputs:
	-------
	iObject (string): name of the polygone object e.g. 'transmon_pad'
	iFiletRay (float): filet radius in mm
	iIndex (int): vertex index

	Outputs:
	--------
	takes existing iObject, applies a filet of radius iFiletRay on vertex iIndex

	this code employs 6 polygone vertices
	0--5
	|  |
	|  4--3
	|     |  
	1-----2
	One would want to filet vertices 1 and 4

	1. Run this function by setting iIndex = 1. This will filet vertex 1
	HFSS then deletes this vertex and replaces it by two vertices left and right
	this shifts the labeling. Previously labeled vertex 4 is now labeled 5.
	0--6
	|  |
	|  5--4
	1     |  
	\2----3
	2. Run this function by setting iIndex = 5
	'''
	vertices = oEditor.GetVertexIDsFromObject(iObject)
	
	oEditor.Fillet(
	[
		"NAME:Selections",
		"Selections:="		, iObject,
		"NewPartsModelFlag:="	, "Model"
	], 
	[
		"NAME:Parameters",
		[
			"NAME:FilletParameters",
			"Edges:="		, [],
			"Vertices:="		, [int(vertices[iIndex])],
			"Radius:="		, "("+str(iFiletRay)+")mm",
			"Setback:="		, "0mm"
		]
	])
	return

def fuse(iObjects):
	'''
	Performs the unions of the elements of iObjects

	Input:
	iObjects (list of strings): list of object names e.g. ['transmon_pad','transmon_pad_extension']

	Returns:
	string: name of the merged object: iObjects[0]
	'''

	for i in range(len(iObjects)-1):
		oEditor.Unite(
		[
			"NAME:Selections",
			"Selections:="		, iObjects[0]+","+iObjects[i+1]
		], 
		[
			"NAME:UniteParameters",
			"KeepOriginals:="	, False
		])
	
	return iObjects[0]

def assignPerfE(iObject, iName):
	'''
	Assigns boundary condition PerfE of name iName to object iObject

	Inputs:
	-------
	iObject: name of the object e.g. transmon_pad
	iName: name of the perfect E e.g. PerfE1

	'''
	oModule.AssignPerfectE(
	[
		"NAME:PerfE_"+iName,
		"Objects:="		, iObject,
		"InfGroundPlane:="	, False
	])
	return

def substract(iObjectTool, iObjectBlanc):
	'''
	suBTracts iObjectTool from iObjectBlanc

	Inputs:
	-------
	iObjectTool (string) : HFSS object name e.g. 'readout_box' 
	iObjectBlanc (string) : HFSS object name e.g. 'ground_plane'

	'''

	oEditor.Subtract(
	[
		"NAME:Selections",
		"Blank Parts:="		, iObjectBlanc,
		"Tool Parts:="		, iObjectTool
	], 
	[
		"NAME:SubtractParameters",
		"KeepOriginals:="	, False
	])
	
	return iObjectTool

	
# main function
def draw(iName, iPoints):
	'''
	Inputs
	------
	iName : object name in HFSS e.g. transmon_pad
	iPoints : list of tuples [(x0, y0), (x1, y1), ..., (x0, y0)] last and first points are the same

	'''

	pointsStr = ["NAME:PolylinePoints"]
	indexsStr = ["NAME:PolylineSegments"]
	for i in range(len(iPoints)):
		pointsStr.append(["NAME:PLPoint", "X:=", "("+str(iPoints[i].x)+")mm", "Y:=", "("+str(iPoints[i].y)+")mm", "Z:=", "0mm"])
		if(i!=len(iPoints)-1):
			indexsStr.append(["NAME:PLSegment", "SegmentType:=", "Line", "StartIndex:=", i, "NoOfPoints:=", 2])

	oEditor.CreatePolyline(
		[
			"NAME:PolylineParameters",
			"IsPolylineCovered:="	, True,
			"IsPolylineClosed:="	, True,
			pointsStr,	
			indexsStr,
			[
				"NAME:PolylineXSection",
				"XSectionType:="	, "None",
				"XSectionOrient:="	, "Auto",
				"XSectionWidth:="	, "0mm",
				"XSectionTopWidth:="	, "0mm",
				"XSectionHeight:="	, "0mm",
				"XSectionNumSegments:="	, "0",
				"XSectionBendType:="	, "Corner"
			]
		],
		[
			"NAME:Attributes",
			"Name:="		, iName,
			"Flags:="		, "",
			"Color:="		, "(143 175 143)",
			"Transparency:="	, 0,
			"PartCoordinateSystem:=", "Global",
			"UDMId:="		, "",
			"MaterialValue:="	, "\"vacuum\"",
			"SurfaceMaterialValue:=", "\"\"",
			"SolveInside:="		, True,
			"IsMaterialEditable:="	, True,
			"UseMaterialAppearance:=", False
		])
		
	return iName



def CreateBondwire(iName, iIn):
	pos, ori = Vector(iIn[POS]), Vector(iIn[ORI])
	bondWidth=0.2
	shiftbond=0.03
	oEditor.CreateBondwire(
		[
			"NAME:BondwireParameters", 
			"WireType:=", "Low", 
			"WireDiameter:=", "0.02mm",
			"NumSides:=", 6,
			"XPadPos:=", "("+str(pos.x-bondWidth/2.*ori.y+shiftbond*ori.x)+")mm",
			"YPadPos:=", "("+str(pos.y+bondWidth/2.*ori.x+shiftbond*ori.y)+")mm",
			"ZPadPos:=", "0mm",
			"XDir:=", ori.y,
			"YDir:=", -ori.x,
			"ZDir:=", 0,
			"Distance:=", str(bondWidth)+"mm",
			"h1:=", "0.1mm",
			"h2:=", "0mm",
			"alpha:=", "80deg",
			"beta:=", "80deg",
			"WhichAxis:=", "Z"
		],
		[
			"NAME:Attributes",
			"Name:="		, iName,
			"Flags:="		, "",
			"Color:="		, "(143 175 143)",
			"Transparency:="	, 0,
			"PartCoordinateSystem:=", "Global",
			"UDMId:="		, "",
			"MaterialValue:="	, "\"vacuum\"",
			"SurfaceMaterialValue:=", "\"\"",
			"SolveInside:="		, True,
			"IsMaterialEditable:="	, True,
			"UseMaterialAppearance:=", False
		])
		
	bondwireObjects.append(iName)
		
	
	




	
class LitExp(object):
	'''
	Class of HFSS variables

	attributes:
	val (float): variable value
	exp (string): variable name

	e.g.

	Lj = LitExp("$Lj", 2)
	Lj+1 is an instance of LitExp, with attributes:
	val = 2+1 = 3
	exp = '$Lj+1'

	'''
	val = 0
	exp = ''
	
	def __init__(self, iName, iDef=0,VarDef=False):
	
		self.val = iDef
		self.exp = str(iName)
		
		# first check if this variable exists, if it does, update its value, 
		# if not, create the variable and assign value.
		VariableTable=oProject.GetVariables()
		
		
			
		if str(iName) in VariableTable:
			oProject.SetVariableValue(str(iName),iDef)
			
		elif VarDef is True:
			#Create a project variable
			#oDesktop.PauseScript(str(iName))
			oProject.ChangeProperty(["NAME:AllTabs",
									["NAME:ProjectVariableTab",
									["NAME:PropServers", "ProjectVariables"],
									["NAME:NewProps",
									["NAME:"+str(iName),
									"PropType:=", "VariableProp",
									"Value:=", iDef,
									"Description:=", "" ]]]])
	
		
		
		
	def __add__(self, other):
		if(type(other) is LitExp):
			ret = self.copy()
			ret.val += other.val
			ret.exp += " + " + other.exp
			return ret
		else:
			ret = self.copy()
			ret.val += other
			ret.exp += " + " + str(other)
			return ret
		return "This is not possible"

	def __sub__(self, other):
		if(type(other) is LitExp):
			ret = self.copy()
			ret.val -= other.val
			ret.exp += " - (" + other.exp + ")"
			return ret
		else:
			ret = self.copy()
			ret.val -= other
			ret.exp += " - (" + str(other) + ")"
			return ret
		return "This is not possible"

	def __mul__(self, other):
		if(type(other) is LitExp):
			ret = self.copy()
			ret.val *= other.val
			ret.exp = "(" + ret.exp + ")*(" + other.exp + ")"
			return ret
		elif(type(other) is Vector):
			return other*self
		else:
			ret = self.copy()
			ret.val *= other
			ret.exp = "(" + ret.exp + ")*(" + str(other) + ")"
			return ret
		return "This is not possible"
	
	def __div__(self, other):
		if(type(other) is LitExp):
			ret = self.copy()
			ret.val /= other.val
			ret.exp = "(" + ret.exp + ")/(" + other.exp + ")"
			return ret
		elif(type(other) is Vector):
			return other/self
		else:
			ret = self.copy()
			ret.val /= other
			ret.exp = "(" + ret.exp + ")/(" + str(other) + ")"
			return ret
		return "This is not possible"


	def __radd__(self, other):
		return LitExp(str(other), other)+self
	def __rsub__(self, other):
		return LitExp(str(other), other)-self
	def __rmul__(self, other):
		return LitExp(str(other), other)*self
	def __rdiv__(self, other):
		return LitExp(str(other), other)/self

	def __neg__(self):
		ret = self.copy()
		ret.val *= -1
		ret.exp = " - (" + ret.exp + ")"
		return ret

	def __lt__(self, other):
		return self.val < other
	def __le__(self, other):
		return self.val <= other
	def __eq__(self, other):
		return self.val == other
		
	def __ne__(self, other):
		return self.val != other
		
	def __gt__(self, other):
		return self.val > other
		
	def __ge__(self, other):
		return self.val >= other
	
	def __abs__(self):
		ret = self.copy()
		ret.val = abs(self.val)
		ret.exp = "abs(" + self.exp + ")"
		return ret

	def __str__(self):
		return self.exp
	def __float__(self):
		return self.val

	def copy(self):
		ret = LitExp(""+self.exp,self.val)
		#ret.val = self.val
		#ret.exp = ""+self.exp
		return ret	

class Vector(object):
	'''
	Class of vectors in the xy plane
	Attributes:
	x (float or LitExp) e.g. x = $witdh+2
	y (float or LitExp)
	'''

	x = 0
	y = 0
	def __init__(self, iX, iY=0):
		if(type(iX) is Vector):
			self.x, self.y = iX.x, iX.y
		elif(type(iX) is list):
			self.x, self.y = iX[0], iX[1]
		else:
			self.x, self.y = iX, iY

	def __add__(self, other):
		if(type(other) is list):
			return Vector(self.x+other[0], self.y+other[1])
		elif(type(other) is Vector):
			return Vector(self.x+other.x, self.y+other.y)
		else:
			return Vector(self.x+other, self.y+other)
		return "This is not possible"

	def __sub__(self, other):
		if(type(other) is list):
			return Vector(self.x-other[0], self.y-other[1])
		elif(type(other) is Vector):
			return Vector(self.x-other.x, self.y-other.y)
		else:
			return Vector(self.x-other, self.y-other)
		return "This is not possible"

	def __mul__(self, other):
		if(type(other) is list):
			return self.x*other[0]+self.y*other[1]
		elif(type(other) is Vector):
			return self.x*other.x+self.y*other.y
		else:
			return Vector(self.x*other, self.y*other)
		return "This is not possible"
		
	def dot(self, other):
		if(type(other) is list):
			return Vector(self.x*other[0], self.y*other[1])
		elif(type(other) is Vector):
			return Vector(self.x*other.x, self.y*other.y)
		else:
			return Vector(self.x*other, self.y*other)
		return "This is not possible"

	def __radd__(self, other):
		return Vector(other)+self

	def __rsub__(self, other):
		return Vector(other)-self

	def __rmul__(self, other):
		return self*other

	def __neg__(self):
		return Vector(-self.x, -self.y)

	def __str__(self):
		return "["+str(self.x)+","+str(self.y)+"]"

	def norm(self):
		return (self.x**2+self.y**2)**0.5

	def abs(self):
		return Vector(abs(self.x), abs(self.y))

	def unit(self):
		norm = self.norm()
		return Vector(self.x/norm, self.y/norm)

	def orth(self):
		return Vector(self.y, -self.x)

	def rot(self, other):
		''' 
		Inputs:
		-------
		other: vector
		
		Returns
		-------
		vector: rotation around z of self by an angle given by other w.r.t to x
		''' 
		unitOther = other.unit()
		return Vector(self*unitOther, self*unitOther.orth())

def drawCapa(iName, iIn, iOut, iLength, iWidth, iSize):
	'''
	Inputs:
	-------
	iName: string name of object
	iIn: (position, direction, track, gap) defines the input port
	iOut: (position, direction, track, gap) defines the output port
	       position and direction are None: this is calculated from
		   other parameters
	iLength: (float) length of pads
	iWidth: (float) width of pads
	iSize: (float) spacing between pads (see drawing)
	
	Outputs:
	--------
	retIn: same as iIn, with flipped vector
	retOut: calculated output port to match all input dimensions
	
		     iSize
		  +--+  +--+
		  |  |  |  |
		+-+  |  |  +-+
	iIn |    |	|    | iOut
		+-+  |  |  +-+
		  |  |  |  |
		  +--+  +--+
	'''

	pos, ori = Vector(iIn[POS]), Vector(iIn[ORI])
	_, _, inTrack, inGap = iIn
	_, _, outTrack, outGap = iOut
	retIn = [pos, -ori, inTrack, inGap]
	retOut = [pos+(inGap+outGap+iSize+2*iWidth)*ori, ori, outTrack, outGap]

	points = [pos + Vector(inGap+iWidth, 0).rot(ori)]
	points.append(points[-1] + Vector(0, -iLength/2).rot(ori))
	points.append(points[-1] + Vector(-iWidth, 0).rot(ori))
	points.append(points[-1] + Vector(0, iLength/2-inTrack/2).rot(ori))
	points.append(points[-1] + Vector(-inGap, 0).rot(ori))
	points.append(points[-1] + Vector(0, inTrack).rot(ori))
	points.append(points[-1] + Vector(inGap, 0).rot(ori))
	points.append(points[-1] + Vector(0, iLength/2-inTrack/2).rot(ori))
	points.append(points[-1] + Vector(iWidth, 0).rot(ori))
	points.append(points[0])
	trackObjects.append(draw(iName+"_track1", points))

	points = [pos + Vector(inGap+iWidth+iSize, 0).rot(ori)]
	points.append(points[-1] + Vector(0, -iLength/2).rot(ori))
	points.append(points[-1] + Vector(+iWidth, 0).rot(ori))
	points.append(points[-1] + Vector(0, iLength/2-outTrack/2).rot(ori))
	points.append(points[-1] + Vector(+outGap, 0).rot(ori))
	points.append(points[-1] + Vector(0, outTrack).rot(ori))
	points.append(points[-1] + Vector(-outGap, 0).rot(ori))
	points.append(points[-1] + Vector(0, iLength/2-outTrack/2).rot(ori))
	points.append(points[-1] + Vector(-iWidth, 0).rot(ori))
	points.append(points[0])
	trackObjects.append(draw(iName+"_track2", points))

	points = [pos]
	points.append(points[-1] + Vector(0, iLength/2+inGap).rot(ori))
	points.append(points[-1] + Vector(inGap+iWidth+iSize/2, 0).rot(ori))
	if(outGap-inGap!=0):
		points.append(points[-1] + Vector(0, outGap-inGap).rot(ori))
	points.append(points[-1] + Vector(outGap+iWidth+iSize/2, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iLength-2*outGap).rot(ori))
	points.append(points[-1] + Vector(-(outGap+iWidth+iSize*0.5), 0).rot(ori))
	if(outGap-inGap!=0):
		points.append(points[-1] + Vector(0, outGap-inGap).rot(ori))
	points.append(points[-1] + Vector(-(inGap+iWidth+iSize*0.5), 0).rot(ori))
	points.append(points[0])
	gapObjects.append(draw(iName+"_gap", points))

	draw(iName+"_mesh", points)
	assignMeshLength(iName+"_mesh",0.05)

	
	return [retIn, retOut]

def drawCombCapa(iName, iIn, iOut, iLength, iWidth, iSize, iN):
	'''
	Draws an interdigitated capacitor
	
	Inputs:
	-------
	see drawCapa, additional inputs:
	iN: (int) the total number of fingers is 4*iN+1
	2*N on the left, 2*N+1 on the right 
	iLength: (float) length of fingers
	iWidth: (float) width of base pad, and width of fingers
	
	Outputs:
	--------
	see drawCapa
	'''

	pos, ori = Vector(iIn[POS]), Vector(iIn[ORI])
	_, _, inTrack, inGap = iIn
	_, _, outTrack, outGap = iOut
	retIn = [pos, -ori, inTrack, inGap]
	retOut = [pos+(inGap+outGap+iSize+iLength+2*iWidth)*ori, ori, outTrack, outGap]

	points = [pos + Vector(inGap+iWidth, 0).rot(ori)]
	for i in range(iN):
		points.append(points[-1] + Vector(0, iSize+iWidth*0.5).rot(ori))
		points.append(points[-1] + Vector(iLength, 0).rot(ori))
		points.append(points[-1] + Vector(0, iWidth).rot(ori))
		points.append(points[-1] + Vector(-iLength, 0).rot(ori))
		points.append(points[-1] + Vector(0, iSize+iWidth*0.5).rot(ori))
	points.append(points[-1] + Vector(0, iWidth*0.5).rot(ori))
	points.append(points[-1] + Vector(-iWidth, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iWidth*0.5-iN*(2*iSize+2*iWidth)+inTrack*0.5).rot(ori))
	points.append(points[-1] + Vector(-inGap, 0).rot(ori))
	points.append(points[-1] + Vector(0, -inTrack).rot(ori))
	points.append(points[-1] + Vector(inGap, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iWidth*0.5-iN*(2*iSize+2*iWidth)+inTrack*0.5).rot(ori))
	points.append(points[-1] + Vector(iWidth, 0).rot(ori))
	points.append(points[-1] + Vector(0, iWidth*0.5).rot(ori))
	for i in range(iN):
		points.append(points[-1] + Vector(0, iSize+iWidth*0.5).rot(ori))
		points.append(points[-1] + Vector(iLength, 0).rot(ori))
		points.append(points[-1] + Vector(0, iWidth).rot(ori))
		points.append(points[-1] + Vector(-iLength, 0).rot(ori))
		if(i==iN-1):
			points.append(points[0])
		else:
			points.append(points[-1] + Vector(0, iSize+iWidth*0.5).rot(ori))
	trackObjects.append(draw(iName+"_track1", points))

	points = [pos + Vector(inGap+iWidth+iSize, 0).rot(ori)]
	for i in range(iN):
		points.append(points[-1] + Vector(0, iWidth*0.5).rot(ori))
		points.append(points[-1] + Vector(iLength, 0).rot(ori))
		points.append(points[-1] + Vector(0, iWidth+2*iSize).rot(ori))
		points.append(points[-1] + Vector(-iLength, 0).rot(ori))
		points.append(points[-1] + Vector(0, iWidth*0.5).rot(ori))
	points.append(points[-1] + Vector(0, iWidth*0.5).rot(ori))
	points.append(points[-1] + Vector(iWidth+iLength, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iWidth*0.5-iN*(2*iSize+2*iWidth)+outTrack*0.5).rot(ori))
	points.append(points[-1] + Vector(outGap, 0).rot(ori))
	points.append(points[-1] + Vector(0, -outTrack).rot(ori))
	points.append(points[-1] + Vector(-outGap, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iWidth*0.5-iN*(2*iSize+2*iWidth)+outTrack*0.5).rot(ori))
	points.append(points[-1] + Vector(-iWidth-iLength, 0).rot(ori))
	points.append(points[-1] + Vector(0, iWidth*0.5).rot(ori))
	for i in range(iN):
		points.append(points[-1] + Vector(0, iWidth*0.5).rot(ori))
		points.append(points[-1] + Vector(iLength, 0).rot(ori))
		points.append(points[-1] + Vector(0, iWidth+2*iSize).rot(ori))
		points.append(points[-1] + Vector(-iLength, 0).rot(ori))
		if(i==iN-1):
			points.append(points[0])
		else:
			points.append(points[-1] + Vector(0, iWidth*0.5).rot(ori))
	trackObjects.append(draw(iName+"_track2", points))

	points = [pos]
	points.append(points[-1] + Vector(0, -(inGap+iWidth*0.5+iN*(2*iSize+2*iWidth))).rot(ori))
	points.append(points[-1] + Vector(inGap+iWidth+iSize*0.5, 0).rot(ori))
	if(outGap-inGap!=0):
		points.append(points[-1] + Vector(0, -outGap+inGap).rot(ori))
	points.append(points[-1] + Vector(outGap+iWidth+iLength+iSize*0.5, 0).rot(ori))
	points.append(points[-1] + Vector(0, 2*outGap+2*iN*(2*iSize+2*iWidth)+iWidth).rot(ori))
	points.append(points[-1] + Vector(-(outGap+iWidth+iLength+iSize*0.5), 0).rot(ori))
	if(outGap-inGap!=0):
		points.append(points[-1] + Vector(0, -outGap+inGap).rot(ori))
	points.append(points[-1] + Vector(-(inGap+iWidth+iSize*0.5), 0).rot(ori))
	points.append(points[0])
	gapObjects.append(draw(iName+"_gap", points))

	draw(iName+"_mesh", points)
	assignMeshLength(iName+"_mesh",0.05)
	
	return [retIn, retOut]

def drawConnector(iName, iIn, iOut, iBoundLength, iSlope=1, iLineTest=False):
	'''
	Draws a CPW connector for inputs and outputs.
	
	Inputs:
	-------
	
	iBoundLength: (float) corresponds to dimension a in the drawing
	iSlope (float): 1 is 45 degrees to adpat lengths between bounding pad and iout 
	iLineTest (Bool): unclear, keep False
	
	    ground plane
		+------+
		|       \
		|        \
		|  	+-a+\|
	iIn	|   |Cond|  iOut
		|   +--+/|
		|        /
		|       /
		+------+
		
	Outputs:
	--------
	returns iIn and recalculated iOut
	'''

	pos, ori = Vector(iIn[POS]), Vector(iIn[ORI])
	_, _, inTrack, inGap = iIn
	_, _, outTrack, outGap = iOut
	adaptDist = abs(outTrack/2-inTrack/2)/iSlope
	retIn = [pos, -ori, inTrack, inGap]
	retOut = [pos+(adaptDist+inGap+iBoundLength)*ori, ori, outTrack, outGap]

	points = [pos + Vector(inGap, inTrack/2.).rot(ori)]
	points.append(points[-1] + Vector(iBoundLength+iLineTest*inGap, 0).rot(ori))
	points.append(points[-1] + Vector(adaptDist, outTrack/2-inTrack/2).rot(ori))
	points.append(points[-1] + Vector(0, -outTrack).rot(ori))
	points.append(points[-1] + Vector(-adaptDist, outTrack/2-inTrack/2).rot(ori))
	points.append(points[-1] + Vector(-iBoundLength-iLineTest*inGap, 0).rot(ori))
	points.append(points[0])
	trackObjects.append(draw(iName+"_track1", points))

	points = [pos + Vector(inGap/2, inGap+inTrack/2.).rot(ori)]
	points.append(points[-1] + Vector(inGap/2+iBoundLength, 0).rot(ori))
	points.append(points[-1] + Vector(adaptDist, (outGap-inGap)+(outTrack-inTrack)*0.5).rot(ori))
	points.append(points[-1] + Vector(0, -2*outGap-outTrack).rot(ori))
	points.append(points[-1] + Vector(-adaptDist, (outGap-inGap)+(outTrack-inTrack)*0.5).rot(ori))
	points.append(points[-1] + Vector(-(inGap/2+iBoundLength), 0).rot(ori))
	points.append(points[0])
	gapObjects.append(draw(iName+"_gap", points))
	
	if(iLineTest==False):
		points = [pos + Vector(inGap/2, 0.5*inTrack).rot(ori)]
		points.append(points[-1] + Vector(inGap/2, 0).rot(ori))
		points.append(points[-1] + Vector(0, -inTrack).rot(ori))
		points.append(points[-1] + Vector(-inGap/2, 0).rot(ori))
		points.append(points[0])
		draw(iName+"_track2", points)
		
		assigne50Ohm(iName+"_track2", pos+Vector(inGap/2., 0).rot(ori), pos+Vector(inGap, 0).rot(ori))
		
	
	# if(iLineTest==False):
		# points = [pos + Vector(0, inGap+0.5*inTrack).rot(ori)]
		# points.append(points[-1] + Vector(inGap/2, 0).rot(ori))
		# points.append(points[-1] + Vector(0, -2*inGap-inTrack).rot(ori))
		# points.append(points[-1] + Vector(-inGap/2, 0).rot(ori))
		# points.append(points[0])
		# trackObjects.append(draw(iName+"_track2", points))

	return [retIn, retOut]

def drawAdaptor(iName, iIn, iOut, iSlope=1):
	'''
	Draws an adaptor between two ports.
	Given input port iIn, and slope of line iSlope, calculates iOut, and draws adpator.
	
	Inputs:
	-------
	iName:
	iIn: tuple, input port
	iOut: tuple, output port, usually Nones
	iSlope: slope of line to connect iIn and iOut. If iSlope=1, 45 degrees.
	
	Returns:
	--------
	reversed iIn and calculated iOut
	'''

	pos, ori = Vector(iIn[POS]), Vector(iIn[ORI])
	_, _, inTrack, inGap = iIn
	_, _, outTrack, outGap = iOut
	retIn = [pos, -ori, inTrack, inGap]
	retOut = [pos+abs(outTrack/2-inTrack/2)*ori, ori, outTrack, outGap]
	adaptDist = abs(outTrack/2-inTrack/2)/iSlope

	points = [pos + Vector(0, inTrack/2).rot(ori)]
	points.append(points[-1] + Vector(adaptDist, outTrack/2-inTrack/2).rot(ori))
	points.append(points[-1] + Vector(0, -outTrack).rot(ori))
	points.append(points[-1] + Vector(-adaptDist, outTrack/2-inTrack/2).rot(ori))
	points.append(points[0])
	trackObjects.append(draw(iName+"_track", points))

	points = [pos + Vector(0, inGap+0.5*inTrack).rot(ori)]
	points.append(points[-1] + Vector(adaptDist, (outGap-inGap)+(outTrack-inTrack)*0.5).rot(ori))
	points.append(points[-1] + Vector(0, -2*outGap-outTrack).rot(ori))
	points.append(points[-1] + Vector(-adaptDist, (outGap-inGap)+(outTrack-inTrack)*0.5).rot(ori))
	points.append(points[0])
	gapObjects.append(draw(iName+"_gap", points))

	return [retIn, retOut]

def drawTriJunction(iName, iIn, iOut1, iOut3):
	'''
	Draws a T object
	.: metal
	
		 iOut1	
		+-+-+-+
		| |.| | 
		+-+.| |
	iIn |...| |
		+-+.| |
		| |.| |
		+-+-+-+
		 iOut2  
		  
	'''

	pos, ori = Vector(iIn[POS]), Vector(iIn[ORI])
	_, _, inTrack, inGap = iIn
	_, _, out1Track, out1Gap = iOut1
	_, _, out3Track, out3Gap = iOut3
	retIn = [pos, -ori, inTrack, inGap]
	retOut1 = [pos+Vector(inGap+inTrack/2, inGap+inTrack/2).rot(ori), Vector(0, 1).rot(ori), out1Track, out1Gap]
	retOut3 = [pos+Vector(inGap+inTrack/2, -(inGap+inTrack/2)).rot(ori), Vector(0, -1).rot(ori), out3Track, out3Gap]
	
	points = [pos + Vector(0, inGap+inTrack/2).rot(ori)]
	points.append(points[-1] + Vector(2*inGap+inTrack, 0).rot(ori))
	points.append(points[-1] + Vector(0, -(2*inGap+inTrack)).rot(ori))
	points.append(points[-1] + Vector(-(2*inGap+inTrack), 0).rot(ori))
	points.append(points[0])
	gapObjects.append(draw(iName+"_gap", points))
	
	points = [pos + Vector(0, inTrack/2).rot(ori)]
	points.append(points[-1] + Vector(inGap, 0).rot(ori))
	points.append(points[-1] + Vector(0, inGap).rot(ori))
	points.append(points[-1] + Vector(inTrack, 0).rot(ori))
	points.append(points[-1] + Vector(0, -2*inGap-inTrack).rot(ori))
	points.append(points[-1] + Vector(-inTrack, 0).rot(ori))
	points.append(points[-1] + Vector(0, inGap).rot(ori))
	points.append(points[-1] + Vector(-inGap, 0).rot(ori))
	points.append(points[0])
	trackObjects.append(draw(iName+"_track", points))
	
	return [retIn, retOut3, retOut1]
	
def drawQuadJunction(iName, iIn, iOut1, iOut2, iOut3):
	'''
	Draws a cross object
	.: metal
	
		 iOut1	
		+-+-+-+
		| |.| |
		+-+.+-+ 
	iIn |.....| iOut2
		+-+.+-+ 
		| |.| |
		+-+-+-+
		 iOut3
		  
	'''

	pos, ori = Vector(iIn[POS]), Vector(iIn[ORI])
	_, _, inTrack, inGap = iIn
	_, _, out1Track, out1Gap = iOut1
	_, _, out2Track, out2Gap = iOut2
	_, _, out3Track, out3Gap = iOut3
	retIn = [pos, -ori, inTrack, inGap]
	retOut1 = [pos+Vector(inGap+inTrack/2, inGap+inTrack/2).rot(ori), Vector(0, 1).rot(ori), out1Track, out1Gap]
	retOut2 = [pos+Vector(2*inGap+inTrack, 0).rot(ori), Vector(1, 0).rot(ori), out2Track, out2Gap]
	retOut3 = [pos+Vector(inGap+inTrack/2, -(inGap+inTrack/2)).rot(ori), Vector(0, -1).rot(ori), out3Track, out3Gap]
	
	points = [pos + Vector(0, inGap+inTrack/2).rot(ori)]
	points.append(points[-1] + Vector(2*inGap+inTrack, 0).rot(ori))
	points.append(points[-1] + Vector(0, -(2*inGap+inTrack)).rot(ori))
	points.append(points[-1] + Vector(-(2*inGap+inTrack), 0).rot(ori))
	points.append(points[0])
	gapObjects.append(draw(iName+"_gap", points))
	
	points = [pos + Vector(0, inTrack/2).rot(ori)]
	points.append(points[-1] + Vector(inGap, 0).rot(ori))
	points.append(points[-1] + Vector(0, inGap).rot(ori))
	points.append(points[-1] + Vector(inTrack, 0).rot(ori))
	points.append(points[-1] + Vector(0, -inGap).rot(ori))
	points.append(points[-1] + Vector(inGap, 0).rot(ori))
	points.append(points[-1] + Vector(0, -inTrack).rot(ori))
	points.append(points[-1] + Vector(-inGap, 0).rot(ori))
	points.append(points[-1] + Vector(0, -inGap).rot(ori))
	points.append(points[-1] + Vector(-inTrack, 0).rot(ori))
	points.append(points[-1] + Vector(0, inGap).rot(ori))
	points.append(points[-1] + Vector(-inGap, 0).rot(ori))
	points.append(points[0])
	trackObjects.append(draw(iName+"_track", points))
	
	return [retIn, retOut3, retOut2, retOut1]
	
def drawEnd(iName, iIn):
	'''
	Draws a square cap on iIn
	
	'''

	pos, ori = Vector(iIn[POS]), Vector(iIn[ORI])
	_, _, inTrack, inGap = iIn
	retIn = [pos, -ori, inTrack, inGap]
	
	points = [pos + Vector(0, inGap+inTrack/2).rot(ori)]
	points.append(points[-1] + Vector(2*inGap+inTrack, 0).rot(ori))
	points.append(points[-1] + Vector(0, -(2*inGap+inTrack)).rot(ori))
	points.append(points[-1] + Vector(-(2*inGap+inTrack), 0).rot(ori))
	points.append(points[0])
	gapObjects.append(draw(iName+"_gap", points))
	
	points = [pos + Vector(0, inTrack/2).rot(ori)]
	points.append(points[-1] + Vector(inGap+inTrack, 0).rot(ori))
	points.append(points[-1] + Vector(0, -inTrack).rot(ori))
	points.append(points[-1] + Vector(-inGap-inTrack, 0).rot(ori))
	points.append(points[0])
	trackObjects.append(draw(iName+"_track", points))
	
	return retIn

	
def drawStraightCable(iName, iIn, iOut):
	inPos, inOri = Vector(iIn[POS]), Vector(iIn[ORI])
	outPos, outOri = Vector(iOut[POS]), Vector(iOut[ORI])
	
	_, _, track, gap = iIn

	points = [inPos -inOri.orth().abs()*0.5*track]
	points.append(outPos-outOri.orth().abs()*0.5*track)
	points.append(outPos+outOri.orth().abs()*0.5*track)
	points.append(inPos+inOri.orth().abs()*0.5*track)
	points.append(points[0])
	trackObjects.append(draw(iName+"_track", points))

	points = [inPos -inOri.orth().abs()*(0.5*track+gap)]
	points.append(outPos-outOri.orth().abs()*(0.5*track+gap))
	points.append(outPos+outOri.orth().abs()*(0.5*track+gap))
	points.append(inPos+inOri.orth().abs()*(0.5*track+gap))
	points.append(points[0])
	gapObjects.append(draw(iName+"_gap", points))
	
	CreateBondwire(iName+"_bondwire", [(inPos+outPos)*0.5,inOri, track, gap])
	
	
def drawElbowCable(iName, iIn, iOut,iMaxFilet=0.5, extrabonds=None):

	inPos, inOri = Vector(iIn[POS]), Vector(iIn[ORI])
	outPos, outOri = Vector(iOut[POS]), Vector(iOut[ORI])

	_, _, track, gap = iIn

	inLength = (outPos-inPos).rot(inOri).abs().x
	outLength = (outPos-inPos).rot(inOri).abs().y

	points = [inPos -outOri*0.5*track]
	points.append(points[-1] + (inLength-0.5*track)*inOri)
	points.append(points[-1] + -(outLength-0.5*track)*outOri)
	points.append(points[-1] + track*inOri)
	points.append(points[-1] + (outLength+0.5*track)*outOri)
	points.append(points[-1] + -(inLength+0.5*track)*inOri)
	points.append(points[0])
	trackObjects.append(draw(iName+"_track", points))

	points = [inPos + -(0.5*track+gap)*outOri]
	points.append(points[-1] + (inLength-(0.5*track+gap))*inOri)
	points.append(points[-1] + -(outLength-(0.5*track+gap))*outOri)
	points.append(points[-1] + (track+2*gap)*inOri)
	points.append(points[-1] + (outLength+(0.5*track+gap))*outOri)
	points.append(points[-1] + -(inLength+(0.5*track+gap))*inOri)
	points.append(points[0])
	gapObjects.append(draw(iName+"_gap", points))

	filetRay = min(iMaxFilet, min(inLength, outLength)-0.5*track-gap-0.02)
	if(iMaxFilet != None and filetRay >0):
		filet(iName+"_track", filetRay+gap, 1)
		filet(iName+"_track", filetRay+gap+track, 5)
		filet(iName+"_gap", filetRay, 1)
		filet(iName+"_gap", filetRay+2*gap+track, 5)
		
	CreateBondwire(iName+"_bondwire", iIn)
	if extrabonds is not None:
		for i in range(1,extrabonds+1):
			iIn_wire=[inPos+i/extrabonds*(inLength-filetRay)*inOri,inOri,track,gap]
			CreateBondwire(iName+"_bondwireIn"+str(i), iIn_wire)
					
			iOut_wire=[outPos+i/extrabonds*(outLength-filetRay)*outOri,outOri,track,gap]
			CreateBondwire(iName+"_bondwireOut"+str(i), iOut_wire)
	
def drawCable(iName, iIn, iOut, iMaxFilet=0.5,extrabonds=None):
	'''
	Draws a CPW transmission line between iIn and iOut
	
	if iIn and iOut are facing eachother, and offset,
	draws a transmission line with two elbows half-way in between.
	
	if iIn and iOut are perpendicular,
	draws a transmission line with one elbow.
	
	if iIn and iOut do not have the same track/ gap size, this function calls
	drawAdaptor before iOut.
	
	N.B: do not separate the two ports by their track or gap size.
	
	Inputs:
    -------
	iName: (string) base-name of object, draws 'iName_adaptor' etc
	iIn: (tuple) input port
	iOut: (tuple) output port
	iMaxFilet: (float), maximum filet radius
	
	'''

	inPos, inOri = Vector(iIn[POS]), Vector(iIn[ORI])
	outPos, outOri = Vector(iOut[POS]), Vector(iOut[ORI])
	_, _, track, gap = iIn
	
	
	if(iIn[TRACK] != iOut[TRACK] or iIn[GAP] != iOut[GAP]):
		if(iIn[TRACK]+iIn[GAP] > iOut[TRACK] + iOut[GAP]):
			_, newOut = drawAdaptor(iName+"_adaptor", iOut, [None, None, iIn[TRACK], iIn[GAP]])
			outPos = newOut[POS]
			outOri = newOut[ORI]
		else:
			_, newIn = drawAdaptor(iName+"_adaptor", iIn, [None, None, iOut[TRACK], iOut[GAP]])
			inPos = newIn[POS]
			inOri = newIn[ORI]
			track = newIn[TRACK]
			gap = newIn[GAP]

	if(inOri*outOri != 0):#if in/out in the same direction
		
		middle = 0.5*(inPos+outPos)

		if abs(((inPos-outPos)*inOri.orth()) == 0.): #if in/out aligned
			
			drawStraightCable(iName, iIn, iOut)
			
			
		else:	#if in/out not aligned
		
			if(inOri.y==0.):
				middleTop = [middle, [0,1], track, gap]
				middleBottom = [middle, [0,-1], track, gap]
				if(inPos.y>outPos.y):
					drawElbowCable(iName+"_in", [inPos, inOri, track, gap], middleTop, iMaxFilet,extrabonds=extrabonds)
					drawElbowCable(iName+"_out", middleBottom, [outPos, outOri, track, gap], iMaxFilet,extrabonds=extrabonds)
					return
				else:
					drawElbowCable(iName+"_in", [inPos, inOri, track, gap], middleBottom, iMaxFilet,extrabonds=extrabonds)
					drawElbowCable(iName+"_out", middleTop, [outPos, outOri, track, gap], iMaxFilet,extrabonds=extrabonds)
					return
			else:
				middleRight = [middle, [1,0], track, gap]
				middleLeft = [middle, [-1,0], track, gap]
				if(inPos.x>outPos.x):
					drawElbowCable(iName+"_in", [inPos, inOri, track, gap], middleRight, iMaxFilet,extrabonds=extrabonds)
					drawElbowCable(iName+"_out", middleLeft, [outPos, outOri, track, gap], iMaxFilet,extrabonds=extrabonds)
					return
				else:
					drawElbowCable(iName+"_in", [inPos, inOri, track, gap], middleLeft, iMaxFilet,extrabonds=extrabonds)
					drawElbowCable(iName+"_out", middleRight, [outPos, outOri, track, gap], iMaxFilet,extrabonds=extrabonds)
					return
					
	else:#if in/out in different directions

		drawElbowCable(iName, iIn, iOut,iMaxFilet,extrabonds=extrabonds)

	#CreateBondwire(iName+"_bondwire", iIn)
	#CreateBondwire(iName+"_bondwire", iOut)
	
def drawWave(iName, iIn, iOut, iNbWaves, iRightShift, iLeftShift, iGndSize):
	'''
	Draws a meander 
	
	Inputs:
	-------
	iIn: input port, imposes track and gap size
	iOut: output port, None, calculated to match other inputs
	iNbWaves: (int) if iNbWaves = 1, draws one period. If iNbWaves=1+N, draws 1+(N/2) periods
			   e.g iNbWaves = 3, draws one period + 2 half periods, so 2 periods
	iRightShift: (float) 
	iLeftShift: (float) for now, needs to be equal to iRightShift
	iRightShift + iLeftShift is the oscillation amplitude
	iGndSize: (float) hald the oscillation period
	
	Outputs:
	--------
	'''

	pos, ori = Vector(iIn[POS]), Vector(iIn[ORI])
	_, _, track, gap = iIn
	period = 2*gap+track+iGndSize/2
	border = 2*gap+track+iGndSize

	retOut = [pos + (2*border+iNbWaves*period)*ori, ori, track, gap]
	retIn = [pos, -ori, track, gap]

	for i in range(iNbWaves+2):
		j = ((i%2)-0.5)*2
		waveIn = [pos + Vector(0, j*(iRightShift*(j%2)+iLeftShift*(1-j%2))).rot(ori) + (border+(i-1)*period)*ori, ori, track, gap]
		waveOut = [pos + Vector(0, -j*(iLeftShift*(j%2)+iRightShift*(1-j%2))).rot(ori) + (border+i*period)*ori, -ori, track, gap]

		if(i==0):
			waveIn=iIn
		if(i==iNbWaves+1):
			waveOut=[pos + (2*border+iNbWaves*period)*ori, -ori, track, gap]

		drawCable(iName+"_waveN"+str(i), waveIn, waveOut)

	return [retIn, retOut]
	
	
def ShiftPortStraigth(iName,iIn,iStraight):
	pos, ori = Vector(iIn[POS]), Vector(iIn[ORI])
	_, _, inTrack, inGap = iIn
	iIn_2 = [pos+iStraight*ori, ori, inTrack, inGap]
	drawStraightCable(iName+'_straight1',iIn,iIn_2)
	return iIn_2
	
	
def drawWave2(iName, iIn, iOut, iNbWaves, iShift=0,iOnesideShift=0,Lres=None, iStraight=None):
	'''
	Draws a meander 
	
	Inputs:
	-------
	iIn: input port, imposes track and gap size
	iOut: output port, None, calculated to match other inputs
	iNbWaves: (int) if iNbWaves = 1, draws one period. If iNbWaves=1+N, draws 1+(N/2) periods
			   e.g iNbWaves = 3, draws one period + 2 half periods, so 2 periods
	iShift: Overshoot on the meanders
	
	Outputs:
	--------
	'''

	if iStraight is not None:
		
		iIn=ShiftPortStraigth(iName+'_straight1',iIn,iStraight)
		iOut=ShiftPortStraigth(iName+'_straight2',iOut,iStraight)
	
	pos_in, ori_in = Vector(iIn[POS]), Vector(iIn[ORI])

	pos_out, ori_out = Vector(iOut[POS]), Vector(iOut[ORI])

	_, _, track, gap = iIn
	
	test=(pos_out-pos_in).rot(ori_in.abs())
	if Lres is not None:
		if iStraight is not None:
			if test.y >= 0:
				iShift=(Lres-abs(test.x)-abs(test.y)*(2*iNbWaves+1))/(4*iNbWaves)
			elif test.y < 0:
				iShift=-(Lres-abs(test.x)-abs(test.y)*(2*iNbWaves+1))/(4*iNbWaves)
		else:
			if test.y >= 0:
				iShift=(Lres-2*iStraight-abs(test.x)-abs(test.y)*(2*iNbWaves+1))/(4*iNbWaves)
			elif test.y < 0:
				iShift=-(Lres-2*iStraight-abs(test.x)-abs(test.y)*(2*iNbWaves+1))/(4*iNbWaves)


	for i in range(2*iNbWaves+1):
		if i==0:
			waveIn = iIn
			waveOut = [pos_in + (pos_out-pos_in).dot(Vector((i+1)/(2.*iNbWaves+1.),((i+1)%2)).rot(ori_in.abs()))+Vector(0,(2*((i+1)%2)-1)*iShift+iOnesideShift).rot(ori_in.abs()),
					   ori_out, track, gap]
					   #+Vector(0,(2*((i+1)%2)-1)*iShift).rot(ori_in)
		elif i==2*iNbWaves:
			#oDesktop.PauseScript(str("test2"))
			waveIn=waveOut
			waveIn[ORI]=ori_in
			waveOut = iOut 	#[pos_in + (pos_out-pos_in)*(Vector(i/(2*n+1),(i+1)%2).rot(ori_in))+Vector(0,(2*((i+1)%2)-1)*extra).rot(ori_in)),
							# ori_in, track, gap]
		else:
			#oDesktop.PauseScript(str("test"))
			waveIn=waveOut
			waveIn[ORI]=ori_in

			waveOut = [pos_in + (pos_out-pos_in).dot(Vector((i+1)/(2.*iNbWaves+1.),((i+1)%2)).rot(ori_in.abs()))+Vector(0,(2*((i+1)%2)-1)*iShift+iOnesideShift).rot(ori_in.abs()),
						ori_out, track, gap]			
		#+Vector(0,(2*((i+1)%2)-1)*2*iShift).rot(ori_in)
		#waveOut[ORI]=-waveOut[ORI]
		drawCable(iName+"_waveN"+str(i),waveIn, waveOut)
		#oDesktop.PauseScript(str(ori_in))

	#CreateBondwire(iName+"_bondwire", iOut)
		
	return [iIn, iOut]	


def drawJSJunc(iName, iIn, iOut, iSize, iWidth, iLength, iInduct=0.1):
	'''
	Draws a Joseph's Son Junction. 
	
	Draws a rectangle, here called "junction", 
	with boundary condition :lumped RLC, C=R=0, L=iInduct in nH
	Draws needed adaptors on each side
	
	Inputs:
	-------
	iName:
	iIn: (tuple) input port
	iOut: (tuple) output port - None, ignored and recalculated
	iSize: (float) length of junction
	iWidth: (float) width of junction
	iLength: (float) distance between iIn and iOut, including
			 the adaptor length 
	iInduct: (float in nH)
	
	Outputs:
	--------
	
	'''

	pos, ori = Vector(iIn[POS]), Vector(iIn[ORI])
	_, _, track, gap = iIn
	adaptDist = abs(iWidth/2-track/2)
	retIn = [pos, -ori, track, gap]
	retOut = [pos+abs(2*adaptDist+2*iLength+iSize)*ori, ori, track, gap]	
	
	points = [pos + Vector(0, 0.5*track).rot(ori)]
	points.append(points[-1] + Vector(adaptDist, iWidth/2-track/2).rot(ori))
	points.append(points[-1] + Vector(iLength, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iWidth).rot(ori))
	points.append(points[-1] + Vector(-iLength, 0).rot(ori))
	points.append(points[-1] + Vector(-adaptDist, iWidth/2-track/2).rot(ori))
	points.append(points[0])
	trackObjects.append(draw(iName+"_track1", points))
	
	points = [pos + Vector(2*adaptDist+2*iLength+iSize, 0.5*track).rot(ori)]
	points.append(points[-1] + Vector(-adaptDist, iWidth/2-track/2).rot(ori))
	points.append(points[-1] + Vector(-iLength, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iWidth).rot(ori))
	points.append(points[-1] + Vector(iLength, 0).rot(ori))
	points.append(points[-1] + Vector(adaptDist, iWidth/2-track/2).rot(ori))
	points.append(points[0])
	trackObjects.append(draw(iName+"_track2", points))

	points = [pos + Vector(0, 0.5*track+gap).rot(ori)]
	points.append(points[-1] + Vector(2*adaptDist+2*iLength+iSize, 0).rot(ori))
	points.append(points[-1] + Vector(0, -2*gap-track).rot(ori))
	points.append(points[-1] + Vector(-2*adaptDist-2*iLength-iSize, 0).rot(ori))
	points.append(points[0])
	gapObjects.append(draw(iName+"_gap", points))
	
	draw(iName+"_mesh", points)
	assignMeshLength(iName+"_mesh",0.01)
	
	points = [pos + Vector(adaptDist+iLength, +0.5*iWidth).rot(ori)]
	points.append(points[-1] + Vector(iSize, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iWidth).rot(ori))
	points.append(points[-1] + Vector(-iSize, 0).rot(ori))
	points.append(points[0])
	draw(iName+"_junc", points)
	assigneInduc(iName+"_junc", pos + Vector(adaptDist+iLength, 0).rot(ori), pos + Vector(adaptDist+iLength+iSize, 0).rot(ori), iInduct)
	
	return [retIn, retOut]

def drawSQUID(iName, iIn, iOut, iSize, iWidth, iLength, iSpace, iInduct=0.1):
	'''
	Draws a SQUID: two rectangles in parallel, each with boundary condition lumped RLC, L=iInduct
	
	Inputs:
	-------
	see drawJSJunc, additional inputs:
	iSpace: (float) distance between the two junctions
	
	Outputs:
	--------
	retIn: iIn with reversed direction
	retOut: calculated output port
	retFluxUp: port for top fast flux bias
	retFluxDown: port for bottom fast flux bias
	'''

	pos, ori = Vector(iIn[POS]), Vector(iIn[ORI])
	_, _, track, gap = iIn
	adaptDist = abs((2*iWidth+iSpace)/2-track/2)
	retIn = [pos, -ori, track, gap]
	retOut = [pos+abs(iLength)*ori, ori, track, gap]
	#retFluxDown = [pos+Vector(0.5*iLength, -gap-track/2).rot(ori), Vector(0, 1).rot(ori), None, None]
	#retFluxDown = [pos+abs(iLength)*ori, ori, track, gap]#[pos+0.5*abs(iLength)*ori+abs(gap+track/2)*ori.orth(), Vector(0, 1).rot(ori), None, None]
	#retFluxUp = [pos+abs(iLength)*ori, ori, track, gap] #[pos+0.5*abs(iLength)*ori-abs(gap+track/2)*ori.orth(), Vector(0, -1).rot(ori), None, None]
	
	retFluxDown = [pos+Vector(0.5*iLength, -gap-track/2).rot(ori), Vector(0, 1).rot(ori), track, gap]
	retFluxUp = [pos+Vector(0.5*iLength, gap+track/2).rot(ori), Vector(0, -1).rot(ori), track, gap]
	
	
	
	points = [pos + Vector(0, 0.5*track).rot(ori)]
	points.append(points[-1] + Vector(adaptDist, (2*iWidth+iSpace)/2-track/2).rot(ori))
	points.append(points[-1] + Vector(iLength/2-adaptDist-iSize/2, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iWidth).rot(ori))
	points.append(points[-1] + Vector(-iLength/2+adaptDist+iSize/2, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iSpace).rot(ori))
	points.append(points[-1] + Vector(iLength/2-adaptDist-iSize/2, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iWidth).rot(ori))
	points.append(points[-1] + Vector(-iLength/2+adaptDist+iSize/2, 0).rot(ori))
	points.append(points[-1] + Vector(-adaptDist, (2*iWidth+iSpace)/2-track/2).rot(ori))
	points.append(points[0])
	trackObjects.append(draw(iName+"_track1", points))
	
	points = [pos + Vector(iLength, 0.5*track).rot(ori)]
	points.append(points[-1] + Vector(-adaptDist, (2*iWidth+iSpace)/2-track/2).rot(ori))
	points.append(points[-1] + Vector(-iLength/2+adaptDist+iSize/2, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iWidth).rot(ori))
	points.append(points[-1] + Vector(iLength/2-adaptDist-iSize/2, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iSpace).rot(ori))
	points.append(points[-1] + Vector(-iLength/2+adaptDist+iSize/2, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iWidth).rot(ori))
	points.append(points[-1] + Vector(iLength/2-adaptDist-iSize/2, 0).rot(ori))
	points.append(points[-1] + Vector(adaptDist, (2*iWidth+iSpace)/2-track/2).rot(ori))
	points.append(points[0])
	trackObjects.append(draw(iName+"_track2", points))

	points = [pos + Vector(0, 0.5*track+gap).rot(ori)]
	points.append(points[-1] + Vector(iLength, 0).rot(ori))
	points.append(points[-1] + Vector(0, -2*gap-track).rot(ori))
	points.append(points[-1] + Vector(-iLength, 0).rot(ori))
	points.append(points[0])
	gapObjects.append(draw(iName+"_gap", points))
	
	draw(iName+"_mesh", points)
	assignMeshLength(iName+"_mesh",0.025)
	
	points = [pos + Vector(iLength/2-iSize/2, iWidth+0.5*iSpace).rot(ori)]
	points.append(points[-1] + Vector(iSize, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iWidth).rot(ori))
	points.append(points[-1] + Vector(-iSize, 0).rot(ori))
	points.append(points[0])
	draw(iName+"_junc1", points)
	assigneInduc(iName+"_junc1", pos + Vector(iLength/2-iSize/2, 0.5*iWidth+0.5*iSpace).rot(ori), pos + Vector(iLength/2+iSize/2, 0.5*iWidth+0.5*iSpace).rot(ori), iInduct)
	
	points = [pos + Vector(iLength/2-iSize/2, -0.5*iSpace).rot(ori)]
	points.append(points[-1] + Vector(iSize, 0).rot(ori))
	points.append(points[-1] + Vector(0, -iWidth).rot(ori))
	points.append(points[-1] + Vector(-iSize, 0).rot(ori))
	points.append(points[0])
	draw(iName+"_junc2", points)
	assigneInduc(iName+"_junc2", pos + Vector(iLength/2-iSize/2, -0.5*iWidth-0.5*iSpace).rot(ori), pos + Vector(iLength/2+iSize/2, -0.5*iWidth-0.5*iSpace).rot(ori), iInduct)
	
	return [retIn, retOut, retFluxUp, retFluxDown]


	
	
def drawFlux(iName, iIn, iSize, iMargin, iGndGap=None):
	'''
	Draws a flux bias line: loop joining center conductor and ground plane
	(not ground on both sides)
	
	Inputs:
	-------
	iName:
	iIn:
	iSize:
	iMargin:
	iGndGap:
	
	Outputs:
	--------
	retIn: iIn with reversed direction
	'''

	# pos, ori = Vector(FluxLine_In[POS]), Vector(FluxLine_In[ORI])
# FluxLine_In = [pos, ori, track/4., gap/4.]
# FluxLine_In=drawFlux("fluxline", FluxLine_In, 0.02, 0.02)
# # drawFlux(iName, iIn, iSize, iMargin, iGndGap=None):

# pos, ori = Vector(FluxLine_In[POS]), Vector(FluxLine_In[ORI])
# _, _, inTrack, inGap = FluxLine_In
# FluxLine_Out = [pos, ori, track, gap]

# FluxLine_In, FluxLine_Out = drawAdaptor("fluxline_adapt",FluxLine_In,FluxLine_Out)
# #drawAdaptor(iName, iIn, iOut, iSlope=1):
	
	
	
	
	pos, ori = Vector(iIn[POS]), Vector(iIn[ORI])
	_, _, track, gap = iIn
	pos = pos-(iMargin+iSize+2*track+gap)*ori
	retIn = [pos, -ori, track, gap]
	
	points = [pos + Vector(0, 0.5*track).rot(ori)]
	points.append(points[-1] + Vector(iMargin, 0).rot(ori))
	points.append(points[-1] + Vector(0, gap).rot(ori))
	points.append(points[-1] + Vector(iSize+2*track, 0).rot(ori))
	points.append(points[-1] + Vector(0, -2*track-2*gap).rot(ori))
	points.append(points[-1] + Vector(-2*track-iSize-iMargin, 0).rot(ori))
	points.append(points[-1] + Vector(0, track).rot(ori))
	points.append(points[-1] + Vector(iMargin+track+iSize, 0).rot(ori))
	points.append(points[-1] + Vector(0, 2*gap).rot(ori))
	points.append(points[-1] + Vector(-iSize, 0).rot(ori))
	points.append(points[-1] + Vector(0, -gap).rot(ori))
	points.append(points[-1] + Vector(-iMargin-track, 0).rot(ori))
	points.append(points[0])
	trackObjects.append(draw(iName+"_track1", points))
	
	
	points = [pos + Vector(0, 0.5*track+gap).rot(ori)]
	points.append(points[-1] + Vector(iMargin-gap, 0).rot(ori))
	points.append(points[-1] + Vector(0, gap).rot(ori))
	points.append(points[-1] + Vector(iSize+2*track+2*gap, 0).rot(ori))
	points.append(points[-1] + Vector(0, -2*track-4*gap).rot(ori))
	points.append(points[-1] + Vector(-2*track-iSize-iMargin-gap, 0).rot(ori))
	# points.append(points[-1] + Vector(0, track).rot(ori))
	# points.append(points[-1] + Vector(iMargin+track+iSize, 0).rot(ori))
	# points.append(points[-1] + Vector(0, 2*gap).rot(ori))
	# points.append(points[-1] + Vector(-iSize, 0).rot(ori))
	# points.append(points[-1] + Vector(0, -gap).rot(ori))
	# points.append(points[-1] + Vector(-iMargin-track, 0).rot(ori))
	points.append(points[0])
	gapObjects.append(draw(iName+"_gap1", points))
	
	
	return retIn

	
def drawBox(iName,center_vec,Width,Length,ori):
	draw(iName,[center_vec+Vector(Width/2., Length/2.).rot(ori),
				center_vec+Vector(-Width/2., Length/2.).rot(ori),
				center_vec+Vector(-Width/2., -Length/2.).rot(ori),
				center_vec+Vector(Width/2., -Length/2.).rot(ori),
				center_vec+Vector(Width/2., Length/2.).rot(ori)] )
	return iName
	
def drawIBMTansmon( iName,
					center_vec,
					cutout_Length,
					cutout_Width,
					pad_spacing,
					pad_Length,
					pad_Width,
					Jwidth,
					track,
					gap,
					Jinduc,
					ori,
					nport=1):

	gapObjects.append(drawBox(iName+"_cutout",center_vec,cutout_Width,cutout_Length,ori))
	trackObjects.append(drawBox(iName+"transmon_pad1",center_vec+Vector(pad_spacing/2.,0).rot(ori)+Vector(pad_Width/2.,0).rot(ori),pad_Width,pad_Length,ori))
	trackObjects.append(drawBox(iName+"transmon_pad2",center_vec-Vector(pad_spacing/2.,0).rot(ori)-Vector(pad_Width/2.,0).rot(ori),pad_Width,pad_Length,ori))
	
	drawBox(iName+"_transmon_mesh",center_vec,cutout_Width,cutout_Length,ori)
	assignMeshLength(iName+"_transmon_mesh",0.05)

	track_J=Jwidth*4.
	iJSin=[center_vec-Vector(pad_spacing/2.,0).rot(ori),ori,track_J,track_J]
	iJSout=[center_vec+Vector(pad_spacing/2.,0).rot(ori),ori,track_J,track_J]
	
	drawJSJunc(iName+"JJ", iJSin, iJSout, Jwidth, Jwidth, pad_spacing/2-track_J/2, iInduct=Jinduc)
				#(iName, iIn, iOut, iSize, iWidth, iLength, iInduct=0.1)
	#oDesktop.PauseScript(str(-chip_width))
	

	Vec_shift=Vector(0,pad_Length/2.)-Vector(0,track/2.)-Vector(0,gap)
	if nport==1:
		return [center_vec+Vector(cutout_Width/2.,0).rot(ori),ori,track,gap]
	elif nport==2:
		return ([center_vec+Vector(cutout_Width/2.,0).rot(ori),ori,track,gap],
				[center_vec-Vector(cutout_Width/2.,0).rot(ori),-ori,track,gap])
	elif nport==3:
		return ([center_vec+Vector(cutout_Width/2.,0).rot(ori)+Vec_shift.rot(ori),ori,track,gap],
				[center_vec+Vector(cutout_Width/2.,0).rot(ori)-Vec_shift.rot(ori),ori,track,gap],
				[center_vec+Vector(-cutout_Width/2.,0).rot(ori),-ori,track,gap])
	elif nport==4:
		return ([center_vec+Vector(cutout_Width/2.,0).rot(ori)+Vec_shift.rot(ori),ori,track,gap],
				[center_vec+Vector(cutout_Width/2.,0).rot(ori)-Vec_shift.rot(ori),ori,track,gap],
				[center_vec-Vector(cutout_Width/2.,0).rot(ori)+Vec_shift.rot(ori),-ori,track,gap],
				[center_vec-Vector(cutout_Width/2.,0).rot(ori)-Vec_shift.rot(ori),-ori,track,gap])
	elif nport==5:
		return ([center_vec+Vector( cutout_Width/2.,0).rot(ori), ori,track,gap],
				[center_vec+Vector(-cutout_Width/2.,0).rot(ori),-ori,track,gap],
				[center_vec+Vector( pad_spacing/2.+pad_Width+track/2.+gap,cutout_Length/2.).rot(ori), ori.orth(),track,gap])
		
def drawHalfCapa(iName, iIn, iLength, iWidth, iGap):
	'''
	Inputs:
	-------
	iName: string name of object
	iIn: (position, direction, track, gap) defines the input port
	iOut: (position, direction, track, gap) defines the output port
	       position and direction are None: this is calculated from
		   other parameters
	iLength: (float) length of pads
	iWidth: (float) width of pads
	
	Outputs:
	--------
	retIn: same as iIn, with flipped vector
	retOut: calculated output port to match all input dimensions
	
		igap iWidth
		     +--+  
		     |  |  
		+----+  | iLength
	iIn |       |	
		+----+  | 
		     |  | 
		     +--+ 
	'''

	pos, ori = Vector(iIn[POS]), -Vector(iIn[ORI])
	_, _, inTrack, _ = iIn
	
	points = [pos + Vector(iGap+iWidth, 0).rot(ori)]
	points.append(points[-1] + Vector(0, -iLength/2).rot(ori))
	points.append(points[-1] + Vector(-iWidth, 0).rot(ori))
	points.append(points[-1] + Vector(0, iLength/2-inTrack/2).rot(ori))
	points.append(points[-1] + Vector(-iGap, 0).rot(ori))
	points.append(points[-1] + Vector(0, inTrack).rot(ori))
	points.append(points[-1] + Vector(iGap, 0).rot(ori))
	points.append(points[-1] + Vector(0, iLength/2-inTrack/2).rot(ori))
	points.append(points[-1] + Vector(iWidth, 0).rot(ori))
	points.append(points[0])
	trackObjects.append(draw(iName+"_TransmonCapa", points))	
	
	CreateBondwire(iName+"_bondwire", iIn)


#def drawJSJunc(iName, iIn, iOut, iSize, iWidth, iLength, iInduct=0.1):
def createbox(iName,x0,y0,z0,DX,DY,DZ,material):
	oEditor.CreateBox(["NAME:BoxParameters",
						"CoordinateSystemID:=", -1,
						"XPosition:=", "("+str(x0)+")mm", 
						"YPosition:=", "("+str(y0)+")mm",
						"ZPosition:=", "("+str(z0)+")mm",
						"XSize:=", "("+str(DX)+")mm",
						"YSize:=", "("+str(DY)+")mm",
						"ZSize:=", "("+str(DZ)+")mm"],
						["NAME:Attributes",
						"Name:=", iName,
						"MaterialName:=",material,
						"PartCoordinateSystem:=", "Global",
						"Transparency:=", 0.75]
						)
	return iName


	
###########################################################
################# DRAWING STARTS HERE #####################
###########################################################



########## Define All Project Varaible ############


chip_width=LitExp("$chip_width", 8.67,VarDef=True)
chip_length=LitExp("$chip_length", 8.12,VarDef=True)
chip_thickness=LitExp("$chip_thickness", 0.28,VarDef=True)
pcb_thickness=LitExp("$pcb_thickness", 0.32,VarDef=True)


oEditor.SetModelUnits(["NAME:Units Parameter","Units:=", "mm","Rescale:=",False])







#oDesktop.PauseScript(str(-chip_width))
	
######### PCB, Chip and Ground Plane ###########
# chip dimensions (mm)
#chipLength, chipWidth = 8.12, 8.67
createbox("pcb","0mm","0mm",-chip_thickness,chip_width,chip_length,-pcb_thickness,"Rogers TMM 10i (tm) no loss")
createbox("chip","0mm","0mm","0mm",chip_width,chip_length,-chip_thickness,"silicon")
createbox("box","0mm","0mm",-chip_thickness-pcb_thickness,chip_width,chip_length,4*chip_thickness,"vacuum")

draw("Ground_plane",[Vector(0, 0),Vector(chip_width, 0),Vector(chip_width, chip_length),Vector(0, chip_length),Vector(0, 0)] )
assignPerfE(["Ground_plane"], "Ground_plane_PE")
assignPerfE(["box"], "box_PE")
#oDesktop.PauseScript(str(-chip_width))


######## Connectors #######
# posistions of connectors: distance (mm) from left or bottom edge 
con1, con2, con3, con4, con5 = 4.06, 3.01, 7.537, 3.01, 6.36
pcbTrack, pcbGap, boundLength, boundSlope = 0.3, 0.2, 0.2, 0.5
track, gap = 0.084, 0.05
track, gap = 0.042, 0.025


FluxIn, FluxOut = [[0, con1], [1, 0], pcbTrack, pcbGap], [None, None, track, gap]
WasteIn, WasteOut = [[con2, chip_length], [0, -1], pcbTrack, pcbGap], [None, None, track, gap]
QubitIn, QubitOut = [[con3, chip_length], [0, -1], pcbTrack, pcbGap], [None, None, track, gap]
BufferIn, BufferOut = [[con4, 0], [0, 1], pcbTrack, pcbGap], [None, None, track, gap]
memoryIn, memoryOut = [[con5, 0], [0, 1], pcbTrack, pcbGap], [None, None, track, gap]

BufferIn, BufferOut = drawConnector("Buffer_con", BufferIn, memoryOut, boundLength, boundSlope)
BufferOut=ShiftPortStraigth('Buffer_con_straight',BufferOut,0.2)

WasteIn, WasteOut = drawConnector("Waste_con", WasteIn, WasteOut, boundLength, boundSlope)
WasteOut=ShiftPortStraigth('Waste_con_straight',WasteOut,0.2)

QubitIn, QubitOut = drawConnector("Qubit_con", QubitIn, QubitOut, boundLength, boundSlope)
FluxIn, FluxOut = drawConnector("Flux_con", FluxIn, FluxOut, boundLength, boundSlope)
memoryIn, memoryOut = drawConnector("memory_con", memoryIn, memoryOut, boundLength, boundSlope)





Lj=LitExp("$Lj", 12,VarDef=True)
T_cutout_Length=LitExp("$T_cutout_Length", 1,VarDef=True)
T_cutout_Width=LitExp("$T_cutout_Width", 0.65,VarDef=True)
T_pad_Width=LitExp("$T_pad_Width", 0.1,VarDef=True)
T_pad_Spacing=LitExp("$T_pad_Spacing", 0.06,VarDef=True)
T_pad_Length=LitExp("$T_pad_Length", 0.75,VarDef=True)
T_offset_width=LitExp("$T_offset_width", 1.0,VarDef=True)
T_offset_length=LitExp("$T_offset_length", 0.0,VarDef=True)


TransmonIn_W,TransmonIn_B,TransmonIn_C=drawIBMTansmon(
																	"transmon",
																	Vector(chip_width/2+T_offset_width,chip_length/2+T_offset_length),		#center_vec
																	T_cutout_Length,				#cutout_Length
																	T_cutout_Width,				#cutout_Width
																	T_pad_Spacing,				#pad_Spacing
																	T_pad_Length,				#pad_Length
																	T_pad_Width,				#pad_Width
																	0.005,				#Jwidth
																	track,				#track
																	gap,				#gap
																	Lj,					#Jinduc
																	Vector(0.,1.),		#ori
																	5					#nport
																	)
																	

drawCable("Control_cable", QubitOut, TransmonIn_C,extrabonds=3)

																	
TW_capa_width=LitExp("$TW_capa_width", 2*track,VarDef=True)
TB_capa_width=LitExp("$TB_capa_width", 2*track,VarDef=True)
TC_capa_width=LitExp("$TC_capa_width", track+0.01,VarDef=True)

TW_capa_length=LitExp("$TW_capa_length", 0.05,VarDef=True)
TB_capa_length=LitExp("$TB_capa_length", 0.05,VarDef=True)
TC_capa_length=LitExp("$TC_capa_length", 0.005,VarDef=True)
												
drawHalfCapa("TransmonWasteCapa",TransmonIn_W,TW_capa_width,0.05,TW_capa_length)
drawHalfCapa("TransmonBufferCapa",TransmonIn_B,TB_capa_width,0.05,TB_capa_length)
drawHalfCapa("TransmonControlCapa",TransmonIn_C,TC_capa_width,0.05,TC_capa_length)


WasteCapaLength, WasteCapaSize = LitExp("$WasteCapa_length", 0.1,VarDef=True), LitExp("$WasteCapa_size", 0.02,VarDef=True)
BufferCapaLength, BufferCapaSize = LitExp("$BufferCapa_length", 0.1,VarDef=True), LitExp("$BufferCapa_size", 0.02,VarDef=True)

WasteCapaIn, WasteCapaOut = WasteOut, [None, None, track, gap]
BufferCapaIn, BufferCapaOut = BufferOut, [None, None, track, gap]

WasteCapaIn, WasteCapaOut = drawCapa("WasteCapa", WasteCapaIn, WasteCapaOut, WasteCapaLength, track/2, WasteCapaSize)
BufferCapaIn, BufferCapaOut = drawCapa("BufferCapa", BufferCapaIn, BufferCapaOut, BufferCapaLength, track/2, BufferCapaSize)

W_resLength=LitExp("$W_resLength", 9,VarDef=True)
W_resShift=LitExp("$W_resShift", 0.,VarDef=True)

B_resLength=LitExp("$B_resLength", 9.5,VarDef=True)
B_resShift=LitExp("$B_resShift", 0.,VarDef=True)




#drawWave2("BufferMeander", BufferCapaOut, TransmonIn_B, 3, 0.35,B_resShift, Lres=B_resLength, iStraight=0.3)
#Number of meanders, Shift, total length (shift is ignored if not None)


pos1, ori = Vector(BufferCapaOut[POS]), Vector(BufferCapaOut[ORI])
pos2 = Vector(TransmonIn_B[POS])
_, _, inTrack, inGap = BufferCapaOut
BufferSQUID_In = [0.5*(pos1+pos2), ori, inTrack, inGap]

LSQUID=LitExp("$LSQUID", 0.1,VarDef=True)
BufferSQUID_In,BufferSQUID_Out,FluxLine_In,_=drawSQUID("BufferSQUID", BufferSQUID_In,BufferSQUID_In,0.01,0.01,0.05,0.01,iInduct=LSQUID)
#(iName, iIn, iOut, iSize, iWidth, iLength, iSpace, iInduct=0.1)

pos, ori = Vector(FluxLine_In[POS]), Vector(FluxLine_In[ORI])
FluxLine_In = [pos, ori, track/4., gap/4.]
FluxLine_In=drawFlux("fluxline", FluxLine_In, 0.02, 0.02)
# drawFlux(iName, iIn, iSize, iMargin, iGndGap=None):

pos, ori = Vector(FluxLine_In[POS]), Vector(FluxLine_In[ORI])
_, _, inTrack, inGap = FluxLine_In
FluxLine_Out = [pos, ori, track, gap]

FluxLine_In, FluxLine_Out = drawAdaptor("fluxline_adapt",FluxLine_In,FluxLine_Out)
#drawAdaptor(iName, iIn, iOut, iSlope=1):

drawCable("Flux_cable", FluxLine_Out, memoryOut, extrabonds=3)

drawWave2("BufferMeander1", BufferCapaOut,BufferSQUID_In, 1, 0.35,B_resShift, Lres=B_resLength/2., iStraight=0.3)
drawWave2("BufferMeander2",  TransmonIn_B, BufferSQUID_Out, 1, 0.35,B_resShift, Lres=B_resLength/2., iStraight=0.3)

drawWave2("WasteMeander", WasteCapaOut, TransmonIn_W, 3, 0.35, W_resShift,Lres=W_resLength, iStraight=0.3) 	








######## Connectors #######
# posistions of connectors: distance (mm) from left or bottom edge 
# con1, con2, con3, con4, con5 = 4.06, 3.01, 7.537, 3.01, 6.36

# chipIn, chipOut = [[0, con1], [1, 0], pcbTrack, pcbGap], [None, None, track, gap]
# tomoIn, tomoOut = [[con2, chip_length], [0, -1], pcbTrack, pcbGap], [None, None, track, gap]
# qubitIn, qubitOut = [[con3, chip_length], [0, -1], pcbTrack, pcbGap], [None, None, track, gap]
# fluxIn, fluxOut = [[con4, 0], [0, 1], pcbTrack, pcbGap], [None, None, track, gap]
# memoryIn, memoryOut = [[con5, 0], [0, 1], pcbTrack, pcbGap], [None, None, track, gap]

# chipIn, chipOut = drawConnector("chip_con", chipIn, chipOut, boundLength, boundSlope)
# tomoIn, tomoOut = drawConnector("tomo_con", tomoIn, tomoOut, boundLength, boundSlope)
# qubitIn, qubitOut = drawConnector("qubit_con", qubitIn, qubitOut, boundLength, boundSlope)
# fluxIn, fluxOut = drawConnector("flux_con", fluxIn, fluxOut, boundLength, boundSlope)
# memoryIn, memoryOut = drawConnector("memory_con", memoryIn, memoryOut, boundLength, boundSlope)



# ######## Drawing of the Transmon #######
# # comment: this transmon is non-standard, it resembles in Xmon, and is
# # formed by two T junctions

# crossGap = LitExp("$cross_gap", 0.02)
# crossSize = LitExp("$cross_size", 0.3)

# # ports
# cross1In = [[memoryIn[POS].x-crossGap-track/2, chipIn[POS].y], [1, 0], track, crossGap]

# cross1OutRight, cross1OutLeft = [None, None, track, crossGap], [None, None, track, crossGap]
# cross1OutUp, cross1OutRight, cross1OutLeft = drawTriJunction("cross1", cross1In, cross1OutRight, cross1OutLeft)

# crossJSJuncIn, crossJSJuncOut = [cross1OutRight[POS], [0, 1], track, crossGap], [None, None, track, crossGap]
# crossJSJuncIn, crossJSJuncOut = drawJSJunc("memory_jsjunc", crossJSJuncIn, crossJSJuncOut, 0.02, 0.01, 0.02, 8.2)

# cross2In = [crossJSJuncOut[POS]+Vector(crossGap+track/2, crossGap+track/2), [-1, 0], track, crossGap]
# cross2OutLeft, cross2OutRight = [None, None, track, crossGap], [None, None, track, crossGap]
# cross2OutDown, cross2OutLeft, cross2OutRight = drawTriJunction("cross2", cross2In, cross2OutLeft, cross2OutRight)

# ####### BUFFER RESONATOR ##############
# # Dessin de la piste du buffer (buffer)

# chipCapaPos, chipWaveLength = LitExp("$chip_capa_pos", 0.6), LitExp("$chip_wave_length", 0.52)
# chipCapaLength, chipCapaSize = LitExp("$chip_capa_length", 0.13), LitExp("$chip_capa_size", 0.05)
# chipCapa2Length, chipCapa2Size = LitExp("$chip_capa2_length", 0.16), LitExp("$chip_capa2_size", 0.05)

# chipCapa1In, chipCapa1Out = [chipOut[POS]+Vector(chipCapaPos, 0), [1, 0], track, gap], [None, None, track, gap]
# chipCapa1In, chipCapa1Out = drawCapa("chip_capa1", chipCapa1In, chipCapa1Out, chipCapaLength, track/2, chipCapaSize)

# chipCapa2In, chipCapa2Out = [cross1OutUp[POS]+Vector(-crossSize, 0), [-1, 0], track, crossGap], [None, None, track, gap]
# chipCapa2In, chipCapa2Out = drawCapa("chip_capa2", chipCapa2In, chipCapa2Out, chipCapa2Length, track/2, chipCapa2Size)


# chipSquideIn, chipSquideOut = [chipCapa1Out[POS]*0.5+chipCapa2Out[POS]*0.5+Vector(-0.052, 0), [1, 0], track, gap], [None, None, track, gap]
# chipSquideIn, chipSquideOut, fluxFluxIn, _ = drawSquide("chip_squide", chipSquideIn, chipSquideOut, 0.02, 0.01, 0.02, 0.02, 2.62)


# chipWave1In, chipWave1Out = [chipCapa1Out[POS], [1, 0], track, gap], [None, None, track, gap]
# chipWave1In, chipWave1Out = drawWave("chip1_wave", chipWave1In, chipWave1Out, 1, chipWaveLength, chipWaveLength, 0.4)

# chipWave2In, chipWave2Out = [chipCapa2Out[POS], [-1, 0], track, gap], [None, None, track, gap]
# chipWave2In, chipWave2Out = drawWave("chip2_wave", chipWave2In, chipWave2Out, 1, chipWaveLength, chipWaveLength, 0.4)

# drawCable("chip_cable1", chipOut, chipCapa1In)
# drawCable("chip_cable2", chipWave1Out, chipSquideIn)
# drawCable("chip_cable3", chipSquideOut, chipWave2Out)
# drawCable("chip_cable4", chipCapa2In, cross1OutUp)

# ####### DIRECT INPUT LINE FOR THE TRANSMON #####
# #Dessin du la piste du transmon (transmon-out)

# qubitCapaLength, qubitCapaSize = LitExp("$qubit_capa_length", 0.11), LitExp("$qubit_capa_size", 0.31)

# qubitCapaIn, qubitCapaOut = [cross2OutDown[POS]+Vector(crossSize, 0), [1, 0], track, crossGap], [None, None, track, gap]
# qubitCapaIn, qubitCapaOut = drawCapa("qubit_capa", qubitCapaIn, qubitCapaOut, qubitCapaLength, track, qubitCapaSize)


# qubitDoubleJuncOut = [qubitCapaOut[POS]+Vector(0.5, 0.4), [0, 1], track, gap]
# qubitDoubleJuncIn = [qubitDoubleJuncOut[POS], [0, -1], track, gap]

# drawCable("qubit_cable1", qubitDoubleJuncOut, qubitOut)
# drawCable("qubit_cable2", qubitCapaOut, qubitDoubleJuncIn)
# drawCable("qubit_cable3", qubitCapaIn, cross2OutDown)

# ######## READOUT RESONATOR #############
# #Dessin de la piste de la tomographie (read-out)

# tomoDecal, tomoCapaPos = LitExp("$tomo_decal", 2.32), LitExp("$tomo_capa_pos", 0.2)
# tomoCapaLength, tomoCapaSize = LitExp("$tomo_capa_length", 0.13), LitExp("$tomo_capa_size", 0.05)
# tomoCapa2Length, tomoCapa2Size = LitExp("$tomo_capa2_length", 0.16), LitExp("$tomo_capa2_size", 0.05)

# tomoCapa1In, tomoCapa1Out = [tomoOut[POS]+Vector(0, -tomoCapaPos), [0, -1], track, gap], [None, None, track, gap]
# tomoCapa1In, tomoCapa1Out = drawCapa("tomo_capa1", tomoCapa1In, tomoCapa1Out, tomoCapaLength, track/2, tomoCapaSize)


# tomoCapa2In, tomoCapa2Out = [cross2OutRight[POS]+Vector(0, crossSize), [0, 1], track, crossGap], [None, None, track, gap]
# tomoCapa2In, tomoCapa2Out = drawCapa("tomo_capa2", tomoCapa2In, tomoCapa2Out, tomoCapa2Length, track/2, tomoCapa2Size)


# tomoDoubleJunc1In = [[tomoCapa1Out[POS].x, tomoCapa1Out[POS].y*0.70+tomoCapa2Out[POS].y*0.3]+Vector(-tomoDecal, 0), [0, 1], track, gap]
# tomoDoubleJunc1Out = [[tomoCapa1Out[POS].x, tomoCapa1Out[POS].y*0.70+tomoCapa2Out[POS].y*0.3]+Vector(-tomoDecal, 0), [0, -1], track, gap]

# tomoDoubleJunc2In = [[tomoCapa1Out[POS].x, tomoCapa1Out[POS].y*0.5+tomoCapa2Out[POS].y*0.5], [-1, 0], track, gap]
# tomoDoubleJunc2Out = [[tomoCapa1Out[POS].x, tomoCapa1Out[POS].y*0.5+tomoCapa2Out[POS].y*0.5], [1, 0], track, gap]

# drawCable("tomo_cable1", tomoOut, tomoCapa1In)
# drawCable("tomo_cable2", tomoCapa1Out, tomoDoubleJunc1In)
# drawCable("tomo_cable3", tomoDoubleJunc1Out, tomoDoubleJunc2In)
# drawCable("tomo_cable4", tomoDoubleJunc2Out, tomoCapa2Out)
# drawCable("tomo_cable5", tomoCapa2In, cross2OutRight)


# ######## MEMORY RESONATOR ###############
# #Dessin de la piste de la memoire (memory)

# memoryCapa1Pos, memoryWaveLength = LitExp("$memory_capa_pos", 0.2), LitExp("$memory_wave_length", 1.7)
# memoryCapaLength, memoryCapaSize = LitExp("$memory_capa_length", 0.05), LitExp("$memory_capa_size", 0.20)
# memoryCapa2Length, memoryCapa2Size = LitExp("$memory_capa2_length", 0.16), LitExp("$memory_capa2_size", 0.05)

# memoryCapa1In, memoryCapa1Out = [memoryOut[POS]+Vector(0, memoryCapa1Pos), [0, 1], track, gap], [None, None, track, gap]
# memoryCapa1In, memoryCapa1Out = drawCapa("memory_capa1", memoryCapa1In, memoryCapa1Out, memoryCapaLength, track/2, memoryCapaSize)

# memoryCapa2In, memoryCapa2Out = [cross1OutLeft[POS]+Vector(0, -crossSize), [0, -1], track, crossGap], [None, None, track, gap]
# memoryCapa2In, memoryCapa2Out = drawCapa("memory_capa2", memoryCapa2In, memoryCapa2Out, memoryCapa2Length, track/2, memoryCapa2Size)

# memoryWaveIn, memoryWaveOut = [memoryCapa1Out[POS], [0, 1], track, gap], [None, None, track, gap]
# memoryWaveIn, memoryWaveOut = drawWave("memory_wave", memoryWaveIn, memoryWaveOut, 2, memoryWaveLength, memoryWaveLength, 0.35)

# drawCable("memory_cable1", memoryOut, memoryCapa1In)
# drawCable("memory_cable2", memoryWaveOut, memoryCapa2Out)
# drawCable("memory_cable3", memoryCapa2In, cross1OutLeft)

# ######### FAST FLUX LINE ################
# #Dessin de la ligne d'application du flux (flux)

# fluxDecal = LitExp("$flux_decal", 0.0)

# fluxFluxIn[TRACK], fluxFluxIn[GAP] = 0.006, 0.022
# fluxFluxIn[POS] = fluxFluxIn[POS] + Vector(fluxDecal, -fluxFluxIn[TRACK])
# fluxFluxIn[ORI] = -fluxFluxIn[ORI]

# fluxDoubleJuncOut = [fluxFluxIn[POS]+Vector(0, -0.2), [0, 1], 0.006, 0.022]
# fluxDoubleJuncIn = [fluxFluxIn[POS]+Vector(0, -0.2), [0, -1], 0.006, 0.022]

# drawCable("flux_cable1", fluxOut, fluxDoubleJuncIn)
# drawCable("flux_cable2", fluxDoubleJuncOut, fluxFluxIn)

######## UNITE - SUBSTRACT - PerfE #####
#Finalisation globale du dessin

trackObject = fuse(trackObjects)
gapObject = fuse(gapObjects)

substract(gapObject, "Ground_plane")
assignPerfE([trackObject], "track")
assignPerfE(bondwireObjects, "bondwire")