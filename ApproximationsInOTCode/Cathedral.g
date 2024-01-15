Read("ReadInSTL.g");;
Read("IsInside.g");;
Read("ApproximateOT.g");;
Read("outer_hull.g");;
Read("triangulation.gd");;
Read("triangulation.gi");;

eps:=1.*10^-6;
# using a approximation of the Aachen Cathedral
data:=ReadSTL("Dom_Aachen");;
Dom:=TriangularComplexByVerticesInFacesNC(data[1]);;
Dom_points:=data[2];;
t:=Dom;
points:=Dom_points;

for i in [1..10] do
	step:=11-i;;
	data_tri:=ApproximateInTetroctahedrille(t,points/step,false);;
	DrawSTLScratchVOF(data_tri[1],Concatenation("Dom_OT_",String(step)),data_tri[2]);;
	PrintTo(Concatenation("Dom_OT_",String(step)),data_tri);;
od;