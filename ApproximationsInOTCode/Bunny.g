Read("ReadInSTL.g");;
Read("IsInside.g");;
Read("ApproximateOT.g");;
Read("outer_hull.g");;
Read("triangulation.gd");;
Read("triangulation.gi");;

eps:=1.*10^-6;
# using a 10% approximation of the Standford Bunny
data:=ReadSTL("Stanford_Bunny_1000");;
t:=SimplicialSurfaceByVerticesInFacesNC(data[1]);;
points:=data[2];;

for i in [1..10] do
	step:=11-i;;
	data_tri:=ApproximateInTetroctahedrille(t,points/step,false);;
	DrawSTLScratchVOF(data_tri[1],Concatenation("Bunny_OT_",String(step)),data_tri[2]);;
	PrintTo(Concatenation("Bunny_OT_",String(step)),data_tri);;
od;