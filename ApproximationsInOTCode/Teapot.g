Read("ReadInSTL.g");;
Read("IsInside.g");;
Read("ApproximateOT.g");;
Read("outer_hull.g");;
Read("triangulation.gd");;
Read("triangulation.gi");;

eps:=1.*10^-6;
# using a approximation of the Utah Teapot
data:=ReadSTL("Utah_teapot_solid");;
t:=TriangularComplexByVerticesInFacesNC(data[1]);;
points:=data[2];;



for i in [2..5] do
	data_tri:=ApproximateInTetroctahedrille(t,points*i,false);;
	DrawSTLScratchVOF(data_tri[1],Concatenation("Teapot_OT_times_",String(i)),data_tri[2]);;
	PrintTo(Concatenation("Teapot_OT_times_",String(i)),data_tri);;
od;

#for i in [1..10] do
#	step:=11-i;;
#	data_tri:=ApproximateInTetroctahedrille(t,points/step,false);;
#	DrawSTLScratchVOF(data_tri[1],Concatenation("Teapot_OT_",String(step)),data_tri[2]);;
#	PrintTo(Concatenation("Teapot_OT_",String(step)),data_tri);;
#od;