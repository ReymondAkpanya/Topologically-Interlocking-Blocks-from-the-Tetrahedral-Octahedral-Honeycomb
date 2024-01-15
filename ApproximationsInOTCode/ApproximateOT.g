e1:=[1.,0.,0.];
	e2:=[0.,1.,0.];
	e3:=[0.,0.,1.];
	v1:=e2+e3;
	v2:=e1+e3;
	v3:=e1+e2;


Tetroctahedrille:=function(l,m,n)
	local e1,e2,e3,v1,v2,v3,oct,tet1,tet2,tet3,tet4,tets,res_octs,res_tets,i,j,k;
	e1:=[1.,0.,0.];
	e2:=[0.,1.,0.];
	e3:=[0.,0.,1.];
	v1:=e2+e3;
	v2:=e1+e3;
	v3:=e1+e2;
	oct:=[v1,v2,v3,v1+v3,v1+v2,v2+v3];
	tet1:=[[0.,0.,0.],v1,v2,v3];
	tet2:=[v2,v3,v2+v3,v2+v3-v1];
	tet3:=[v1,v2,v1+v2,v1+v2-v3];
	tet4:=[v2,v1+v2,2*v2,v2+v3];
	tets:=[tet1,tet2,tet3,tet4];

	res_tets:=[];
	res_octs:=[];
	for i in [-l..l] do
		for j in [-m..m] do
			for k in [-n..n] do
				res_tets:=Concatenation(res_tets,tets+i*v1+j*v2+k*v3);
				res_octs:=Concatenation(res_octs,[oct]+i*v1+j*v2+k*v3);
			od;
		od; 
	od;
	return [res_tets,res_octs];
end;;


## new approach

## check for each vertex which translational cell it lies in
WhichCell:=function(point)
	local x_trans,y_trans,z_trans;	

	x_trans:=Int(point[1]/2);
	y_trans:=Int(point[2]/2);
	z_trans:=Int(point[3]/2);
	if x_trans<0 then
		x_trans:=x_trans;
	fi;
	if y_trans<0 then
		y_trans:=y_trans;
	fi;
	if z_trans<0 then
		z_trans:=z_trans;
	fi;
	return x_trans*[-1,1,1]+y_trans*[1,-1,1]+z_trans*[1,1,-1];
end;;



StrPercent:=function(f,digits)
	local pow,n,a,b;
	f:=Float(f);
	pow:=10^digits;
	n:=Int(f*pow*100);
	a:=QuoInt(n,pow);
	b:=n-a*pow;
	return Concatenation(String(a),".",String(b),"%");
end;

ApproximateInTetroctahedrille:=function(t,points,only_oct)

	check_points:=[];
	for f in Faces(t) do
		vof:=points{VerticesOfFace(t,f)};
		length1:=MyNorm(vof[2]-vof[1]);
		length2:=MyNorm(vof[3]-vof[1]);
		length3:=MyNorm(vof[3]-vof[2]);
		parts1:=Int(length1/Sqrt(2.));
		parts2:=Int(length2/Sqrt(2.));
		parts3:=Int(length3/Sqrt(2.));
		parts:=Maximum([parts1,parts2,parts3,1]);
		for i in [1..parts] do
			for j in [i..parts] do
				Add(check_points,Float((parts-i-j)/parts)*vof[1]+Float(i/parts)*vof[2]+Float(j/parts)*vof[3]);
			od;
		od;
	od;
	check_points:=Concatenation(check_points,points);;
	cells:=List(check_points,p->WhichCell(p));;


	cells:=Set(cells);
	cells_region:=[];
	# we want to check regions around it
	for i in [-1..1] do
		for j in [-1..1] do
			for k in [-1..1] do
				cells_region:=Concatenation(cells_region,cells+[i,j,k]);
			od;
		od;
	od;
	cells_region:=Set(cells_region);;
	#Print("Check ",cells_region, "cells!\n");
	if only_oct then
		coordinates_assemblies:=[];
		faces_assemblies:=[];
		honey:=Tetroctahedrille(0,0,0);
		# go over all possible cells and check if centers lie inside the t
		for c in [1..Size(cells_region)] do
			Print("\r",StrPercent(c/Size(cells_region),3));
			c:=cells_region[c];
			cell:=honey+v1*c[1]+v2*c[2]+v3*c[3];
			oct:=cell[2][1];
			if IsInside(t,points,Sum(oct)/6,eps) then
				faces_assemblies:=Concatenation(faces_assemblies,VerticesOfFaces(Octahedron())+Size(coordinates_assemblies));
				coordinates_assemblies:=Concatenation(coordinates_assemblies,oct);
			fi;
		od;
	else
		coordinates_assemblies:=[];
		faces_assemblies:=[];
		honey:=Tetroctahedrille(0,0,0);
		# go over all possible cells and check if centers lie inside the t
		for c in [1..Size(cells_region)] do
			Print("\r",StrPercent(c/Size(cells_region),3));
			c:=cells_region[c];
			cell:=honey+v1*c[1]+v2*c[2]+v3*c[3];
			for tet in cell[1] do
				if IsInside(t,points,Sum(tet)/4,eps) then
					faces_assemblies:=Concatenation(faces_assemblies,VerticesOfFaces(Tetrahedron())+Size(coordinates_assemblies));
					coordinates_assemblies:=Concatenation(coordinates_assemblies,tet);
				fi;
			od;
			oct:=cell[2][1];
			if IsInside(t,points,Sum(oct)/6,eps) then
				faces_assemblies:=Concatenation(faces_assemblies,VerticesOfFaces(Octahedron())+Size(coordinates_assemblies));
				coordinates_assemblies:=Concatenation(coordinates_assemblies,oct);
			fi;
		od;
	fi;
	return [faces_assemblies,coordinates_assemblies];
end;;
	

ApproximateInTetroctahedrilleSimple:=function(t,points,only_oct)

	check_points:=[];
	for f in Faces(t) do
		vof:=points{VerticesOfFace(t,f)};
		length1:=MyNorm(vof[2]-vof[1]);
		length2:=MyNorm(vof[3]-vof[1]);
		length3:=MyNorm(vof[3]-vof[2]);
		parts1:=Int(length1/Sqrt(2.));
		parts2:=Int(length2/Sqrt(2.));
		parts3:=Int(length3/Sqrt(2.));
		parts:=Maximum([parts1,parts2,parts3,1]);
		for i in [1..parts] do
			for j in [i..parts] do
				Add(check_points,Float((parts-i-j)/parts)*vof[1]+Float(i/parts)*vof[2]+Float(j/parts)*vof[3]);
			od;
		od;
	od;
	check_points:=Concatenation(check_points,points);;
	cells:=List(check_points,p->WhichCell(p));;
	cells_region:=Set(cells);;
	#Print("Check ",cells_region, "cells!\n");
	if only_oct then
		coordinates_assemblies:=[];
		faces_assemblies:=[];
		honey:=Tetroctahedrille(0,0,0);
		# go over all possible cells and check if centers lie inside the t
		for c in [1..Size(cells_region)] do
			Print("\r",StrPercent(c/Size(cells_region),3));
			c:=cells_region[c];
			cell:=honey+v1*c[1]+v2*c[2]+v3*c[3];
			oct:=cell[2][1];
			faces_assemblies:=Concatenation(faces_assemblies,VerticesOfFaces(Octahedron())+Size(coordinates_assemblies));
			coordinates_assemblies:=Concatenation(coordinates_assemblies,oct);
		od;
	else
		coordinates_assemblies:=[];
		faces_assemblies:=[];
		honey:=Tetroctahedrille(0,0,0);
		# go over all possible cells and check if centers lie inside the t
		for c in [1..Size(cells_region)] do
			Print("\r",StrPercent(c/Size(cells_region),3));
			c:=cells_region[c];
			cell:=honey+v1*c[1]+v2*c[2]+v3*c[3];
			for tet in cell[1] do
					faces_assemblies:=Concatenation(faces_assemblies,VerticesOfFaces(Tetrahedron())+Size(coordinates_assemblies));
					coordinates_assemblies:=Concatenation(coordinates_assemblies,tet);
			od;
			oct:=cell[2][1];
				faces_assemblies:=Concatenation(faces_assemblies,VerticesOfFaces(Octahedron())+Size(coordinates_assemblies));
				coordinates_assemblies:=Concatenation(coordinates_assemblies,oct);
		od;
	fi;
	return [faces_assemblies,coordinates_assemblies];
end;;

