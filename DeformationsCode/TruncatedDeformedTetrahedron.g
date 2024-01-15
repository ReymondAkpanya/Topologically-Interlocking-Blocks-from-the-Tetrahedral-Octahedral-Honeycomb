



	#pi:=3.1415926535;
	#lambda:=function(t) t:=Float(t); t:=t*2*pi; return 3./(4*pi)*(2./3*(t-Sin(t))); end;
	#f:=function(t) 	if Float(t)>0. then return 2.71828^(1/-Float(t)); else return 0; fi; end;
	#g:=function(t) t:=Float(t); return f(t)/(f(t)+f(1-t)); end;
	#lambda:=function(t) return l*g(t/l); end;
	#lambda:=function(t) return Sqrt(1.*t); end;
	#l:=0.999;
	#only allow even number of steps
	#steps:=8;



LoadPackage("GapIC");

# Create Deformed and truncated Tetrahedron
# Truncation with factor l
# lambda gives deformation with approximation steps given by steps (only even)
# red faces indicate non-contact faces
# grey faces indicate contact faces
# Example Runs:
# DeformedTetrahedron(0.8,t->t^2,10);;
# DeformedTetrahedron(0.99,t->Sqrt(1.*t),10);;
# You can create an assembly of this block with 
# DeformedTetrahedronAssembly(last[1],last[2],2,2);
DeformedTetrahedron:=function(l,lambda,steps)
	local coordinates,faces,faces_middle_up,faces_middle_down,last_b,s,i,a,b,pr,non_interlocked_faces,list,FaceByVerticesOfFace,TriangleArea,Dot,MyNorm,CrossProduct;
	if steps mod 2 = 1 then
		Print("Only even number of steps!\n");
		return fail;
	fi;	


	FaceByVerticesOfFace:=function(s,vof)
		local vofs, i;
		vofs:=VerticesOfFaces(s);
		for i in [1..Size(vofs)] do
			if Set(vofs[i])=Set(vof) then
				return i;
			fi;
		od;
	end;;
	Dot:=function(a,b)
		return a[1]*b[1]+a[2]*b[2]+a[3]*b[3];
	end;;

	MyNorm:=function(a)
		return Sqrt(Dot(a,a));
	end;

	CrossProduct := function(v1, v2)
	    local x1, y1, z1, x2, y2, z2;
	    
	    x1 := v1[1];
	    y1 := v1[2];
	    z1 := v1[3];
	    
	    x2 := v2[1];
	    y2 := v2[2];
	    z2 := v2[3];
	    
	    return [y1 * z2 - z1 * y2, z1 * x2 - x1 * z2, x1 * y2 - y1 * x2];
	end;

	TriangleArea := function(triangle)
	    local p1, p2, p3, v1, v2, cross_product, area;
	    if Length(triangle) <> 3 then
	        Error("Input must be a list of three 3D points.");
	    fi;
	    
	    p1 := triangle[1];
	    p2 := triangle[2];
	    p3 := triangle[3];
	    
	    v1 := p2 - p1;
	    v2 := p3 - p1;
	    
	    cross_product := CrossProduct(v1, v2); # Vector cross product
	    area := 0.5 * MyNorm(cross_product);
	    
	    return area;
	end;


	# give coordinates for first two layers
	coordinates:=[

	[0,0,0],
	[0,-1,0],
	[0,-2,0],
	[1,-2,0],

	[2,-2,0],
	[2,-1,0],
	[2,0,0],
	[1,0,0],

	[lambda(l/steps),lambda(l/steps),(l/steps)],
	[lambda(l/steps),-lambda(l/steps),(l/steps)],
	[lambda(l/steps),-2+lambda(l/steps),(l/steps)],
	[lambda(l/steps),-2-lambda(l/steps),(l/steps)],

	[2-lambda(l/steps),-2-lambda(l/steps),(l/steps)],
	[2-lambda(l/steps),-2+lambda(l/steps),(l/steps)],
	[2-lambda(l/steps),-lambda(l/steps),(l/steps)],
	[2-lambda(l/steps),lambda(l/steps),(l/steps)],



	[-lambda(l/steps),-lambda(l/steps),-(l/steps)],
	[-lambda(l/steps),-2+lambda(l/steps),-(l/steps)],
	[lambda(l/steps),-2+lambda(l/steps),-(l/steps)],
	[2-lambda(l/steps),-2+lambda(l/steps),-(l/steps)],

	[2+lambda(l/steps),-2+lambda(l/steps),-(l/steps)],
	[2+lambda(l/steps),-lambda(l/steps),-(l/steps)],
	[2-lambda(l/steps),-lambda(l/steps),-(l/steps)],
	[lambda(l/steps),-lambda(l/steps),-(l/steps)]

	];
	# give faces by relabeling the vertices as they were first described
	faces_middle_up:=[
	[8,16,7],
	[8,9,16],
	[8,9,1],

	[1,9,10],
	[1,2,10],
	[2,10,11],
	[2,3,11],
	[11,12,3],

	[3,4,12],
	[4,12,13],
	[4,13,5],

	[13,14,5],
	[5,6,14],
	[14,15,6],
	[6,7,15],
	[15,16,7]
	];;

	faces_middle_down:=[
	[1,2,17],
	[17,18,2],
	[2,3,18],

	[18,19,3],
	[3,4,19],
	[19,20,4],
	[4,5,20],
	[20,21,5],

	[5,6,21],
	[21,22,6],
	[6,7,22],

	[22,23,7],
	[7,8,23],
	[23,24,8],
	[1,8,24],
	[24,17,1]

	];;
	last_b:=8;
	faces:=Concatenation(faces_middle_up,faces_middle_down);
	non_interlocked_faces:=[[1,9,10],[11,12,3],[24,17,1],[22,23,7],[15,16,7],[13,14,5],[20,21,5],[18,19,3]];;
	# now add layer by layer (the number of steps is the number of layers)
	for i in [2..steps] do
		# new coordinates
		if i mod 2=0 then
			# this is the even case
			# introduce  more coordinates and faces here
			# coordinates
			a:=last_b;
			b:=Size(coordinates);
			last_b:=b;
			coordinates:=Concatenation(coordinates,[
	[lambda(i*l/steps),lambda(i*l/steps),(i*l/steps)],
	[lambda(i*l/steps),-lambda(i*l/steps),(i*l/steps)],
	[lambda(i*l/steps),-1,(i*l/steps)],
	[lambda(i*l/steps),-2+lambda(i*l/steps),(i*l/steps)],
	
	[lambda(i*l/steps),-2-lambda(i*l/steps),(i*l/steps)],
	[1,-2-lambda(i*l/steps),(i*l/steps)],
	[2-lambda(i*l/steps),-2-lambda(i*l/steps),(i*l/steps)],
	[2-lambda(i*l/steps),-2+lambda(i*l/steps),(i*l/steps)],

	[2-lambda(i*l/steps),-1,(i*l/steps)],
	[2-lambda(i*l/steps),-lambda(i*l/steps),(i*l/steps)],
	[2-lambda(i*l/steps),lambda(i*l/steps),(i*l/steps)],
	[1,lambda(i*l/steps),(i*l/steps)],


	[-lambda(i*l/steps),-lambda(i*l/steps),-(i*l/steps)],
	[-lambda(i*l/steps),-1,-(i*l/steps)],
	[-lambda(i*l/steps),-2+lambda(i*l/steps),-(i*l/steps)],
	[lambda(i*l/steps),-2+lambda(i*l/steps),-(i*l/steps)],

	[1,-2+lambda(i*l/steps),-(i*l/steps)],
	[2-lambda(i*l/steps),-2+lambda(i*l/steps),-(i*l/steps)],
	[2+lambda(i*l/steps),-2+lambda(i*l/steps),-(i*l/steps)],
	[2+lambda(i*l/steps),-1,-(i*l/steps)],

	[2+lambda(i*l/steps),-lambda(i*l/steps),-(i*l/steps)],
	[2-lambda(i*l/steps),-lambda(i*l/steps),-(i*l/steps)],
	[1,-lambda(i*l/steps),-(i*l/steps)],
	[lambda(i*l/steps),-lambda(i*l/steps),-(i*l/steps)]
			
			]);
			# new faces
			faces:=Concatenation(faces,[
			[1+a,1+b,12+b],
			[1+a,8+a,12+b],
			[12+b,11+b,8+a],
			
			[1+b,2+b,1+a],
			[1+a,2+a,2+b],
			[2+b,3+b,2+a],
			[2+a,3+a,3+b],
			[3+b,4+b,3+a],
			[3+a,4+a,4+b],
			[4+b,5+b,4+a],
			
			[4+a,5+b,6+b],
			[4+a,5+a,6+b],
			[6+b,7+b,5+a],
			
			[7+b,8+b,5+a],
			[5+a,6+a,8+b],
			[8+b,9+b,6+a],
			[6+a,7+a,9+b],
			[9+b,10+b,7+a],
			[7+a,8+a,10+b],
			[10+b,11+b,8+a]
			]);
			non_interlocked_faces:=Concatenation(non_interlocked_faces, [
					[1+b,2+b,1+a],
					[1+a,2+a,2+b],
					[3+a,4+a,4+b],
					[4+b,5+b,4+a],
					[7+b,8+b,5+a],
					[5+a,6+a,8+b],
					[7+a,8+a,10+b],
					[10+b,11+b,8+a]
			]);
			a:=a;
			b:=b+4;
			faces:=Concatenation(faces,[
			[9+a,9+b,10+b],
			[9+a,10+a,10+b],
			[10+b,11+b,10+a],

			[11+b,12+b,10+a],
			[10+a,11+a,12+b],
			[12+b,13+b,11+a],
			[11+a,12+a,13+b],
			[13+b,14+b,12+a],
			[12+a,13+a,14+b],
			[14+b,15+b,13+a],

			[15+b,16+b,13+a],
			[13+a,14+a,16+b],
			[16+b,17+b,14+a],

			[17+b,18+b,14+a],
			[14+a,15+a,18+b],
			[18+b,19+b,15+a],
			[15+a,16+a,19+b],
			[19+b,20+b,16+a],
			[16+a,9+a,20+b],
			[9+b,20+b,9+a],
			]);
			non_interlocked_faces:=Concatenation(non_interlocked_faces,[
					[11+b,12+b,10+a],
					[10+a,11+a,12+b],
					[12+a,13+a,14+b],
					[14+b,15+b,13+a],
					[17+b,18+b,14+a],
					[14+a,15+a,18+b],
					[16+a,9+a,20+b],
					[9+b,20+b,9+a]
			]);
		else
			# coordinates
			b:=last_b;
			a:=Size(coordinates);
			last_b:=a;
			coordinates:=Concatenation(coordinates,[
			[lambda(i*l/steps),lambda(i*l/steps),(i*l/steps)],
	[lambda(i*l/steps),-lambda(i*l/steps),(i*l/steps)],
	[lambda(i*l/steps),-2+lambda(i*l/steps),(i*l/steps)],
	[lambda(i*l/steps),-2-lambda(i*l/steps),(i*l/steps)],

	[2-lambda(i*l/steps),-2-lambda(i*l/steps),(i*l/steps)],
	[2-lambda(i*l/steps),-2+lambda(i*l/steps),(i*l/steps)],
	[2-lambda(i*l/steps),-lambda(i*l/steps),(i*l/steps)],
	[2-lambda(i*l/steps),lambda(i*l/steps),(i*l/steps)],



	[-lambda(i*l/steps),-lambda(i*l/steps),-(i*l/steps)],

	[-lambda(i*l/steps),-2+lambda(i*l/steps),-(i*l/steps)],
	[lambda(i*l/steps),-2+lambda(i*l/steps),-(i*l/steps)],
	[2-lambda(i*l/steps),-2+lambda(i*l/steps),-(i*l/steps)],
	[2+lambda(i*l/steps),-2+lambda(i*l/steps),-(i*l/steps)],

	[2+lambda(i*l/steps),-lambda(i*l/steps),-(i*l/steps)],
	[2-lambda(i*l/steps),-lambda(i*l/steps),-(i*l/steps)],
	[lambda(i*l/steps),-lambda(i*l/steps),-(i*l/steps)]
	]);
		# new faces
			# new faces
			faces:=Concatenation(faces,[
			[1+a,1+b,12+b],
			[1+a,8+a,12+b],
			[12+b,11+b,8+a],
			
			[1+b,2+b,1+a],
			[1+a,2+a,2+b],
			[2+b,3+b,2+a],
			[2+a,3+a,3+b],
			[3+b,4+b,3+a],
			[3+a,4+a,4+b],
			[4+b,5+b,4+a],
			
			[4+a,5+b,6+b],
			[4+a,5+a,6+b],
			[6+b,7+b,5+a],
			
			[7+b,8+b,5+a],
			[5+a,6+a,8+b],
			[8+b,9+b,6+a],
			[6+a,7+a,9+b],
			[9+b,10+b,7+a],
			[7+a,8+a,10+b],
			[10+b,11+b,8+a]
			]);
			non_interlocked_faces:=Concatenation(non_interlocked_faces, [
					[1+b,2+b,1+a],
					[1+a,2+a,2+b],
					[3+a,4+a,4+b],
					[4+b,5+b,4+a],
					[7+b,8+b,5+a],
					[5+a,6+a,8+b],
					[7+a,8+a,10+b],
					[10+b,11+b,8+a]
			]);
			a:=a;
			b:=b+4;
			faces:=Concatenation(faces,[
			[9+a,9+b,10+b],
			[9+a,10+a,10+b],
			[10+b,11+b,10+a],

			[11+b,12+b,10+a],
			[10+a,11+a,12+b],
			[12+b,13+b,11+a],
			[11+a,12+a,13+b],
			[13+b,14+b,12+a],
			[12+a,13+a,14+b],
			[14+b,15+b,13+a],

			[15+b,16+b,13+a],
			[13+a,14+a,16+b],
			[16+b,17+b,14+a],

			[17+b,18+b,14+a],
			[14+a,15+a,18+b],
			[18+b,19+b,15+a],
			[15+a,16+a,19+b],
			[19+b,20+b,16+a],
			[16+a,9+a,20+b],
			[9+b,20+b,9+a]
			]);
			non_interlocked_faces:=Concatenation(non_interlocked_faces,[
					[11+b,12+b,10+a],
					[10+a,11+a,12+b],
					[12+a,13+a,14+b],
					[14+b,15+b,13+a],
					[17+b,18+b,14+a],
					[14+a,15+a,18+b],
					[16+a,9+a,20+b],
					[9+b,20+b,9+a]
			]);
		fi;
	od;;
	coordinates:=coordinates*[[1.,0,0],[0,1.,0],[0,0,Sqrt(2.)]];


	# Print the ratio of boundary faces (non-contact contact ) faces
	s:=TriangularComplexByVerticesInFaces(faces);;
	contact_area:=0.;
	non_contact_area:=0.;
	list:=List(non_interlocked_faces,i->FaceByVerticesOfFace(s,i));;
	for i in Faces(s) do 
		if i in list then
			non_contact_area:=non_contact_area+TriangleArea(coordinates{VerticesOfFace(s,i)});
		else
			contact_area:=contact_area+TriangleArea(coordinates{VerticesOfFace(s,i)});
		fi;
	od;
	Print("Non contact area: ",non_contact_area,"\n");
	Print("contact area: ",contact_area,"\n");

	# add top and bottom faces
	faces:=Concatenation(faces,[[18,21,27],[18,20,21],[18,19,20],[17,18,28],[18,27,28],[26,24,27],[21,24,27],[26,24,25],[22,23,24],[21,22,24],[ 7, 13, 16 ], [ 6, 7, 16 ], [ 5, 6, 16 ], [ 14, 15, 16 ], [ 13, 14, 16 ], [ 10, 12, 13 ], [ 7, 10, 13 ], 
  [ 10, 11, 12 ], [ 8, 9, 10 ], [ 7, 8, 10 ] ]+[40,40,40]*(steps-1)/2
	);;
	non_interlocked_faces:=Concatenation(non_interlocked_faces,[[18,21,27],[18,20,21],[18,19,20],[17,18,28],[18,27,28],[26,24,27],[21,24,27],[26,24,25],[22,23,24],[21,22,24],[ 7, 13, 16 ], [ 6, 7, 16 ], [ 5, 6, 16 ], [ 14, 15, 16 ], [ 13, 14, 16 ], [ 10, 12, 13 ], [ 7, 10, 13 ], 
  [ 10, 11, 12 ], [ 8, 9, 10 ], [ 7, 8, 10 ] ]+[40,40,40]*(steps-1)/2);;



	
	s:=TriangularComplexByVerticesInFaces(faces);;
	pr:=SetVertexCoordinates3D(s,coordinates);;  
	list:=List(non_interlocked_faces,i->FaceByVerticesOfFace(s,i));;
	for i in Faces(s) do
		if i in list then
			SetFaceColour(s,i,"0xEA050D",pr);;
		else
			SetFaceColour(s,i,"0x757573",pr);;
		fi;
	od;
	#pr:= DeactivateVertices(s,pr);;
	#Deactivate Edges Labels
	#pr:= DeactivateEdges(s,pr);;
	DrawComplexToJavaScript(s,"test",pr);;
	return [s,coordinates];
end;;

#### draw interlocking in x and y direction
DeformedTetrahedronAssembly:=function(s,coordinates,x,y)
	local vof,new_vof,rot90,points,points4,translations,new_points,s4,pr4,i,j;
	# join 4 copies together for translation building block

	vof:=  VerticesOfFaces(s);;
	new_vof:=Concatenation(vof,vof+NumberOfVertices(s),vof+NumberOfVertices(s)*2,vof+NumberOfVertices(s)*3);;
	s:=SimplicialSurfaceByVerticesInFaces(new_vof);;
	rot90:=[[0,1,0],[-1,0,0],[0,0,1]];;
	points:=coordinates;
	points4:=Concatenation(points,points*rot90,points*rot90^2,points*rot90^3);



	translations:=[[4,0,0],[0,4,0]];
	new_vof:=Concatenation(List([1..x*y], x->VerticesOfFaces(s)+NumberOfVertices(s)*(x-1)));;
	new_points:=[];
	for i in [0..(x-1)] do
		for j in [0..(y-1)] do
			new_points:=Concatenation(new_points,List(points4,point->point+i*translations[1]+j*translations[2]));
		od;
	od;
	s4:=SimplicialSurfaceByVerticesInFaces(new_vof);;
	points4:=new_points;
	pr4:=SetVertexCoordinates3D(s4,points4);;
	DrawComplexToJavaScript(s4,"DeformedTetAssembly",pr4);
end;;
