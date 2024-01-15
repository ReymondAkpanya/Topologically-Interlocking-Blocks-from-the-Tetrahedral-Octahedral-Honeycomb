
lambda:=function(t) return Concatenation("(Math.max(",String(t),",a)-a)/(1-a)"); end;
lambda:=function(t) return Concatenation("Math.pow(",String(t),",a)"); end;
lambda:=function(t) return Concatenation("Math.asinh(",String(t),"*(a))/Math.asinh(a)"); end;

DeformedTetrahedronString:=function(l,lambda,steps)
	local coordinates,faces,faces_middle_up,faces_middle_down,last_b,s,i,a,b,pr,non_interlocked_faces,list,FaceByVerticesOfFace,params;
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

	# give coordinates for first two layers
	coordinates:=[

	["0", "0", "0"],
    ["0", "-1", "0"],
    ["0", "-2", "0"],
    ["1", "-2", "0"],

    ["2", "-2", "0"],
    ["2", "-1", "0"],
    ["2", "0", "0"],
    ["1", "0", "0"],

	[lambda(l/steps),lambda(l/steps),String(l*Sqrt(2.)/steps)],
	[lambda(l/steps),Concatenation("-",lambda(l/steps)),String(l*Sqrt(2.)/steps)],
	[lambda(l/steps),Concatenation("-2+",lambda(l/steps)),String(l*Sqrt(2.)/steps)],
	[lambda(l/steps),Concatenation("-2-",lambda(l/steps)),String(l*Sqrt(2.)/steps)],

	[Concatenation("2-",lambda(l/steps)),Concatenation("-2-",lambda(l/steps)),String(l*Sqrt(2.)/steps)],
	[Concatenation("2-",lambda(l/steps)),Concatenation("-2+",lambda(l/steps)),String(l*Sqrt(2.)/steps)],
	[Concatenation("2-",lambda(l/steps)),Concatenation("-",lambda(l/steps)),String(l*Sqrt(2.)/steps)],
	[Concatenation("2-",lambda(l/steps)),lambda(l/steps),String(l*Sqrt(2.)/steps)],



	[Concatenation("-",lambda(l/steps)),Concatenation("-",lambda(l/steps)),String(-l*Sqrt(2.)/steps)],
	[Concatenation("-",lambda(l/steps)),Concatenation("-2+",lambda(l/steps)),String(-l*Sqrt(2.)/steps)],
	[lambda(l/steps),Concatenation("-2+",lambda(l/steps)),String(-l*Sqrt(2.)/steps)],
	[Concatenation("2-",lambda(l/steps)),Concatenation("-2+",lambda(l/steps)),String(-l*Sqrt(2.)/steps)],

	[Concatenation("2+",lambda(l/steps)),Concatenation("-2+",lambda(l/steps)),String(-l*Sqrt(2.)/steps)],
	[Concatenation("2+",lambda(l/steps)),Concatenation("-",lambda(l/steps)),String(-l*Sqrt(2.)/steps)],
	[Concatenation("2-",lambda(l/steps)),Concatenation("-",lambda(l/steps)),String(-l*Sqrt(2.)/steps)],
	[lambda(l/steps),Concatenation("-",lambda(l/steps)),String(-l*Sqrt(2.)/steps)]

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
	[lambda(i*l/steps),lambda(i*l/steps),String(i*l*Sqrt(2.)/steps)],
	[lambda(i*l/steps),Concatenation("-",lambda(i*l/steps)),String(i*l*Sqrt(2.)/steps)],
	[lambda(i*l/steps),"-1",String(i*l*Sqrt(2.)/steps)],
	[lambda(i*l/steps),Concatenation("-2+",lambda(i*l/steps)),String(i*l*Sqrt(2.)/steps)],
	
	[lambda(i*l/steps),Concatenation("-2-",lambda(i*l/steps)),String(i*l*Sqrt(2.)/steps)],
	["1",Concatenation("-2-",lambda(i*l/steps)),String(i*l*Sqrt(2.)/steps)],
	[Concatenation("2-",lambda(i*l/steps)),Concatenation("-2-",lambda(i*l/steps)),String(i*l*Sqrt(2.)/steps)],
	[Concatenation("2-",lambda(i*l/steps)),Concatenation("-2+",lambda(i*l/steps)),String(i*l*Sqrt(2.)/steps)],

	[Concatenation("2-",lambda(i*l/steps)),"-1",String(i*l*Sqrt(2.)/steps)],
	[Concatenation("2-",lambda(i*l/steps)),Concatenation("-",lambda(i*l/steps)),String(i*l*Sqrt(2.)/steps)],
	[Concatenation("2-",lambda(i*l/steps)),lambda(i*l/steps),String(i*l*Sqrt(2.)/steps)],
	["1",lambda(i*l/steps),String(i*l*Sqrt(2.)/steps)],


	[Concatenation("-",lambda(i*l/steps)),Concatenation("-",lambda(i*l/steps)),String(-i*l*Sqrt(2.)/steps)],
	[Concatenation("-",lambda(i*l/steps)),"-1",String(-i*l*Sqrt(2.)/steps)],
	[Concatenation("-",lambda(i*l/steps)),Concatenation("-2+",lambda(i*l/steps)),String(-i*l*Sqrt(2.)/steps)],
	[lambda(i*l/steps),Concatenation("-2+",lambda(i*l/steps)),String(-i*l*Sqrt(2.)/steps)],

	["1",Concatenation("-2+",lambda(i*l/steps)),String(-i*l*Sqrt(2.)/steps)],
	[Concatenation("2-",lambda(i*l/steps)),Concatenation("-2+",lambda(i*l/steps)),String(-i*l*Sqrt(2.)/steps)],
	[Concatenation("2+",lambda(i*l/steps)),Concatenation("-2+",lambda(i*l/steps)),String(-i*l*Sqrt(2.)/steps)],
	[Concatenation("2+",lambda(i*l/steps)),"-1",String(-i*l*Sqrt(2.)/steps)],
	
	[Concatenation("2+",lambda(i*l/steps)),Concatenation("-",lambda(i*l/steps)),String(-i*l*Sqrt(2.)/steps)],
	[Concatenation("2-",lambda(i*l/steps)),Concatenation("-",lambda(i*l/steps)),String(-i*l*Sqrt(2.)/steps)],
	["1",Concatenation("-",lambda(i*l/steps)),String(-i*l*Sqrt(2.)/steps)],
	[lambda(i*l/steps),Concatenation("-",lambda(i*l/steps)),String(-i*l*Sqrt(2.)/steps)]
			
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
			[lambda(i*l/steps),lambda(i*l/steps),String(Sqrt(2.)*i*l/steps)],
	[lambda(i*l/steps),Concatenation("-",lambda(i*l/steps)),String(Sqrt(2.)*i*l/steps)],
	[lambda(i*l/steps),Concatenation("-2+",lambda(i*l/steps)),String(Sqrt(2.)*i*l/steps)],
	[lambda(i*l/steps),Concatenation("-2-",lambda(i*l/steps)),String(Sqrt(2.)*i*l/steps)],

	[Concatenation("2-",lambda(i*l/steps)),Concatenation("-2-",lambda(i*l/steps)),String(Sqrt(2.)*i*l/steps)],
	[Concatenation("2-",lambda(i*l/steps)),Concatenation("-2+",lambda(i*l/steps)),String(Sqrt(2.)*i*l/steps)],
	[Concatenation("2-",lambda(i*l/steps)),Concatenation("-",lambda(i*l/steps)),String(Sqrt(2.)*i*l/steps)],
	[Concatenation("2-",lambda(i*l/steps)),lambda(i*l/steps),String(Sqrt(2.)*i*l/steps)],



	[Concatenation("-",lambda(i*l/steps)),Concatenation("-",lambda(i*l/steps)),String(-Sqrt(2.)*i*l/steps)],

	[Concatenation("-",lambda(i*l/steps)),Concatenation("-2+",lambda(i*l/steps)),String(-Sqrt(2.)*i*l/steps)],
	[lambda(i*l/steps),Concatenation("-2+",lambda(i*l/steps)),String(-Sqrt(2.)*i*l/steps)],
	[Concatenation("2-",lambda(i*l/steps)),Concatenation("-2+",lambda(i*l/steps)),String(-Sqrt(2.)*i*l/steps)],
	[Concatenation("2+",lambda(i*l/steps)),Concatenation("-2+",lambda(i*l/steps)),String(-Sqrt(2.)*i*l/steps)],

	[Concatenation("2+",lambda(i*l/steps)),Concatenation("-",lambda(i*l/steps)),String(-Sqrt(2.)*i*l/steps)],
	[Concatenation("2-",lambda(i*l/steps)),Concatenation("-",lambda(i*l/steps)),String(-Sqrt(2.)*i*l/steps)],
	[lambda(i*l/steps),Concatenation("-",lambda(i*l/steps)),String(-Sqrt(2.)*i*l/steps)]
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

	faces:=Concatenation(faces,[[18,21,27],[18,20,21],[18,19,20],[17,18,28],[18,27,28],[26,24,27],[21,24,27],[26,24,25],[22,23,24],[21,22,24],[ 7, 13, 16 ], [ 6, 7, 16 ], [ 5, 6, 16 ], [ 14, 15, 16 ], [ 13, 14, 16 ], [ 10, 12, 13 ], [ 7, 10, 13 ], 
  [ 10, 11, 12 ], [ 8, 9, 10 ], [ 7, 8, 10 ] ]+[40,40,40]*(steps-1)/2
	);;
	non_interlocked_faces:=Concatenation(non_interlocked_faces,[[18,21,27],[18,20,21],[18,19,20],[17,18,28],[18,27,28],[26,24,27],[21,24,27],[26,24,25],[22,23,24],[21,22,24],[ 7, 13, 16 ], [ 6, 7, 16 ], [ 5, 6, 16 ], [ 14, 15, 16 ], [ 13, 14, 16 ], [ 10, 12, 13 ], [ 7, 10, 13 ], 
  [ 10, 11, 12 ], [ 8, 9, 10 ], [ 7, 8, 10 ] ]+[40,40,40]*(steps-1)/2);;



	
	s:=TriangularComplexByVerticesInFaces(faces);;
	#coordinates:=coordinates*[[1,0,0],[0,1,0],[0,0,Sqrt(2.)]];
	pr := SetVertexCoordinatesParameterized(s, coordinates, rec());
	params := [["a", 0, [0,10]]];
	pr := SetVertexParameters(s, params, pr);
	#pr:=SetVertexCoordinates3D(s,coordinates);;  
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
	return [s,coordinates,pr];
end;;



#### draw interlocking in x and y direction
DeformedTetrahedronAssemblyString:=function(s,coordinates,x,y)
	local vof,new_vof,rot90,points,points4,translations,new_points,s4,pr4,i,j,p;
	# join 4 copies together for translation building block

	vof:=  VerticesOfFaces(s);;
	new_vof:=Concatenation(vof,vof+NumberOfVertices(s),vof+NumberOfVertices(s)*2,vof+NumberOfVertices(s)*3);;
	s:=SimplicialSurfaceByVerticesInFaces(new_vof);;
	rot90:=[[0,1,0],[-1,0,0],[0,0,1]];;
	points:=coordinates;

	#rot90^0
	points4:=StructuralCopy(points);
	#rot90^1
	for p in points do
		Add(points4,[Concatenation("-(",p[2],")"),p[1],p[3]]);
	od;
	#rot90^2
	for p in points do
		Add(points4,[Concatenation("-(",p[1],")"),Concatenation("-(",p[2],")"),p[3]]);
	od;
	#rot90^3
	for p in points do
		Add(points4,[p[2],Concatenation("-(",p[1],")"),p[3]]);
	od;




	translations:=[[4,0,0],[0,4,0]];
	new_vof:=Concatenation(List([1..x*y], x->VerticesOfFaces(s)+NumberOfVertices(s)*(x-1)));;
	new_points:=[];
	for i in [0..(x-1)] do
		for j in [0..(y-1)] do
			for p in points4 do
				Add(new_points,[Concatenation(p[1],"+",String(i),"*4"),Concatenation(p[2],"+",String(j),"*4"),p[3]]);
			od;
		od;
	od;
	s4:=SimplicialSurfaceByVerticesInFaces(new_vof);;
	points4:=new_points;
	pr4 := SetVertexCoordinatesParameterized(s4, points4, rec());
	params := [["a", 0, [0,10]]];
	pr4 := SetVertexParameters(s4, params, pr4);
	#pr4:=SetVertexCoordinates3D(s4,points4);;
	DrawComplexToJavaScript(s4,"DeformedTetAssembly",pr4);
end;;


#### draw interlocking in x and y direction
DeformedTetrahedronAssemblyStringRecord:=function(s,coordinates,pr,x,y)
	local vof,new_vof,rot90,points,points4,translations,new_points,s4,pr4,i,j,p;
	# join 4 copies together for translation building block

	vof:=  VerticesOfFaces(s);;
	new_vof:=Concatenation(vof,vof+NumberOfVertices(s),vof+NumberOfVertices(s)*2,vof+NumberOfVertices(s)*3);;
	s:=SimplicialSurfaceByVerticesInFaces(new_vof);;
	rot90:=[[0,1,0],[-1,0,0],[0,0,1]];;
	points:=coordinates;

	#rot90^0
	points4:=StructuralCopy(points);
	#rot90^1
	for p in points do
		Add(points4,[Concatenation("-(",p[2],")"),p[1],p[3]]);
	od;
	#rot90^2
	for p in points do
		Add(points4,[Concatenation("-(",p[1],")"),Concatenation("-(",p[2],")"),p[3]]);
	od;
	#rot90^3
	for p in points do
		Add(points4,[p[2],Concatenation("-(",p[1],")"),p[3]]);
	od;




	translations:=[[4,0,0],[0,4,0]];
	new_vof:=Concatenation(List([1..x*y], x->VerticesOfFaces(s)+NumberOfVertices(s)*(x-1)));;
	new_points:=[];
	for i in [0..(x-1)] do
		for j in [0..(y-1)] do
			for p in points4 do
				Add(new_points,[Concatenation(p[1],"+",String(i),"*4"),Concatenation(p[2],"+",String(j),"*4"),p[3]]);
			od;
		od;
	od;
	s4:=SimplicialSurfaceByVerticesInFaces(new_vof);;
	points4:=new_points;
	pr4 := SetVertexCoordinatesParameterized(s4, points4, rec());
	params := [["a", 0, [0,10]]];
	pr4 := SetVertexParameters(s4, params, pr4);
	pr4.faceColours:=Concatenation(List([1..x*y*4],i->pr.faceColours));
	#pr4:=SetVertexCoordinates3D(s4,points4);;
	DrawComplexToJavaScript(s4,"DeformedTetAssembly",pr4);
end;;
