# Dot Product
Dot:=function(a,b)
	return a[1]*b[1]+a[2]*b[2]+a[3]*b[3];
end;;

MyNorm:=function(a)
	return Sqrt(Dot(a,a));
end;

MyRoundVector:=function(n,eps)
	local i;
	n:=1.0*n;
	for i in [1..3] do 
		if Sqrt(n[i]^2)<eps then
			n[i]:=0.;
		fi;
	od;
	return n;
end;

VectorAngle:=function(a,b)
	return Acos(Dot(a,b)/(MyNorm(a)*MyNorm(b)));
end;;

VectorAnglePlane:=function(a,b,n)
	return Atan2(Determinant(1.*[a,b,n]),Dot(a,b));
end;;

Rotate_Vector_Axis:=function(v,u,alpha)
	local r_alpha_u;
	u:=u/MyNorm(u);
	r_alpha_u:= Cos(alpha)*IdentityMat(3) + Sin(alpha)*[[0,-u[3],u[2]],[u[3],0,-u[1]],[-u[2],u[1],0]] + (1- Cos(alpha)) *(TransposedMat([u])*[u]);
	return (v*TransposedMat(r_alpha_u));
end;;



# Given an edge e, from a triangular complex s, with vertexCoordinates vC and a face f with MyNormal vector n
# calculate the fan [(f,n),(f2,n2),...] in direction of MyNormal
CalculateFan:=function(s,e,vC,f,MyNormal)
	local FacesOfEdge_e, ThirdPoint, VoE, a, p, n, t, vec, face; 
	# for all faces with edge e in s arrange them first
	FacesOfEdge_e:=FacesOfEdge(s,e);
	ThirdPoint:=[];
	VoE:=VerticesOfEdge(s,e);
	for face in FacesOfEdge_e do
		Add(ThirdPoint,[face,DifferenceLists(VerticesOfFace(s,face),VoE)[1]]);
	od;
	a:=vC[VoE[1]];
	n:=vC[VoE[2]]-vC[VoE[1]];
	n:=n/MyNorm(n);
	for t in ThirdPoint do
		p:=vC[t[2]];
		t[3]:=p-a-Dot((p-a),n)*n;
		# if we consider face f save the vector
		if t[1]=f then
			vec:=t[3];
		fi;
	od;
	# now compute the angle of all derived vectors with vec and rotate MyNormal vector n according to the angle
	for t in ThirdPoint do
		if Determinant(1.*[vec,n,MyNormal]) < 0. then
			t[4]:=VectorAnglePlane(vec,t[3],n);
			t[5]:=Rotate_Vector_Axis(MyNormal,n,t[4]);
		else
			t[4]:=VectorAnglePlane(vec,t[3],-n);
			t[5]:=Rotate_Vector_Axis(MyNormal,-n,t[4]);
		fi;
	od;
	SortBy(ThirdPoint,x->x[4]);
	return ThirdPoint;
end;;

UpwardContinuation:=function(s,e,vC,f,MyNormal)
	local Fan, index, i;
	Fan:=CalculateFan(s,e,vC,f,MyNormal);
	for i in [1..Size(Fan)] do 
		if Fan[i][1]=f then
			index:=(i mod (Size(Fan))) + 1;
			return [Fan[index][1],-Fan[index][5]];
		fi;
	od;
	return fail;
end;;


# Algorithm from Paper by M. Attene Title: As-exact-as-possible repair of unprintable stl files
# Input: a triangular complex s with coordinates vC and an outer face f with normal vector n
# Output: s, Restricted Complex to outer faces, Outer Faces, and correctly oriented normal vectors
OuterHull:=function(s,vC,f,n)
	local B, e, OuterTriangles, InnerTriangles, NormalVectors, b, edge, t, t_new, eps;
	eps := 1.0/(10^(12));
	B:=[];
	for e in EdgesOfFace(s,f) do
		Add(B,[f,e]);
	od;
	OuterTriangles:=[f];
	NormalVectors:=[];
	NormalVectors[f]:=n;
	InnerTriangles:=[];
	while not IsEmpty(B) do 
		t:=Remove(B);
		t_new:=UpwardContinuation(s,t[2],vC,t[1],NormalVectors[t[1]]);
		if not t_new[1] in OuterTriangles then

			Add(OuterTriangles,t_new[1]);
			NormalVectors[t_new[1]]:=MyRoundVector(t_new[2], eps);
			for edge in EdgesOfFace(s,t_new[1]) do
				if edge <> t[2] then
					Add(B,[t_new[1],edge]);
				fi;
			od;
		fi;
	od;
	return [s,SubcomplexByFaces(s,OuterTriangles),OuterTriangles,NormalVectors];
end;;

DrawSTLScratch:=function(t,fileName,vC)
	local f,normals,normal,x,y,ccoords;
	normals:=[];
	for f in Faces(t) do
		ccoords:=vC{VerticesOfFace(t,f)};
		x := ccoords[2]-ccoords[1];
        y := ccoords[3]-ccoords[1];
        normal := Crossproduct(x,y);
        normal := normal / Sqrt(normal*normal);
        normals[f] := normal;
    od;
    DrawSTLwithNormals(t,fileName,vC,normals,[]);
end;;

# kleine aenderung
DrawSTLwithNormals:=function(s,fileName, vC,normals,visualize_normal_list)
        local file, filesepr, name, output, x,y, i, j, k, r,l, coords, Copy_Coords, normal, eps, VoF, f,
        		middle, v2, v3, new_coords, new_normal, perms, perm,n;

        #######################################################################################################################################################################
        #   INPUTS
        ##
        # This method takes a string and a list l in the coordinate format (l=[face1,face2,,face3,....]). The faces are also lists in the format 
        # face1 = [vertex1,vertex2,vertex3,normal,vertexNumbs]. The vertices and normal then are lists as well: vertex1 = [x,y,z] where x,y,z are floats
        # and vertexNumbs = [v1,v2,v3] where v1,v2,v3 are integers corresponding to the vertex indices in the simplicial face. Since those change quite a lot when fixing self-
        # intersections, vertexNumbs is !depricated!
        #
        ##
        #   METHOD
        ##
        # The method iself saves a STL file corresponding to the object with name specified by fileName
        #
        #######################################################################################################################################################################

        eps := 1.0/(10^(12));

        filesepr := SplitString(fileName, ".");
        name := filesepr[1];
        # test file name
        file := Filename( DirectoryCurrent(), Concatenation(name,".stl") );
        output := OutputTextFile( file, false ); # override other files
            
        if output = fail then
            Error(Concatenation("File ", String(file), " can't be opened.") );
        fi;        

        AppendTo(output, Concatenation("solid ", name, "\n"));
        
        for f in Faces(s) do 
            
                VoF:=VerticesOfFace(s,f);
                # get coords of vertices
                coords := [vC[VoF[1]],vC[VoF[2]],vC[VoF[3]]];
                normal := normals[f];
                
                # write normal

                AppendTo(output, "\tfacet normal ");
                for j in [1..3] do
                    AppendTo(output, Concatenation(String(normal[j])," "));
                od;
                AppendTo(output, "\n");
                    
                # write vertex coords according to the right hand rule
                perms:=[[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]];
                for perm in perms do
                	n:=Crossproduct(coords[perm[2]]-coords[perm[1]],coords[perm[3]]-coords[perm[1]]);
                	n:=-n/Norm2(n);
                	if FlVEq(n,normal,eps) then

                		break;
                	fi;
                od;

                AppendTo(output, "\t\touter loop\n");
                for j in [1..3] do
                    AppendTo(output,"\t\t\tvertex ");
                       
                    for k in [1..3] do
                        AppendTo(output, Concatenation(String(coords[perm[j]][k])," "));
                    od;
                    AppendTo(output,"\n");
                od;
                AppendTo(output, "\t\tendloop\n");
                AppendTo(output,"\tendfacet\n");
                
                # visualize normal of face
                if f in visualize_normal_list then
                middle := 1./3 * coords[1] + 1./3 * coords[2] + 1./3 * coords[3];
                v2 := middle + 0.02 * (coords[1]-coords[2])/Norm2(coords[1]-coords[2]);
                v3 := middle + 0.1 * normal;
                
                new_coords := [middle,v2,v3];
                new_normal := Crossproduct(middle-v2, middle - v3);
                new_normal := new_normal / Norm2(new_normal);
                
                AppendTo(output, "\tfacet normal ");
                
                for j in [1..3] do
                    AppendTo(output, Concatenation(String(new_normal[j])," "));
                od;
                AppendTo(output, "\n");
                
                 # write vertex coords
                AppendTo(output, "\t\touter loop\n");
                for j in [1..3] do
                    AppendTo(output,"\t\t\tvertex ");
                       
                    for k in [1..3] do
                        AppendTo(output, Concatenation(String(new_coords[j][k])," "));
                    od;
                    AppendTo(output,"\n");
                od;
                AppendTo(output, "\t\tendloop\n");
                AppendTo(output,"\tendfacet\n");
                fi;
        od;
        AppendTo(output, Concatenation("endsolid ", name));
        Print("\n Saved file");
        CloseStream(output);
        return;
    end;;


# kleine aenderung
DrawSurfaceToObj:=function(s,fileName, vC)
        local file, filesepr, name, output, x,y, i, j, k, r,l, coords, Copy_Coords, normal, eps, VoF, f,
        		middle, v2, v3, new_coords, new_normal,v;

        
        

        filesepr := SplitString(fileName, ".");
        name := filesepr[1];
        # test file name
        file := Filename( DirectoryCurrent(), Concatenation(name,".obj") );
        output := OutputTextFile( file, false ); # override other files
            
        if output = fail then
            Error(Concatenation("File ", String(file), " can't be opened.") );
        fi;        

        #AppendTo(output, Concatenation("solid ", name, "/n"));

        for v in vC do
        	AppendTo(output,"v ");
        	for i in [1..3] do
        		AppendTo(output,v[i]);
        		if i<3 then
        			AppendTo(output," ");
        		fi;
        	od;
        	AppendTo(output,"\n");
        od;
        
        for f in VerticesOfFaces(s) do 
        	AppendTo(output,"f ");
        	for i in [1..3] do
        		AppendTo(output,f[i]);
        		if i<3 then
        			AppendTo(output," ");
        		fi;
        	od;
        	AppendTo(output,"\n");            
        od;
        Print("\n Saved file");
        CloseStream(output);
        return;
end;;




DrawSTLScratchVOF:=function(vof,fileName,vC)
	local f,normals,normal,x,y,ccoords,i;
	normals:=[];
	for i in [1..Size(vof)] do
		ccoords:=vC{vof[i]};
		x := ccoords[2]-ccoords[1];
        y := ccoords[3]-ccoords[1];
        normal := Crossproduct(x,y);
        normal := normal / Sqrt(normal*normal);
        normals[i] := normal;
    od;
    DrawSTLwithNormalsVOF(vof,fileName,vC,normals,[]);
end;;

# kleine aenderung
DrawSTLwithNormalsVOF:=function(vof,fileName, vC,normals,visualize_normal_list)
        local file, filesepr, name, output, x,y, i, j, k, r,l, coords, Copy_Coords, normal, eps, VoF, f,
        		middle, v2, v3, new_coords, new_normal, perms, perm,n;

        #######################################################################################################################################################################
        #   INPUTS
        ##
        # This method takes a string and a list l in the coordinate format (l=[face1,face2,,face3,....]). The faces are also lists in the format 
        # face1 = [vertex1,vertex2,vertex3,normal,vertexNumbs]. The vertices and normal then are lists as well: vertex1 = [x,y,z] where x,y,z are floats
        # and vertexNumbs = [v1,v2,v3] where v1,v2,v3 are integers corresponding to the vertex indices in the simplicial face. Since those change quite a lot when fixing self-
        # intersections, vertexNumbs is !depricated!
        #
        ##
        #   METHOD
        ##
        # The method iself saves a STL file corresponding to the object with name specified by fileName
        #
        #######################################################################################################################################################################

        eps := 1.0/(10^(12));

        filesepr := SplitString(fileName, ".");
        name := filesepr[1];
        # test file name
        file := Filename( DirectoryCurrent(), Concatenation(name,".stl") );
        output := OutputTextFile( file, false ); # override other files
            
        if output = fail then
            Error(Concatenation("File ", String(file), " can't be opened.") );
        fi;        

        AppendTo(output, Concatenation("solid ", name, "\n"));
        
        for f in [1..Size(vof)] do 
            
                VoF:=vof[f];
                # get coords of vertices
                coords := [vC[VoF[1]],vC[VoF[2]],vC[VoF[3]]];
                normal := normals[f];
                
                # write normal

                AppendTo(output, "\tfacet normal ");
                for j in [1..3] do
                    AppendTo(output, Concatenation(String(normal[j])," "));
                od;
                AppendTo(output, "\n");
                    
                # write vertex coords according to the right hand rule
                perms:=[[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]];
                for perm in perms do
                	n:=Crossproduct(coords[perm[2]]-coords[perm[1]],coords[perm[3]]-coords[perm[1]]);
                	n:=-n/Norm2(n);
                	if FlVEq(n,normal,eps) then

                		break;
                	fi;
                od;

                AppendTo(output, "\t\touter loop\n");
                for j in [1..3] do
                    AppendTo(output,"\t\t\tvertex ");
                       
                    for k in [1..3] do
                        AppendTo(output, Concatenation(String(coords[perm[j]][k])," "));
                    od;
                    AppendTo(output,"\n");
                od;
                AppendTo(output, "\t\tendloop\n");
                AppendTo(output,"\tendfacet\n");
                
                # visualize normal of face
                if f in visualize_normal_list then
                middle := 1./3 * coords[1] + 1./3 * coords[2] + 1./3 * coords[3];
                v2 := middle + 0.02 * (coords[1]-coords[2])/Norm2(coords[1]-coords[2]);
                v3 := middle + 0.1 * normal;
                
                new_coords := [middle,v2,v3];
                new_normal := Crossproduct(middle-v2, middle - v3);
                new_normal := new_normal / Norm2(new_normal);
                
                AppendTo(output, "\tfacet normal ");
                
                for j in [1..3] do
                    AppendTo(output, Concatenation(String(new_normal[j])," "));
                od;
                AppendTo(output, "\n");
                
                 # write vertex coords
                AppendTo(output, "\t\touter loop\n");
                for j in [1..3] do
                    AppendTo(output,"\t\t\tvertex ");
                       
                    for k in [1..3] do
                        AppendTo(output, Concatenation(String(new_coords[j][k])," "));
                    od;
                    AppendTo(output,"\n");
                od;
                AppendTo(output, "\t\tendloop\n");
                AppendTo(output,"\tendfacet\n");
                fi;
        od;
        AppendTo(output, Concatenation("endsolid ", name));
        Print("\n Saved file");
        CloseStream(output);
        return;
    end;;



