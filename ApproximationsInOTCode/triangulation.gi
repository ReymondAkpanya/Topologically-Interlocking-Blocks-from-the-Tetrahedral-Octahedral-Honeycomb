InstallMethod( Signum,
    [IsFloat],
    function(x)
        if x > 0. then
            return 1.;
        fi;
        if x = 0. then
            return 0.;
        else 
            return -1.;
        fi;
    end
);

InstallMethod( Norm2,
    [IsList],
    function(x)
        # 3-D euklid. norm
        return Sqrt(x[1]^2+x[2]^2+x[3]^2);
    end
);

InstallMethod( SurfaceTriangle,
    [IsList,IsList,IsList],
    function(x,y,z)
    	local a,b,c;
        # surface of triangle
        a := Norm2(x-y);
        b := Norm2(x-z);
        c := Norm2(y-z);
        
        if FlEq((a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c),0.,1./(10^9)) then
        	# bc of numerical errors, expression can be negative instead of 0
        	return 0.;
        else
        	return 0.25 * CubeRoot( (a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c) )^(3./2.);
        fi;
    end
);

InstallMethod( FlEq,
    [IsFloat,IsFloat,IsFloat],
    function(x,y,epsilon)
        # check if x is equal to y with a numerical error margin, equivalent to being both FlLeq and FlGeq
        if x > y-epsilon then
            if x < y + epsilon then
                return true;
            fi;
        fi;
        return false;
    end
);

InstallMethod( FlVEq,
    [IsList,IsList,IsFloat],
    function(x,y,epsilon)
        # check if List x is equal to List y with a numerical error margin, equivalent to being both FlLeq and FlGeq
        if x[1] > y[1]-epsilon and x[1] < y[1] + epsilon then
            if x[2] > y[2]-epsilon and x[2] < y[2] + epsilon then
                if x[3] > y[3]-epsilon and x[3] < y[3] + epsilon then
                    return true;
                fi;
            fi;
        fi;
        return false;
    end
);

InstallMethod( NumericalPosition,
    [IsList,IsList,IsFloat],
    function(list,entry,epsilon)
        local i, n;

        n := Length(list);
        i := 1;
        
        while i <= n do
            if FlVEq(list[i],entry, epsilon) then
                return i;
            fi;
            i := i + 1;
        od;
        
        return fail;
    end
);

InstallMethod( FlLeq,
    [IsFloat,IsFloat,IsFloat],
    function(x,y,epsilon)
        # check if x is less or equal to y with a numerical error margin
        if x <= y+epsilon then
            return true;
        fi;
        return false;
    end
);
InstallMethod( FlGeq,
    [IsFloat,IsFloat,IsFloat],
    function(x,y,epsilon)
        # check if x is greater or equal to y with a numerical error margin
        if y <= x+epsilon then
            return true;
        fi;
        return false;
    end
);

InstallMethod( Crossproduct,
    [IsList, IsList],
    function(x,y)
        return [x[2]*y[3]-x[3]*y[2],x[3]*y[1]-x[1]*y[3],x[1]*y[2]-x[2]*y[1]];
    end
);

InstallMethod( NumericalUniqueListOfLists,
    [IsList, IsFloat],
    function(list, eps)
        local n, I, i, unique;
        
        n := Length(list);
        I := [2..n];
        unique := [];
        unique[1] := list[1];
        
        for i in I do
            if ForAll(unique, x-> not FlVEq(x,list[i],eps)) then
                unique[Length(unique)+1] := list[i];
            fi;
        od;
        return unique;
    end
);



InstallMethod( PointsInOneLine,
    [IsList, IsList, IsList, IsFloat],
    function(v1, v2, v3, eps)
        local area;
        
        area := SurfaceTriangle(v1,v2,v3);
        
        return AbsoluteValue(area)<eps;
    end
);


InstallMethod( CheckForMistakes,
   [IsList, IsList, IsFloat],
   function(verts,I,eps)
   	local check_verts, r;
   	
   	for r in I do
        	check_verts := NumericalUniqueListOfLists([verts[r][1],verts[r][2],verts[r][3]],eps);
        	#if Length(check_verts) < 3 then
        	if Length(check_verts) < 3 or PointsInOneLine(check_verts[1],check_verts[2],check_verts[3],1./(10^11)) then
                    return [r,true];
        	fi;
        od;
   	return [0,false];
   end
);

InstallMethod( ProjectLineOnTriangle,
    [IsList, IsList, IsFloat, IsList,IsList, IsFloat],
    function(x,l,d,normal,points, eps)
        local lambda, p, orthog, coords, o, o_true, i;
        if not FlEq(d - x*normal,0.,eps) then
            # the point is not on the triangle plane already
            lambda := (d - x*normal) / (l*normal);
            
        elif FlEq(d - x*normal,0.,eps) and not FlEq(l*normal,0.,eps) then
            # the edge is not orthogonal to the triangle and the point already on the triangle plane
            lambda := 0.;
        else
            # then the edge is orthogonal to the triangle and the point already on the triangle plane
            orthog := points[1];
            coords := points[2];
            o := [(coords[1]*orthog[1] - x*orthog[1]) / (l*orthog[1]), (coords[2]*orthog[2] - x*orthog[2]) / (l*orthog[2]), (coords[3]*orthog[3] - x*orthog[3]) / (l*orthog[3])];
            o_true := [];
            for i in [1..3] do
                if FlGeq(o[i],0.,eps) then
                    o_true[i] := o[i];
                else
                    o_true[i] := Maximum(o)+ 1.;
                fi;
            od;
            lambda := Minimum(o_true);
            if FlGeq(lambda,1.,eps) then
                # if the face is fully inside the triangle then we always need to triangulate
                lambda := 1.;
            fi;
        fi;
        p := x+lambda*l;
        return [p,lambda];
    end
);

InstallMethod( RemoveDuplicateFaces,
    [IsList,IsBool,IsFloat],
    function(Coords,write,eps)
        local r, l, Copy_Coords, coords, test_coords,w, t;
        Copy_Coords := [];
        t := 1;
        w := 0;

        for r in [1..Length(Coords)] do
            if IsBound(Coords[r]) and Coords[r] <> [] then
                Copy_Coords[t] := ShallowCopy(Coords[r]);
                t := t + 1;
            fi;
        od;
        
        for r in [1..Length(Copy_Coords)] do
            if IsBound(Copy_Coords[r]) and Copy_Coords[r] <> [] then
                coords := [Copy_Coords[r][1],Copy_Coords[r][2],Copy_Coords[r][3]];
                for l in [1..Length(Copy_Coords)] do
                    if IsBound(Copy_Coords[l]) and Copy_Coords[l] <> [] and l <> r then
                        test_coords := [Copy_Coords[l][1],Copy_Coords[l][2],Copy_Coords[l][3]];
                        if Length(NumericalUniqueListOfLists([coords[1],coords[2],coords[3],test_coords[1],test_coords[2],test_coords[3]],eps)) < 4 then
                            Copy_Coords[l] := [];
                            w := w + 1;
                        fi;
                    fi;
                od;
            fi;
        od;
        if write then
		Print("\n","\n",Length(Coords), " starting length\n");
		Print(Length(Coords)-t+1, " empty faces\n");
		Print(w, " duplicate faces\n");
	fi;
        return Copy_Coords;
    end
);

InstallMethod( PlaneEquation,
    [IsList,IsList,IsList],
    function(point,normal,v1)
        return point*normal - v1*normal;
    end
);

InstallMethod( InPlane,
    [IsList,IsList,IsList,IsFloat],
    function(point,normal,v1,eps)
        return FlEq(PlaneEquation(point,normal,v1),0.,eps);
    end
);

InstallMethod( OnEdges,
    [IsList,IsList,IsFloat],
    function(verts,x,eps)
        local a1,a2,a3, o, orthog, normal;
        
        a1 := verts[1];
        a2 := verts[2];
        a3 := verts[3];

        normal := Crossproduct(a2-a1,a3-a1)/ Norm2 (Crossproduct(a2-a1,a3-a1));
        orthog := [Crossproduct(a2-a1,normal),Crossproduct(a3-a2,normal),Crossproduct(a1-a3,normal)];
                 
        orthog[1] := orthog[1] / Sqrt(orthog[1]*orthog[1]);
        orthog[2] := orthog[2] / Sqrt(orthog[2]*orthog[2]);
        orthog[3] := orthog[3] / Sqrt(orthog[3]*orthog[3]);
                        
        # must be right of all the planes, need to have the orthogonal vectors point to the inside of the triangle
        orthog[1] := orthog[1] * ((orthog[1]*a3-orthog[1]*a1) / AbsoluteValue(orthog[1]*a3-orthog[1]*a1));
        orthog[2] := orthog[2] * ((orthog[2]*a1-orthog[2]*a2) / AbsoluteValue(orthog[2]*a1-orthog[2]*a2));
        orthog[3] := orthog[3] * ((orthog[3]*a2-orthog[3]*a3) / AbsoluteValue(orthog[3]*a2-orthog[3]*a3));
                

        o := [orthog[1]*x-orthog[1]*a1, orthog[2]*x-orthog[2]*a2, orthog[3]*x-orthog[3]*a3];
        
        return [FlEq(o[1],0.,eps),FlEq(o[2],0.,eps),FlEq(o[3],0.,eps)];        
    end
);

InstallMethod( PointInTriangle,
    [IsList, IsList, IsList, IsList, IsFloat],
    function(a1, a2, a3, x, eps)
        local orthog, normal, o, m;
        # checks if point x (already on triangle plane) is inside the triangle
        normal := Crossproduct(a2-a1,a3-a1);
        orthog := [Crossproduct(a2-a1,normal),Crossproduct(a3-a2,normal),Crossproduct(a1-a3,normal)];
                 
        orthog[1] := orthog[1] / Sqrt(orthog[1]*orthog[1]);
        orthog[2] := orthog[2] / Sqrt(orthog[2]*orthog[2]);
        orthog[3] := orthog[3] / Sqrt(orthog[3]*orthog[3]);
                        
        # must be right of all the planes, need to have the orthogonal vectors point to the inside of the triangle
        orthog[1] := orthog[1] * ((orthog[1]*a3-orthog[1]*a1) / AbsoluteValue(orthog[1]*a3-orthog[1]*a1));
        orthog[2] := orthog[2] * ((orthog[2]*a1-orthog[2]*a2) / AbsoluteValue(orthog[2]*a1-orthog[2]*a2));
        orthog[3] := orthog[3] * ((orthog[3]*a2-orthog[3]*a3) / AbsoluteValue(orthog[3]*a2-orthog[3]*a3));
                

        o := [orthog[1]*x-orthog[1]*a1, orthog[2]*x-orthog[2]*a2, orthog[3]*x-orthog[3]*a3];
        
        m := Minimum(o[1],o[2],o[3]);
        return FlGeq(m,0.,eps);
    end
);

InstallMethod( TriangleFullyInside,
    [IsList, IsList, IsFloat],
    function(target, current, eps)
        local c_normal;
        c_normal := Crossproduct(current[1]-current[2],current[1]-current[3]);
        if InPlane(target[1],c_normal,current[1],eps) and InPlane(target[2],c_normal,current[1],eps) and InPlane(target[3],c_normal,current[1],eps) then
            return PointInTriangle(current[1],current[2],current[3],target[1],eps) and PointInTriangle(current[1],current[2],current[3],target[2],eps) and PointInTriangle(current[1],current[2],current[3],target[3],eps);
        else
            return false;
        fi;
    end
);


InstallMethod( PushLineToEdges,
    [IsList, IsList, IsList, IsList, IsList, IsFloat],
    function(a1, a2, a3, x, l, epsilon)
        local inters_at_edge, normal, o1, o2, o3, alpha_1, alpha_2, alpha_3, alpha, alpha2, abs_alpha, temp, val_1, val_2, po1, al1, al2, po2, p1, p2;
        

        inters_at_edge := false;
        
        normal := Crossproduct(a2-a1,a3-a1);
        
        o1 := Crossproduct(a1-a2,normal);
        o2 := Crossproduct(a2-a3,normal);
        o3 := Crossproduct(a3-a1,normal);
        if (o1*l) = 0. then
            alpha_1 := 10000000000000.;
        else
            alpha_1 := -( (o1*x) - (a1*o1) )/ (o1*l);
        fi;
        
        if (o2*l) = 0. then
            alpha_2 := 10000000000000.;
        else
            alpha_2 := - ( (o2*x) - (a2*o2) ) / (o2*l);
        fi;
        
        if (o3*l) = 0. then
            alpha_3 := 10000000000000.;
        else
            alpha_3 := - ( (o3*x) - (a3*o3) ) / (o3*l);
        fi;

        alpha2 := [alpha_1, alpha_2, alpha_3];
        alpha := [alpha_1, alpha_2, alpha_3];
        abs_alpha := [AbsoluteValue(alpha_1),AbsoluteValue(alpha_2),AbsoluteValue(alpha_3)];
        
        
        
        
        po1 := PositionMinimum(abs_alpha);
        al1 := alpha[po1];
        
        Remove(alpha, po1);
        Remove(abs_alpha, po1);
        
        
        po2 := PositionMinimum(abs_alpha);
        if not PointInTriangle(a1,a2,a3,x + alpha[po2]*l,epsilon) then
            temp := [1,2];
            Remove(temp,po2);
            po2 := temp[1];
        fi;
        al2 := alpha[po2];
        
        Remove(alpha,po2);
        Remove(abs_alpha, po2);
        
        if FlEq(al1,al2,epsilon) then
            # then the line pushes to a vertex
    
            al2 := alpha[1];
            
        fi;   
        
        
        p1 := x + al1*l;
        p2 := x + al2*l;
        if FlEq(al1,0.,epsilon) or FlEq(al2,0.,epsilon) then
            inters_at_edge := true;
        fi;
        
        return [p1,p2,inters_at_edge, po1, po2, al1, al2, alpha2];
    end
);


InstallMethod( SplitAdjacent,
    [IsTriangularComplex and IsNotEdgeRamified, IsList, IsList, IsList, IsList, IsList],
    function(surf, params, edge_v1, edge_v2, split_v, info)
        local I,j, J1, J2, face1, face2, face3, s1, s2, s3, surf1, surf2, surf3, e3, copy_verts, split_v1, split_v2, eps, l;
        
        # find faces to split
        l := Length(split_v);
        eps := info[1];

        if l = 1 then
            split_v := split_v[1];
        else
            # split edge in three
            if Norm2(split_v[1]-edge_v1) < Norm2(split_v[2]-edge_v1) then
                split_v1 := split_v[1];
                split_v2 := split_v[2];
            else
                split_v1 := split_v[2];
                split_v2 := split_v[1];
            fi;
        fi;        

        I := ShallowCopy(Faces(surf));
        for j in I do 
       	    if IsBound(params[1][j]) then
		    if Length(NumericalUniqueListOfLists([params[1][j][1],params[1][j][2],params[1][j][3],edge_v1],eps)) < 4 and 
		       Length(NumericalUniqueListOfLists([params[1][j][1],params[1][j][2],params[1][j][3],edge_v2],eps)) < 4 and 
		       Length(NumericalUniqueListOfLists([params[1][j][1],params[1][j][2],params[1][j][3]],eps)) = 3 then
		        # face adjacent to split up edge
		        copy_verts := [params[1][j][1],params[1][j][2],params[1][j][3]];
		        Remove(copy_verts,NumericalPosition(copy_verts,edge_v1,eps));
		        Remove(copy_verts,NumericalPosition(copy_verts,edge_v2,eps));

		        e3 := copy_verts[1];
		        
		        
		        
		        if l = 1 then
		            # split edge in two
		            
		            face1 :=[split_v,
		                     edge_v1,
		                     e3];
		            face2 :=[split_v,
		                     edge_v2,
		                     e3];
		            
		            surf1 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
		            surf2 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );

			    
		            J1 := ShallowCopy(Faces(surf));
		            
		            
		            s1 := DisjointUnion(surf, surf1);
		            surf := s1[1];
		            
		            s2 := DisjointUnion(surf, surf2);
		            surf := s2[1];
		            
		            J2 := ShallowCopy(Faces(surf));
		            
		            SubtractSet(J2,J1);

		            surf := RemoveFaces(surf,[j]);
		            
		            Unbind\[\](params[1],j);
			    Unbind\[\](params[3],j);
			    

		            params[1][J2[1]] := [face1[1], face1[2], face1[3], Crossproduct(face1[2]-face1[1],face1[3]-face1[1])/Norm2(Crossproduct(face1[2]-face1[1],face1[3]-face1[1])), VerticesOfFace(surf,J2[1])];
		            params[1][J2[2]] := [face2[1], face2[2], face2[3], Crossproduct(face2[2]-face2[1],face2[3]-face2[1])/Norm2(Crossproduct(face2[2]-face2[1],face2[3]-face2[1])), VerticesOfFace(surf,J2[2])];
		            
		            # remove faces that have been replaced
		            Unbind\[\](params[2],Position(params[2],j));
		            Add(params[2], J2[1],J2[1]);
		            Add(params[2], J2[2], J2[2]);
		            
		            params[3][J2[1]] := [J2[1]];
			    params[3][J2[2]] := [J2[2]];
		        else
		            # split edge in three

		            face1 :=[split_v1,
		                     e3,
		                     edge_v1];
		            face2 :=[split_v2,
		                     e3,
		                     split_v1];
		            face3 :=[edge_v2,
		                     e3,
		                     split_v2];

		            surf1 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
		            surf2 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
		            surf3 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );


		            J1 := ShallowCopy(Faces(surf));
		            
		            
		            s1 := DisjointUnion(surf, surf1);
		            surf := s1[1];
		            
		            s2 := DisjointUnion(surf, surf2);
		            surf := s2[1];
		            
		            s3 := DisjointUnion(surf, surf3);
		            surf := s3[1];
	    
		            J2 := ShallowCopy(Faces(surf));
		            
		            SubtractSet(J2,J1);
		            
		            surf := RemoveFaces(surf,[j]);

		            Unbind\[\](params[1],j);
			    Unbind\[\](params[3],j);
		            params[1][J2[1]] := [face1[1], face1[2], face1[3], Crossproduct(face1[2]-face1[1],face1[3]-face1[1])/Norm2(Crossproduct(face1[2]-face1[1],face1[3]-face1[1])), VerticesOfFace(surf,J2[1])];
		            params[1][J2[2]] := [face2[1], face2[2], face2[3], Crossproduct(face2[2]-face2[1],face2[3]-face2[1])/Norm2(Crossproduct(face2[2]-face2[1],face2[3]-face2[1])), VerticesOfFace(surf,J2[2])];
		            params[1][J2[3]] := [face3[1], face3[2], face3[3], Crossproduct(face3[2]-face3[1],face3[3]-face3[1])/Norm2(Crossproduct(face3[2]-face3[1],face3[3]-face3[1])), VerticesOfFace(surf,J2[3])];
		            
		            # remove faces that have been replaced
		            Unbind\[\](params[2],Position(params[2],j));
		            Add(params[2], J2[1],J2[1]);
		            Add(params[2], J2[2], J2[2]);
		            Add(params[2], J2[3], J2[3]);
		            
		            params[3][J2[1]] := [J2[1]];
			    params[3][J2[2]] := [J2[2]];
			    params[3][J2[3]] := [J2[3]];
		        fi;
		    fi;
            fi;
        od;
        
        return [surf,params];
    end
);


InstallMethod( Retriangulate,
    [IsTriangularComplex and IsNotEdgeRamified, IsRecord, IsList, IsList, IsList, IsList],
    function(surf, reco, params, inters, res, li)
        local triang, surf1, surf2, surf3, surf4, surf5, verts, bad_v, v1, v2, v3, v4, v5, v6, verticesPos1, verticesPos2, verticesPos3, verticesPos4, verticesPos5, vertex, J1, J2, s1, s2, s3, s4, s5, dist, area, split_adj1, split_adj2, split_adj3, temp, j, k, eps, p1, p2;
        # need to split up the triangle in 5 parts here

        j := li[1];
        k := li[3];
        eps := li[2];
        
        vertex := 0;
        
        triang := ShallowCopy(params[1][j]);
        verts := [triang[1],triang[2],triang[3]];
        

        area := [SurfaceTriangle(inters[1],inters[2],triang[1]),SurfaceTriangle(inters[1],inters[2],triang[2]),SurfaceTriangle(inters[1],inters[2],triang[3])];
        bad_v := Position(area, Minimum(area));
        v1 := 0.;
        v2 := 0.;
        v3 := 0.;
        v4 := 0.;
        v5 := 0.;
        
        # connect to only other ones 
                
        
        
        
        v3 := verts[bad_v]; # not connected to both of the intersection vertices
        Remove(verts,bad_v);
        v2 := verts[2];
        v1 := verts[1];
        
        surf1 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
        surf2 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
        surf3 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
        surf4 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
        surf5 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
        
        # need to check geometrical arrangement here
        dist :=  [Norm2(inters[1]-v3),Norm2(inters[2]-v3)];
        if dist[1] < dist[2] then
        	v4 := inters[1];
        	v5 := inters[2];
        else
        	v4 := inters[2];
        	v5 := inters[1];
        fi;
        
        
        verticesPos1 :=[v1,
                        v4,
                        v5];
        verticesPos2 :=[v1,
                        v3,
                        v4];
                        
        verticesPos3 :=[v1,
                        v5,
                        v2];
        
        verticesPos4 :=[v4,
                        v5,
                        v2];
        
        verticesPos5 :=[v4,
                        v3,
                        v2];
                        
        
                 
        J1 := ShallowCopy(Faces(surf));
        
        
        s1 := DisjointUnion(surf, surf1);
        surf := s1[1];
        
        s2 := DisjointUnion(surf, surf2);
        surf := s2[1];
        
        s3 := DisjointUnion(surf, surf3);
        surf := s3[1];
        
        s4 := DisjointUnion(surf, surf4);
        surf := s4[1];
         
        s5 := DisjointUnion(surf, surf5);
        surf := s5[1];
        
        J2 := ShallowCopy(Faces(surf));
        
        SubtractSet(J2,J1);
        
        surf := RemoveFaces(surf,[j]);
        
        Unbind\[\](params[1],j);
	Unbind\[\](params[3],j);
        params[1][J2[1]] := [v1, v4, v5, Crossproduct(v1-v4,v1-v5)/Norm2(Crossproduct(v1-v4,v1-v5)), VerticesOfFace(surf,J2[1])];
        params[1][J2[2]] := [v1, v3, v4, Crossproduct(v1-v3,v1-v4)/Norm2(Crossproduct(v1-v3,v1-v4)), VerticesOfFace(surf,J2[2])];
        params[1][J2[3]] := [v1, v5, v2, Crossproduct(v1-v5,v1-v2)/Norm2(Crossproduct(v1-v5,v1-v2)), VerticesOfFace(surf,J2[3])];
        params[1][J2[4]] := [v4, v5, v2, Crossproduct(v4-v5,v4-v2)/Norm2(Crossproduct(v4-v5,v4-v2)), VerticesOfFace(surf,J2[4])];
        params[1][J2[5]] := [v4, v3, v2, Crossproduct(v4-v3,v4-v2)/Norm2(Crossproduct(v4-v3,v4-v2)), VerticesOfFace(surf,J2[5])];
        
         # remove faces that have been replaced
        Unbind\[\](params[2],Position(params[2],j));
        Add(params[2], J2[1], J2[1]);
        Add(params[2], J2[2], J2[2]);
        Add(params[2], J2[3], J2[3]);
        Add(params[2], J2[4], J2[4]);
        Add(params[2], J2[5], J2[5]);
        
        params[3][J2[1]] := [J2[1]];
	params[3][J2[2]] := [J2[2]];
	params[3][J2[3]] := [J2[3]];
	params[3][J2[4]] := [J2[4]];
	params[3][J2[5]] := [J2[5]];
        
        IsNotEdgeRamified(surf);
        IsTriangularComplex(surf);

        return [surf,params];
        
        
        
    end
);

InstallMethod( RetriangulateEasy,
    [IsTriangularComplex and IsNotEdgeRamified, IsRecord, IsList, IsList, IsList, IsList],
    function(surf, reco, params, inters, res, li)
        local current, surf1, surf2, surf3, surf4, surf5, surf6, surf7, surf8, verts, triang, vPos1, vPos2, vPos3, in1_on_edge, in2_on_edge, inters_verts, vertex, inters_vertex, same, other, bad_v, i_verts, int_which_v, which_one, along_edges, which_edges, v1, v2, v3, v4, v5, v6, v7, v8, orthog, alpha1, alpha2, verticesPos1, verticesPos2, verticesPos3, verticesPos4, verticesPos5, face1, face2, face3, face4, face5, face6, face7, face8, res_temp, not_split, which_v1, near, J1, J2, s1, s2, s3, s4, s5, s6, s7, s8, split_adj1, split_adj2, split_adj3, o1, o2, m1, m2, j, l, k, one_or_zero, test, eps, case, len, area, c_normal, dist, l1, leng1;
        j := li[1];
        k := li[3];
        eps := li[2];
	
	current := params[1][j];
        c_normal := current[4];
        
        # need to split up the triangle in 2-6 parts here
        # as the intersection pushes to an vertex
        # depending on if the intersection points lie on vertices or edges

        alpha1 := AbsoluteValue(res[6]);
        alpha2 := AbsoluteValue(res[7]);
        inters_vertex := false;
        
        

	# find out on which of the edges (if at all) the intersection points lie
	# depending on this we need to divide the triangle differently
	
	orthog := [Crossproduct(current[2]-current[1],c_normal),Crossproduct(current[3]-current[2],c_normal),Crossproduct(current[1]-current[3],c_normal)];
                        
        orthog[1] := orthog[1] / Sqrt(orthog[1]*orthog[1]);
        orthog[2] := orthog[2] / Sqrt(orthog[2]*orthog[2]);
        orthog[3] := orthog[3] / Sqrt(orthog[3]*orthog[3]);
	orthog[1] := orthog[1] * ((orthog[1]*current[3]-orthog[1]*current[1]) / AbsoluteValue(orthog[1]*current[3]-orthog[1]*current[1]));
        orthog[2] := orthog[2] * ((orthog[2]*current[1]-orthog[2]*current[2]) / AbsoluteValue(orthog[2]*current[1]-orthog[2]*current[2]));
        orthog[3] := orthog[3] * ((orthog[3]*current[2]-orthog[3]*current[3]) / AbsoluteValue(orthog[3]*current[2]-orthog[3]*current[3]));
                

        o1 := [AbsoluteValue(orthog[1]*inters[1]-orthog[1]*current[1]), AbsoluteValue(orthog[2]*inters[1]-orthog[2]*current[2]), AbsoluteValue(orthog[3]*inters[1]-orthog[3]*current[3])];
        o2 := [AbsoluteValue(orthog[1]*inters[2]-orthog[1]*current[1]), AbsoluteValue(orthog[2]*inters[2]-orthog[2]*current[2]), AbsoluteValue(orthog[3]*inters[2]-orthog[3]*current[3])];
        
        in1_on_edge := FlEq(o1[1],0.,eps) or FlEq(o1[2],0.,eps) or FlEq(o1[3],0.,eps);
        in2_on_edge := FlEq(o2[1],0.,eps) or FlEq(o2[2],0.,eps) or FlEq(o2[3],0.,eps);
        
        case := 3;
        
        if in1_on_edge or in2_on_edge then
            # then one intersection point is at a edge or vertex
            case := 2;
        fi;
        if in1_on_edge and in2_on_edge then
            # then both of the intersection points are at an edge or a vertex
            # minimum of alpha has to be 0 then
            case := 1;
        fi;
        
        
       	
        if case = 1 then
            
            if Length(NumericalUniqueListOfLists([inters[1],current[1],current[2],current[3]],eps)) < 4 then
            	v4 := inters[2];
            	other := inters[1];
            	inters_vertex := true;
            elif Length(NumericalUniqueListOfLists([inters[2],current[1],current[2],current[3]],eps)) < 4 then
            	v4 := inters[1];
            	other := inters[2];
            	inters_vertex := true;
            else
            	v4 := inters[1];
            	v5 := inters[2];
            fi;
            
            
            
            # test if the intersection is just along an edge
            along_edges := [FlEq(o1[1],0.,eps) and FlEq(o2[1],0.,eps), FlEq(o1[2],0.,eps) and FlEq(o2[2],0.,eps), FlEq(o1[3],0.,eps) and FlEq(o2[3],0.,eps)];
            
            
            
            
            if inters_vertex then
            	    if along_edges[1] or along_edges[2] or along_edges[3] then
            	    	# the intersection points both lie on the same edge
            	    	
            		# choose points of edge on which intersection lies
            		if along_edges[1] then
            			l1 := 1;
            		elif along_edges[2] then
            			l1 := 2;
            		else
            			l1 := 3;
            		fi;
            		leng1 := Length(NumericalUniqueListOfLists([inters[1],inters[2],current[l1 mod 3 + 1],current[l1]],eps));
		        
		        
		        
		        if leng1 = 4 then
		            # split edge in three parts
		            IsNotEdgeRamified(surf);
		            IsTriangularComplex(surf);

		            split_adj1 := SplitAdjacent(surf,params,current[l1 mod 3 + 1],current[l1],[inters[1],inters[2]],[eps,j,li[3]]);
		            surf := split_adj1[1];
		            params := split_adj1[2];

		        elif leng1 = 3 then
		            # split edge in two parts
		            IsNotEdgeRamified(surf);
		            IsTriangularComplex(surf);
		            
		            if Length(NumericalUniqueListOfLists([inters[1],current[l1 mod 3 + 1],current[l1]],eps)) < 3 then
		            
		            	
		                split_adj1 := SplitAdjacent(surf,params,current[l1 mod 3 + 1],current[l1],[inters[2]],[eps,j,li[3]]);
		                
		                surf := split_adj1[1];
		                params := split_adj1[2];
		            else
		                split_adj1 := SplitAdjacent(surf,params,current[l1 mod 3 + 1],current[l1],[inters[1]],[eps,j,li[3]]);
		                surf := split_adj1[1];
		                params := split_adj1[2];
		            fi;
		            
		        fi;
            	    else	
		    	    # split triangle into two parts easily
		    	    
		    	    surf1 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
			    surf2 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );

			    same := [Norm2(other-current[1]),Norm2(other-current[2]),Norm2(other-current[3])];
			    which_one := Position(same, Minimum(same));
			    
			    v1 := current[which_one];
			    
			    l := [1,2,3];
			    Remove(l,Position(l,which_one));
			    v2 := current[l[1]];
			    v3 := current[l[2]];
			    
			    
			    verticesPos1 :=[v1,
				            v3,
				            v4 ];
			    verticesPos2 :=[v1,
				            v2,
				            v4 ];
			    
			    

			    J1 := ShallowCopy(Faces(surf));
			    
			    
			    s1 := DisjointUnion(surf, surf1);
			    surf := s1[1];
			    
			    s2 := DisjointUnion(surf, surf2);
			    surf := s2[1];
			    
			    J2 := ShallowCopy(Faces(surf));
			    
			    SubtractSet(J2,J1);
			    
			    surf := RemoveFaces(surf,[j]);
			    
			    Unbind\[\](params[1],j);
			    Unbind\[\](params[3],j);
			    params[1][J2[1]] := [v1, v3, v4, Crossproduct(v1-v3,v1-v4)/Norm2(Crossproduct(v1-v3,v1-v4)), VerticesOfFace(surf,J2[1])];
			    params[1][J2[2]] := [v1, v2, v4, Crossproduct(v1-v2,v1-v4)/Norm2(Crossproduct(v1-v2,v1-v4)), VerticesOfFace(surf,J2[2])];
			    
			    # remove faces that have been replaced
			    Unbind\[\](params[2],Position(params[2],j));
			    Add(params[2], J2[1],J2[1]);
			    Add(params[2], J2[2], J2[2]);
			    
			    params[3][J2[1]] := [J2[1]];
			    params[3][J2[2]] := [J2[2]];

			    IsNotEdgeRamified(surf);
			    IsTriangularComplex(surf);

			    split_adj1 := SplitAdjacent(surf,params,v2,v3,[v4],[eps]);
			    surf := split_adj1[1];
			    params := split_adj1[2];
		    fi;
		    
            else
            	# split triangle in three parts
            	surf1 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
		surf2 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
		surf3 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
		
		v4 := inters[1];
		v5 := inters[2];
		
		
                
                which_edges := [FlEq(o1[1],0.,eps) or FlEq(o2[1],0.,eps), FlEq(o1[2],0.,eps) or FlEq(o2[2],0.,eps), FlEq(o1[3],0.,eps) or FlEq(o2[3],0.,eps)];
  
                  
                if not which_edges[1] then
                	v1 := current[3];
                	verts := [current[1],current[2]];
                	
                	if FlEq(o1[2],0.,eps) then
                		int_which_v := [2,1];
                	else
                		int_which_v := [1,2];
                	fi;
                	
                elif not which_edges[2] then
                	v1 := current[1];
                	verts := [current[2],current[3]];
                	if FlEq(o1[3],0.,eps) then
                		int_which_v := [2,1];
                	else
                		int_which_v := [1,2];
                	fi;
                else
                	v1 := current[2];
                	verts := [current[1],current[3]];
                	if FlEq(o1[2],0.,eps) then
                		int_which_v := [2,1];
                	else
                		int_which_v := [1,2];
                	fi;
                fi;
                
                v2 := verts[int_which_v[1]];
		v3 := verts[int_which_v[2]];
                
		area := [Minimum(SurfaceTriangle(v2,v4,v5),SurfaceTriangle(v2,v3,v5)),Minimum(SurfaceTriangle(v3,v4,v5),SurfaceTriangle(v2,v3,v4))];
		
		
		
		if area[1] > area [2] then
			vPos1 :=[v1,
		                 v4,
		                 v5];
			vPos2 :=[v2,
				 v4,
				 v5];
			vPos3 :=[v2,
				 v3,
				 v5];
		else
			vPos1 :=[v1,
		                 v4,
		                 v5];
			vPos2 :=[v3,
				 v4,
				 v5];
			vPos3 :=[v2,
				 v3,
				 v4];
		
		fi;
		
		           
		
		J1 := ShallowCopy(Faces(surf));
		    
		    
		s1 := DisjointUnion(surf, surf1);
		surf := s1[1];
		    
		s2 := DisjointUnion(surf, surf2);
		surf := s2[1];
		
		s3 := DisjointUnion(surf, surf3);
		surf := s3[1];
		    
		J2 := ShallowCopy(Faces(surf));
		    
		SubtractSet(J2,J1);
		    
		surf := RemoveFaces(surf,[j]);
		    
		Unbind\[\](params[1],j);
		Unbind\[\](params[3],j);
		params[1][J2[1]] := [vPos1[1], vPos1[2], vPos1[3], Crossproduct(vPos1[1]-vPos1[2],vPos1[1]-vPos1[3])/Norm2(Crossproduct(vPos1[1]-vPos1[2],vPos1[1]-vPos1[3])), VerticesOfFace(surf,J2[1])];
		params[1][J2[2]] := [vPos2[1], vPos2[2], vPos2[3], Crossproduct(vPos2[1]-vPos2[2],vPos2[1]-vPos2[3])/Norm2(Crossproduct(vPos2[1]-vPos2[2],vPos2[1]-vPos2[3])), VerticesOfFace(surf,J2[2])];
		params[1][J2[3]] := [vPos3[1], vPos3[2], vPos3[3], Crossproduct(vPos3[1]-vPos3[2],vPos3[1]-vPos3[3])/Norm2(Crossproduct(vPos3[1]-vPos3[2],vPos3[1]-vPos3[3])), VerticesOfFace(surf,J2[3])];
		    
		# remove faces that have been replaced
		Unbind\[\](params[2],Position(params[2],j));
	        Add(params[2], J2[1], J2[1]);
		Add(params[2], J2[2], J2[2]);
		Add(params[2], J2[3], J2[3]);
		
		params[3][J2[1]] := [J2[1]];
		params[3][J2[2]] := [J2[2]];
		params[3][J2[3]] := [J2[3]];
		
		IsNotEdgeRamified(surf);
		IsTriangularComplex(surf);
		
		
		split_adj1 := SplitAdjacent(surf,params,v1,verts[int_which_v[1]],[v4],[eps]);
		surf := split_adj1[1];
		params := split_adj1[2];
		
		IsNotEdgeRamified(surf);
		IsTriangularComplex(surf);
		
		split_adj2 := SplitAdjacent(surf,params,v1,verts[int_which_v[2]],[v5],[eps]);
		surf := split_adj2[1];
		params := split_adj2[2];
		
		IsNotEdgeRamified(surf);
		IsTriangularComplex(surf);
            fi;
            

        elif case = 2 then
        
            
            # test first which of the intersection points lies in a vertex or edge
            verts := [current[1], current[2], current[3]];

            inters_verts := [Length(NumericalUniqueListOfLists([verts[1], verts[2], verts[3],inters[1]],eps)) < 4, 
                             Length(NumericalUniqueListOfLists([verts[1], verts[2], verts[3],inters[2]],eps)) < 4];
            
            
 	    
            if inters_verts[1] or inters_verts[2] then
                # a intersection point is a vertex, thus we can split the triangle in 3 parts
                if inters_verts[1] then
                    v1 := inters[1];
                    v2 := inters[2];
                    Remove(verts, NumericalPosition(verts, inters[1],eps));
                else
                    v1 := inters[2];
                    v2 := inters[1];
                    Remove(verts, NumericalPosition(verts, inters[2],eps));
                fi;
                
                v3 := verts[1];
                v4 := verts[2];

                surf1 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
                surf2 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
                surf3 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );

                face1 :=[v1,
                         v2,
                         v3 ];
                face2 :=[v1,
                         v2,
                         v4 ];

                face3 :=[v2,
                         v3,
                         v4 ];


                J1 := ShallowCopy(Faces(surf));
            
            
                s1 := DisjointUnion(surf, surf1);
                surf := s1[1];
                
                s2 := DisjointUnion(surf, surf2);
                surf := s2[1];

                s3 := DisjointUnion(surf, surf3);
                surf := s3[1];
                
                J2 := ShallowCopy(Faces(surf));
                
                SubtractSet(J2,J1);
                
                surf := RemoveFaces(surf,[j]);
                
                Unbind\[\](params[1],j);
		Unbind\[\](params[3],j);
		
                params[1][J2[1]] := [face1[1], face1[2], face1[3], Crossproduct(face1[2]-face1[1],face1[3]-face1[1])/Norm2(Crossproduct(face1[2]-face1[1],face1[3]-face1[1])), VerticesOfFace(surf,J2[1])];
                params[1][J2[2]] := [face2[1], face2[2], face2[3], Crossproduct(face2[2]-face2[1],face2[3]-face2[1])/Norm2(Crossproduct(face2[2]-face2[1],face2[3]-face2[1])), VerticesOfFace(surf,J2[2])];
                params[1][J2[3]] := [face3[1], face3[2], face3[3], Crossproduct(face3[2]-face3[1],face3[3]-face3[1])/Norm2(Crossproduct(face3[2]-face3[1],face3[3]-face3[1])), VerticesOfFace(surf,J2[3])];
                
                # remove faces that have been replaced
                Unbind\[\](params[2],Position(params[2],j));
                Add(params[2], J2[1],J2[1]);
                Add(params[2], J2[2], J2[2]);
                Add(params[2], J2[3], J2[3]);
                
                params[3][J2[1]] := [J2[1]];
		params[3][J2[2]] := [J2[2]];
		params[3][J2[3]] := [J2[3]];

            else
                # no intersection point is a vertex, thus we need to split the triangle in 4 parts
                # Find out which vertex lies on an edge first
                
                orthog := [Crossproduct(current[2]-current[1],c_normal),Crossproduct(current[3]-current[2],c_normal),Crossproduct(current[1]-current[3],c_normal)];

                m1 := Minimum(AbsoluteValue(o1[1]),AbsoluteValue(o1[2]),AbsoluteValue(o1[3]));
                m2 := Minimum(AbsoluteValue(o2[1]),AbsoluteValue(o2[2]),AbsoluteValue(o2[3]));

                if FlEq(m1,0.,eps) then
                    # then inters[2] is inside the vertex
                    v2 := inters[1];
                    v4 := inters[2];
                    if FlEq(o1[1],0.,eps) then
                    	v1 := current[1];
                    	v3 := current[2];
                    	v5 := current[3];
                    elif FlEq(o1[2],0.,eps) then
                    	v1 := current[2];
                    	v3 := current[3];
                    	v5 := current[1];
                    else
                    	v1 := current[1];
                    	v3 := current[3];
                    	v5 := current[2];
                    fi;
                    
                elif FlEq(m2,0.,eps) then
                    v2 := inters[2];
                    v4 := inters[1];
                    
                    if FlEq(o2[1],0.,eps) then
                    	v1 := current[1];
                    	v3 := current[2];
                    	v5 := current[3];
                    elif FlEq(o2[2],0.,eps) then
                    	v1 := current[2];
                    	v3 := current[3];
                    	v5 := current[1];
                    else
                    	v1 := current[1];
                    	v3 := current[3];
                    	v5 := current[2];
                    fi;
                else
                    # can not happen geometrically
                    ThrowError();
                fi;


                surf1 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
                surf2 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
                surf3 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
                surf4 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );

                face1 :=[v1,
                         v2,
                         v4 ];
                face2 :=[v1,
                         v4,
                         v5 ];

                face3 :=[v2,
                         v3,
                         v4 ];
                face4 :=[v3,
                         v4,
                         v5 ];
                         
		

                J1 := ShallowCopy(Faces(surf));
            
            
                s1 := DisjointUnion(surf, surf1);
                surf := s1[1];
                
                s2 := DisjointUnion(surf, surf2);
                surf := s2[1];

                s3 := DisjointUnion(surf, surf3);
                surf := s3[1];

                s4 := DisjointUnion(surf, surf4);
                surf := s4[1];
                
                J2 := ShallowCopy(Faces(surf));
                
                SubtractSet(J2,J1);
                
                surf := RemoveFaces(surf,[j]);

                Unbind\[\](params[1],j);
		Unbind\[\](params[3],j);
		
                params[1][J2[1]] := [face1[1], face1[2], face1[3], Crossproduct(face1[2]-face1[1],face1[3]-face1[1])/Norm2(Crossproduct(face1[2]-face1[1],face1[3]-face1[1])), VerticesOfFace(surf,J2[1])];
                params[1][J2[2]] := [face2[1], face2[2], face2[3], Crossproduct(face2[2]-face2[1],face2[3]-face2[1])/Norm2(Crossproduct(face2[2]-face2[1],face2[3]-face2[1])), VerticesOfFace(surf,J2[2])];
                params[1][J2[3]] := [face3[1], face3[2], face3[3], Crossproduct(face3[2]-face3[1],face3[3]-face3[1])/Norm2(Crossproduct(face3[2]-face3[1],face3[3]-face3[1])), VerticesOfFace(surf,J2[3])];
                params[1][J2[4]] := [face4[1], face4[2], face4[3], Crossproduct(face4[2]-face4[1],face4[3]-face4[1])/Norm2(Crossproduct(face4[2]-face4[1],face4[3]-face4[1])), VerticesOfFace(surf,J2[4])];
                
                # remove faces that have been replaced
                Unbind\[\](params[2],Position(params[2],j));
                Add(params[2], J2[1],J2[1]);
                Add(params[2], J2[2], J2[2]);
                Add(params[2], J2[3], J2[3]);
                Add(params[2], J2[4], J2[4]);
                
                params[3][J2[1]] := [J2[1]];
		params[3][J2[2]] := [J2[2]];
		params[3][J2[3]] := [J2[3]];
		params[3][J2[4]] := [J2[4]];

                IsNotEdgeRamified(surf);
                IsTriangularComplex(surf);

                split_adj1 := SplitAdjacent(surf,params,v1,v3,[v2],[eps]);
                surf := split_adj1[1];
                params := split_adj1[2];
            fi;
            
        else
        	j:= li[1];
		eps := li[2];
		
		vertex := 0;
        
		triang := ShallowCopy(params[1][j]);
		verts := [triang[1],triang[2],triang[3]];
		

		area := [SurfaceTriangle(inters[1],inters[2],triang[1]),SurfaceTriangle(inters[1],inters[2],triang[2]),SurfaceTriangle(inters[1],inters[2],triang[3])];
		bad_v := Position(area, Minimum(area));
		v1 := 0.;
		v2 := 0.;
		v3 := 0.;
		v4 := 0.;
		v5 := 0.;
		
		# connect to only other ones 
		        
		
		
		
		v3 := verts[bad_v]; # not connected to both of the intersection vertices
		Remove(verts,bad_v);
		v2 := verts[2];
		v1 := verts[1];
		
		surf1 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
		surf2 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
		surf3 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
		surf4 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
		surf5 := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
		
		# need to check geometrical arrangement here
		dist :=  [Norm2(inters[1]-v3),Norm2(inters[2]-v3)];
		if dist[1] < dist[2] then
			v4 := inters[1];
			v5 := inters[2];
		else
			v4 := inters[2];
			v5 := inters[1];
		fi;
		
		
		verticesPos1 :=[v1,
		                v4,
		                v5];
		verticesPos2 :=[v1,
		                v3,
		                v4];
		                
		verticesPos3 :=[v1,
		                v5,
		                v2];
		
		verticesPos4 :=[v4,
		                v5,
		                v2];
		
		verticesPos5 :=[v4,
		                v3,
		                v2];
		         
		J1 := ShallowCopy(Faces(surf));
		
		
		s1 := DisjointUnion(surf, surf1);
		surf := s1[1];
		
		s2 := DisjointUnion(surf, surf2);
		surf := s2[1];
		
		s3 := DisjointUnion(surf, surf3);
		surf := s3[1];
		
		s4 := DisjointUnion(surf, surf4);
		surf := s4[1];
		 
		s5 := DisjointUnion(surf, surf5);
		surf := s5[1];
		
		J2 := ShallowCopy(Faces(surf));
		
		SubtractSet(J2,J1);
		
		surf := RemoveFaces(surf,[j]);
		
		Unbind\[\](params[1],j);
		Unbind\[\](params[3],j);
		
		params[1][J2[1]] := [v1, v4, v5, Crossproduct(v1-v4,v1-v5)/Norm2(Crossproduct(v1-v4,v1-v5)), VerticesOfFace(surf,J2[1])];
		params[1][J2[2]] := [v1, v3, v4, Crossproduct(v1-v3,v1-v4)/Norm2(Crossproduct(v1-v3,v1-v4)), VerticesOfFace(surf,J2[2])];
		params[1][J2[3]] := [v1, v5, v2, Crossproduct(v1-v5,v1-v2)/Norm2(Crossproduct(v1-v5,v1-v2)), VerticesOfFace(surf,J2[3])];
		params[1][J2[4]] := [v4, v5, v2, Crossproduct(v4-v5,v4-v2)/Norm2(Crossproduct(v4-v5,v4-v2)), VerticesOfFace(surf,J2[4])];
		params[1][J2[5]] := [v4, v3, v2, Crossproduct(v4-v3,v4-v2)/Norm2(Crossproduct(v4-v3,v4-v2)), VerticesOfFace(surf,J2[5])];
		
		 # remove faces that have been replaced
		Unbind\[\](params[2],Position(params[2],j));
		Add(params[2], J2[1], J2[1]);
		Add(params[2], J2[2], J2[2]);
		Add(params[2], J2[3], J2[3]);
		Add(params[2], J2[4], J2[4]);
		Add(params[2], J2[5], J2[5]);
		
		params[3][J2[1]] := [J2[1]];
		params[3][J2[2]] := [J2[2]];
		params[3][J2[3]] := [J2[3]];
		params[3][J2[4]] := [J2[4]];
		params[3][J2[5]] := [J2[5]];
		
		IsNotEdgeRamified(surf);
		IsTriangularComplex(surf);

        fi;

        return [surf,params];
        
    end
);

InstallMethod( FixSelfIntersections,
    [IsTriangularComplex and IsNotEdgeRamified, IsRecord, IsList, IsList, IsList, IsFloat],
    function(surface, printRecord, params, inters, info, epsilon)
        local current, target, c_normal, t_normal, j, k, d,d2, point_in_plane,  point_in_plane2, split_adj1, split_adj2, res1, res1_temp, res2, res2_temp, s1, s2, leng2, leng1, l1, c_inside, t_inside, l2, c_on_edges, t_on_edges, c_on_edge_both, t_on_edge_both, c_on_edge1, c_on_edge2, c_on_edge3, t_on_edge1, t_on_edge2, t_on_edge3, checks_c, checks_t,inters_c_vertices ,inters_t_vertices;
        
        
        j := info[1];
        k := info[2];
        
        
        current := params[1][j];
        target := params[1][k];
        
        c_normal := current[4];
        t_normal := target[4];

        t_inside := TriangleFullyInside(target,current,epsilon);
        c_inside := TriangleFullyInside(current,target,epsilon);
	if (j = 5119 and k = 5116) or (k = 5119 and j = 5116) then
		#ThrowError();
	fi;
        # if one triangle is fully inside the other one we just delete it
        if t_inside then
            # target inside current
            surface := RemoveFaces(surface,[k]);
            Unbind\[\](params[1],k);
	    Unbind\[\](params[3],k);
            Unbind\[\](params[2],Position(params[2],k));
            
            return [surface,params];
        elif c_inside then
            surface := RemoveFaces(surface,[j]);
            Unbind\[\](params[1],j);
	    Unbind\[\](params[3],j);
            Unbind\[\](params[2],Position(params[2],j));
            return [surface,params];
        fi;
        
        # check if the intersection is one edge of a triangle
        inters_c_vertices := Length(NumericalUniqueListOfLists([current[1],current[2],current[3],inters[1],inters[2]],epsilon));
        inters_t_vertices := Length(NumericalUniqueListOfLists([target[1],target[2],target[3],inters[1],inters[2]],epsilon));

        c_on_edge1 := OnEdges(current,inters[1],epsilon);
        c_on_edge2 := OnEdges(current,inters[2],epsilon);
        c_on_edge_both := [c_on_edge1[1] and c_on_edge2[1],c_on_edge1[2] and c_on_edge2[2],c_on_edge1[3] and c_on_edge2[3]];
        c_on_edges := c_on_edge_both[1] or c_on_edge_both[2] or c_on_edge_both[3];
        

        t_on_edge1 := OnEdges(target,inters[1],epsilon);
        t_on_edge2 := OnEdges(target,inters[2],epsilon);
        t_on_edge_both := [t_on_edge1[1] and t_on_edge2[1],t_on_edge1[2] and t_on_edge2[2],t_on_edge1[3] and t_on_edge2[3]];
        t_on_edges := t_on_edge_both[1] or t_on_edge_both[2] or t_on_edge_both[3];

        
        checks_c := inters_c_vertices  <=  3 or c_on_edges;
        checks_t := inters_t_vertices  <=  3 or t_on_edges;

        d := current[1]*c_normal;
        d2 := target[1]*t_normal;
        
        
        
        
        
        res1 := PushLineToEdges(current[1], current[2], current[3], inters[1], inters[1]-inters[2], epsilon);
        res1_temp := PushLineToEdges(current[1], current[2], current[3], inters[2], inters[2]-inters[1], epsilon);
        res1[3] := res1[3] or res1_temp[3];
        
	
	

        if (not info[4][1]) and (not checks_c) then
            if res1[3] or inters_c_vertices <= 4 then
                s1 := RetriangulateEasy(surface, printRecord, params,inters, res1, [info[1], epsilon, info[2]]);
            else
                s1 := Retriangulate(surface, printRecord, params,inters, res1, [info[1], epsilon, info[2]]);
            fi;
            params := s1[2];
            surface := ShallowCopy(s1[1]);
        
        # ignore intersection if it is an edge    
        elif not inters_c_vertices <=  3 then
            # the edge of the intersection is perpendicular to our face
            # so either all three or only two of the target vertices lie in the current triangle plane
            point_in_plane := [FlEq(target[1]*c_normal,d,epsilon),FlEq(target[2]*c_normal,d,epsilon),FlEq(target[3]*c_normal,d,epsilon)];
            
            if not (point_in_plane[1] and point_in_plane[2] and point_in_plane[3]) then
                # the intersection is along an edge
                # test if intersection spans a full edge
                l1 := info[4][2];
                leng1 := Length(NumericalUniqueListOfLists([inters[1],inters[2],current[l1 mod 3 + 1],current[l1]],epsilon));
                
                if leng1 = 4 then
                    # split edge in three parts
                    IsNotEdgeRamified(surface);
                    IsTriangularComplex(surface);

                    split_adj1 := SplitAdjacent(surface,params,current[l1 mod 3 + 1],current[l1],[inters[1],inters[2]],[epsilon,j,k]);
                    surface := split_adj1[1];
                    params := split_adj1[2];

                elif leng1 = 3 then
                    # split edge in two parts
                    IsNotEdgeRamified(surface);
                    IsTriangularComplex(surface);
                       
                    if Length(NumericalUniqueListOfLists([inters[1],current[l1 mod 3 + 1],current[l1]],epsilon)) < 3 then
                        split_adj1 := SplitAdjacent(surface,params,current[l1 mod 3 + 1],current[l1],[inters[2]],[epsilon,j,k]);
                        
                        surface := split_adj1[1];
                        params := split_adj1[2];
                    else
                        split_adj1 := SplitAdjacent(surface,params,current[l1 mod 3 + 1],current[l1],[inters[1]],[epsilon,j,k]);
                        surface := split_adj1[1];
                        params := split_adj1[2];
                    fi;

                fi;
            else
                if res1[3] or inters_c_vertices <= 4 then
                    s1 := RetriangulateEasy(surface, printRecord, params,inters, res1, [info[1], epsilon, info[2]]);
                else
                    s1 := Retriangulate(surface, printRecord, params,inters, res1, [info[1], epsilon, info[2]]);
                fi;
                params := s1[2];
                surface := ShallowCopy(s1[1]);
            fi;
                    
        fi;
        
        res2 := PushLineToEdges(target[1], target[2], target[3], inters[2], inters[2]-inters[1], epsilon);
        res2_temp := PushLineToEdges(target[1], target[2], target[3], inters[1], inters[1]-inters[2], epsilon);
        res2[3] := res2[3] or res2_temp[3];
        
        IsNotEdgeRamified(surface);
        IsSimplicialSurface(surface);
        IsTriangularComplex(surface);

        
       
            		
        if not info[3][1] and (not checks_t) then
            if res2[3] or inters_t_vertices <= 4 then
                s2 :=  RetriangulateEasy(surface, printRecord, params,inters, res2, [info[2], epsilon, info[1]]);
            else
                s2 :=  Retriangulate(surface, printRecord, params,inters, res2, [info[2], epsilon, info[1]]);
            fi;
            surface := s2[1]; 
            params := s2[2];
        elif not inters_t_vertices  <=  3 then
            # the edge of the intersection is perpendicular to our face
            # so either all three or only two of the target vertices lie in the current triangle plane

            point_in_plane2 := [FlEq(current[1]*t_normal,d2,epsilon),FlEq(current[2]*t_normal,d2,epsilon),FlEq(current[3]*t_normal,d2,epsilon)];
	    
            if not (point_in_plane2[1] and point_in_plane2[2] and point_in_plane2[3]) then
                # so intersection is along an edge
                # test if intersection spans a full edge
                l2 := info[3][2];
                leng2 := Length(NumericalUniqueListOfLists([inters[1],inters[2],target[l2 mod 3 + 1],target[l2]],epsilon));
                if leng2 = 4 then
                    # split edge in three parts
                    IsNotEdgeRamified(surface);
                    IsTriangularComplex(surface);

                    split_adj2 := SplitAdjacent(surface,params,target[l2 mod 3 + 1],target[l2],[inters[1],inters[2]],[epsilon,j,k]);
                    surface := split_adj2[1];
                    params := split_adj2[2];

                elif leng2 = 3 then
                    # split edge in two parts
                    IsNotEdgeRamified(surface);
                    IsTriangularComplex(surface);
                    
                    if Length(NumericalUniqueListOfLists([inters[1],target[l2 mod 3 + 1],target[l2]],epsilon)) < 3 then
                        split_adj2 := SplitAdjacent(surface,params,target[l2 mod 3 + 1],target[l2],[inters[2]],[epsilon,j,k]);
                        surface := split_adj2[1];
                        params := split_adj2[2];
                    else
                        split_adj2 := SplitAdjacent(surface,params,target[l2 mod 3 + 1],target[l2],[inters[1]],[epsilon,j,k]);
                        surface := split_adj2[1];
                        params := split_adj2[2];
                    fi;

                fi;
            else
                if res2[3] or inters_t_vertices <= 4 then
                	s2 :=  RetriangulateEasy(surface, printRecord, params,inters, res2, [info[2], epsilon, info[1]]);
                else
                    s2 :=  Retriangulate(surface, printRecord, params,inters, res2, [info[2], epsilon, info[1]]);
                fi;
                surface := s2[1]; 
                params := s2[2];
            fi;
                    
        fi;
        
        return [surface,params];
    end
);

InstallMethod( TestTrianglesConstellation,
    [IsList,IsInt],
    function(vertices, i)
        local temp, pRecord, surf, surf2, VofEdges, EofFaces, a, numb, j;
        # good way to test errors, only working for 2 now
        numb := [2..i];
        surf := SimplicialSurfaceByVerticesInFaces( [[1,2,3]] );
        
        for j in numb do
            # need to do annoying calculations to not have the indices of the vertices shift up, TEST
            
            a := (j-1)*3;
            VofEdges := [];
            EofFaces := [];
            VofEdges[1+a] := [1+(j-1)*3,2+(j-1)*3];
            VofEdges[2+a] := [2+(j-1)*3,3+(j-1)*3];
            VofEdges[3+a] := [1+(j-1)*3,3+(j-1)*3];
            EofFaces[j] := [1+(j-1)*3,2+(j-1)*3,3+(j-1)*3];
            
            surf2 := SimplicialSurfaceByDownwardIncidence( VofEdges, EofFaces );
            temp := DisjointUnion(surf, surf2);
            surf := temp[1];
        od;

        pRecord := SetVertexCoordinates3D(surf, vertices, rec());
        DrawSurfaceToSTL(surf,"test_constellation.txt",pRecord, false);
    end    
);

InstallMethod( EdgeOnFace,
    [IsList, IsList, IsFloat, IsList,  IsFloat],
    function(c_coords, c_normal, d1, t_coords, eps)
        local edge_param, min, l;
        
        for l in [1..3] do
            edge_param := t_coords[l mod 3 + 1]-t_coords[l];
            # test if an edge is perpendicuclar to the face
            if FlEq(AbsoluteValue(edge_param * c_normal),0., eps) then
                # test if a vertex of this edge is on the plane
                min := Minimum(AbsoluteValue(t_coords[l mod 3 + 1]*c_normal-d1),AbsoluteValue(t_coords[l]*c_normal-d1));
                # if one is, the edge is on the plane as well
                if FlEq(min,0.,eps) then
                    return [true,l];
                fi;
            fi;
        od;
        return [false,0];
    end    
);

InstallMethod( CheckSelfIntersections,
    "for a triagular complex without edge ramifications, a filename and a record",
    [IsTriangularComplex and IsNotEdgeRamified, IsRecord, IsList, IsString, IsInt, IsList],
    function(surface, printRecord, Coords, fileName, printSteps, Steps)
        local old_surface, Compared, ints, FaceIntersections, ccoords, t, error, orthog, orthog2, normal, c_normal, c_edges, t_normal, vertsOfFace, edge, edge_param, edge_param2, possb_intersec, failed, ints_points, not_on_edges, not_on_edge, step, steps, name, info, numb_int, one, c_coords, t_coords, check_verts, vOedge, vOedge2, dist, res, res2, lambda, lambda2, val_1, val_2, A, p, p2, x, y, v, d1, d2, min, I, o, o1, o2, i, j, k, l, r, diff, different, step_info_length, maxSteps, print_start, print_stop, stop_step, eps;
            ints := false;
            different := true;            

            ints_points := [];
            Compared := [];
            
            I := ShallowCopy(Faces(surface));
            
            for t in I do
            	Compared[t] := [t];
            od;
            eps := 10.0/(10^(7));
            steps := 0;
		

            # parameters for checking surfaces in certain step ranges
            maxSteps := Steps[1];

            step_info_length := Length(Steps);
            print_start := 0;
            print_stop := 0;
            stop_step := 0;

            if step_info_length = 2 then
                print_start := Steps[2];
            elif step_info_length = 3 then
                print_start := Steps[2];
                print_stop := Steps[3];
            elif step_info_length = 4 then
                print_start := Steps[2];
                print_stop := Steps[3];
                stop_step := Steps[4];
            fi;

            possb_intersec := 0;
            
            
            for j in I do      
                for k in I do
                    # set coords again since the coordinate matrix may change every iteration
                    different := true;

                    if j in I then
                    	if Coords[j] <> [] then
		                c_normal := Coords[j][4];
		                d1 := Coords[j][1]*c_normal;
		            
		                c_coords := [Coords[j][1],Coords[j][2],Coords[j][3]];
                        fi;
                    fi;
                    
                    step := [];
                    t_coords := ShallowCopy([Coords[k][1],Coords[k][2],Coords[k][3]]);
                    # cant intersect if incident

                    if j in I and c_coords <> t_coords and (not j in Compared[k]) and (not k in Compared[j]) then
                        t_normal := Coords[k][4];
                        d2 := Coords[k][1]*t_normal;

                        orthog := [Crossproduct(c_coords[2]-c_coords[1],c_normal),Crossproduct(c_coords[3]-c_coords[2],c_normal),Crossproduct(c_coords[1]-c_coords[3],c_normal)];
                        
                        orthog[1] := orthog[1] / Sqrt(orthog[1]*orthog[1]);
                        orthog[2] := orthog[2] / Sqrt(orthog[2]*orthog[2]);
                        orthog[3] := orthog[3] / Sqrt(orthog[3]*orthog[3]);
                        
                        # must be right of planes, need to have the orthogonal vectors point to the inside of the triangle
                        orthog[1] := orthog[1] * ((orthog[1]*c_coords[3]-orthog[1]*c_coords[1]) / AbsoluteValue(orthog[1]*c_coords[3]-orthog[1]*c_coords[1]));
                        orthog[2] := orthog[2] * ((orthog[2]*c_coords[1]-orthog[2]*c_coords[2]) / AbsoluteValue(orthog[2]*c_coords[1]-orthog[2]*c_coords[2]));
                        orthog[3] := orthog[3] * ((orthog[3]*c_coords[2]-orthog[3]*c_coords[3]) / AbsoluteValue(orthog[3]*c_coords[2]-orthog[3]*c_coords[3]));
                        # check if triangle k intersects with j
                        numb_int := 0;
                        ints_points := [];
                        
                        for l in [1..3] do

                               
                                vOedge := [Coords[k][l],Coords[k][l mod 3 + 1]];
                                
                                # if two faces share a edge than this edge cant be an intersection
                                if Length(NumericalUniqueListOfLists([c_coords[1],c_coords[2],c_coords[3],vOedge[1],vOedge[2]],eps)) >  3 then
                                    edge_param := Coords[k][l mod 3 + 1]-Coords[k][l];
                                    
                                    
                                    # find point that intersects with plane
                                    res := ProjectLineOnTriangle(Coords[k][l],edge_param, d1, c_normal,[orthog,c_coords],eps);
                                    
                                    p := res[1];
                                    lambda := res[2];
                                        
                                    
                                    
                                    o := [orthog[1]*p-orthog[1]*c_coords[1], orthog[2]*p-orthog[2]*c_coords[2], orthog[3]*p-orthog[3]*c_coords[3]];
                                    
                                    # check if point inside triangle
                                     
                                    if FlGeq(Minimum(o[1],o[2],o[3]),0.,eps) and FlGeq(lambda,0.,eps) and FlLeq(lambda,1.,eps) then
                                    
                                        
                                        ints := true;
                                        possb_intersec := possb_intersec + 1;
                                        numb_int := numb_int + 1;
                                        ints_points[numb_int] := p;
                                       	
                                    fi;
                                        
                                    
                                fi;
                             
                        od;
                        
                        if Length(ints_points) > 0 then
                            ints_points := NumericalUniqueListOfLists(ints_points,eps);
                            diff := numb_int - Length(ints_points);
                            possb_intersec := possb_intersec - diff;
                        fi;                        
                        numb_int := Length(ints_points);
			
                
                        if numb_int = 2 then
                            ints_points := NumericalUniqueListOfLists(ints_points,eps);
                            if Length(ints_points) = 2 then
                                # if there is no second point, there is no intersection
                                IsNotEdgeRamified(surface);
                                IsSimplicialSurface(surface);
                                IsTriangularComplex(surface);
                                old_surface := surface;
                                
                                info := [j,k,EdgeOnFace(c_coords, c_normal, d1, t_coords, eps), EdgeOnFace(t_coords, t_normal, d2, c_coords, eps)];
                              
                                step := FixSelfIntersections(surface, printRecord, [Coords,I,Compared], ints_points, info, eps);
                                
                                
                                surface := step[1];
                                Coords := step[2][1];
                                I := step[2][2];
                                Compared := step[2][3];
                                
                                if Length(Faces(surface)) = Length(Faces(old_surface)) then
		                        # wrongly counted
		                        possb_intersec := possb_intersec - numb_int;	
                               	fi;
                               	
                                error := CheckForMistakes(Coords,I,eps);
                                
                                if error[2] then
                                	ThrowError();
                                fi;
    
                                steps := steps + 1;
                                
                                
                                
                                IsNotEdgeRamified(surface);
                                IsSimplicialSurface(surface);
                                IsTriangularComplex(surface);
                                
                                Print(Concatenation("\n",String(steps)));
                                
                                
                                if printSteps <> 0 then
                                    if steps mod printSteps = 0 and steps >= print_start and (steps < print_stop or print_stop = 0) then
                                        name := Concatenation("fails/",fileName, "_step_", String(steps));
                                        Print("\t", possb_intersec," intersections found, ", Length(Faces(surface)), " faces in model\n");
                                        if steps > 900 then
                                        	#eps := eps + 1.0/(10^(9));
                                        fi;
                                        Print("Epsilon is equal to ", eps, "\n");
                                        DrawSTL(name, Coords, true);
                                    fi;
                                fi;

                                if steps = stop_step then
                                    ThrowError();
                                fi;

                                if steps = maxSteps then
                                    name := Concatenation("fails/", fileName, "_step_", String(steps));
                                    Print("\t", possb_intersec," intersections found, ", Length(Faces(surface)), " faces in model\n");
                                    DrawSTL(name, Coords, true);
                                    return [ints,surface, Coords];
                                    
                                fi;
                            else 
                                # wrongly counted
                                possb_intersec := possb_intersec - numb_int;
                                Compared[j][Length(Compared[j])+1] := k;
                            fi;
                            
                        elif numb_int = 1 then
                            # need to find second point that parameterises the intersection
                            
                            d2 := Coords[k][1]*t_normal;
                            
                            orthog2 := [Crossproduct(t_coords[2]-t_coords[1],t_normal),Crossproduct(t_coords[3]-t_coords[2],t_normal),Crossproduct(t_coords[1]-t_coords[3],t_normal)];
                            
                            orthog2[1] := orthog2[1] / Sqrt(orthog2[1]*orthog2[1]);
                            orthog2[2] := orthog2[2] / Sqrt(orthog2[2]*orthog2[2]);
                            orthog2[3] := orthog2[3] / Sqrt(orthog2[3]*orthog2[3]);
                            
                            # correct orientation of normals
                            orthog2[1] := orthog2[1] * ((orthog2[1]*t_coords[3]-orthog2[1]*t_coords[1]) / AbsoluteValue(orthog2[1]*t_coords[3]-orthog2[1]*t_coords[1]));
                            orthog2[2] := orthog2[2] * ((orthog2[2]*t_coords[1]-orthog2[2]*t_coords[2]) / AbsoluteValue(orthog2[2]*t_coords[1]-orthog2[2]*t_coords[2]));
                            orthog2[3] := orthog2[3] * ((orthog2[3]*t_coords[2]-orthog2[3]*t_coords[3]) / AbsoluteValue(orthog2[3]*t_coords[2]-orthog2[3]*t_coords[3]));
                            
                            
                            for l in [1..3] do
                                vOedge2 := [Coords[j][l],Coords[j][l mod 3 + 1]];
                                
                                # if faces are connected by a vertex, and current edge is incident to that vertex,
                                # that edge can't be used to detect intersections
                                
                                if Length(NumericalUniqueListOfLists([t_coords[1],t_coords[2],t_coords[3],vOedge2[1],vOedge2[2]],eps)) >  3 then

                                    edge_param2 := Coords[j][l mod 3 + 1]-Coords[j][l];
                                    
                                    res2 := ProjectLineOnTriangle(Coords[j][l],edge_param2, d2, t_normal,[orthog2,t_coords],eps);
                                    p2 := res2[1];
                                    lambda2 := res2[2];
                                            
                                    o2 := [orthog2[1]*p2-orthog2[1]*t_coords[1], orthog2[2]*p2-orthog2[2]*t_coords[2], orthog2[3]*p2-orthog2[3]*t_coords[3]];
                                            
                                    # check if point inside triangle
                                                                                    
                                    if FlGeq(Minimum(o2[1],o2[2],o2[3]),0.,eps) and FlGeq(lambda2,0.,eps) and FlLeq(lambda2,1.,eps) then
                                    
                                        possb_intersec := possb_intersec + 1;
                                        numb_int := numb_int + 1;
                                        ints_points[numb_int] := p2;
                                    fi;
                                fi;
                            od;
                            
                            # need to to remove (numerical) duplicate entries that are possible in certain geometrical constellations
                            ints_points := NumericalUniqueListOfLists(ints_points,eps);
                            
                            if Length(ints_points) = 2 then
                                # if there is not second point, there is no intersection
                                IsNotEdgeRamified(surface);
                                IsSimplicialSurface(surface);
                                IsTriangularComplex(surface);
                                old_surface := surface;
                                
                                info := [j,k,EdgeOnFace(c_coords, c_normal, d1, t_coords, eps), EdgeOnFace(t_coords, t_normal, d2, c_coords, eps)];
                                
                                step := FixSelfIntersections(surface, printRecord, [Coords,I,Compared], ints_points, info, eps);
                                
                                surface := step[1];
                                Coords := step[2][1];
                                I := step[2][2];
                                Compared := step[2][3];

				if Length(Faces(surface)) = Length(Faces(old_surface)) then
		                        # wrongly counted
		                        possb_intersec := possb_intersec - numb_int;	
                               	fi;
                                
                                error := CheckForMistakes(Coords,I,eps);
                                
                                if error[2] then
                                	ThrowError();
                                fi;
                                steps := steps + 1;
                                
                                
                                
                                IsNotEdgeRamified(surface);
                                IsSimplicialSurface(surface);
                                IsTriangularComplex(surface);
                                Print(Concatenation("\n",String(steps)));
                                if printSteps <> 0 then
                                    if steps mod printSteps = 0 and steps >= print_start and (steps < print_stop or print_stop = 0) then
                                        name := Concatenation("fails/",fileName, "_step_", String(steps));
                                        Print("\t", possb_intersec," intersections found, ", Length(Faces(surface)), " faces in model\n");
                                        if steps > 900 then
                                        	#eps := eps + 1.0/(10^(9));
                                        fi;
                                        Print("Epsilon is equal to ", eps, "\n");
                                        DrawSTL(name, Coords, true);
                                    fi;
                                fi;

                                if steps = stop_step then
                                    ThrowError();
                                fi;


                                if steps = maxSteps then
                                    name := Concatenation("fails/", fileName, "_step_", String(steps));
                                    Print("\t", possb_intersec," intersections found, ", Length(Faces(surface)), " faces in model\n");
                                    DrawSTL(name, Coords, true);
                                    return [ints,surface, Coords];
                                    
                                fi;
                                
                            else 
                                # wrongly counted
                                possb_intersec := possb_intersec - numb_int;
                                Compared[j][Length(Compared[j])+1] := k;
                            fi;
                        else
                        	Compared[j][Length(Compared[j])+1] := k;
                        
                        fi;
                        
                        #only check one triangle since the loop goes through all of them
                        
                    fi;
                od;
            od;
            
            Print("\n", possb_intersec," intersections found in ", String(steps), " iterations");
            
            return [ints,surface, Coords, steps, I];
    end
);


InstallMethod( CountSelfIntersections,
    "for a triagular complex without edge ramifications, a filename and a record",
    [IsTriangularComplex and IsNotEdgeRamified, IsRecord, IsList, IsString,],
    function(surface, printRecord, Coords, fileName)
        local ints, FaceIntersections, ccoords, orthog, orthog2, normal, target, checks_c, c_on_edge1, c_on_edge2, c_on_edges, c_on_edge_both, inters_c_vertices, I, different, t, d, info, Compared, leng1, point_in_plane, current, c_normal, c_edges, t_normal, vertsOfFace, edge, edge_param, edge_param2, possb_intersec, failed, ints_points, not_on_edges, not_on_edge, step, steps, name, numb_int, one, c_coords, t_coords, vOedge, vOedge2, dist, res, res2, lambda, lambda2, val_1, val_2, A, p, p2, x, y, v, d1, d2, min, o, o1, o2, i, j, k, l, intersec_faces, diff, eps;
            # This method counts the number of self-intersections and writes it in the console output
            ints := false;
            different := true;            

            ints_points := [];
            Compared := [];
            
            I := ShallowCopy(Faces(surface));
            
            for t in I do
            	Compared[t] := [t];
            od;
            eps := 10.0/(10^(7));
            steps := 0;
		

            possb_intersec := 0;
            
            
            for j in I do      
                for k in I do
                    # set coords again since the coordinate matrix may change every iteration
                    different := true;
		    
                    if j in I then
                    	if Coords[j] <> [] then
		                c_normal := Coords[j][4];
		                d1 := Coords[j][1]*c_normal;
		            
		                c_coords := [Coords[j][1],Coords[j][2],Coords[j][3]];
                        fi;
                    fi;
                    
                    step := [];
                    t_coords := ShallowCopy([Coords[k][1],Coords[k][2],Coords[k][3]]);
                    # cant intersect if incident

		    current := c_coords;
		    target := t_coords;
		    
                    if j in I and c_coords <> t_coords and (not j in Compared[k]) and (not k in Compared[j]) then
                    	
                    	Compared[j][Length(Compared[j])+1] := k;
                    	
                        t_normal := Coords[k][4];
                        d2 := Coords[k][1]*t_normal;

                        orthog := [Crossproduct(c_coords[2]-c_coords[1],c_normal),Crossproduct(c_coords[3]-c_coords[2],c_normal),Crossproduct(c_coords[1]-c_coords[3],c_normal)];
                        
                        orthog[1] := orthog[1] / Sqrt(orthog[1]*orthog[1]);
                        orthog[2] := orthog[2] / Sqrt(orthog[2]*orthog[2]);
                        orthog[3] := orthog[3] / Sqrt(orthog[3]*orthog[3]);
                        
                        # must be right of planes, need to have the orthogonal vectors point to the inside of the triangle
                        orthog[1] := orthog[1] * ((orthog[1]*c_coords[3]-orthog[1]*c_coords[1]) / AbsoluteValue(orthog[1]*c_coords[3]-orthog[1]*c_coords[1]));
                        orthog[2] := orthog[2] * ((orthog[2]*c_coords[1]-orthog[2]*c_coords[2]) / AbsoluteValue(orthog[2]*c_coords[1]-orthog[2]*c_coords[2]));
                        orthog[3] := orthog[3] * ((orthog[3]*c_coords[2]-orthog[3]*c_coords[3]) / AbsoluteValue(orthog[3]*c_coords[2]-orthog[3]*c_coords[3]));
                        # check if triangle k intersects with j
                        numb_int := 0;
                        ints_points := [];
                        
                        for l in [1..3] do

                               
                                vOedge := [Coords[k][l],Coords[k][l mod 3 + 1]];
                                
                                # if two faces share a edge than this edge cant be an intersection
                                if Length(NumericalUniqueListOfLists([c_coords[1],c_coords[2],c_coords[3],vOedge[1],vOedge[2]],eps)) >  3 then
                                    edge_param := Coords[k][l mod 3 + 1]-Coords[k][l];
                                    
                                    
                                    # find point that intersects with plane
                                    res := ProjectLineOnTriangle(Coords[k][l],edge_param, d1, c_normal,[orthog,c_coords],eps);
                                    
                                    p := res[1];
                                    lambda := res[2];
                                        
                                    
                                    
                                    o := [orthog[1]*p-orthog[1]*c_coords[1], orthog[2]*p-orthog[2]*c_coords[2], orthog[3]*p-orthog[3]*c_coords[3]];
                                    
                                    # check if point inside triangle
                                     
                                    if FlGeq(Minimum(o[1],o[2],o[3]),0.,eps) and FlGeq(lambda,0.,eps) and FlLeq(lambda,1.,eps) then
                                    
                                        
                                        ints := true;
                                        possb_intersec := possb_intersec + 1;
                                        numb_int := numb_int + 1;
                                        ints_points[numb_int] := p;
                                       	
                                    fi;
                                        
                                    
                                fi;
                             
                        od;
                        
                        if Length(ints_points) > 0 then
                            ints_points := NumericalUniqueListOfLists(ints_points,eps);
                            diff := numb_int - Length(ints_points);
                            possb_intersec := possb_intersec - diff;
                        fi;                        
                        numb_int := Length(ints_points);
			
                
                        if numb_int = 2 then
                            ints_points := NumericalUniqueListOfLists(ints_points,eps);
                            if Length(ints_points) = 2 then
                                # if there is no second point, there is no intersection
                                IsNotEdgeRamified(surface);
                                IsSimplicialSurface(surface);
                                IsTriangularComplex(surface);
                                
                                info := [j,k,EdgeOnFace(c_coords, c_normal, d1, t_coords, eps), EdgeOnFace(t_coords, t_normal, d2, c_coords, eps)];
                              	d := current[1]*c_normal;
                              	
                              	c_on_edge1 := OnEdges(current,ints_points[1],eps);
				c_on_edge2 := OnEdges(current,ints_points[2],eps);
				c_on_edge_both := [c_on_edge1[1] and c_on_edge2[1],c_on_edge1[2] and c_on_edge2[2],c_on_edge1[3] and c_on_edge2[3]];
				c_on_edges := c_on_edge_both[1] or c_on_edge_both[2] or c_on_edge_both[3];
				
        			inters_c_vertices := Length(NumericalUniqueListOfLists([current[1],current[2],current[3],ints_points[1],ints_points[2]],eps));
        			
				checks_c := inters_c_vertices  <=  3 or c_on_edges;
        			
                                if checks_c or info[4][1] then
			        	# the edge of the intersection is perpendicular to our face
			        	# so either all three or only two of the target vertices lie in the current triangle plane
			        	point_in_plane := [FlEq(target[1]*c_normal,d,eps),FlEq(target[2]*c_normal,d,eps),FlEq(target[3]*c_normal,d,eps)];
			    
			        	if not (point_in_plane[1] and point_in_plane[2] and point_in_plane[3]) then
						# the intersection is along an edge
						# test if intersection spans a full edge
						l1 := info[4][2];
						leng1 := Length(NumericalUniqueListOfLists([ints_points[1],ints_points[2],current[l1 mod 3 + 1],current[l1]],eps));
				
						if leng1 = 2 then
				 	 		# not an intersection
				 	 		possb_intersec := possb_intersec - numb_int;	
				 	 	fi;
					fi;
                                fi;
                                
                                
    
                                steps := steps + 1;
                                
                                
                                
                                IsNotEdgeRamified(surface);
                                IsSimplicialSurface(surface);
                                IsTriangularComplex(surface);
                                

                            else 
                                # wrongly counted
                                possb_intersec := possb_intersec - numb_int;
                                Compared[j][Length(Compared[j])+1] := k;
                            fi;
                            
                        elif numb_int = 1 then
                            # need to find second point that parameterises the intersection
                            
                            d2 := Coords[k][1]*t_normal;
                            
                            orthog2 := [Crossproduct(t_coords[2]-t_coords[1],t_normal),Crossproduct(t_coords[3]-t_coords[2],t_normal),Crossproduct(t_coords[1]-t_coords[3],t_normal)];
                            
                            orthog2[1] := orthog2[1] / Sqrt(orthog2[1]*orthog2[1]);
                            orthog2[2] := orthog2[2] / Sqrt(orthog2[2]*orthog2[2]);
                            orthog2[3] := orthog2[3] / Sqrt(orthog2[3]*orthog2[3]);
                            
                            # correct orientation of normals
                            orthog2[1] := orthog2[1] * ((orthog2[1]*t_coords[3]-orthog2[1]*t_coords[1]) / AbsoluteValue(orthog2[1]*t_coords[3]-orthog2[1]*t_coords[1]));
                            orthog2[2] := orthog2[2] * ((orthog2[2]*t_coords[1]-orthog2[2]*t_coords[2]) / AbsoluteValue(orthog2[2]*t_coords[1]-orthog2[2]*t_coords[2]));
                            orthog2[3] := orthog2[3] * ((orthog2[3]*t_coords[2]-orthog2[3]*t_coords[3]) / AbsoluteValue(orthog2[3]*t_coords[2]-orthog2[3]*t_coords[3]));
                            
                            
                            for l in [1..3] do
                                vOedge2 := [Coords[j][l],Coords[j][l mod 3 + 1]];
                                
                                # if faces are connected by a vertex, and current edge is incident to that vertex,
                                # that edge can't be used to detect intersections
                                
                                if Length(NumericalUniqueListOfLists([t_coords[1],t_coords[2],t_coords[3],vOedge2[1],vOedge2[2]],eps)) >  3 then

                                    edge_param2 := Coords[j][l mod 3 + 1]-Coords[j][l];
                                    
                                    res2 := ProjectLineOnTriangle(Coords[j][l],edge_param2, d2, t_normal,[orthog2,t_coords],eps);
                                    p2 := res2[1];
                                    lambda2 := res2[2];
                                            
                                    o2 := [orthog2[1]*p2-orthog2[1]*t_coords[1], orthog2[2]*p2-orthog2[2]*t_coords[2], orthog2[3]*p2-orthog2[3]*t_coords[3]];
                                            
                                    # check if point inside triangle
                                                                                    
                                    if FlGeq(Minimum(o2[1],o2[2],o2[3]),0.,eps) and FlGeq(lambda2,0.,eps) and FlLeq(lambda2,1.,eps) then
                                    
                                        possb_intersec := possb_intersec + 1;
                                        numb_int := numb_int + 1;
                                        ints_points[numb_int] := p2;
                                    fi;
                                fi;
                            od;
                            
                            # need to to remove (numerical) duplicate entries that are possible in certain geometrical constellations
                            ints_points := NumericalUniqueListOfLists(ints_points,eps);
                            
                            if Length(ints_points) = 2 then
                                # if there is not second point, there is no intersection
                                IsNotEdgeRamified(surface);
                                IsSimplicialSurface(surface);
                                IsTriangularComplex(surface);
                                
                                info := [j,k,EdgeOnFace(c_coords, c_normal, d1, t_coords, eps), EdgeOnFace(t_coords, t_normal, d2, c_coords, eps)];
                                d := current[1]*c_normal;
                                
                                c_on_edge1 := OnEdges(current,ints_points[1],eps);
				c_on_edge2 := OnEdges(current,ints_points[2],eps);
				c_on_edge_both := [c_on_edge1[1] and c_on_edge2[1],c_on_edge1[2] and c_on_edge2[2],c_on_edge1[3] and c_on_edge2[3]];
				c_on_edges := c_on_edge_both[1] or c_on_edge_both[2] or c_on_edge_both[3];
				
        			inters_c_vertices := Length(NumericalUniqueListOfLists([current[1],current[2],current[3],ints_points[1],ints_points[2]],eps));
        			
				checks_c := inters_c_vertices  <=  3 or c_on_edges;
        			
                                if checks_c or info[4][1] then
			        	# the edge of the intersection is perpendicular to our face
			        	# so either all three or only two of the target vertices lie in the current triangle plane
			        	point_in_plane := [FlEq(target[1]*c_normal,d,eps),FlEq(target[2]*c_normal,d,eps),FlEq(target[3]*c_normal,d,eps)];
			    
			        	if not (point_in_plane[1] and point_in_plane[2] and point_in_plane[3]) then
						# the intersection is along an edge
						# test if intersection spans a full edge
						l1 := info[4][2];
						leng1 := Length(NumericalUniqueListOfLists([ints_points[1],ints_points[2],current[l1 mod 3 + 1],current[l1]],eps));
				
						if leng1 = 2 then
				 	 		# not an intersection
				 	 		possb_intersec := possb_intersec - numb_int;		
				 	 	fi;
					fi;
                                fi;
                               	
                                steps := steps + 1;
                                
                                
                                
                                IsNotEdgeRamified(surface);
                                IsSimplicialSurface(surface);
                                IsTriangularComplex(surface);
                                
                                
                            else 
                                # wrongly counted
                                possb_intersec := possb_intersec - numb_int;
                                Compared[j][Length(Compared[j])+1] := k;
                            fi;
                        else
                        	Compared[j][Length(Compared[j])+1] := k;
                        
                        fi;
                        
                        #only check one triangle since the loop goes through all of them
                        
                    fi;
                od;
            od;
            
            Print("\n", possb_intersec," intersection points found");
            
            return [ints,surface, Coords, steps, I];
    end
);



InstallMethod(DrawSTL,
    [IsString, IsList, IsBool],
    function(fileName, Coords, write)
        local file, filesepr, name, output, x,y, i, j, k, r,l, coords, Copy_Coords, normal, eps;

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

        Coords := RemoveDuplicateFaces(Coords,write,eps);
        
        for i in [1..Length(Coords)] do 
            if IsBound(Coords[i]) and Coords[i] <> [] then
            
                # get coords of vertices
                coords := [Coords[i][1],Coords[i][2],Coords[i][3]];
                x := Coords[i][2]-Coords[i][1];
                y := Coords[i][3]-Coords[i][1];
                normal := Crossproduct(x,y);
                normal := normal / Sqrt(normal*normal);
                
                # write normal
                AppendTo(output, "\tfacet normal ");
                for j in [1..3] do
                    AppendTo(output, Concatenation(String(normal[j])," "));
                od;
                AppendTo(output, "\n");
                    
                # write vertex coords
                AppendTo(output, "\t\touter loop\n");
                for j in [1..3] do
                    AppendTo(output,"\t\t\tvertex ");
                       
                    for k in [1..3] do
                        AppendTo(output, Concatenation(String(coords[j][k])," "));
                    od;
                    AppendTo(output,"\n");
                od;
                AppendTo(output, "\t\tendloop\n");
            AppendTo(output,"\tendfacet\n");
            fi;
        od;
        AppendTo(output, Concatenation("endsolid ", name));
        if write then
        	Print("\n Saved file");
        fi;
        CloseStream(output);
        return;
    end
);

InstallMethod( DrawSurfaceToSTLCalculate,
    "for a polygonal complex without edge ramifications, a filename and a record",
    [IsTriangularComplex and IsNotEdgeRamified, IsString, IsRecord, IsBool, IsInt, IsList],
    function(surface, fileName, printRecord, test, printSteps, Steps)
        local file, filesepr, name, output, Coords, intersections, ccoords, new, old_surface, new_surface, new_coords, coords, i, j, k, vertsOfFace, coord1, coord2, coord3, x, y, normal, intersects;
        
        #######################################################################################################################################################################
        #   INPUTS
        ##
        # A triangular or more commonly simplicial surface that is to be analysed for intersections, a name for the created STL files, a printRecord with the coordinate rmation
        # of the vertices, a bool test if the surface should be tested on self-intersections, a int printSteps at which interval intermediate results are to be saved as STL files,
        # a list Steps containg information about at which amount of steps to stop the fixing of self-intersections 
        ##
        #   METHOD
        ##
        # The method iself gets the coordinate information saved in the printRecord, converts it into a list format and then calls the main method to triangulate the surface.
        #
        #######################################################################################################################################################################
    
            Coords := [];
            old_surface := surface;
            # get data
            for i in Faces(surface) do
                    vertsOfFace := VerticesAsList(PerimeterPathOfFace(surface,i));
                    # get coords of current vertices
                    ccoords := [GetVertexCoordinates3DNC(surface,vertsOfFace[1], printRecord), GetVertexCoordinates3DNC(surface, vertsOfFace[2], printRecord),GetVertexCoordinates3DNC(surface, vertsOfFace[3], printRecord)];
                    # get normal 
                    x := ccoords[2]-ccoords[1];
                    y := ccoords[3]-ccoords[1];
                    normal := Crossproduct(x,y);
                    normal := normal / Sqrt(normal*normal);
                    ccoords[4] := normal;
                    ccoords[5] := vertsOfFace;
                    
                    Coords[i] := ccoords;
            od;
            
            
            # check for self intersections
            
            if test then
            	DrawSTL(Concatenation("base_models/", fileName,"_unchanged"), Coords, true);
            	
            	Print("\n",fileName, " number of intersections:");
            	CountSelfIntersections(surface, printRecord, Coords, fileName);
            	
                intersections := CheckSelfIntersections(surface, printRecord, Coords, fileName, printSteps, Steps);
                intersects := intersections[1];
                # CheckSelfIntersections(intersections[2], printRecord, intersections[3], fileName);
                surface := intersections[2];
                Coords := intersections[3];
                if intersects then
                    Print("\nSelf-Intersections present, produced ", Length(Faces(surface)) - Length(Faces(old_surface)), " new faces");
                fi;
                if Steps[1] = 0 or intersections[4] < Steps[1] then
                    DrawSTL(Concatenation(fileName,"_fixed"), Coords, true);
                fi;
		
		new := SimplicialSurfaceFromCoordinates([Coords,intersections[5]],1./(10^8));
            	new_surface := new[1];
            	new_coords := new[2];
            
            	return [old_surface, new_surface, new_coords];
                
            else
                Print("\n",fileName, " number of intersections:");
                intersections := CountSelfIntersections(surface, printRecord, Coords, fileName);
                intersects := intersections[1];
                # CheckSelfIntersections(intersections[2], printRecord, intersections[3], fileName);
                surface := intersections[2];
                Coords := intersections[3];
                DrawSTL(Concatenation(fileName,""), Coords, false);
                
                return [surface];
            fi;
            
            
    end
);

InstallMethod( DrawSurfaceToSTL,
    "for a polygonal complex without edge ramifications, a filename and a record",
    [IsTriangularComplex and IsNotEdgeRamified, IsString, IsRecord, IsBool],
    function(surface, fileName, printRecord, test)
        return DrawSurfaceToSTLCalculate(surface,fileName,printRecord, test, 5000,[0]);
    end
);


InstallMethod( SimplicialSurfaceFromCoordinates,
    [IsList,IsFloat],
    function(params,eps)
        local Coords, faces, f, i, j, l, pos, VerticesInFaces, VerticesCoords, verts, surf;
        Coords := params[1];
        faces := params[2];
        
        # vertex counter
        i := 1;

        # face counter
        j := 1;
        
        VerticesInFaces := [];
        VerticesCoords := [];
        
        for f in faces do
            verts := [];

            for l in [1,2,3] do 
                pos := NumericalPosition(VerticesCoords,Coords[f][l],eps);
            
                if pos = fail then
                    # vertex coord. is new
                    VerticesCoords[i] := Coords[f][l];
                    
                    verts[l] := i;

                    i := i + 1;
                
                else
                    verts[l] := pos;
                fi;

                
            od;

            VerticesInFaces[j] := verts;
            j := j + 1;
        od;
        
        surf := TriangularComplexByVerticesInFaces(VerticesInFaces);

        return [surf,VerticesCoords];
    end
);

InstallMethod( SimplicialSurfaceFromChangedCoordinates,
    [IsList,IsFloat],
	function(params,eps)
        local Coords, faces, f, i, j, l, pos, old_surf, VerticesInFaces, VerticesCoords, verts, surf;
        Coords := params[1];
        old_surf := params[2];
        faces := ShallowCopy(Faces(old_surf));
        
        VerticesInFaces := [];
        
        
        for f in faces do
            verts := Coords[f][5];

            VerticesInFaces[f] := verts;

        od;
        
        surf := TriangularComplexByVerticesInFaces(VerticesInFaces);

        return [surf];
    end
);



