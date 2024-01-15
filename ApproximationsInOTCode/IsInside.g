IsInside:=function(t,points,point,eps)
	local max_x,p,line,count,face;
	max_x:=0.;
	for p in points do
		if p[1]>max_x then
			max_x:=p[1];
		fi;
	od;


	line:=[point,[max_x+100.,0,0]];
	count:=0;
	for face in VerticesOfFaces(t) do
		if IntersectLineFace(line,points{face},eps) then
			count:=count+1;
			#Print(face,"\n");
		fi;
	od;
	#Print(count,"\n");
	return count mod 2 = 1;
end;;

# returns true when p lies inside the triangle given by a,b,c
# use https://stackoverflow.com/questions/995445/determine-if-a-3d-point-is-within-a-triangle
MyPointInTriangle:=function(a,b,c,p,eps)
local n1,n2,n3;
# if N1=0 or N2=0 or N3=0 check if P lies on the line
#1. the unit normal of triange (A, B, P)  - call it N1
n1:=Crossproduct(p-a,b-a);
if Norm2(n1)<eps then
	if Norm2(p-a)+Norm2(p-b)<=Norm2(a-b)+eps then
		return true;
	else
		return false;
	fi;
else
	n1:=n1/Norm2(n1);
fi;
#2. the unit normal of triangle (B, C, P) - call it N2
n2:=Crossproduct(p-b,c-b);
if Norm2(n2)<eps then
	if Norm2(p-b)+Norm2(p-c)<=Norm2(b-c)+eps then
		return true;
	else
		return false;
	fi;
else
	n2:=n2/Norm2(n2);
fi;
#3.  the unit normal (C,A,P) called N3
n3:=Crossproduct(p-c,a-c);
if Norm2(n3)<eps then
	if Norm2(p-a)+Norm2(p-c)<=Norm2(a-c)+eps then
		return true;
	else
		return false;
	fi;
else
	n3:=n3/Norm2(n3);
fi;

#N1*N2 == 1.0 ?
#N2*N3 == 1.0 ?  
if AbsoluteValue(Dot(n1,n2)-1)<=eps and AbsoluteValue(Dot(n2,n3)-1)<=eps then
	return true;
fi;

return false;

end;;

# Define a function to check if a line intersects a face
IntersectLineFace := function(line, face, eps)
    local d, t, normal, v1, v2, v3, edge1, edge2, intersection_point;

    # Extract the vertices of the face
    v1 := face[1];
    v2 := face[2];
    v3 := face[3];

    # Calculate the normal vector of the face using the cross product of two edges
    edge1 := v2 - v1;
    edge2 := v3 - v1;
    normal := Crossproduct(edge1, edge2);

    # Calculate the signed distance from the origin to the plane of the face
    d := -Dot(normal, v1);

    # Calculate the dot product between the line direction and the face normal
    t := Dot(normal, line[2] - line[1]);

    # Check if the line is parallel to the face
    if AbsoluteValue(t) < eps then
        return false; # The line is parallel to the face
    fi;

    # Calculate the parameter t at which the line intersects the plane of the face
    t := (d + Dot(normal, line[1])) / t;

    # Check if the intersection point is within the line segment
    if t < 0. or t > 1. then
        return false; # The intersection point is outside the line segment
    fi;

    # Calculate the intersection point
    intersection_point := line[1] - t * (line[2] - line[1]);

    # Check if the intersection point is inside the face
    if MyPointInTriangle(v1, v2, v3, intersection_point, eps) then
        return true; # The line intersects the face
    else
        return false; # The intersection point is outside the face
    fi;
end;

