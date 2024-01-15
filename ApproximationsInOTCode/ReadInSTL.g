ReadSTL:=function(fileName)
	# reads a file from the current dir
	local surf, file, name, r, r2, eps, filesepr, endsign, normal, data, normals, points, test, i,j, index, verts, coords, input, Coords;
	eps := 1./10^6;
	filesepr := SplitString(fileName, ".");
        name := filesepr[1];
        file := Filename( DirectoryCurrent(), Concatenation(name,".stl") );
        points := [];
        Coords:=[];
        faces:=[];
        i := 1;
       	# test file name
	if IsReadableFile(file) then
		
        input := InputTextFile(file);
		r := ReadLine(input);
		endsign := SplitString(ShallowCopy(r), " ")[1];
		
		while not endsign = "endsolid" do
			
			
			r := ReadLine(input);
			r := Concatenation(SplitString(ShallowCopy(r),"\n"));
			if r="" then
				continue;
			fi;
			while r[1]=' ' do
				Remove(r,1);
			od;
			r2 := SplitString(ShallowCopy(r), " ");
			endsign := r2[1];
			if not endsign = "endsolid" then
				Coords[i]:=[];
				# TODO:  maybe find way to round less?
				normal := [Float(r2[3]),Float(r2[4]),Float(r2[5])];
				Coords[i][4]:=normal;
				Print(i,"\n");
				r := ReadLine(input);
				
				j := 1;
				verts := [];
				face:=[];
				while j < 4 do
					r := ReadLine(input);
					r := Concatenation(SplitString(ShallowCopy(r),"\n"));
					while r[1]=' ' do
						Remove(r,1);
					od;
					r2 := SplitString(ShallowCopy(r), " ");
					coords := [Float(r2[2]),Float(r2[3]),Float(r2[4])];
					
					#test := ShallowCopy(points);
					#Add(test,coords);
					entry:=NumericalPosition(points,coords,eps);
					if entry = fail then
						Add(points,coords);
						index := Length(points); 
					else
						index := entry;
					fi;
					Add(face,index);
					verts[j] := index;
					Coords[i][j] := coords;
					j := j+1;
				od;
				if Size(Set(face))=3 then
					Add(faces,face);
				fi;
				Coords[i][5] := verts;
				r := ReadLine(input);
				r := ReadLine(input);
				i := i + 1;
			fi;
			
		od;
		
	else
		Print("file does not exist");
		return 0;
	fi;
	
	#faces := [1..Length(Coords)];
	#data := SimplicialSurfaceFromCoordinates([Coords,faces],eps);
	return [faces,points];
end;
