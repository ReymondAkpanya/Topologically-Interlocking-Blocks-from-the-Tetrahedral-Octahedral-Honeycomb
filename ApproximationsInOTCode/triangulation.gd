#############################################################################
##
##  SimplicialSurface package
##
##  Copyright 2012-2019
##    Markus Baumeister, RWTH Aachen University
##    Alice Niemeyer, RWTH Aachen University
##    Jens Brandt, RWTH Aachen University
##
## Licensed under the GPL 3 or later.
##
#############################################################################

#! @Chapter 3D Printing
#! @ChapterLabel STL files and self intersections

#! @BeginGroup Miscanelous methods
#! @Description
#! Many small method that are needed for triangulation
#! Cross Product between two 3-D vectors
#! @Arguments x,y
DeclareOperation( "Crossproduct",              [IsList, IsList] );
#! Euclidean norm of a 3-D vector
#! @Arguments x
DeclareOperation( "Norm2",                     [IsList] );
#! Surface of 3-D triangle
#! @Arguments x
DeclareOperation( "SurfaceTriangle",           [IsList,IsList,IsList] );
#! Signum function for float
#! @Arguments x,y,z
DeclareOperation( "Signum",                    [IsFloat]);
#! Comparisons for floats and vectors of floats (up to an epsilon-error)
#! @Arguments x,y,epsilon
DeclareOperation( "FlEq",                      [IsFloat,IsFloat,IsFloat] );
#! @Arguments x,y,epsilon
DeclareOperation( "FlVEq",                     [IsList,IsList,IsFloat]);
#! @Arguments x,y,epsilon
DeclareOperation( "FlLeq",                     [IsFloat,IsFloat,IsFloat]);
#! @Arguments x,y,epsilon
DeclareOperation( "FlGeq",                     [IsFloat,IsFloat,IsFloat]);
#! Numerical Position of a vector in a list (so up to epsilon norm difference)
#! @Arguments list,entry,epsilon
DeclareOperation( "NumericalPosition",         [IsList,IsList,IsFloat]);
#! Remove all the numerical non-unique vectors in a list (so up to epsilon norm difference)
#! @Returns a new list in which every vector is numerically unique
#! @Arguments list,entry,epsilon
DeclareOperation( "NumericalUniqueListOfLists",[IsList, IsFloat]);
#! Calculates if three points are along one line (up to epsilon)
#! @Returns a boolean
#! @Arguments v1, v2, v3, eps
DeclareOperation( "PointsInOneLine",[IsList, IsList, IsList, IsFloat]);
#! Checks if there are any errors in the calculated triangulation till now
#! @Returns   [index,boolean]
#! @Arguments verts, I, eps
DeclareOperation( "CheckForMistakes",   [IsList, IsList, IsFloat]);
#! @EndGroup

#! @BeginGroup Remove Duplicate Faces
#! @Description
#! A small method that removes any faces whose vertices lie in the same positions
#! @Arguments Coords,eps
DeclareOperation( "RemoveDuplicateFaces",      [IsList,IsBool,IsFloat]);
#! @EndGroup


#! @BeginGroup ProjectLineOnTriangle 
#! @Description
#! Calculates the point where the line x+t*l hits the plane specified by the point d and the normal
#! @Arguments x,l,d,normal
DeclareOperation( "ProjectLineOnTriangle",     [IsList, IsList, IsFloat, IsList,IsList, IsFloat]);
#! @EndGroup

#! @BeginGroup PlaneEquation
#! @Description
#! Calculates the plane equation for a point wrt to a normal and a vector on the plane
#! @Arguments point,normal,v1
DeclareOperation( "PlaneEquation",    [IsList,IsList,IsList]);
#! @Arguments point,normal,v1,eps
DeclareOperation( "InPlane",    [IsList,IsList,IsList,IsFloat]);
#! @EndGroup

DeclareOperation( "OnEdges",    [IsList,IsList,IsFloat]);
#! @BeginGroup PointInTriangle
#! @Description
#! Calculates if the point x is in the triangle spanned by a1, a2 and a3. x is assumed to be in the plane of the triangle
#! The calculation is complicated
#! @Arguments a1, a2, a3, x, eps
DeclareOperation( "PointInTriangle",           [IsList, IsList, IsList, IsList, IsFloat]);
#! @EndGroup

#! @BeginGroup TriangleFullyInside
#! @Description
#! Checks if the triangle target is fully inside the triangle current by using PointInTriangle
#! If so, we will delete it
#! @Arguments target, current, eps
DeclareOperation( "TriangleFullyInside",    [IsList, IsList, IsFloat]);
#! @EndGroup

#! @BeginGroup PushLineToEdges
#! @Description
#! Takes the line x + t*l and calculates the nearest point where it hits the edge of the triangle spanned by a1, a2 and a3. x is assumed to be in the plane of the triangle
#! This is mainly used to calcuate how to divide triangles when intersections occur
#! @Arguments a1, a2, a3, x, l, epsilon
DeclareOperation( "PushLineToEdges",           [IsList, IsList, IsList, IsList, IsList, IsFloat]);
#! @EndGroup

#! @BeginGroup EdgePerpendicularToFace
#! @Description
#! Test if an edge is perpendicuclar to the face
#! @Arguments c_coords, c_normal, d1, t_coords, eps
DeclareOperation( "EdgeOnFace",   [IsList, IsList, IsFloat, IsList, IsFloat]);
#! @EndGroup

#! @BeginGroup CheckSelfIntersections
#! @Description
#! This is the main method that loops through all the faces and checks for and repairs self intersections
#! @Arguments surface, printRecord, Coords, fileName, printSteps, Steps
DeclareOperation( "CheckSelfIntersections",    [IsTriangularComplex and IsNotEdgeRamified, IsRecord, IsList, IsString, IsInt, IsList] );
#! @EndGroup

#! @BeginGroup CountSelfIntersections
#! @Description
#! This version of the main method just checks for intersections and does not fix them
#! @Arguments surface, printRecord, Coords, fileName
DeclareOperation( "CountSelfIntersections",    [IsTriangularComplex and IsNotEdgeRamified, IsRecord, IsList, IsString] );
#! @EndGroup

#! @BeginGroup DrawSTL
#! @Description
#! This takes a list of coordinates of faces and a string to save an stl file of the described surface
#! @Arguments fileName, Coords
DeclareOperation( "DrawSTL",                   [IsString, IsList,IsBool] );
#! @EndGroup

#! @BeginGroup DrawSTL
#! @Description
#! This starts the main method with default configuration
#! @Arguments surface, fileName, printRecord, test
DeclareOperation( "DrawSurfaceToSTL",          [IsTriangularComplex and IsNotEdgeRamified, IsString, IsRecord, IsBool] );
#! @Arguments surface, fileName, printRecord, test, printSteps, Steps
DeclareOperation( "DrawSurfaceToSTLCalculate", [IsTriangularComplex and IsNotEdgeRamified, IsString, IsRecord, IsBool, IsInt, IsList] );
#! @EndGroup


#! @BeginGroup TestTrianglesConstellation
#! @Description
#! Kind of depricated method to test errors
#! @Arguments vertices, i
DeclareOperation( "TestTrianglesConstellation",[IsList,IsInt] );
#! @EndGroup

#! @BeginGroup FixSelfIntersections
#! @Description
#! This method is called after a intersection is detected.
#! It finds out how the geometric constellation of the intersection is and then calles the appropriate method to fix the intersection.
#! @Arguments surface, printRecord, params, inters, info, epsilon
DeclareOperation( "FixSelfIntersections",      [IsTriangularComplex and IsNotEdgeRamified, IsRecord, IsList, IsList, IsList, IsFloat]);
#! @EndGroup

#! @BeginGroup RetriangulateEasy
#! @Description
#! In this case, both of the intersection points lie inside the triangle (so on no edge or vertex)
#! Depending on where the intersection line pushes to, the triangle will be intersected in 2-6 parts
#! @Arguments surf, reco, params, inters, res, li
DeclareOperation( "RetriangulateEasy",         [IsTriangularComplex and IsNotEdgeRamified, IsRecord, IsList, IsList, IsList, IsList]);
#! @EndGroup

#! @BeginGroup Retriangulate
#! @Description
#! In this case, one or both of the intersection points lie on an edge or vertex.
#! We always split up the triangle in four parts.
#! @Arguments surf, reco, params, inters, res, li
DeclareOperation( "Retriangulate",             [IsTriangularComplex and IsNotEdgeRamified, IsRecord, IsList, IsList, IsList, IsList]);
#! @EndGroup

#! @BeginGroup SplitAdjacent
#! @Description
#! This method is used for the clean-up after fixing an intersection by splitting a triangle.
#! It splits triangles adjacent to split triangles as well so the resulting surface is still well-formed
#! @Arguments surf, params, edge_v1, edge_v2, split_v, info
DeclareOperation( "SplitAdjacent",             [IsTriangularComplex and IsNotEdgeRamified, IsList, IsList, IsList, IsList, IsList]);
#! @EndGroup

#! @BeginGroup SimplicialSurfaceFromCoordinates
#! @Description
#! After fixing the intersections in a surface, this method is called to create a simplicial surface from the resulting coordinates
#! @Arguments params,eps
DeclareOperation( "SimplicialSurfaceFromCoordinates", [IsList,IsFloat]);
#! @EndGroup

#! @BeginGroup SimplicialSurfaceFromCoordinates
#! @Description
#! After changing Coordinates of a simplicial surface, the old surface can be updated using this method
#! @Arguments params,eps
DeclareOperation( "SimplicialSurfaceFromChangedCoordinates", [IsList,IsFloat]);
#! @EndGroup






