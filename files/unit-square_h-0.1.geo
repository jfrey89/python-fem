/*********************************************************************
 *
 *  Gmsh tutorial 1
 *
 *  Variables, elementary entities (points, lines, surfaces), physical
 *  entities (points, lines, surfaces)
 *
 *********************************************************************/

h = 1e-1 ;

Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h} ;
Point(3) = {1, 1, 0, h} ;
Point(4) = {0, 1, 0, h} ;

Line(1) = {1,2} ;
Line(2) = {3,2} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;
Line Loop(5) = {4,1,-2,3} ;
Plane Surface(6) = {5} ;
