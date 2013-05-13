z_coord = 0.0;
 
// define a mesh characteristic length (cl) - smaller=more refined mesh
cl = 0.025;
 
// define points to describe your 2D domain
Point(1)={0.0,0.0,z_coord,cl};
Point(2)={1.0,0.0,z_coord,cl};
Point(3)={1.0,1.0,z_coord,cl};
Point(4)={0.0,1.0,z_coord,cl};
 
// make these points into lines
Line(1)={1,2}; //Line(#)={Start Point,End Point}
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
 
Line Loop(6)={1,2,3,4};
 
Plane Surface(7)={6}; //Plane Surface(#)={boundary Line Loop,interior hole Line Loop}
 
Physical Line(8)={1};//bottom boundary
Physical Line(9)={2};//right boundary
Physical Line(10)={3};//top boundary
Physical Line(11)={4};//left boundary
 
//You also have to give a name to the Plane surface
Physical Surface(12)={7};
