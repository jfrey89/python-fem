// define some variables to describe your domain:
x_left = 0.0;
x_right = 9.0;
y_bottom = 0.0;
y_top = 3.0;
 
z_coord = 0.0;
 
// define a mesh characteristic length (cl) - smaller=more refined mesh
cl = 0.1;
 
// define points to describe your 2D domain
Point(1)={x_left,y_bottom,z_coord,cl};
Point(2)={x_right,y_bottom,z_coord,cl};
Point(3)={x_right,y_top,z_coord,cl};
Point(4)={x_left,y_top,z_coord,cl};
 
// make these points into lines
Line(1)={1,2}; //Line(#)={Start Point,End Point}
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Line Loop(5)={1,2,3,4}; //Line Loop(#)={Line Segment 1,...,Line Segment n}
Plane Surface(6)={5};
 
//Name the surfaces to allow assignment of boundary conditions
Physical Line(7)={1};//bottom boundary
Physical Line(8)={2};//right boundary
Physical Line(9)={3};//top boundary
Physical Line(10)={4};//left boundary
 
//You also have to give a name to the Plane surface
Physical Surface(11)={6};
