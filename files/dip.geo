// define some variables to describe your domain:
x_left = 0.0;
x_right = 5.0;
y_bottom = 0.0;
y_top = 2.0;
left_corner = 2.0;
right_corner = 3.0;
dip = -1.0;
 
z_coord = 0.0;
 
// define a mesh characteristic length (cl) - smaller=more refined mesh
cl = 0.15;
 
// define points to describe your 2D domain
Point(1)={x_left,y_bottom,z_coord,cl};
Point(2)={left_corner,y_bottom,z_coord,cl};
Point(3)={left_corner,dip,z_coord,cl};
Point(4)={right_corner,dip,z_coord,cl};
Point(5)={right_corner,y_bottom,z_coord,cl};
Point(6)={x_right,y_bottom,z_coord,cl};
Point(7)={x_right,y_top,z_coord,cl};
Point(8)={x_left,y_top,z_coord,cl};
 
// make these points into lines
Line(1)={1,2}; //Line(#)={Start Point,End Point}
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,5};
Line(5)={5,6};
Line(6)={6,7};
Line(7)={7,8};
Line(8)={8,1};
Line Loop(9)={1,2,3,4,5,6,7,8};
Plane Surface(10)={9};
 
//Name the surfaces to allow assignment of boundary conditions
Physical Line(11)={1};//bottom boundary
Physical Line(12)={2};//right boundary
Physical Line(13)={3};//top boundary
Physical Line(14)={4};//left boundary
Physical Line(15)={5};//left boundary
Physical Line(16)={6};//left boundary
Physical Line(17)={7};//left boundary
Physical Line(18)={8};//left boundary
 
//You also have to give a name to the Plane surface
Physical Surface(19)={10};
