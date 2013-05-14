// define some variables to describe your domain:
x_left = 0.0;
x_right = 5.0;
y_bottom = 0.0;
y_top = 3.0;
 
z_coord = 0.0;
 
// define a mesh characteristic length (cl) - smaller=more refined mesh
cl = 0.05;
 
// define points to describe your 2D domain
Point(1)={x_left,y_bottom,z_coord,cl};
Point(2)={x_right,y_bottom,z_coord,cl};
Point(3)={x_right,y_top,z_coord,cl};
Point(4)={x_left,y_top,z_coord,cl};
Point(5)={x_left+1,y_bottom,z_coord,cl};
Point(6)={x_right-1,y_bottom,z_coord,cl};
Point(7)={x_left+1,y_bottom+1,z_coord,cl};
Point(8)={x_right-1,y_bottom+1,z_coord,cl};
Point(9)={x_left+1,y_top,z_coord,cl};
Point(10)={x_right-1,y_top,z_coord,cl};
Point(11)={x_left+1,y_top-1,z_coord,cl};
Point(12)={x_right-1,y_top-1,z_coord,cl};
 
// make these points into lines
Line(1)={1,5}; //Line(#)={Start Point,End Point}
Line(2)={5,7};
Line(3)={7,8};
Line(4)={8,6};
Line(5)={6,2};
Line(6)={2,3};
Line(7)={3,10};
Line(8)={10,12};
Line(9)={12,11};
Line(10)={11,9};
Line(11)={9,4};
Line(12)={4,1};
 
// Now connect the lines of the outer rectangle into a line loop.
Line Loop(13)={1,2,3,4,5,6,7,8,9,10,11,12}; //Line Loop(#)={Line Segment 1,...,Line Segment n}
 
//Create a plane surface that can be meshed:
Plane Surface(14)={13}; //Plane Surface(#)={boundary Line Loop,interior hole Line Loop}
 
//You also have to give a name to the Plane surface
Physical Surface(15)={14};

