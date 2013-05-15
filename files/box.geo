// define some variables to describe your domain:
x_left = 0.0;
x_right = 5.0;
y_bottom = 0.0;
y_top = 3.0;
 
z_coord = 0.0;
 
// define a mesh characteristic length (cl) - smaller=more refined mesh
cl = 0.15;
 
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

 
// Now connect the lines of the outer rectangle into a line loop.
Line Loop(13)={1,2,3,4}; //Line Loop(#)={Line Segment 1,...,Line Segment n}
 
//Create a plane surface that can be meshed:
Plane Surface(14)={13}; //Plane Surface(#)={boundary Line Loop,interior hole Line Loop}
 
//You also have to give a name to the Plane surface
Physical Surface(15)={14};

