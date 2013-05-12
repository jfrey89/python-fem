z_coord = 0.0;
 
// define a mesh characteristic length (cl) - smaller=more refined mesh
cl = 0.2;
 
// define points to describe your 2D domain
Point(1)={-3.0,1.0,z_coord,cl};
Point(2)={3.0,1.0,z_coord,cl};
Point(3)={3.0,0.5,z_coord,cl};
Point(4)={1.0,0.5,z_coord,cl};
Point(5)={1.0,-0.5,z_coord,cl};
Point(6)={3.0,-0.5,z_coord,cl};
Point(7)={3.0,-1.0,z_coord,cl};
Point(8)={-3.0,-1.0,z_coord,cl};
Point(9)={-3.0,-0.5,z_coord,cl};
Point(10)={-1.0,-0.5,z_coord,cl};
Point(11)={-1.0,0.5,z_coord,cl};
Point(12)={-3.0,0.5,z_coord,cl};
 
// make these points into lines
Line(1)={1,2}; //Line(#)={Start Point,End Point}
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,5};
Line(5)={5,6};
Line(6)={6,7};
Line(7)={7,8};
Line(8)={8,9};
Line(9)={9,10};
Line(10)={10,11};
Line(11)={11,12};
Line(12)={12,1};
 
// now, define a circle in the center of the rectangle
/*x_c = (0.5)*((0.5)*((0.5)*(x_left + x_right) + x_right) + x_right);*/
/*y_c = (0.5)*(y_bottom + y_top);*/
/*r = (0.1)*(y_top-y_bottom);*/
 
//To define a circle, place a point at the center + at points around the radius
/*Point(5)={x_c,y_c,z_coord,cl}; //center point*/
/*Point(6)={x_c+r,y_c,z_coord,cl}; //3 o'clock position*/
/*Point(7)={x_c,y_c+r,z_coord,cl}; //12 o'clock position*/
/*Point(8)={x_c-r,y_c,z_coord,cl}; //etc...*/
/*Point(9)={x_c,y_c-r,z_coord,cl};*/
 
//connect these points with circle segments
/*Circle(5)={6,5,7}; //Circle(#)={arc start, center point, arc end}.*/
/*Circle(6)={7,5,8};*/
/*Circle(7)={8,5,9};*/
/*Circle(8)={9,5,6};*/
 
// Now connect the lines of the outer rectangle into a line loop.
Line Loop(13)={1,2,3,4,5,6,7,8,9,10,11,12}; 
 
//...do the same with the circular hole...
/*Line Loop(10)={5,6,7,8};*/
 
//Create a plane surface that can be meshed:
/*Plane Surface(11)={9,10}; //Plane Surface(#)={boundary Line Loop,interior hole Line Loop}*/
Plane Surface(14)={13}; //Plane Surface(#)={boundary Line Loop,interior hole Line Loop}
 
//Name the surfaces to allow assignment of boundary conditions
/*Physical Line(15)={1};//bottom boundary*/
/*Physical Line(13)={2};//right boundary*/
/*Physical Line(14)={3};//top boundary*/
/*Physical Line(15)={4};//left boundary*/
 
//Name the segments of the circle
/*Physical Line(16)={5};*/
/*Physical Line(17)={6};*/
/*Physical Line(18)={7};*/
/*Physical Line(19)={8};*/
 
//You also have to give a name to the Plane surface
Physical Surface(15)={14};
