// define some variables to describe your domain:
x_left = 0.0;
x_right = 1.0;
y_bottom = 0.0;
y_top = 1.0;

z_coord = 0.0;
 
// define a mesh characteristic length (cl) - smaller=more refined mesh
cl = 0.025;
 
// define points to describe your 2D domain
Point(1)={x_left,y_bottom,z_coord,cl};
Point(2)={(0.375)*(x_left + x_right),y_bottom,z_coord,cl};
Point(3)={(0.375)*(x_left + x_right),(0.375)*(y_bottom + y_top),z_coord,cl};
Point(4)={(0.625)*(x_left + x_right),(0.375)*(y_bottom + y_top),z_coord,cl};
Point(5)={(0.625)*(x_left + x_right),y_bottom,z_coord,cl};
Point(6)={x_right,y_bottom,z_coord,cl};
Point(7)={x_right,y_top,z_coord,cl};
Point(8)={(0.625)*(x_left + x_right),y_top,z_coord,cl};
Point(9)={(0.625)*(x_left + x_right),(0.625)*(y_top + y_bottom),z_coord,cl};
Point(10)={(0.375)*(x_left + x_right),(0.625)*(y_top + y_bottom),z_coord,cl};
Point(11)={(0.375)*(x_left + x_right),y_top,z_coord,cl};
Point(12)={x_left,y_top,z_coord,cl};
 
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
Line Loop(13)={1,2,3,4,5,6,7,8,9,10,11,12};
Plane Surface(14)={13};

//You also have to give a name to the Plane surface
Physical Surface(19)={10};
