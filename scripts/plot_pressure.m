x = load('./../files/x_p.txt');
y = load('./../files/y_p.txt');
z = load('./../files/z_p.txt');
tri = load('./../files/tri_p.txt');
tri = tri + 1;

figure()
trisurf(tri, x, y, z)