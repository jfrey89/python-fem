close all
x_p = load('./../files/x_p.txt');
y_p = load('./../files/y_p.txt');
p = load('./../files/p.txt');
tri = load('./../files/tri.txt');
tri = tri + 1;

figure()
trisurf(tri, x_p, y_p, p)

x = load('./../files/x.txt');
y = load('./../files/y.txt');
u = load('./../files/u.txt');
v = load('./../files/v.txt');
figure()
quiver(x, y, u, v)