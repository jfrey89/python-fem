close all
clear all
clc
x = load('./files/x.txt');
y = load('./files/y.txt');
p = load('./files/p.txt');
u = load('./files/u.txt');
v = load('./files/v.txt');
tri = load('./files/tri.txt');
UVP = load('./files/UVP.txt');
tri = tri + 1;

figure()
trisurf(tri, x, y, p)
figure()
trisurf(tri, x, y, u)
figure()
trisurf(tri, x, y, v)
figure()
quiver(x, y, u, v)