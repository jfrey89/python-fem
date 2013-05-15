close all
clc
x = load('./files/x.txt');
y = load('./files/y.txt');
p_s = load('./files/p.txt');
u_s = load('./files/u.txt');
v_s = load('./files/v.txt');
tri = load('./files/tri.txt');

% figure()
% femMesh = trisurf(tri,x,y,p, 'EdgeColor', 'flat', ...
%   'FaceLighting','phong');
% view(-50,30)
% axis tight
% camlight left
% femTitle = title({'pressure p_{h}(x,y)'});
% femXLabel = xlabel('x-axis');
% femYLabel = ylabel('y-axis');
% femZLabel = zlabel('p_{h}(x,y)');
% set( gca                       , ...
%   'FontName'   , 'Helvetica' );
% set([femXLabel, femYLabel, femZLabel], ...
%   'FontName'   , 'AvanteGarde');
% set([femXLabel, femYLabel, femZLabel]  , ...
%   'FontSize'   , 10          );
% set( femTitle                    , ...
%   'FontSize'   , 12          , ...
%   'FontWeight' , 'bold'      );
% 
% set(gca, ...
%   'Box'         , 'off'     , ...
%   'TickDir'     , 'in'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'XMinorTick'  , 'on'      , ...
%   'YMinorTick'  , 'on'      , ...
%   'ZMinorTick', 'on', ...
%   'XGrid', 'on', ...,
%   'YGrid', 'on', ...
%   'ZGrid', 'on', ...
%   'LineWidth'   , 1         );
% 
% figure()
% femMesh = trisurf(tri,x,y,u, 'EdgeColor', 'flat', ...
%   'FaceLighting','phong');
% view(-50,30)
% axis tight
% camlight left
% femTitle = title({'x velocity: u_{h}(x,y)'});
% femXLabel = xlabel('x-axis');
% femYLabel = ylabel('y-axis');
% femZLabel = zlabel('u_{h}(x,y)');
% set( gca                       , ...
%   'FontName'   , 'Helvetica' );
% set([femXLabel, femYLabel, femZLabel], ...
%   'FontName'   , 'AvanteGarde');
% set([femXLabel, femYLabel, femZLabel]  , ...
%   'FontSize'   , 10          );
% set( femTitle                    , ...
%   'FontSize'   , 12          , ...
%   'FontWeight' , 'bold'      );
% 
% set(gca, ...
%   'Box'         , 'off'     , ...
%   'TickDir'     , 'in'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'XMinorTick'  , 'on'      , ...
%   'YMinorTick'  , 'on'      , ...
%   'ZMinorTick', 'on', ...
%   'XGrid', 'on', ...,
%   'YGrid', 'on', ...
%   'ZGrid', 'on', ...
%   'LineWidth'   , 1         );
% 
% figure()
% femMesh = trisurf(tri,x,y,v, 'EdgeColor', 'flat', ...
%   'FaceLighting','phong');
% view(-50,30)
% axis tight
% camlight left
% femTitle = title({'y velocity: v_{h}(x,y)'});
% femXLabel = xlabel('x-axis');
% femYLabel = ylabel('y-axis');
% femZLabel = zlabel('v_{h}(x,y)');
% set( gca                       , ...
%   'FontName'   , 'Helvetica' );
% set([femXLabel, femYLabel, femZLabel], ...
%   'FontName'   , 'AvanteGarde');
% set([femXLabel, femYLabel, femZLabel]  , ...
%   'FontSize'   , 10          );
% set( femTitle                    , ...
%   'FontSize'   , 12          , ...
%   'FontWeight' , 'bold'      );
% 
% set(gca, ...
%   'Box'         , 'off'     , ...
%   'TickDir'     , 'in'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'XMinorTick'  , 'on'      , ...
%   'YMinorTick'  , 'on'      , ...
%   'ZMinorTick', 'on', ...
%   'XGrid', 'on', ...,
%   'YGrid', 'on', ...
%   'ZGrid', 'on', ...
%   'LineWidth'   , 1         );
% 
% figure()
% quiver(x, y, u, v)
% title('velocity vector field')