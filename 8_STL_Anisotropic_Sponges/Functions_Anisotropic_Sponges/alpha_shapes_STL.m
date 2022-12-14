function TR_alpha=alpha_shapes_STL(Points,filter_size)
% Description: this function calculates the alpha shape of a 3D point cloud
% and returns its triangulation

% Input:  
% 1. Points (Mx3): coordinates of the mesh vertices
% 2. filter_size (1x1): alpha radius, defines the level of detail of the
% alpha shape
% Output: TR_alpha (struct): triangulation of the alpha shape

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

%% scale from µm to nm
Points=Points*1e3;

%% calculate alpha shape
shp = alphaShape(Points,filter_size);

%% extract outer boundary and create a triangulation
[con_alpha,Points_alpha] = boundaryFacets(shp);
TR_alpha=triangulation(con_alpha,Points_alpha);

%% plot alpha shape
figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
trisurf(con_alpha,Points_alpha(:,1),Points_alpha(:,2),Points_alpha(:,3))
title('\alpha-shape')
xlabel('x [nm]')
ylabel('y [nm]')
zlabel('z [nm]')
view(15,35)
axis equal
grid on
drawnow
end