function [Unit_normals,Midpoints,Surface_area,Volume,Euler_characteristic]=...
    calculate_mesh_quantities(Connectivity,Points)
% Description: this function calculates global geometric quantities for a
% given 3D triangulation

% Input:  
% 1. Connectivity (Nx3): connectivity matrix of the mesh
% 2. Points (Mx3): coordinates of the mesh vertices
% Output:
% 1. Unit_Normals (Mx3): normalized face normals for each triangle
% 2. Midpoints (Mx3): coordinates of the midpoints for each triangle
% 3. Surface_Area (1x1): total area of all triangles
% 4. Volume (1x1): enclosed volume of the triangulation
% 5. Euler_Charakteristik (1x1): Euler characteristic, i.e. E=M-0.5*N
% of the triangulation 

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

%% coordinates of all mesh vertices of all triangles
x0=Points(Connectivity(:,1),1);
x1=Points(Connectivity(:,2),1);
x2=Points(Connectivity(:,3),1);

y0=Points(Connectivity(:,1),2);
y1=Points(Connectivity(:,2),2);
y2=Points(Connectivity(:,3),2);

z0=Points(Connectivity(:,1),3);
z1=Points(Connectivity(:,2),3);
z2=Points(Connectivity(:,3),3);

d12= [x0-x2,y0-y2,z0-z2];
d23= [x1-x2,y1-y2,z1-z2];
d13= [x0-x2,y0-y2,z0-z2];

%% lengths of all edges of the triangulation
L12=sqrt(sum(d12.^2,2));
L13=sqrt(sum(d23.^2,2));
L32=sqrt(sum(d13.^2,2));

%% midpoint coordinates of all triangles
Midpoints=(Points(Connectivity(:,3),:).*L12+...
    Points(Connectivity(:,2),:).*L13+...
    Points(Connectivity(:,1),:).*L32)./(L12+L13+L32);

%% calculate geometric quantities
Normals = cross_vectorized(d12,d23);
Areas  = 0.5*sqrt(sum(Normals.^2,2));
Unit_normals =Normals ./(2*Areas);
Surface_area=sum(Areas);

% volume of the triangulation from int(div(F))dV=int(F)dA for
% F=1/3*(x,y,z) --> div(F)=1 and thus V=int(F)dA
Volume = 1/6*sum(Points(Connectivity(:,1),1).*Normals(:,1)+...
                 Points(Connectivity(:,1),2).*Normals(:,2)+...
                 Points(Connectivity(:,1),3).*Normals(:,3));    

% Gauss bonnet
Euler_characteristic=length(Points(:,1))-0.5*length(Connectivity(:,1));
end


