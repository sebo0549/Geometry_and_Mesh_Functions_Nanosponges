function centroid=calculate_centroid(Connectivity,Points)
% Description: this function calculates the centroid coordinates of a
% closed 3D triangulation

% Input:  
% 1. Connectivity (Nx3): connectivity matrix of the mesh
% 2. Points (Mx3): coordinates of the mesh vertices
% Output:
% 1. centroid (1x3): centroid coordinates 

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

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

Normals = cross_vectorized(d12,d23);

% volume of the triangulation from int(div(F))dV=int(F)dA for
% F=1/3*(x,y,z) --> div(F)=1 and thus V=int(F)dA
Volume = 1/6*sum(Points(Connectivity(:,1),1).*Normals(:,1)+Points(Connectivity(:,1),2).*Normals(:,2)+Points(Connectivity(:,1),3).*Normals(:,3));    

a0 = Points(Connectivity(:,1),:);
a1 = Points(Connectivity(:,2),:);
a2 = Points(Connectivity(:,3),:);

% centroid is calculated using gauss theorem:
% S_i= 1/V*int(x*e_i)dV = 1/(2V)*int(0.5*div(x*e_i)^2*e_i)dV
%    = 1/(2*V)*int((x*e_i)^2*(n*e_i))dA 
%    = 1/(2*V)*sum(int((x*e_i)^2*(n_j*e_i))dA) 
% the integral is evaluated using: int(f(x))dA=A/3*(f(x_1)+f(x_2)+f(x_3))
% (exact only for second degree polynomials but this is good enough for the
% centroid) where x_i are the midpoints of each triangle
Sx = 1/(48*Volume)*(sum(Normals(:,1).*((a0(:,1)+a1(:,1)).^2+(a1(:,1)+a2(:,1)).^2 + (a2(:,1)+a0(:,1)).^2)));
Sy = 1/(48*Volume)*(sum(Normals(:,2).*((a0(:,2)+a1(:,2)).^2+(a1(:,2)+a2(:,2)).^2 + (a2(:,2)+a0(:,2)).^2)));
Sz = 1/(48*Volume)*(sum(Normals(:,3).*((a0(:,3)+a1(:,3)).^2+(a1(:,3)+a2(:,3)).^2 + (a2(:,3)+a0(:,3)).^2)));

centroid=[Sx,Sy,Sz];
end


