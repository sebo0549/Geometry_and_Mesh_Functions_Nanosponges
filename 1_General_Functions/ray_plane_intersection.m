function P_intersection=ray_plane_intersection(Normal_vec,V0_plane,P0,P1)
% Description: this function calculates the intersection points of M rays 
% specified by two points P0 and P1 with a plane specified by a point 
% V0_plane and a normal vector Normal_Plane

% Input:  
% 1. Normal_vec (1x3): normal vector of the plane
% 2. V0_plane (1x3): point on the plane
% 3. P0 (Mx3): start point of the ray
% 4. P1 (Mx3): end point of the ray
% Output: 
% 1. P_intersection (Mx3): coordinates of the intersection point

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

if(norm(Normal_vec)~=1)
    Normal_vec=Normal_vec./norm(Normal_vec);
end

a1=V0_plane-P0;
b1=P1-P0;

A=dot(repmat(Normal_vec,length(P0(:,1)),1),a1,2);
B=dot(repmat(Normal_vec,length(P0(:,1)),1),b1,2);

t=A./B;

P_intersection=P0+t.*b1;
end