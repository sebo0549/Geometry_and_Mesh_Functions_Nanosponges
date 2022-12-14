function [Sideness,Signed_Distance]=determine_sideness_of_points(Points,Normal_plane,V0_plane)
% Description: this function calculates the signed distances of M Points 
% in regards to a plane specified by a point V0_plane and a normal vector
% Normal_plane

% Input: 
% 1. Points (Mx3): coordinates of the points
% 2. Normal_plane (1x3): normal vector of the plane
% 3. V0_plane (1x3): point on the plane
% Output: 
% 1. Sideness (Mx1): +1 if distance>0, -1 if distance <0 and 0 if distance=0
% 2. Signed_Distance (Mx1): signed distance of the points to the plane

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

if(norm(Normal_plane)~=1)
    Normal_plane=Normal_plane./norm(Normal_plane);
end

Signed_Distance=dot(repmat(Normal_plane,length(Points(:,1)),1),Points-V0_plane,2);
Sideness=sign(Signed_Distance);
end