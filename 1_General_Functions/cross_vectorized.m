function N=cross_vectorized(a,b)
% Description: vectorized version of cross product

% Input:  
% 1. a (Nx3): list of points
% 2. b (Nx3): list of points
% Output:
% 1. N (Nx3): result of the cross product, i.e. N_i=a_i x b_i 

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

N = [a(:,2).*b(:,3) - a(:,3).*b(:,2)...
     a(:,3).*b(:,1) - a(:,1).*b(:,3)...
     a(:,1).*b(:,2) - a(:,2).*b(:,1)]; 
end
