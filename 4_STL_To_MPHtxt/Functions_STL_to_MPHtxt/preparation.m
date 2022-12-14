function [TR,Points_Scaled]=preparation(filename_input,scaling_factor,...
    ax,ay,az,centralize_true_false)
% Description: this function loads and modifies a specified STL file

% Input:  
% 1. filename_input (string): file name of the STL file
% 2. scaling_factor (1x1): scaling factor, the coordinates are multiplied 
% by it
% 3. ax (1x1): rotation angle around x-axis
% 4. ay (1x1): rotation angle around y-axis
% 5. az (1x1): rotation angle around z-axis
% 6. centralize_true_false (1x1): boolean, if true the triangulation is
% centralized at its centroid
% 
% Output:  
% 1. fileID (1x1): updated file

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

%% read STL file
TR = stlread(sprintf('%s.stl',filename_input));
Points_Scaled=TR.Points*scaling_factor;

%% if wanted centralize the triangulation at the centroid coordinates
if(centralize_true_false)
    Centroid=calculate_centroid(TR.ConnectivityList,Points_Scaled);  
    Points_Scaled(:,1)=Points_Scaled(:,1)-Centroid(1);
    Points_Scaled(:,2)=Points_Scaled(:,2)-Centroid(2);
    Points_Scaled(:,3)=Points_Scaled(:,3)-Centroid(3);
end

%% rotate the geometry around the given angles
R=rot_matrices(ax,ay,az);
for n=1:length(Points_Scaled(:,1))
    Points_Scaled(n,:)=R*Points_Scaled(n,:)';
end

%% calculate some geometric quantites of the triangulation
[Unit_Normals,Midpoints,Surface_Area,Volume,Euler_Characteristic]=...
    calculate_mesh_quantities(TR.ConnectivityList,Points_Scaled);

% print some quantities to the command line to check if everything works as
% expected
fprintf('\nSurface area: %1.8g\n',Surface_Area)
fprintf('Volume: %1.8g\n',Volume)
fprintf('Euler-Characteristic: %1.0f\n',Euler_Characteristic)

%% plot the triangulation and the face normals
figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
hold on
trisurf(TR.ConnectivityList,Points_Scaled(:,1),Points_Scaled(:,2),...
    Points_Scaled(:,3),'Facecolor',[1 0.9 0])
quiver3(Midpoints(:,1),Midpoints(:,2),Midpoints(:,3),...
    Unit_Normals(:,1),Unit_Normals(:,2),Unit_Normals(:,3),'-b')
xlabel('x [nm]')
ylabel('y [nm]')
zlabel('z [nm]')
view(15,35)
axis equal
grid on
end

function R=rot_matrices(ax,ay,az)
% Description: this function returns rotation matrix for a rotation
% around the x,y- and z-axis with specified angles 

% Input: 
% 1. ax (1x1): rotation angles around x-axis
% 2. ay (1x1): rotation angles around y-axis
% 3. az (1x1): rotation angles around z-axis
% Output:
% 1. R (3x3): combined rotation matrix 

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

Rx=[1 0 0; 0 cosd(ax) -sind(ax); 0 sind(ax) cosd(ax)];
Ry=[cosd(ay) 0 sind(ay); 0 1 0; -sind(ay) 0 cosd(ay)];
Rz=[cosd(az) -sind(az) 0; sind(az) cosd(az) 0; 0 0 1];

R=Rz*Ry*Rx;
end