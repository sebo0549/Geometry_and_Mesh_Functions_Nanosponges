function STL_To_MPHtxt
% Description: This function creates a constraint mesh from an STL in an
% mphtxt file, which can be read into COMSOL and used for error-free 
% geometry creation. For this purpose, a physical surface is created from
% each triangle.

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

close all

%% add all folders to path
pathstr = mfilename('fullpath');
[pathstr,~,~] = fileparts( pathstr );
[pathstr,~,~] = fileparts( pathstr );
addpath(genpath(pathstr));

%% parameters
filename='fine_pored_reference_curvature';
scaling_factor=1e3; % from µm to nm
rot_angle_around_xAxis=0;
rot_angle_around_yAxis=0;
rot_angle_around_zAxis=0;
centralize_true_false=true;

%% preparation
[TR,Points_Scaled]=preparation(filename,scaling_factor,...
    rot_angle_around_xAxis,rot_angle_around_yAxis,rot_angle_around_zAxis,...
    centralize_true_false);

%% create and store the .mphtxt mesh file
folder_name_results=sprintf('%s/4_STL_to_MPHtxt/Results_STL_to_MPHtxt',pathstr);
if not(isfolder(folder_name_results))
    mkdir(folder_name_results)
end
fileID = fopen(sprintf('%s/4_STL_to_MPHtxt/Results_STL_to_MPHtxt/%s.mphtxt'...
                       ,pathstr,filename),'w');
fprintf(fileID,'# Created by COMSOL Multiphysics.\n\n');
fprintf(fileID,'# Major & minor version\n');
fprintf(fileID,'0 1\n');
fprintf(fileID,'1 # number of tags\n');
fprintf(fileID,'# Tags\n');
fprintf(fileID,'5 mesh1\n');
fprintf(fileID,'1 # number of types\n');
fprintf(fileID,'# Types\n');
fprintf(fileID,'3 obj\n\n');
fprintf(fileID,'# --------- Object 0 ----------\n\n');
fprintf(fileID,'0 0 1\n');
fprintf(fileID,'4 Mesh # class\n');
fprintf(fileID,'4 # version\n');
fprintf(fileID,'3 # sdim\n');
fprintf(fileID,sprintf('%1.0f # number of mesh vertices\n',length(TR.Points(:,1))));
fprintf(fileID,'0 # lowest mesh vertex index\n\n');
fprintf(fileID,'# Mesh vertex coordinates\n');
fprintf(fileID,'%1.32f %1.32f %1.32f\n',Points_Scaled');
fileID=create_mesh_edge_constrained(fileID,TR.ConnectivityList);
fclose(fileID);
end