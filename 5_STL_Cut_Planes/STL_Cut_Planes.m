function STL_Cut_Planes
% Description: This function calculates the intersection polygons that 
% result from intersections of a 3D triangulation with planes. The
% resulting closed polygons are determined and geometric quantities are
% calculated. The results are saved into a .mat-file.

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

close all

%% add all folders to path
pathstr = mfilename('fullpath');
[pathstr,~,~] = fileparts( pathstr );
[pathstr,~,~] = fileparts( pathstr );
addpath(genpath(pathstr));

warning('off','MATLAB:polyshape:repairedBySimplify')

%% parameters
Nz=100; % number of cut planes
plot_true_false=true;
filename='fine_pored_phase_field';

%% read STL file
TR = stlread(sprintf('%s/2_Inputfiles/STLs/%s.stl',pathstr,filename));
Points=TR.Points;
Points=Points*1e3; % from µm to nm
Connectivity=TR.ConnectivityList;

%% create vector with the z-coordinates of all cut planes
z_vec=linspace(-10,max(Points(:,3))+10,Nz);

%% calculate cut planes
results=calculate_properties_cutplanes(Points,Connectivity,z_vec,plot_true_false);

%% plot results
figure()
subplot(1,2,1)
hold on
plot(z_vec,results.A_ges./results.A_convhull)
grid on
xlabel('z [nm]')
ylabel('A_{Regions}/A_{CH}')

subplot(1,2,2)
hold on
plot(z_vec,results.L_ges./results.L_convhull)
grid on
xlabel('z [nm]')
ylabel('U_{Regions}/U_{CH}')

%% save results
folder_name_results=sprintf('%s/5_STL_cut_planes/Results_STL_cut_planes',pathstr);
if not(isfolder(folder_name_results))
    mkdir(folder_name_results)
end
save(sprintf('%s/5_STL_cut_planes/Results_STL_cut_planes/%s_results',...
    pathstr,filename), '-struct', 'results');
end