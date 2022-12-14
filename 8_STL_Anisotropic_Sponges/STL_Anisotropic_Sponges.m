function STL_Anisotropic_Sponges
% Description: This function performs all the necessary steps to create an
% anisotropic sponge based on a reference sponge. In the first step, the
% alpha shape of the reference sponge is calculated. Then the density is 
% fitted as a function of the threshold value. Finally, a threshold matrix
% Phi is created and exported, which can be used directly for the initialization of the phase field
% function. 

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

close all

%% add all folders to path
pathstr = mfilename( 'fullpath' );
[pathstr,~,~] = fileparts( pathstr );
[pathstr,~,~] = fileparts( pathstr );
addpath(genpath(pathstr));

%% parameters
Ngrid=100;
filter_size=25;
filename='fine_pored_reference';
filename_output_fit_phi='fine_pored';
foldername_densities='Anisotropic_Densities';

%% create results folder
folder_name_results=sprintf('%s/8_STL_Anisotropic_Sponges/Results_Anisotropic_Sponges',pathstr);
if not(isfolder(folder_name_results))
    mkdir(folder_name_results)
end

%% determine and store alpha shape of the specified sponge
Tri_sponge = stlread(sprintf('%s/2_Inputfiles/STLs/%s.stl',pathstr,filename));
Tri_alpha_shape=alpha_shapes_STL(Tri_sponge.Points,filter_size);
stlwrite(Tri_alpha_shape,sprintf('%s/8_STL_Anisotropic_Sponges/Results_Anisotropic_Sponges/%s_alpha_shape.stl'...
    ,pathstr,filename))

%% fit threshold function
directory=sprintf('%s/2_Inputfiles/%s',pathstr,foldername_densities);
results_fit=fit_threshold_function(directory);
save(sprintf('%s/8_STL_Anisotropic_Sponges/Results_Anisotropic_Sponges/results_fit_%s.mat',...
    pathstr,filename_output_fit_phi), '-struct', 'results_fit');

%% calculate threshold matrix 
[X_mat,Y_mat,Z_mat,Phi_mat]=calculate_threshold_phi(Tri_sponge,...
                                        Tri_alpha_shape,Ngrid,results_fit);

%% store the results
writematrix([X_mat(:) Y_mat(:) Z_mat(:) Phi_mat(:)],...
    sprintf('%s/8_STL_Anisotropic_Sponges/Results_Anisotropic_Sponges/results_phi_%s.txt',...
    pathstr,filename_output_fit_phi));
end