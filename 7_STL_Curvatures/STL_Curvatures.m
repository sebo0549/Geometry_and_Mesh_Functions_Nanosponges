function STL_Curvatures
% Description: This function calculates, plots and stores the curvature
% properties of multiple closed 3D triangulations. 

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

close all

%% add all folders to path
pathstr = mfilename('fullpath');
[pathstr,~,~] = fileparts( pathstr );
[pathstr,~,~] = fileparts( pathstr );
addpath(genpath(pathstr));

%% define number of bins and axis limits
kappa1_min=-250; % limits should be a multiple of ten
kappa1_max= 100;
kappa2_min=-100;
kappa2_max= 250;

%% define all  input files
filenames = ["coarse_pored_reference"...
             "fine_pored_reference_curvature"];

%% initialize arrays
Nfiles=length(filenames);
res=cell(Nfiles,1);

%% read mesh files
for i = 1:Nfiles    
    res{i}.TR = stlread(sprintf('%s/2_Inputfiles/STLs/%s.stl',pathstr,filenames(i)));
    res{i}.Points = res{i}.TR.Points;
    res{i}.con = res{i}.TR.ConnectivityList;
end

%% perform curvature calculation
for i = 1:Nfiles
    [res{i}.K_Mean,res{i}.K_Gauss,res{i}.kappa1,res{i}.kappa2,res{i}.K_Vektor,...
    res{i}.Total_mean_curvature,res{i}.Total_gauss_curvature,...
    res{i}.inds_undecidable,res{i}.inds_equal,res{i}.A_mixed,res{i}.Thetas_ges,res{i}.VN]=...
    calculate_curvatures(res{i}.con,res{i}.Points);
end

%% plot results
% plot mean curvature over the surface of the triangulation
figure(1);
set(gcf,'units','normalized','outerposition',[0 0 0.5 1]) 

% deetermine a nice dsitribution for the subplots
subplots_rows_cols=num_subplots(Nfiles);
for i = 1:Nfiles
    subplot(subplots_rows_cols(1),subplots_rows_cols(2),i);
    hold on
    trisurf(res{i}.con,res{i}.Points(:,1)*1e3,res{i}.Points(:,2)*1e3,res{i}.Points(:,3)*1e3,res{i}.K_Mean);
    set(gca,'CLim',[-150,150]);
    shading interp
    axis equal
    view(10,15)
    grid on
    xlabel('x [nm]')
    ylabel('y [nm]')
    zlabel('z [nm]')
end

% plot 2D histograms of the principal curvatures
figure(2);
set(gcf,'units','normalized','outerposition',[0.5 0 0.5 1])  
for i = 1:Nfiles
    N_bins_x=ceil((kappa1_max-kappa1_min)/10)+1;
    N_bins_y=ceil((kappa2_max-kappa2_min)/10)+1;
    kappa_lim = [res{i}.kappa1,res{i}.kappa2];

    % filter region
    kappa_lim = kappa_lim(kappa1_min<=kappa_lim(:,1) & kappa_lim(:,1)<=kappa1_max &...
                          kappa2_min<=kappa_lim(:,2) & kappa_lim(:,2)<=kappa2_max,:);

    [Counts,Centers] = hist3(kappa_lim,'Nbins',[N_bins_x,N_bins_y]);
    Counts=Counts';

    % store counts, centers and N_bins_x,y in res-struct
    res{i}.Counts=Counts;
    res{i}.Centers=Centers;
    res{i}.N_bins_x=N_bins_x;
    res{i}.N_bins_y=N_bins_y;
    
    subplot(subplots_rows_cols(1),subplots_rows_cols(2),i);
    hold on;   
    surf(Centers{1},Centers{2},zeros(size(Counts)),Counts);
    colormap jet;
    shading flat;
    plot3([-1000,1000],[1000,-1000],[0 0],'r','LineWidth',1.5); 
    axis equal
    xlabel('\kappa_1 [µm^{-1}]')
    ylabel('\kappa_2 [µm^{-1}]')
    view(2)
    axis([kappa1_min kappa1_max kappa2_min kappa2_max])
end

%% store results
folder_name_results=sprintf('%s/7_STL_Curvatures/Results_STL_Curvatures',pathstr);
if not(isfolder(folder_name_results))
    mkdir(folder_name_results)
end
for i = 1:Nfiles
    res_save=res{i};
    
    save(sprintf('%s/7_STL_Curvatures/Results_STL_Curvatures/%s_curvatures',...
        pathstr,filenames(i)), '-struct', 'res_save');
end
end