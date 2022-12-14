function STL_Geometric_Quantities
% Description: This function calculates global geometric properties 
% (surface area convex hull, volume convex hull, surface area, volume,
% volume fraction, surface volume fraction, surface fraction,
% total Gauss curvature (Euler), total Gauss curvature (integration),
% total mean curvature (integration), Euler characteristic, pore number)
% for multiple closed 3D triangulations and stores the results in a
% .xls-file

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

close all

%% add all folders to path
pathstr = mfilename('fullpath');
[pathstr,~,~] = fileparts( pathstr );
[pathstr,~,~] = fileparts( pathstr );
addpath(genpath(pathstr));

%% parameters
filename_output='Geometric_Quantities_all_Sponges';
files=dir(sprintf('%s/2_Inputfiles/STLs/*.stl',pathstr));

%% calculate geometric quantities for each sponge
N=length(files);

% initialize all arrays and structs
TR=cell(N,1); k=cell(N,1);
Total_mean_curvature=zeros(N,1); Total_gauss_curvature=zeros(N,1);
Surface_Area_Convex_Hull=zeros(N,1); Volume_Convex_Hull=zeros(N,1); 
Surface_Area_Sponge=zeros(N,1); Volume_Sponge=zeros(N,1);
Euler_Characteristic=zeros(N,1); Volume_fraction=zeros(N,1);
Surface_Volume_fraction=zeros(N,1); Surface_fraction=zeros(N,1);
Number_Pores=zeros(N,1); Total_gauss_curvature_Euler=zeros(N,1);

for n=1:N
    fprintf('%1.0f. sponge: %s\n',n,files(n).name)
    TR{n} = stlread(files(n).name);
    k{n} = convhull(TR{n}.Points);
    
    [~,~,~,~,~,Total_mean_curvature(n),Total_gauss_curvature(n)]=...
        calculate_curvatures(TR{n}.ConnectivityList,TR{n}.Points);
    
    [~,~,Surface_Area_Convex_Hull(n),Volume_Convex_Hull(n)]=...
        calculate_mesh_quantities(k{n},TR{n}.Points);
    
    [~,~,Surface_Area_Sponge(n),Volume_Sponge(n),Euler_Characteristic(n)]=...
        calculate_mesh_quantities(TR{n}.ConnectivityList,TR{n}.Points);
    
    if(Volume_Sponge(n)<0)
        fprintf('%1.0f. Sponge Warning: Normal vectors point inside!!!\n',n);
        Volume_Sponge(n)=-Volume_Sponge(n);
    end
    
    if(mod(Euler_Characteristic(n),1)~=0)
        fprintf('%1.0f. Sponge Error: Euler Characteristic is not a whole number!!! There is more than one connected domain.\n',n);
        Euler_Characteristic(n)=NaN;
    end

    if(mod(Euler_Characteristic(n),2)~=0)
        fprintf('%1.0f. Sponge Error: Euler Characteristic is not an even number!!! Check triangulation.\n',n);
    end
       
    Volume_fraction(n)=Volume_Sponge(n)/Volume_Convex_Hull(n);
    Surface_Volume_fraction(n)=Surface_Area_Sponge(n)/Volume_Convex_Hull(n);
    Surface_fraction(n)=Surface_Area_Sponge(n)/Surface_Area_Convex_Hull(n);
    Number_Pores(n)=1-Euler_Characteristic(n)/2;
    Total_gauss_curvature_Euler(n)=2*pi*Euler_Characteristic(n);
end

%% save results to a .xls-file
if(N>0)
    col1=cell(N,1);
    for n=1:N
        col1{n}=files(n).name;
    end
    
    Table = table(col1,Surface_Area_Convex_Hull,Volume_Convex_Hull...
                       ,Surface_Area_Sponge,Volume_Sponge...
                       ,Volume_fraction,Surface_Volume_fraction...
                       ,Surface_fraction,Total_gauss_curvature_Euler...
                       ,Total_gauss_curvature,Total_mean_curvature...
                       ,Euler_Characteristic,Number_Pores,...
        'VariableNames',{'Dateiname STL-File'...
        'Surface_Area_Convex_Hull' 'Volume_Convex_Hull'...
        'Surface_Area_Sponge' 'Volume_Sponge'...
        'Volume_fraction' 'Surface_Volume_fraction'...
        'Surface_fraction' 'Total Gauss Curvature (Euler)' 'Total Gauss Curvature (Integration)' 'Total mean Curvature (Integration)' 'Euler-characteristics' ...
        'Number of holes'});
    
    folder_name_results=sprintf('%s/6_STL_Geometric_Quantities/Results_STL_Geometric_Quantities',pathstr);
    if not(isfolder(folder_name_results))
        mkdir(folder_name_results)
    end
    writetable(Table,sprintf('%s/6_STL_Geometric_Quantities/Results_STL_Geometric_Quantities/%s_%s.xls',...
        pathstr,datetime('now','Format','yyyyMMdd_hhmmss'),filename_output),'Sheet',1)
end
end