function STL_To_Voxel
% Description: This function translates a 3D triangulation into a 3D grid
% with voxel data. Each voxel inside the triangulation is set to 1 and each
% point outside the triangulation is set to 0. For each cut in z-direction
% a png-file and if desired a txt-file with the data is saved.

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

close all

%% add all folders to path
pathstr = mfilename('fullpath');
[pathstr,~,~] = fileparts( pathstr );
[pathstr,~,~] = fileparts( pathstr );
addpath(genpath(pathstr));

%% parameters
filename='coarse_pored_reference';

plot_true_false=true;
save_txt_files_true_false=true;
Nx=200;
Ny=200;
Nz=200;

%% read file
TR = stlread(sprintf('%s/2_Inputfiles/STLs/%s.stl',pathstr,filename));
Points=TR.Points;
Connectivity=TR.ConnectivityList;

%% Scaling from µm to nm
Points=Points*1e3;

%% create the 3D grid
x_min=min(Points(:,1))-1;
x_max=max(Points(:,1))+1;
y_min=min(Points(:,2))-1;
y_max=max(Points(:,2))+1;
z_min=min(Points(:,3))-1;
z_max=max(Points(:,3))+1;

x_vec=linspace(x_min,x_max,Nx);
y_vec=linspace(y_min,y_max,Ny);
z_vec=linspace(z_min,z_max,Nz);

delta_x=x_vec(2)-x_vec(1);
delta_y=y_vec(2)-y_vec(1);
delta_z=z_vec(2)-z_vec(1);

[X,Y]=meshgrid(x_vec,y_vec);
X=X(:);
Y=Y(:);

%% calculate cut planes
results=calculate_properties_cutplanes(Points,Connectivity,z_vec,false);

%% make a new folder for the results
folder_name_results=sprintf('%s/3_STL_to_Voxel/Results_STL_to_Voxel',pathstr);
if not(isfolder(folder_name_results))
    mkdir(folder_name_results)
end
folder_name=sprintf('%s/3_STL_to_Voxel/Results_STL_to_Voxel/%s',pathstr,filename);
mkdir(folder_name)

%% save parameters to text file
file_name_save=sprintf('%s/3_STL_to_Voxel/Results_STL_to_Voxel/%s/%s.txt',pathstr,filename,'parameter');
fileID = fopen(file_name_save,'w');
fprintf(fileID,'%6.24f\r\n',[delta_x delta_y delta_z]');
fclose(fileID);

if(plot_true_false)
    figure(1)
    set(figure(1),'units','normalized','outerposition',[0 0 1 1])   % Figure maximieren auf ganzen Bildschirm
end

%% determine voxel data for each cut plane
for n=1:Nz
    fprintf('Determine voxel data for the %1.0f. cut plane\n',n);
    IN = false(Nx*Ny,1);
    if(~isempty(results.Ordered_Intersection_Points{n}))
        for k=1:length(results.Ordered_Connectivity_single{n})
            x_test=results.Ordered_Intersection_Points{n}(results.Ordered_Connectivity_single{n}{k}(:,1),1);
            y_test=results.Ordered_Intersection_Points{n}(results.Ordered_Connectivity_single{n}{k}(:,1),2);
            
            x_test_min=min(x_test);
            x_test_max=max(x_test);
            y_test_min=min(y_test);
            y_test_max=max(y_test);
            
            X_test=X;
            Y_test=Y;
            
            ind_test=X_test>=x_test_min & X_test<=x_test_max &...
                Y_test>=y_test_min & Y_test<=y_test_max;
            
            in_test =inpoly2([X(ind_test),Y(ind_test)],[x_test,y_test]);
            
            IN_neu = false(Nx*Ny,1);
            IN_neu(ind_test)=in_test;
            IN=xor(IN_neu,IN);
        end
    end
    IN_mat=reshape(IN,Nx,Ny);
    
    % save png of the current cut plane
    IN_mat_image=uint8(IN_mat);
    IN_mat_image(IN_mat_image==0)=255;
    IN_mat_image(IN_mat_image==1)=0;
    
    file_name_save=sprintf('%s/3_STL_to_Voxel/Results_STL_to_Voxel/%s/%s_%1.0f.png',pathstr,filename,filename,n);
    imwrite(IN_mat_image,file_name_save)
    
    if(save_txt_files_true_false)
        % save text-file for the current cut plane
        [i_in,j_in]=find(IN_mat==1);
        [i_out,j_out]=find(IN_mat==0);
        
        file_name_save=sprintf('%s/3_STL_to_Voxel/Results_STL_to_Voxel/%s/%s_%1.0f.txt',pathstr,filename,filename,n);
        fileID = fopen(file_name_save,'w');
        if(~isempty(i_in))
            fprintf(fileID,'%1.0f %1.0f %1.0f\r\n',[i_in ,j_in ,ones(length(i_in),1)]');
        end
        if(~isempty(i_out))
            fprintf(fileID,'%1.0f %1.0f %1.0f\r\n',[i_out,j_out,zeros(length(i_out),1)]');
        end
        fclose(fileID);
    end
    
    % if wanted plot the results for the current cut plane
    if(plot_true_false)
        clf
        hold on
        Z=ones(Nx*Ny,1)*z_vec(n);
        plot3(X(IN),Y(IN),Z(IN),'.b')
        
        if(~isempty(results.Ordered_Intersection_Points{n}))
            plot3(results.Ordered_Intersection_Points{n}(:,1),...
                results.Ordered_Intersection_Points{n}(:,2),...
                results.Ordered_Intersection_Points{n}(:,3),'.k')
        end
        axis equal
        set(gca,'XLim',[x_vec(1) x_vec(end)],...
            'YLim',[y_vec(1) y_vec(end)],'ZLim',[z_vec(1) z_vec(end)])
        grid on
        xlabel('x [nm]')
        ylabel('y [nm]')
        zlabel('z [nm]')
        view(15,35)
        drawnow
    end
end
end