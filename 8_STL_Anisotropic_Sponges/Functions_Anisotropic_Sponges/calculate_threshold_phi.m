function [X_mat,Y_mat,Z_mat,Phi_mat]=calculate_threshold_phi(TR_ref_sponge,...
                                        TR_alpha_shape,Ngrid,results_fit)
% Description: this function generates the 3D grid of the threshold
% function Phi_mat, to recreate the alpha shape of TR_alpha_shape and 
% anisotropies present in TR_ref_sponge, based on the fit data in
% results_fit

% Input:  
% 1. TR_ref_sponge (struct): triangulation object of the reference sponge
% 2. TR_alpha_shape (struct): triangulation object of the alpha shape
% 3. Ngrid (1x1): number of grid points for the x-, y- and z-direction
% 4. results_fit (struct): struct with results of the fit of the density as
% a function of threshold value
% Output:
% 1. X_mat (Ngrid x Ngrid x Ngrid): 3D grid with all x-coordiantes
% 2. Y_mat (Ngrid x Ngrid x Ngrid): 3D grid with all y-coordiantes
% 3. Z_mat (Ngrid x Ngrid x Ngrid): 3D grid with all z-coordiantes
% 4. Phi_mat (Ngrid x Ngrid x Ngrid): 3D grid with all threshold values

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

%% extract points and connectivities and scale points
Points_ref_sponge=TR_ref_sponge.Points*1e3;
% Points_ref_sponge=Points_ref_sponge*1e3;
Points_ref_sponge(:,3)=Points_ref_sponge(:,3)-min(Points_ref_sponge(:,3));
Connectivity_ref_sponge=TR_ref_sponge.ConnectivityList;

% alpha shaped points are already scaled
Points_alpha_shape=TR_alpha_shape.Points;
Points_alpha_shape(:,3)=Points_alpha_shape(:,3)-min(Points_alpha_shape(:,3));
Connectivity_alpha_shape=TR_alpha_shape.ConnectivityList;

%% define z-coordinates of all cut planes
z_vec=linspace(-10,max([max(Points_ref_sponge(:,3)) max(Points_alpha_shape(:,3))])+10,Ngrid);

%% calculate properties of all cut planes for both triangulations
cut_plane_props_ref_sponge=calculate_properties_cutplanes(Points_ref_sponge,...
    Connectivity_ref_sponge,z_vec,false);
cut_plane_props_alpha_shape=calculate_properties_cutplanes(Points_alpha_shape,...
    Connectivity_alpha_shape,z_vec,false);

%% calculate z dependent threshold value
A_ratio=cut_plane_props_ref_sponge.A_ges./cut_plane_props_alpha_shape.A_ges;
phi_threshold_z=results_fit.zfunc_model(results_fit.p_sol,A_ratio);

%% define grid
x_vec=linspace(min([cut_plane_props_ref_sponge.x_min cut_plane_props_alpha_shape.x_min]),...
    max([cut_plane_props_ref_sponge.x_max cut_plane_props_alpha_shape.x_max]),Ngrid);
y_vec=linspace(min([cut_plane_props_ref_sponge.y_min cut_plane_props_alpha_shape.y_min]),...
    max([cut_plane_props_ref_sponge.y_max cut_plane_props_alpha_shape.y_max]),Ngrid);

[X,Y]=meshgrid(x_vec,y_vec);
X_vec=X(:);
Y_vec=Y(:);

[X_mat,Y_mat,Z_mat]=meshgrid(x_vec,y_vec,z_vec);

%% calculate threshold matrix 
Phi_mat=zeros(size(Z_mat));
for n=1:Ngrid
    % determine which point of the current cut plane is inside and which 
    % point is outside of the alpha shape
    if(~isempty(cut_plane_props_alpha_shape.Poly{n}))
        inside=false(Ngrid^2,1);
        for k=1:cut_plane_props_alpha_shape.Poly{n}.numboundaries
            [xbound,ybound]=boundary(cut_plane_props_alpha_shape.Poly{n},k);
            inside=xor(inside,inpoly2([X_vec,Y_vec],[xbound,ybound]));
        end
    else
        inside=false(Ngrid^2,1);
    end   

    inside_mat=reshape(inside,Ngrid,Ngrid);
    Phi_cut=zeros(Ngrid);
    Phi_cut(inside_mat)=phi_threshold_z(n); % set threshold for all points inside
    Phi_cut(~inside_mat)=-10; % set large negative value for all points outside
    Phi_mat(:,:,n)=Phi_cut; 
end

%% plot results
figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
hold on
plot(z_vec,cut_plane_props_ref_sponge.A_ges,'-r','LineWidth',1.5)
plot(z_vec,cut_plane_props_alpha_shape.A_ges,'-g','LineWidth',1.5)
plot(z_vec,cut_plane_props_ref_sponge.A_convhull,'-b','LineWidth',1.5)
plot(z_vec,cut_plane_props_alpha_shape.A_convhull,'-k','LineWidth',1.5)
legend('A_{ges,sponge}','A_{ges,concave hull}',...
       'A_{convex,sponge}','A_{convex,concave hull}')
xlabel('z [nm]')
ylabel('A')
grid on

subplot(1,2,2)
plot(z_vec,phi_threshold_z,'-b','LineWidth',1.5)
xlabel('z [nm]')
ylabel('\phi_T')
grid on
drawnow
end