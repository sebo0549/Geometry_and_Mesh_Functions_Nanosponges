function results=calculate_properties_cutplanes(Points,Connectivity,z_vec,...
                                                plot_true_false)
% Description: this function calculates the geometric quantities for a
% given set of cut planes through a 3D triangulation

% Input:  
% 1. Points (Mx3): coordinates of the mesh vertices
% 2. Connectivity (Nx3): connectivity matrix of the mesh
% 3. z_vec (Kx1): z-coordinates of the cut planes
% 4. plot_true_false (1x1): boolean true if results should be plotted
% Output:
% 1. results: struct containing the geometric information for each cut
% plane 

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

warning('off','MATLAB:polyshape:repairedBySimplify');

%% determine region bounds
x_min=min(Points(:,1))-1;
x_max=max(Points(:,1))+1;
y_min=min(Points(:,2))-1;
y_max=max(Points(:,2))+1;
z_min=min(z_vec);
z_max=max(z_vec);

%% set normals and plane points
Nz=length(z_vec);
n_plane=[zeros(Nz,2),ones(Nz,1)];
V0_plane=[zeros(Nz,2),z_vec'];

%% calculate polyshapes for all cut planes 
[Ordered_Connectivity_single,Ordered_Intersection_Points]=...
    determine_cut_planes_mesh(Points,Connectivity,n_plane,V0_plane);

if(plot_true_false)
    figure(1)
    set(figure(1),'units','normalized','outerposition',[0 0 1 1])   % Figure maximieren auf ganzen Bildschirm
end

%% initialize arrays
L_ges=zeros(Nz,1);
L_convhull=zeros(Nz,1);
A_convhull=zeros(Nz,1);
A_ges=zeros(Nz,1);
num_regions=zeros(Nz,1);
Poly_save=cell(Nz,1);

%% perform calculations for each cut plane
for n=1:Nz
    if(~isempty(Ordered_Intersection_Points{n}))
        for k=1:length(Ordered_Connectivity_single{n})
            inds=Ordered_Connectivity_single{n}{k}(:,1);

            % handle nested shapes
            if(k==1)
                polyges=polyshape(Ordered_Intersection_Points{n}(inds,1),Ordered_Intersection_Points{n}(inds,2));
            else
                pol=polyshape(Ordered_Intersection_Points{n}(inds,1),Ordered_Intersection_Points{n}(inds,2));
                polyges = xor(pol,polyges);
            end
        end
        L_ges(n)=polyges.perimeter;
        A_ges(n)=area(polyges);
        num_regions(n)=polyges.NumRegions;

        Poly_save{n}=polyges;
        
        % determine convex hull of the polyshapes
        [inds,A_convhull(n)] = convhull(Ordered_Intersection_Points{n}(:,1),Ordered_Intersection_Points{n}(:,2));
        con_conv=[inds(1:end-1),circshift(inds(1:end-1),-1)];

        L_convhull(n)=sum(sqrt((Ordered_Intersection_Points{n}(con_conv(:,1),1)-Ordered_Intersection_Points{n}(con_conv(:,2),1)).^2+...
            (Ordered_Intersection_Points{n}(con_conv(:,1),2)-Ordered_Intersection_Points{n}(con_conv(:,2),2)).^2));

        if(plot_true_false)
            clf
            M = makehgtform('translate',[0 0 z_vec(n)]);
            t=hgtransform('Matrix',M);
            plot(polyges,'Parent',t,'FaceColor','y')

            axis equal
            set(gca,'XLim',[x_min x_max]...
                ,'YLim',[y_min y_max]...
                ,'ZLim',[z_min z_max])
            grid on
            view(15,35)
            drawnow
        end
    else
        L_convhull(n)=0;
        A_convhull(n)=0;
        L_ges(n)=0;
        A_ges(n)=0;
        num_regions(n)=0;
        Poly_save{n}=[];
    end
end

%% store results in a single struct
results.x_min=x_min;
results.x_max=x_max;
results.y_min=y_min;
results.y_max=y_max;
results.z_min=z_min;
results.z_max=z_max;
results.Points=Points;
results.Connectivity=Connectivity;
results.A_ges=A_ges;
results.A_convhull=A_convhull;
results.num_regions=num_regions;
results.L_ges=L_ges;
results.L_convhull=L_convhull;
results.z_vec=z_vec;
results.z_start=min(Points(:,3))-z_min;
results.z_end=z_max-max(Points(:,3));
results.Poly=Poly_save;
results.Ordered_Connectivity_single=Ordered_Connectivity_single;
results.Ordered_Intersection_Points=Ordered_Intersection_Points;
end