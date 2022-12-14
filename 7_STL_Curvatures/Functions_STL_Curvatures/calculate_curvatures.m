function [K_Mean,K_Gauss,kappa1,kappa2,K_Vector,Total_mean_curvature,...
    Total_gauss_curvature,inds_undecidable,inds_equal,A_mixed,Thetas_ges,VN]...
    =calculate_curvatures(varargin)
% Description: this function calculates the curvature of a 3D triangulation
% based on: https://link.springer.com/chapter/10.1007/978-3-662-05105-4_2
% DOI: 10.1007/978-3-662-05105-4_2

% Input:  
% 1. Connectivity (Nx3): connectivity matrix of the mesh
% 2. Points (Mx3): coordinates of the mesh vertices
% 3. string determining the type of vertex normal calculation: 
%    VertexNormals_Faces' (default) or 'VertexNormals_Midpoints' 
% 4. threshold_angle (1x1): threshold in rad (default 0) that determines
%    which deviations between the directions of the normal vectors and the 
%    curvature vectors are tolerated while determining the sign of the mean
%    curvature
% Output:
% 1. K_Mean (Mx1):  mean vertex curvature
% 2. K_Gauss (Mx1): gauss vertex curvature
% 3. kappa1 (Mx1): larger principal vertex curvature
% 4. kappa2 (Mx1): smaller principal vertex curvature
% 5. K_Vektor (Mx3): vertex curvature vector
% 6. Total_mean_curvature (1x1): summed up mean curvature
% 7. Total_gauss_curvature (1x1): summed up gauss curvature
% 8. inds_undecidable (Px1): vector with indices for the vertices whose 
% normal direction was undecidable
% 9. inds_equal (Jx1): vector with indices of the vertices where 
% K_Mean.^2-K_Gauss is smaller than zero
% 10. A_mixed (Mx1): area for each vertex which was used for the curvature
% calculation
% 11. Thetas_ges (Mx1): summed angles for each vertex of the attached
% 1-ring
% 12. VN (Mx3): normalized vertex normal vectors

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

%% Processing of the input variables
if(length(varargin)<2)
    disp('Not enough input arguments')
    return
else
    con=varargin{1};
    Points=varargin{2};
    N_Points=length(Points(:,1));
    tr = triangulation(con,Points);
    threshold_angle=0;
    vertex_midpoints=false;
    A_mixed=zeros(N_Points,1);
    Thetas_ges=zeros(N_Points,1);
    inds_triangles_per_vertex = vertexAttachments(tr);

    if(length(varargin)==2) % set default normals but avoid calculating twice
        VN = vertexNormal(tr);
    end
end

if(length(varargin)>=3) % set user specified way to calculate normals
    if(isequal(varargin{3},'VertexNormals_Faces'))
        VN = vertexNormal(tr);
    elseif(isequal(varargin{3},'VertexNormals_Midpoints'))
        VN=zeros(N_Points,3);
        vertex_midpoints=true;
        face_center = incenter(tr);
        face_normal = faceNormal(tr);
    else 
        disp('Unsupported input. Should be VertexNormals_Faces or VertexNormals_Midpoints. I continue with VertexNormals_Faces.')
        VN = vertexNormal(tr);
    end

    if(length(varargin)==4) % set user specified angle threshold
        threshold_angle=varargin{4};
    end
end

%% calculation of the curvature vector for each vertex
K_vertex=zeros(N_Points,3);
for n=1:N_Points
    attached_tris=con(inds_triangles_per_vertex{n},:);
    inds_other_vertices=attached_tris~=n;

    if(vertex_midpoints) % if true calculate midpoint weighted vertex normals
        w_normals=1./sqrt(sum((face_center(inds_triangles_per_vertex{n},:)-Points(n,:)).^2,2));
        VN(n,:)=sum(face_normal(inds_triangles_per_vertex{n},:).*w_normals);
    end

    % loop over all attached triangles of the current vertex
    for k=1:length(attached_tris(:,1))
        other_vertices=attached_tris(k,inds_other_vertices(k,:)); 
        
        % edge vectors of the current attached triangle
        v1 = Points(n,:) - Points(other_vertices(1),:); 
        v2 = Points(n,:) - Points(other_vertices(2),:);
        v3 = Points(other_vertices(2),:) - Points(other_vertices(1),:);
        
        gamma = atan2(norm(cross(v1,v2)),dot(v1,v2)); % angle at x for act. T
        beta  = atan2(norm(cross(v1,v3)),dot(v1,v3));
        alpha = pi - beta - gamma;

        Thetas_ges(n) = Thetas_ges(n) + gamma; 

        K_vertex(n,:) = K_vertex(n,:) + v1/tan(alpha) + v2/tan(beta); % curvature vector

        if(gamma >= pi/2 || beta >= pi/2 || alpha >= pi/2) % T is obtuse
            A = 0.5*norm(cross(v1,v2)); % area triangle
            if(gamma >= pi/2) % T is obtuse at x
                A_mixed(n) = A_mixed(n) + A/2;
            else
                A_mixed(n) = A_mixed(n) + A/4;
            end
        else % voronoi area
            A_mixed(n) = A_mixed(n) + (1/8)*(norm(v1)^2/tan(alpha)+norm(v2)^2/tan(beta));
        end
    end
end
K_Vector= K_vertex./(2*A_mixed);

% determine sign of the mean curvature
VN = VN./sqrt(sum(VN.^2,2)); % normalize normals
angle_normals = atan2(sqrt(sum(cross(K_Vector,VN).^2,2)),dot(K_Vector,VN,2));
in_or_out=-sign(angle_normals-pi/2);
inds_undecidable = abs(pi/2-abs(angle_normals)) < threshold_angle;   % exclude small angles

% include small curvature vectors
length_K=sqrt(sum(K_Vector.^2,2));
mean_K=mean(length_K);
inds_undecidable(length_K<mean_K*1e-10)=false;               

% calculate mean and gauss curvature
K_Mean  = 0.5*length_K.*in_or_out;
K_Gauss = (2*pi-Thetas_ges)./A_mixed;

% calculate principal curvatures --> this fails for points where K_mean^2~K_Gauss
% --> set kappa_1=kappa_2 and return the indices
delta_x = K_Mean.^2-K_Gauss;
inds_equal = delta_x<0;
delta_x(inds_equal) = 0;
sqrt_delta_x=sqrt(delta_x);

% kappa1<kappa2
kappa1 = K_Mean - sqrt_delta_x;
kappa2 = K_Mean + sqrt_delta_x;

% set undecidable values equal to NaN
kappa1(inds_undecidable)=NaN;
kappa2(inds_undecidable)=NaN;
kappa1(inds_equal)=NaN;
kappa2(inds_equal)=NaN;
K_Mean(inds_undecidable)=NaN;
K_Gauss(inds_undecidable)=NaN;
A_mixed(inds_undecidable)=NaN;
Thetas_ges(inds_undecidable)=NaN;
K_Vector(inds_undecidable,:)=NaN;
VN(inds_undecidable,:)=NaN;

% calculate total curvatures as integral over the surface
Total_gauss_curvature = sum(A_mixed.*K_Gauss,'omitnan'); 
Total_mean_curvature  = sum(A_mixed.*K_Mean,'omitnan'); 
end