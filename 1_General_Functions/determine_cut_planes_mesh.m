function [Ordered_Connectivity_single,Ordered_Intersection_Points]=...
    determine_cut_planes_mesh(Points,Connectivity,Normals_plane,V0_plane)
% Description: this function returns the cut line of all planes specified
% by the points V0_plane and the normal vectors n_plane and the 
% triangulation of the geometry

% Input: 
% 1. Points (Mx3): coordinates of the mesh vertices
% 2. Connectivity (Nx3): connectivity matrix of the mesh
% 3. Normals_plane (Kx3): normal vectors of the cut planes
% 4. V0_plane (Kx3): point on each plane
% Output: 
% 1. Ordered_Connectivity_single {K}x{T}: cell arrays with the ordered connectivities of each cut loop
% {j} for each cut plane {i} i.e. Ordered_Connectivity_single{i}{j}=[v1,v2; v2 v3;...]
% 2. Ordered_Intersection_Points {K}: cell array with the coordinates of all intersection points of the {i} cut
% plane

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

Num_planes=length(Normals_plane(:,1));
Ordered_Connectivity_single=cell(Num_planes,1);
Ordered_Intersection_Points=cell(Num_planes,1);
for p=1:Num_planes
    fprintf('Determine polyshapes for the %1.0f. cut plane\n',p);
    Points_shifted=Points;
    
    % calculate the sideness of all vertices
    Sideness=determine_sideness_of_points(Points_shifted,Normals_plane(p,:),V0_plane(p,:));
    
    % pertubate the single vertices until the sideness can be clearly
    % determined
    if(any(Sideness==0))
        random_edge_length=norm(Points(Connectivity(1,1),:)-Points(Connectivity(1,2),:))*1e-6;
        
        while(any(Sideness==0))
            Points_shifted=Points;
            Points_shifted(Sideness==0,:)=Points_shifted(Sideness==0,:)+random_edge_length*rand(sum(Sideness==0),3);
            
            Sideness=Determine_Sideness_of_Points_regarding_Plane(Points_shifted,Normals_plane(p,:),V0_plane(p,:));
        end
    end
    % calculate the sideness of all vertices
    Sideness_v1=Sideness(Connectivity(:,1));
    Sideness_v2=Sideness(Connectivity(:,2));
    Sideness_v3=Sideness(Connectivity(:,3));
    
    % determine the sideness of all vertices regarding the specified plane
    Sideness_Triangles=Sideness_v1+Sideness_v2+Sideness_v3;
    
    % indices of all sliced triangles of the current connectivity
    inds_triangles_sliced=Sideness_Triangles==1 | Sideness_Triangles==-1;
    
    % combine the sideness in a matrix --> each row contains either two
    % times +1 or two times -1 --> the sliced edges can be determined due to
    % the combination of +1 and -1
    Sideness_per_vertex_sliced_triangles=[Sideness_v1(inds_triangles_sliced)...
        Sideness_v2(inds_triangles_sliced) Sideness_v3(inds_triangles_sliced)];
    
    sum_rows_minus_one=sum(Sideness_per_vertex_sliced_triangles,2)==-1;
    Sideness_per_vertex_sliced_triangles(sum_rows_minus_one,:)=...
        -Sideness_per_vertex_sliced_triangles(sum_rows_minus_one,:);
    
    Con_Sliced_Tris=Connectivity(inds_triangles_sliced,:);
    
    if(~isempty(Con_Sliced_Tris))
        [~,inds_other_sideness]=max(Sideness_per_vertex_sliced_triangles==-1,[],2);
        
        % initialization
        inds_all_triangles=[];
        inds_all_triangles(:,1)=1:length(Con_Sliced_Tris(:,1));
        vec_inds_triangles=1:3;
        inds_triangles_tested=zeros(length(Con_Sliced_Tris(:,1)),1);
        inds_triangles_tested(1)=1;
        inds_B_C=vec_inds_triangles(vec_inds_triangles~=inds_other_sideness(1));
        
        current_con=[Con_Sliced_Tris(1,inds_other_sideness(1))...
                     Con_Sliced_Tris(1,inds_B_C(1))];
        
        num_loop=1;
        ind_next_triangle=1;
        j=2;
        ind_start=1;
        Ordered_Connectivity_single{p}{num_loop}=[1,2];    
        for k=1:length(Con_Sliced_Tris(:,1))
            % find all triangles that contain the current edge (always two triangles)
            triangles_with_current_edge=sum((current_con(end,1)==Con_Sliced_Tris)...
                                           +(current_con(end,2)==Con_Sliced_Tris),2)==2;
            
            % find the next triangle with the current edge that is not
            % equal to current one
            ind_next_triangle=inds_all_triangles(triangles_with_current_edge...
                                      & inds_all_triangles~=ind_next_triangle);
            
            % find the first vertex of the next line segment --> this
            % is the vertex that is not part of the current edge
            ind_next_vertex_1=vec_inds_triangles(Con_Sliced_Tris(ind_next_triangle,:)...
                    ~=current_con(end,1) & Con_Sliced_Tris(ind_next_triangle,:)...
                    ~=current_con(end,2));
            
            % determine the sideness of the first vertex
            sideness_vertex_1=Sideness_per_vertex_sliced_triangles(ind_next_triangle,...
                                                                   ind_next_vertex_1);
            
            % find all vertices of the next triangle with a opposite
            % sideness regarding vertex 1 --> one or two vertices can
            % be returned
            inds_next_vertex_2=vec_inds_triangles(Sideness_per_vertex_sliced_triangles(ind_next_triangle,:)...
                             ~=sideness_vertex_1);
            
            % always use the first vertex as second vertex of the next
            % edge element
            ind_next_vertex_2=inds_next_vertex_2(1);
            
            % update the current connectivity
            current_con(j,:)=[Con_Sliced_Tris(ind_next_triangle,ind_next_vertex_1)...
                              Con_Sliced_Tris(ind_next_triangle,ind_next_vertex_2)];
            
            % track the connectivity of the intersection points
            Ordered_Connectivity_single{p}{num_loop}(end+1,:)=[j,j+1];
            
            % test if the next triangle was already tested
            if(ismember(ind_next_triangle,inds_triangles_tested))
                % now there are two possibilities --> 1. triangle was
                % tested and we have already tested all other triangles
                % 2. there are still other triangles we have to test
                % because the plane slices our domain in multiple
                % unconnected polygons
                
                ind_next_triangle_other_dom=find(~ismember(inds_all_triangles,...
                                            inds_triangles_tested),1);
                
                % 2. case
                if(~isempty(ind_next_triangle_other_dom))
                    % start again with the next unused triangle of the
                    % remaining unchecked domain
                    ind_next_triangle=ind_next_triangle_other_dom(1);
                    inds_B_C=vec_inds_triangles(vec_inds_triangles...
                                 ~=inds_other_sideness(ind_next_triangle));
                             
                    current_con(j+1,:)=[Con_Sliced_Tris(ind_next_triangle,...
                                        inds_other_sideness(ind_next_triangle))...
                                        Con_Sliced_Tris(ind_next_triangle,inds_B_C(1))];
                    
                    Ordered_Connectivity_single{p}{num_loop}(end,:)  =[j,ind_start];
                    ind_start=j+1;
                    j=j+1;
                    num_loop=num_loop+1;
                    
                    Ordered_Connectivity_single{p}{num_loop}(1,:)=[j,j+1];
                else
                    % update the connectivity list
                    Ordered_Connectivity_single{p}{num_loop}(end,:)=[j,ind_start];
                    ind_start=j+1;
                end
            end
            
            % update the list with all tested triangles
            inds_triangles_tested(k+1)=ind_next_triangle;
            j=j+1;
        end
        % get all start and end points of the intersecting line
        % elements
        P1=Points_shifted(current_con(:,1),:);
        P2=Points_shifted(current_con(:,2),:);
        
        % calculate intersection points
        Intersection_Points=ray_plane_intersection(Normals_plane(p,:),...
                                                   V0_plane(p,:),P1,P2);        
        Ordered_Intersection_Points{p}=Intersection_Points;
    else
        Ordered_Connectivity_single{p}=[];
        Ordered_Intersection_Points{p}=[];
    end
end
end