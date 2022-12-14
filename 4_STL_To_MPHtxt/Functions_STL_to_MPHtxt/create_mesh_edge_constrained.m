function fileID=create_mesh_edge_constrained(fileID,ConnectivityList)
% Description: this function creates the constraint entries for the 
% triangles and edges in a .mphtxt mesh file

% Input:  
% 1. fileID (1x1): file to store the data in
% 2. ConnectivityList (Mx3): connectivity list of the triangulation
% Output:  
% 1. fileID (1x1): updated file

% Author: Sebastian Bohm (sebastian.bohm@tu-ilmenau.de)
% Date: 08-12-2022

edges=unique(sort([ConnectivityList(:,1),ConnectivityList(:,2);...
    ConnectivityList(:,1),ConnectivityList(:,3);...
    ConnectivityList(:,2),ConnectivityList(:,3)],2),'rows');

fprintf(fileID,'\n2 # number of element types\n\n');

fprintf(fileID,'# Type #0\n\n');
fprintf(fileID,'3 edg # type name\n\n\n');
fprintf(fileID,'2 # number of vertices per element\n');
fprintf(fileID,sprintf('%1.0f # number of elements\n',length(edges(:,1))));
fprintf(fileID,'# Elements\n');
fprintf(fileID,'%1.0f %1.0f\n',edges'-1);
fprintf(fileID,sprintf('\n%1.0f # number of geometric entity indices\n',...
    length(edges(:,1))));
fprintf(fileID,'# Geometric entity indices\n');
fprintf(fileID,'%1.0f\n',ones(length(edges(:,1)),1));

fprintf(fileID,'\n# Type #1\n\n');
fprintf(fileID,'3 tri # type name\n\n\n');
fprintf(fileID,'3 # number of vertices per element\n');
fprintf(fileID,sprintf('%1.0f # number of elements\n',length(ConnectivityList(:,1))));
fprintf(fileID,'# Elements\n');
fprintf(fileID,'%1.0f %1.0f %1.0f\n',ConnectivityList'-1);
fprintf(fileID,sprintf('\n%1.0f # number of geometric entity indices\n',...
    length(ConnectivityList(:,1))));
fprintf(fileID,'# Geometric entity indices\n');
fprintf(fileID,'%1.0f\n',0:length(ConnectivityList(:,1))-2);
fprintf(fileID,'%1.0f\r',length(ConnectivityList(:,1))-1);
end