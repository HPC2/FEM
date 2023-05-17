clear 
close 



nodes = readmatrix('../Problem/rectangle_3x3.co', 'FileType', 'text');
Element = readmatrix('../Problem/rectangle_3x3.el', 'FileType', 'text');
BD = readmatrix('../Problem/rectangle_3x3.bd', 'FileType', 'text');
%
nodes(:,1) = nodes(:,1)/max(nodes(:,1));
nodes(:,2) = nodes(:,2)/max(nodes(:,2));

plot_mesh(nodes,Element,BD)



figure(2)

nodes_ref = readmatrix('../Problem/rectangle_3x3_refined.co', 'FileType', 'text');
Element_ref = readmatrix('../Problem/rectangle_3x3_refined.el', 'FileType', 'text');
BD_ref = readmatrix('../Problem/rectangle_3x3_refined.bd', 'FileType', 'text');

%
nodes_ref(:,1) = nodes_ref(:,1)/max(nodes_ref(:,1));
nodes_ref(:,2) = nodes_ref(:,2)/max(nodes_ref(:,2));


plot_mesh(nodes_ref,Element_ref,BD_ref)







