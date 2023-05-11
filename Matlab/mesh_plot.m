clear
close

nodes = readmatrix('../Problem/rectangle_3x3.co', 'FileType', 'text');
Element = readmatrix('../Problem/rectangle_3x3.el', 'FileType', 'text');
BD = readmatrix('../Problem/rectangle_3x3.bd', 'FileType', 'text');

plot_mesh(nodes, Element, BD);

nodes = readmatrix('../Problem/rectangle_3x3_refined.co', 'FileType', 'text');
Element = readmatrix('../Problem/rectangle_3x3_refined.el', 'FileType', 'text');
BD = readmatrix('../Problem/rectangle_3x3_refined.bd', 'FileType', 'text');

plot_mesh(nodes, Element, BD);