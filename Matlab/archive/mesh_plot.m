clear
close

str = 'rectangle_1x1_refined2times';

nodes = readmatrix(strcat('../Problem/',str,'.co'), 'FileType', 'text');
Element = readmatrix(strcat('../Problem/',str,'.el'), 'FileType', 'text');
BD = readmatrix(strcat('../Problem/',str,'.bd'), 'FileType', 'text');

plot_mesh(nodes, Element, BD);

%%

nodes = readmatrix('../Problem/rectangle_3x3_refined.co', 'FileType', 'text');
Element = readmatrix('../Problem/rectangle_3x3_refined.el', 'FileType', 'text');
BD = readmatrix('../Problem/rectangle_3x3_refined.bd', 'FileType', 'text');

plot_mesh(nodes, Element, BD);