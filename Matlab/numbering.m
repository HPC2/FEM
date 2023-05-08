clc; clear; close all;

elements = readmatrix('../Problem/rectangle_3x3_refined.el', 'FileType','text');
rows = 3;
cols = 3;
nof_nodes_per_process = 3*size(elements,1);

node_idx = 1:3;
m_idx    = 4:6;
aff_idx  = 7;

node_ptr = ones(rows*cols,1); % which process is at which node(s)?, startet in C bei 0
for j=1:rows*cols
    nodeslist{j} = zeros(nof_nodes_per_process,1);
end

for i=1:size(elements,1)
    affiliation = elements(i,aff_idx) + 1;
    nodeslist{affiliation}(node_ptr(affiliation)+[0,1,2]) = [elements(i,node_idx)]';
    node_ptr(affiliation) = node_ptr(affiliation) + 3;
end

for j=1:rows*cols
    % quicksort/mergesort und dann doppelte rauswerfen -> O(n*log(n))
    nodeslist{j} = unique(nodeslist{j}(1:node_ptr(j)-1));
end