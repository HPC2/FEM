clc; clear; close all;

elements = readmatrix('../Problem/rectangle_3x3_refined.el', 'FileType','text');
rows = 3;
cols = 3;
nof_nodes = 100 * cols*(rows-1) + rows*(cols-1); %ncoords in C

node_idx = 1:3;
%m_idx    = 4:6;
aff_idx  = 7;

%node_ptr = ones(rows*cols,1); % which process is at which node(s)?, startet in C bei 0
for j=1:rows*cols
    nodeslist{j} = zeros(nof_nodes,1);
end

for i=1:size(elements,1)
    affiliation = elements(i,aff_idx) + 1; %index shift!
    nodeslist{affiliation}(elements(i,node_idx) + 1) = 1; %attention, node_idx is [1,2,3]; 3 lines in C!
                                                          % no index shift
                                                          % needed in C
    %nodeslist{affiliation}(node_ptr(affiliation)+[0,1,2]) = [elements(i,node_idx)]';
    %node_ptr(affiliation) = node_ptr(affiliation) + 3;
end

for j=1:rows*cols
    % unique nodes has (2^f + 1)^2 entries where f is the refinement step
    % quicksort/mergesort und dann doppelte rauswerfen -> O(n*log(n))
    %nodeslist{j} = unique(nodeslist{j}(1:node_ptr(j)-1));
    % jetzt kein qicksort mehr, sondern für alle Prozesse linear durch
    % logical array gehen und wenn an der Stelle eine 1 steht, gehört der
    % Knoten zum Prozessor
    nodeslist{j} = find(nodeslist{j}) -1 ; %index shift!!!
end