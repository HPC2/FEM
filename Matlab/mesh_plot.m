clear

nodes = readmatrix('../Problem/rectangle_3x3_refined.co', 'FileType', 'text');
Element = readmatrix('../Problem/rectangle_3x3_refined.el', 'FileType', 'text');
BD = readmatrix('../Problem/rectangle_3x3_refined.bd', 'FileType', 'text');

% nodes = readmatrix('../Problem/rectangle_3x3.co', 'FileType', 'text');
% Element = readmatrix('../Problem/rectangle_3x3.el', 'FileType', 'text');
% BD = readmatrix('../Problem/rectangle_3x3.bd', 'FileType', 'text');


edges = Element(:, 4:6);

Element = Element(:,1:3) + 1;

midpoints1 = (nodes(Element(:, 1), :) + nodes(Element(:, 2), :))/2;
midpoints2 = (nodes(Element(:, 2), :) + nodes(Element(:, 3), :))/2;
midpoints3 = (nodes(Element(:, 3), :) + nodes(Element(:, 1), :))/2;

BD(:,1:2) = BD(:,1:2) + 1;

Element2 = Element';
nodes2 = nodes';

figure;

pdeplot(nodes2, Element2, 'ElementLabels', "on", 'NodeLabels', "on");
hold on 

text(midpoints1(:, 1), midpoints1(:, 2), num2str(edges(:, 1)));
text(midpoints2(:, 1), midpoints2(:, 2), num2str(edges(:, 2)));
text(midpoints3(:, 1), midpoints3(:, 2), num2str(edges(:, 3)));

BD_mat = nodes(BD(:,1),:);
BD_mat(end+1,:) = BD_mat(1,:);
plot(BD_mat(:,1),BD_mat(:,2),'r','LineWidth',2)