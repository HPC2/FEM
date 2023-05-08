clear 
close 

nodes = readmatrix('rectangle_10x10_co.txt');
Element = readmatrix('rectangle_10x10_el.txt');
BD = readmatrix('rectangle_10x10_bd.txt');

Element = Element(:,1:3) + 1;
nodes = nodes/10;

BD(:,1:2) = BD(:,1:2) + 1;

Element2 = Element';
nodes2 = nodes';

pdeplot(nodes2, Element2,ElementLabels="on",NodeLabels="on")
hold on 

BD_mat = nodes(BD(:,1),:);
BD_mat(end+1,:) = BD_mat(1,:);
plot(BD_mat(:,1),BD_mat(:,2),'r','LineWidth',2)