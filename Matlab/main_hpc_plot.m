clear
close all


[nodes, elements, boundaries] = hpc_utils.load_problem("rectangle_2x2_global_0ref");
%interfaces = hpc_utils.load_interfaces("rectangle_3x3");

tiledlayout(1,1, 'Padding', 'none', 'TileSpacing', 'compact'); 
nexttile

hpc_plot.plot_skeleton(elements, nodes);

hpc_plot.mark_nodes(nodes);
hpc_plot.label_nodes(elements, nodes, "color", "r", "label", "N");
hpc_plot.label_edges(elements, nodes,"label","M");

hpc_plot.label_elements(elements, nodes);

%hpc_plot.plot_interfaces(interfaces, nodes);
%hpc_plot.label_interfaces(interfaces, nodes);

hpc_plot.relim(nodes, 0.2, 0.2);
axis equal