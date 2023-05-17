clear
close all


[nodes, elements, boundaries] = hpc_utils.load_problem("rectangle_3x3");
interfaces = hpc_utils.load_interfaces("rectangle_3x3");

figure();

hpc_plot.plot_skeleton(elements, nodes);

hpc_plot.mark_nodes(nodes);
hpc_plot.label_nodes(elements, nodes, "x_offset", -0.1, "y_offset", 0.1, "color", "r", "label", "N_g");
hpc_plot.label_edges(elements, nodes,"label","M_g");

hpc_plot.label_elements(elements, nodes);

hpc_plot.plot_interfaces(interfaces, nodes);
hpc_plot.label_interfaces(interfaces, nodes);

hpc_plot.relim(nodes, 0.2, 0.2);
