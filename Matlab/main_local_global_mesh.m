clear
close all

% Set vars
n_rows = 3;
n_cols = 3;
n_refs = 3;
process_id_under_investigation = 6;

% Build file name
fname = "rectangle_%ix%i_global_%iref";
fname_local_no_ref = sprintf(fname, 1, 1, 0);
fname_local_ref = sprintf(fname, 1, 1, n_refs);
fname_global_no_ref = sprintf(fname, n_rows, n_cols, 0);
fname_global_ref = sprintf(fname, n_rows, n_cols, n_refs);
fname_l2g = sprintf("l2g_%ix%i_%iref", n_rows, n_cols,n_refs);

% Load local mesh without refinements
[; ...
    nodes_local_no_ref, ...
    elements_local_no_ref, ...
    boundaries_local_no_ref, ...
    ] = hpc_utils.load_problem(fname_local_no_ref);

% Load local mesh with refinements
[; ...
    nodes_local_ref, ...
    elements_local_ref, ...
    boundaries_local_ref, ...
    ] = hpc_utils.load_problem(fname_local_ref);

% Load global mesh without refinements
[; ...
    nodes_global_no_ref, ...
    elements_global_no_ref, ...
    boundaries_global_no_ref, ...
    ] = hpc_utils.load_problem(fname_global_no_ref);

% Load global mesh with refinements
[; ...
    nodes_global_ref, ...
    elements_global_ref, ...
    boundaries_global_ref, ...
    ] = hpc_utils.load_problem(fname_global_ref);

% Load local to global numbering
l2g = hpc_utils.load_l2g(fname_l2g);

% Calc offset for local mesh
local_offset = [
    mod(process_id_under_investigation,n_cols);
    fix(process_id_under_investigation/n_cols);
];
nodes_local_ref_off = nodes_local_ref + local_offset';
nodes_local_no_ref_off = nodes_local_no_ref + local_offset';

% Plot
figure();

hpc_plot.plot_skeleton(elements_global_ref, nodes_global_ref, "color", "#999999");
hpc_plot.plot_skeleton(elements_local_no_ref, nodes_local_no_ref_off, "linewidth", 2);

hpc_plot.mark_nodes(nodes_global_ref);

hpc_plot.label_nodes(elements_global_ref, nodes_global_ref, "x_offset", 0.01, "y_offset", 0.01, "color", "r", "label", "N_G");
hpc_plot.label_nodes(elements_local_ref, nodes_local_ref_off, "x_offset", 0.01, "y_offset", -0.01, "color", "b", "label", "N_L");

hpc_plot.label_l2g(l2g(process_id_under_investigation+1,:),nodes_global_ref);

axis equal;
hpc_plot.relim(nodes_local_ref_off, 0.2, 0.2);

