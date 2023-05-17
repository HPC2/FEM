function label_nodes_from_nodes(nodes, format_str)
was_holded = ishold;
hold on

if ~exist('format_str', 'var')
    format_str = "N_N%i";
end

x_offset = 0.05;
y_offset = -0.05;

for node_id = 0:size(nodes, 1) - 1
    x = nodes(node_id+1, 1) + x_offset;
    y = nodes(node_id+1, 2) + y_offset;
    text(x, y, sprintf(format_str, node_id), "Color", "red");
end

if ~was_holded
    hold off
end
end
