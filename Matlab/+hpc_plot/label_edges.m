function label_edges(elements, nodes, varargin)
p = inputParser;
addRequired(p, "elements");
addRequired(p, "nodes");
addParameter(p, 'x_offset', 0);
addParameter(p, 'y_offset', 0);
addParameter(p, 'label', "M");
addParameter(p, "color", "magenta");
parse(p, elements, nodes, varargin{:});

x_offset = p.Results.x_offset;
y_offset = p.Results.y_offset;
format_str = p.Results.label + "%i";
color = p.Results.color;

was_holded = ishold;
hold on

n_edges = max(elements(:, 4:6), [], "all") + 1;
was_labeled = false(n_edges);

for i = 1:size(elements, 1)
    node_id_1 = elements(i, 1);
    node_id_2 = elements(i, 2);
    node_id_3 = elements(i, 3);
    edge_id_1 = elements(i, 4);
    edge_id_2 = elements(i, 5);
    edge_id_3 = elements(i, 6);


    if ~was_labeled(edge_id_1+1)
        x = mean(nodes([node_id_1 + 1, node_id_2 + 1], 1)) + x_offset;
        y = mean(nodes([node_id_1 + 1, node_id_2 + 1], 2)) + y_offset;
        text(x, y, sprintf(format_str, edge_id_1), "Color", color, "Clipping", "on");
        was_labeled(edge_id_1+1) = true;
    end

    if ~was_labeled(edge_id_2+1)
        x = mean(nodes([node_id_2 + 1, node_id_3 + 1], 1)) + x_offset;
        y = mean(nodes([node_id_2 + 1, node_id_3 + 1], 2)) + y_offset;
        text(x, y, sprintf(format_str, edge_id_2), "Color", color, "Clipping", "on");
        was_labeled(edge_id_2+1) = true;
    end

    if ~was_labeled(edge_id_3+1)
        x = mean(nodes([node_id_3 + 1, node_id_1 + 1], 1)) + x_offset;
        y = mean(nodes([node_id_3 + 1, node_id_1 + 1], 2)) + y_offset;
        text(x, y, sprintf(format_str, edge_id_3), "Color", color, "Clipping", "on");
        was_labeled(edge_id_3+1) = true;
    end
end

if ~was_holded
    hold off
end

end
