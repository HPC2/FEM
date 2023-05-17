function label_nodes(elements, nodes, varargin)
p = inputParser;
addRequired(p, "elements");
addRequired(p, "nodes");
addParameter(p, 'x_offset', 0);
addParameter(p, 'y_offset', 0);
addParameter(p, 'label', "N");
addParameter(p, "color", "red");
parse(p, elements, nodes, varargin{:});

x_offset = p.Results.x_offset;
y_offset = p.Results.y_offset;
format_str = p.Results.label + "%i";
color = p.Results.color;

was_holded = ishold;
hold on

was_labeled = false(size(nodes, 1));

for i = 1:size(elements, 1)
    for j = 1:3
        node_id = elements(i, j);
        if ~was_labeled(node_id+1)
            x = nodes(node_id+1, 1) + x_offset;
            y = nodes(node_id+1, 2) + y_offset;
            text(x, y, sprintf(format_str, node_id), "Color", color);
            was_labeled(node_id+1) = true;
        end
    end
end

if ~was_holded
    hold off
end

end
