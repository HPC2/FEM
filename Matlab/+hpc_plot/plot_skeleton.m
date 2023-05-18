function plot_skeleton(elements, nodes, varargin)
p = inputParser;
addRequired(p, "elements");
addRequired(p, "nodes");
addParameter(p, "color", "k");
addParameter(p, "linewidth",1);
parse(p, elements, nodes, varargin{:});

color = p.Results.color;
linewidth = p.Results.linewidth;

was_holded = ishold;
hold on
for i = 1:size(elements, 1)
    node_id_1 = elements(i, 1);
    node_id_2 = elements(i, 2);
    node_id_3 = elements(i, 3);
    x = [; ...
        nodes(node_id_1+1, 1); ...
        nodes(node_id_2+1, 1); ...
        nodes(node_id_3+1, 1); ...
        nodes(node_id_1+1, 1); ...
        ];
    y = [; ...
        nodes(node_id_1+1, 2); ...
        nodes(node_id_2+1, 2); ...
        nodes(node_id_3+1, 2); ...
        nodes(node_id_1+1, 2); ...
        ];
    plot(x, y, "Color",color, "LineWidth",linewidth);
end
if ~was_holded
    hold off
end
end
