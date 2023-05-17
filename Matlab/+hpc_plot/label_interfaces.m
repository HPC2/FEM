function label_interfaces(interfaces, nodes, varargin)
p = inputParser;
addRequired(p, "interfaces");
addRequired(p, "nodes");
addParameter(p, 'x_offset', -0.1);
addParameter(p, 'y_offset', +0.1);
addParameter(p, 'label', "Interface");
addParameter(p, "color", "k");
parse(p, interfaces, nodes, varargin{:});

x_offset = p.Results.x_offset;
y_offset = p.Results.y_offset;
format_str = p.Results.label + "\n(%i->%i)\n(l_p=%i,r_p=%i)\n(c=%i)";
color = p.Results.color;

was_holded = ishold;
hold on

for i = 1:size(interfaces, 1)
    node_id_1 = interfaces(i, 1);
    node_id_2 = interfaces(i, 2);
    p_left = interfaces(i, 3);
    p_right = interfaces(i, 4);
    c = interfaces(i, 5);

    x = mean([nodes(node_id_1+1, 1), nodes(node_id_2+1, 1)]) + x_offset;
    y = mean([nodes(node_id_1+1, 2), nodes(node_id_2+1, 2)]) + y_offset;
    text_str = sprintf(format_str, node_id_1, node_id_2, p_left, p_right, c);
    text(x, y, text_str, "Color", color);

end

if ~was_holded
    hold off
end

end