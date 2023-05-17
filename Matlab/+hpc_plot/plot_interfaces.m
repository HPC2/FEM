function plot_interfaces(interfaces, nodes, varargin)
p = inputParser;
addRequired(p, "interfaces");
addRequired(p, "nodes");
addParameter(p, 'x_offset', -0.05);
addParameter(p, 'y_offset', -0.05);
addParameter(p, 'shrink', 0.2);
addParameter(p, "colors", ["green", "magenta", "red", "blue"]);
parse(p, interfaces, nodes, varargin{:});

offset = [p.Results.x_offset, p.Results.y_offset];
shrink = p.Results.shrink;
colors = p.Results.colors;

was_holded = ishold;
hold on


for i = 1:size(interfaces, 1)
    node_id_1 = interfaces(i, 1);
    node_id_2 = interfaces(i, 2);

    xy_1 = nodes(node_id_1+1, :);
    xy_2 = nodes(node_id_2+1, :);
    dir = xy_2 - xy_1;
    dir_ortho = [dir(2),+dir(1)];

    xy_1_shrinked = xy_1 + shrink * dir + dir_ortho .* offset;
    xy_2_shrinked = xy_1 + (1-shrink) * dir + dir_ortho .* offset;
    
    x = [xy_1_shrinked(1),xy_2_shrinked(1)];
    y = [xy_1_shrinked(2),xy_2_shrinked(2)];

    color = colors(interfaces(i,5) + 1);

    plot(x,y,"Color",color);
end

if ~was_holded
    hold off
end

end
