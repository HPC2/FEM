function label_l2g(l2g,nodes)
[n,m] = size(l2g);
for i=1:n
    for j=1:m
        local_node_id = j-1;
        global_node_id = l2g(i,j);
        x = nodes(global_node_id+1,1)-0.06;
        y = nodes(global_node_id+1,2);
        text(x,y,sprintf("G%i\nL%i",global_node_id,local_node_id), "Color","magenta","Clipping", "on");
    end
end
end

