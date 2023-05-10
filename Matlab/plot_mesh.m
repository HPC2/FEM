function plot_mesh(nodes,Element,BD)

    edges = Element(:, 4:6);
    Element = Element(:,1:3) + 1;
    BD(:,1:2) = BD(:,1:2) + 1;

    pdeplot(nodes', Element',ElementLabels="on",NodeLabels="on")
    hold on  
    BD_mat = nodes(BD(:,1),:);
    BD_mat(end+1,:) = BD_mat(1,:);
    plot(BD_mat(:,1),BD_mat(:,2),'r','LineWidth',2)
    
    
    
    midpoints1 = (nodes(Element(:, 1), :) + nodes(Element(:, 2), :))/2;
    midpoints2 = (nodes(Element(:, 2), :) + nodes(Element(:, 3), :))/2;
    midpoints3 = (nodes(Element(:, 3), :) + nodes(Element(:, 1), :))/2;
    
    text(midpoints1(:, 1), midpoints1(:, 2), num2str(edges(:, 1)));
    text(midpoints2(:, 1), midpoints2(:, 2), num2str(edges(:, 2)));
    text(midpoints3(:, 1), midpoints3(:, 2), num2str(edges(:, 3)));
    



end