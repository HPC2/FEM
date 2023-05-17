function label_elements(elements, nodes)
was_holded = ishold;
hold on

x_offset = 0;
y_offset = 0;

was_labeled = false(size(elements, 1));

for i = 1:size(elements, 1)
    for j = 1:3
        if ~was_labeled(i)
            x = mean(nodes(elements(i, 1:3)+1, 1)) + x_offset;
            y = mean(nodes(elements(i, 1:3)+1, 2)) + y_offset;
            text(x, y, sprintf("E%i\nAff%i", i-1, elements(i, 7)), "Color", "green");
            was_labeled(i) = true;
        end
    end
end

if ~was_holded
    hold off
end

end
