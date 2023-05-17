function relim(nodes, x_margin, y_margin)
xlim([min(nodes(:, 1)) - x_margin, max(nodes(:, 1)) + x_margin])
ylim([min(nodes(:, 2)) - y_margin, max(nodes(:, 2)) + y_margin])
end
