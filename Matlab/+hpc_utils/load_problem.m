function [nodes, elements, boundaries] = load_problem(problem)
path = sprintf('../Problem/%s',problem);

nodes       = readmatrix(sprintf('%s.co',path), 'FileType', 'text');
elements    = readmatrix(sprintf('%s.el',path), 'FileType', 'text');
boundaries  = readmatrix(sprintf('%s.bd',path), 'FileType', 'text');
end

