function [interfaces] = load_interfaces(problem)
path = sprintf('../Problem/%s',problem);

interfaces = readmatrix(sprintf('%s.coupl.if',path), 'FileType', 'text');
end

