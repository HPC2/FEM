function l2g = load_l2g(fname)
path = sprintf('../Problem/%s',fname);

l2g = readmatrix(sprintf('%s.l2g',path), 'FileType', 'text');
end

