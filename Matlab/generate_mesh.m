% Def
n_rows = 4;
n_cols = 1;

% Calc number of ...
n_nodes = (n_rows + 1) * (n_cols+1);
n_elems = 2 * n_rows * n_cols;
n_edges_h = (n_rows + 1) * n_cols;
n_edges_v = (n_cols + 1) * n_rows;

% coords
fprintf("N\tX\tY\n")
for i = 0 : n_nodes - 1
    x = mod(i, n_cols + 1); % In C: mod --> i%(n_cols+1)
    y = fix(i / n_rows); % In C: integer div --> i/n_cols
    fprintf("%i\t%i\t%i\n", i, x, y);
end

% elements
fprintf("E\tN1\tN2\tN3\tM1\tM2\tM3\n")
for j = 0 : n_elems - 1
    alpha = fix(j / (2 * n_cols)); % In C: integer div --> j/(2*n_cols)
    beta = fix(j / 2); % In C: integer div --> j/2
    if mod(j,2) == 0 % In C: mod --> j%2==0
        N1 = alpha + beta;
        N2 = alpha + beta + n_cols + 2;
        N3 = alpha + beta + n_cols + 1;
        M1 = beta + n_edges_h + n_edges_v;
        M2 = beta + n_cols;
        M3 = alpha + beta + n_edges_h;
    else
        N1 = alpha + beta;
        N2 = alpha + beta + 1;
        N3 = alpha + beta + n_cols + 2;
        M1 = beta;
        M2 = alpha + beta + n_edges_h + 1;
        M3 = beta + n_edges_h + n_edges_v;
    end
    fprintf("%i\t%i\t%i\t%i\t%i\t%i\t%i\n", j, N1, N2, N3, M1, M2, M3);
end

