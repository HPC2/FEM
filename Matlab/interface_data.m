n_cols = 2;
n_rows = 4;

fprintf("Source\tTarget\tPleft\tPright\tCorrespondingMidpoint\n")
for i=0:(n_rows-1)*n_cols-1
    beta = fix(i/n_cols);
    source = beta + n_cols + i + 1;
    target = beta + n_cols + i + 2;
    left_p = n_cols + i;
    right_p = i;
    corresponding_midpoint = n_cols + i;
    fprintf("%i\t%i\t%i\t%i\t%i\n",source,target,left_p,right_p,corresponding_midpoint)
end
for i=0:(n_cols-1)*n_rows-1
    alpha = mod(i,n_rows);
    beta = fix(i/n_rows);
    gamma = mod(i,n_cols-1);
    delta = fix(i/(n_cols-1));
    source = (n_cols+1)*alpha + beta +1;
    target = (n_cols+1)*alpha + beta +1 + n_cols + 1;
    left_p =  n_cols*alpha + beta;
    right_p = n_cols*alpha + beta + 1;
    corresponding_midpoint = gamma + delta * (n_cols + 1) + n_cols*(n_rows+1) + 1;
    fprintf("%i\t%i\t%i\t%i\t%i\n",source,target,left_p,right_p,corresponding_midpoint)
end
