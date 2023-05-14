clear 
close all

n_rows = 4;
n_cols = 3;


fprintf("Source\tTarget\tPleft\tPrigh\tMidP\tColor\n")

source_offset        = n_cols + 1;
target_offset        = n_cols + 2;
left_process_offset  = n_cols;
right_process_offset = 0;
midpoint_offset      = n_cols;

for i=0 : n_cols * (n_rows - 1) - 1

    alpha    = mod(i,n_cols);
    beta     = fix(i/n_cols);   

    source   = beta + i + source_offset;
    target   = beta + i + target_offset;
    left_p   =        i + left_process_offset;
    right_p  =        i + right_process_offset;
    midpoint =        i + midpoint_offset;
    color    = mod(alpha + mod(beta,2),2);
    
    fprintf("%i\t%i\t%i\t%i\t%i\t%i\n",source,target,left_p,right_p,midpoint,color);
end

source_offset        = 1;
target_offset        = n_cols + 2;
left_process_offset  = 0;
right_process_offset = 1;
midpoint_offset      = n_cols*(n_rows+1) + 1;

for i=0 : (n_cols - 1) * n_rows - 1
    
    alpha    = mod(i,n_cols-1); % Obacht: -1
    beta     = fix(i/(n_cols-1)); % Obacht: -1

    source   = 2 * beta + i + source_offset;
    target   = 2 * beta + i + target_offset;
    left_p   =     beta + i + left_process_offset;
    right_p  =     beta + i + right_process_offset;
    midpoint = 2 * beta + i + midpoint_offset;
    color    = mod(alpha + mod(beta,2) + 1,2); % Obacht: +1

    fprintf("%i\t%i\t%i\t%i\t%i\t%i\n",source,target,left_p,right_p,midpoint,color);
end
