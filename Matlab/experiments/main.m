clear;clc;close all;

res = Result("test.csv");

% Example: Unique settings for 12x1 Domain, cg only, parallel only
experiment_catalog = res.experiment_catalog();
rf = rowfilter(experiment_catalog);
experiment_catalog_12x1_cg_para = experiment_catalog(...
    rf.n_rows==12 &...
    rf.n_cols==1 &...
    rf.solver=="cg" &...
    rf.n_processors > 1 ...
    ,:);


% Iterate over example setting
for row_idx=1:height(experiment_catalog)
    experiment = experiment_catalog(row_idx,:);
    setting_as_cell = table2cell(experiment);
    data = res.find_experiments(setting_as_cell{:});
    %...
end
