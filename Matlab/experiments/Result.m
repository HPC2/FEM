classdef Result < handle

    properties
        data
    end
    
    methods
        function obj = Result(path)
            obj.data = readtable(path);
            
        end
        
        function partial_data = find_experiments(obj,n_rows,n_cols,n_refs,solver,n_processors,varargin)
            rf = rowfilter(obj.data);
            partial_data = obj.data(...
                rf.n_rows == n_rows &...
                rf.n_cols == n_cols &...
                rf.n_refs == n_refs &...
                rf.solver == solver &...
                rf.n_processors == n_processors...
                ,:);
        end

        function settings = experiment_catalog(obj)
            all_settings = obj.data(:,["n_rows","n_cols","n_refs","solver","n_processors"]);
            [settings,~,affiliation] = unique(all_settings,"rows");
            settings.count=histcounts(affiliation,"BinMethod","integers")';
        end
    end
end

