classdef Result < handle

    properties
        data
    end
    
    methods
        function obj = Result(path)
            obj.data = readtable(path);
            
        end
        
        function partial_data = find_experiments(obj,n_rows,n_cols,n_refs,solver,n_processors,varargin)
            mask = (obj.data.n_rows==n_rows) &...
                (obj.data.n_cols==n_cols) &...
                (obj.data.n_refs==n_refs)&...
                (obj.data.n_processors==n_processors)&...
                ismember(obj.data.solver,solver);
            partial_data=obj.data(mask,:);
        end

        function settings = experiment_catalog(obj)
            all_settings = obj.data(:,["n_rows","n_cols","n_refs","solver","n_processors"]);
            [settings,~,affiliation] = unique(all_settings,"rows");
            settings.count=histcounts(affiliation,"BinMethod","integers")';
        end
    end
end

