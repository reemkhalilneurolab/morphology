function main()
    disp(['PROCESSING STARTED @@@@:   ' , datestr(now,'HH:MM:SS.FFF')]);
    flag = config();
    load_switch = 1;
        if flag == 1 
            if load_switch == 1
                %load swc files into cell arrays
                swc_vector   = CreateTree.import_swc_files();
                %covert swc files into trees
                trees_vector = CreateTree.build_trees(swc_vector);
            end

        %---------------Descriptors-------------------
         tortuosity.get_tortuosity();
         flux.get_flux();
         wiring.get_wiring();
         branchingPattern.get_branching_pattern();       
         leafIndex.get_leaf_index();             
         energy.get_energy();
         Tmd_sholl.get_tmd_dendrogram();
         Tmd_classical.get_tmd_dendrogram(); 
          taperRate.get_taper_rate();
%         ----------------------------------------------- 
        disp(['PROCESSING ENDED @@@@:   ' , datestr(now,'HH:MM:SS.FFF')]);
        print_step_functions();
        PrintTree.print_all_trees(); 
        combination_of_distances.combine_distances(1);
        cluster_classes.cluster();
        metric.metriclearn();
        metric.vectorize_unsupervised();

        end 

end

