function print_step_functions()
load(fullfile(pwd, CreateTree.save_directory('vect_energy.mat')));
load(fullfile(pwd, CreateTree.save_directory('vect_flux.mat')));
load(fullfile(pwd, CreateTree.save_directory('vect_leaf_index.mat')));
load(fullfile(pwd, CreateTree.save_directory('vect_taper_rate.mat')));
load(fullfile(pwd, CreateTree.save_directory('vect_wiring.mat')));
load(fullfile(pwd, CreateTree.save_directory('vect_tortuosity.mat')));
load(fullfile(pwd, CreateTree.save_directory('vect_branching_pattern.mat')));

  vects = {vect_energy vect_flux vect_leaf_index vect_branching_pattern vect_wiring vect_tortuosity  vect_taper_rate };

for i=1:length(vects)
    if i == 1     
    descr_name = 'energy';
    elseif i==2
    descr_name = 'flux';  
    elseif i==3
    descr_name = 'leaf_index';
    elseif i==4
    descr_name = 'branching_pattern';
    elseif i==5
    descr_name = 'wiring';
    elseif i==6
    descr_name = 'tortuosity';
    elseif i==7
    descr_name = 'taper_rate';
 
    end 
    
    
PrintTree.print_step_function(vects{i},descr_name);
    
end 



end 

