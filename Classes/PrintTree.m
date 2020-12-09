classdef PrintTree
 
    methods(Static)
        
        function print_one_tree(tree, vect_of_nodes, com )
            %--------------------------------------
            %prints the tree in the order specified by the mat structure
            %tested with swc file and worked like a charm
            %-----------------------------------
            predicate_branch = tree.is_branching_point();
            predicate_leaf = tree.is_leaf();
            predicate=predicate_branch + predicate_leaf;
            all_path=tree.get_paths_from_leafs_to_root();
            mat=all_path{3};
            
            lnwdth =1
            matrix=[];
                for i=1:length(vect_of_nodes)
                    matrix =[matrix ;vect_of_nodes(i).id 1 vect_of_nodes(i).coordinates(1) vect_of_nodes(i).coordinates(2) vect_of_nodes(i).coordinates(3) 1 vect_of_nodes(i).parent];
                end 
            assignin('base','matrix',matrix);
            save('MyMatrix.txt', 'matrix', '-ascii', '-double', '-tabs')

                for i=1:length(vect_of_nodes)-1
                    current_node_parent=vect_of_nodes(i).parent;
                    current_node_x=vect_of_nodes(i).coordinates(1);
                    current_node_y=vect_of_nodes(i).coordinates(2);
                    current_node_z=vect_of_nodes(i).coordinates(3);
                    
                    if current_node_parent==-1
                    current_node_parent=1;
                    end

                    for j=1:length(vect_of_nodes)
                        if vect_of_nodes(j).id==current_node_parent
                            parent_x = vect_of_nodes(j).coordinates(1) ;
                            parent_y = vect_of_nodes(j).coordinates(2) ;
                            parent_z = vect_of_nodes(j).coordinates(3);
                            break;
                        end                          
                    end 
                    
                    a=0;
                    b=0;
                    c=0;
                    lnwdth =0.5;
                    mrkr='none';
%                     mrkr='.';

                    if current_node_x == 0 && current_node_y == 0 && current_node_z == 0 
                        a=1;
                        b=0;
                        c=0;
                        lnwdth =5;
                        mrkr='+';
                    end 

                    line([parent_x,current_node_x] ,[parent_y,current_node_y] ,[parent_z,current_node_z],'Marker',mrkr,'LineWidth',lnwdth,'Color',[a b c]);
%                     if (i<=983 && i >= 979) || i ==936
%                         text(current_node_x,current_node_y,num2str(i));
%                     end

%                     for j =1:length(mat)
%                         if mat(j) ==i
%                          t=text(current_node_x,current_node_y,num2str(i));
%                          t(1).Color = 'red';
%                          t(1).FontSize = 7;
%                         end 
%                             
%                     end 
%                     if predicate(i) ==1
%                        t=text(current_node_x,current_node_y,num2str(i));
%                        t(1).Color = 'red';
%                        t(1).FontSize = 7;
%                     end
    
                end





            %     figure;

            hold on;
            x=com(1);
            y=com(2);
            z=com(3);

            scatter3(com(1),com(2),com(3),'MarkerEdgeColor','k', 'MarkerFaceColor',[0 .75 .75])
            hold off;

        end


        function print_all_trees()   
            close all;
            load(fullfile(pwd, CreateTree.save_directory('trees_vector.mat')));
            load(fullfile(pwd, CreateTree.save_directory('swc_vector.mat')));
            vect=[1:length(trees_vector)];
            matrix_all_48=[];
                for i=1:length(vect)
                    figure(vect(i));
                    mytree1=trees_vector(vect(i));
                     mytree=normalize_tree(mytree1);
                    PrintTree.print_one_tree(mytree, mytree.get_nodesVector(),mytree.get_center_of_mass());
                    num_nodes=mytree.get_num_nodes();
                    num_bif=mytree.get_num_bif_nodes();
                    num_leaf=mytree.get_num_leaf_nodes();
                    name=swc_vector{4,i};
                    name = name(1:end-4);

                    %                 [max_eigen_value average_tortuosity] =Spectral_descriptors.tortuosity_vector(mytree);
                    %                 v= [i num_nodes num_bif num_leaf max_eigen_value average_tortuosity ];
                    v= [i num_nodes num_bif num_leaf];
                    matrix_all_48=[matrix_all_48 ;v];


                    %                 load('C:\Users\rkhalil\Documents\GitHub\PARS-TMD\test_cases\spectral_radii.mat');
                    %                 load('C:\Users\rkhalil\Documents\GitHub\PARS-TMD\test_cases\energetic_volume.mat');
%                     fname = 'D:\GitHub\Morphology_Descriptors\save';
                    fname = fullfile(pwd,'\save');
                    filename=sprintf('%d_%s_%d_%d_%d_.jpg',i,name,num_nodes,num_bif,num_leaf);
%                     saveas(gca, fullfile(fname, filename));
                    %Get Current Figure (GCF) & Set image size before saving image
                    width = 9;  % cm 
                    height = 7; % cm
                    set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])
                    %Set the resolution of 1000dpi and save the plot in TIFF format 
%                     print -dpng -r600 fullfile(fname, filename)
                    print('-dtiff','-r600',fullfile(fname, filename))
                   
                    close();
                end
            assignin('base','matrix_all_48',matrix_all_48);
        end 

        function print_step_function(cellArray,descr_name)
            close all;
            new_dir = strcat(pwd,'\save\',descr_name);
            mkdir(new_dir);           
%             dir = fullfile(pwd,new_dir);
            fname = fullfile(new_dir);
            
            load(fullfile(pwd, CreateTree.save_directory('trees_vector.mat')));
            load(fullfile(pwd, CreateTree.save_directory('swc_vector.mat')));
            vect=[1:length(trees_vector)];
                for i=1:length(vect)
                    figure(vect(i));
                    mytree=trees_vector(vect(i));
                    num_nodes=mytree.get_num_nodes();
                    num_bif=mytree.get_num_bif_nodes();
                    num_leaf=mytree.get_num_leaf_nodes();
                    name=swc_vector{4,i};
                    name = name(1:end-4);
                    x=cellArray{i}(:,1);
                    y=cellArray{i}(:,2);
                    stairs(x,y,'LineWidth',1, 'Color' ,'r','Marker','o','MarkerFaceColor','c') 
                    h = stairs(x,y) ;
                    h(1).LineWidth = 1;
                    h(1).Color = 'r';
                    h(1).Marker = 'o';
                    h(1).MarkerSize = 4;
                    h(1).MarkerFaceColor = 'c';
                    x0=10;
                    y0=10;
                    width=10;
                    height=10;
                    set(gcf,'units','inches','position',[x0,y0,width,height])
                    pos=get(gca,'position');  % retrieve the current values
                      pos(3)=1*pos(3);        % try reducing width 10%
                     set(gca,'position',pos);  % write the new values
   
%                     fname = fullfile(pwd,'\save');
                    filename=sprintf('%d_step_%s_%s_%d_%d_%d_.jpg',i,descr_name,name,num_nodes,num_bif,num_leaf);
                    saveas(gca, fullfile(fname, filename));
                    close();
                end
%             assignin('base','matrix_all_48',matrix_all_48);
%             
%         stairs(x,y,'LineWidth',2, 'Color' ,'r','Marker','d','MarkerFaceColor','c') 
%         fname = fullfile(pwd,'\save');
%         filename=sprintf('%d_%s_%d_%d_%d_.jpg',i,name,num_nodes,num_bif,num_leaf);
%         saveas(gca, fullfile(fname, filename));
%         close();
            
        end 
    
   
    end
end

