classdef CreateTree    

    methods(Static)
        
        function directory = save_directory(name)
             directory = fullfile('\save',name);

        end
      
        function swc_vector=import_swc_files() 
            [swc_files,Paths, set_name] = CreateTree.recdir('.\data', '*.swc'); %laptop
            swc_vector = CreateTree.swc_file_to_vector( swc_files,Paths, set_name );         
        end
        
        function [files,paths,set_name] = recdir(mydir,pattern)
        % Recursively list filenames in all subdirectories. Pattern uses DIR syntax.
        % Example: [N,P] = recdir(pwd, '*.txt')
        files = {};
        paths = {};
        set_name = {};
        flag  = 0;
        nestdir(mydir,'', flag)
            function nestdir(P,sname,flag)
                
                S = dir(fullfile(P,pattern));
                files = [files,{S.name}];
                paths = [paths,repmat({P},1,numel(S))];
                set_name = [set_name,repmat({sname},1,numel(S))];
                S = dir(fullfile(P,'*'));
                C = setdiff({S([S.isdir]).name},{'..','.'});
                for k = 1:numel(C) 
                nestdir(fullfile(P,C{k}),C{k},1)
                end
%                 for k = 1:numel(C) 
%                     if flag ==0
%                         if contains(C{k},'_1') 
%                             nestdir(fullfile(P,C{k}),C{k},1)
%                         end
%                     else
%                         nestdir(fullfile(P,C{k}),C{k}, 1)
%                     end
%                 end
            end
        end      
        
  
        function swc_vector = swc_file_to_vector(swcFiles,Paths, set_name)
            swc_vector = cell( 1 , length(swcFiles) );   
            % Load the files into a cell array
            for i=1:length(swcFiles)
                %Removing extra characters from swc files downloaded from neuromorpho.org 
                %swc files comments start with # character.
                filecontent = fileread(swcFiles{i});
                newcontent = regexprep(filecontent, '^[#:].*$', '', 'lineanchors', 'dotexceptnewline');
                % 7 columns from the swc file
                numcols =7;
                %   %f converts floating-point values in the swc file to text in all 7 columns
                fmt = repmat('%f', 1, numcols);    
                % data is a cell array 
                data= textscan(newcontent, fmt, 'CollectOutput', 1);
                %mtx is an array 
                %data{1} will have all the numeric columns
                mtx= [data{1}];
    %             swc_vector{i}=mtx;
                swc_vector{1,i}=mtx;
                swc_vector{2,i}=swcFiles{i};
                swc_vector{3,i}=Paths{i};
                swc_vector{4,i}=set_name{i};
                
            end
        end
        
        function swcFiles = files_in_this_folder( path )
            % read SWC files from a directory
            swcFiles = dir(path); 
        end        
        

        function trees_vector=build_trees(swc_vector)
            [trees_vector swc_vector]=  CreateTree.swc_vector_to_tree(swc_vector,'tree');
            assignin('base','trees_vector',trees_vector);  
            assignin('base','swc_vector',swc_vector);  
            save(fullfile( pwd , CreateTree.save_directory('trees_vector')),'trees_vector' );
            save(fullfile( pwd , CreateTree.save_directory('swc_vector')),'swc_vector' );
        end 
        

        function [trees_vector swc_vector] =  swc_vector_to_tree(swc_vector,tree_type)
            trees_vector =[];            
            switch tree_type          
               
                
                case 'tree'

                    for i=1:size(swc_vector,2)
                        disp(i);
                        tree=Tree();        
                        mtx= swc_vector{1,i};
                        tree_mtx = CreateTree.remove_rows(mtx,2);
                        [r c]=size(tree_mtx);
                            for j=1:r
                               tree.addNode( Node(tree_mtx(j,:)) ); 
                            end 
                        trees_vector =[trees_vector tree];  
                        swc_vector{5,i}=tree;
                        swc_vector{6,i}='whole tree';
                    end 

                case 'basal'
                     for i=1:size(swc_vector,2)
                        disp(i);
                        tree=Tree();        
                        tree_mtx= swc_vector{1,i};
                        [r c]=size(tree_mtx);
                            for j=1:r
                                if (tree_mtx(j,2)~=2 )
                                if (tree_mtx(j,2)==1 || tree_mtx(j,2)==3) 
                                  tree.addNode( Node(tree_mtx(j,:)) );  
                                end 
                                end
                            end 
                        trees_vector =[trees_vector tree];
                     end 
                    
                case 'apical'
                     for i=1:size(swc_vector,2)
                        disp(i);
                        tree=Tree();        
                        tree_mtx= swc_vector{1,i};
                        [r c]=size(tree_mtx);
                            for j=1:r
                                if (tree_mtx(j,2)~=2 )
                                if (tree_mtx(j,2)==1 || tree_mtx(j,2)==4)
                                    tree.addNode( Node(tree_mtx(j,:)) ); 
                                end
                                end
                            end 
                        trees_vector =[trees_vector tree]; 
                    end 
                otherwise
                    disp('Not a valid tree type')
            end               
        end 
          
       
       function mtx = remove_rows(tree_mtx,col_val)
           mtx = tree_mtx(tree_mtx(:,2)~=col_val,:); % delete all rows for a specific structure ..example axons
           mtx(1,8)=-1;
           for i=2:length(mtx)              
             mtx(i,8)=find(mtx(:,1)==mtx(i,7));           
           end 
           
           mtx(1,9)=1;
           for i=2:length(mtx)              
             mtx(i,9)=find(mtx(:,1)==mtx(i,1));           
           end 
           
            mtx(:,1) = mtx(:,9);% replace column one with the new values from column 9
            mtx(:,7) = mtx(:,8);
            mtx(:,8:9) = [];
           
           
       end 
            
    end
end

