classdef optimization_functions 
   
methods(Static)
    function opt
        tStart = tic; 
        %In this file we will read distance matrices comong from various descriptor
        %and build a combinations of them to get one distance matrix to rule them
        %all.
        load_switch =1;
%         load swc vector
        if load_switch ==1
            load(fullfile(pwd, CreateTree.save_directory('swc_vector.mat')));
            %First we read the distance matrices:
            load(fullfile(pwd, CreateTree.save_directory('spectral_dist_matrix.mat')));
%             load(fullfile(pwd, CreateTree.save_directory('taper_rate_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('tortuosity_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('branching_dist_matrix.mat')));
%             load(fullfile(pwd, CreateTree.save_directory('energy_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('flux_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('south_flux_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('north_flux_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('lateral_flux_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('eigen_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('tmd_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('branching_angle_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('leaf_index_dist_matrix.mat')));
        end 
        spectral   = spectral_dist_matrix;
%         taper      = taper_rate_dist_matrix;
        tortuosity = tortuosity_dist_matrix;        
        branching  = branching_dist_matrix;
%         energy     = energy_dist_matrix;
        flux       = flux_dist_matrix;
        SFlux      = south_flux_dist_matrix;
        NFlux      = north_flux_dist_matrix;
        LFlux      = lateral_flux_dist_matrix;
        eigen      = eigen_dist_matrix;
        tmd      = tmd_dist_matrix;
        branching_angle = branching_angle_dist_matrix;
        leaf_index =leaf_index_dist_matrix
        
        %optionally we can normalize distnace matrices by dividing each of them by
        %the maximal value therein. Then they will be in the same scale for the
        %further processing:
        spectral   = spectral./max(max(spectral));
%         taper      = taper./max(max(taper));  
        tortuosity = tortuosity./max(max(tortuosity));
        branching = branching./max(max(branching));
%         energy     = energy./max(max(energy));
        flux       = flux./max(max(flux));        
        SFlux = SFlux./max(max(SFlux));
        NFlux = NFlux./max(max(NFlux));
        LFlux = LFlux./max(max(LFlux));
        eigen = eigen./max(max(eigen)); 
        tmd = tmd./max(max(tmd)); 
        branching_angle = branching_angle_dist_matrix;
        leaf_index =leaf_index_dist_matrix;

        %Now we create weights with which they should be used. 
%           dist_mat = {spectral taper tortuosity branching energy flux SFlux NFlux LFlux eigen tmd branching_angle leaf_index};
          dist_mat = {spectral tortuosity branching flux SFlux NFlux LFlux eigen tmd branching_angle leaf_index};
%         
        [classes,~,idx]=unique(swc_vector(4,:)');
        
        for j=1:length(dist_mat) % descriptors 
                mat=[];
                for i=1:length(classes) %neurons classes
                    % get all the indices of the class i 
                    k = find(idx==i);
                    %starting index
                    s = min(k);
                    %ending index
                    e = max(k); 
                    mat = [ mat; [i s e]];
                    max_distance_within_i_class(j,i) = max( max(dist_mat{j}(s:e,s:e)));
                end 
                % max of all the internal distances across class for each
                % descriptor
                internal(1,j) = sum(max_distance_within_i_class(j,:));                
                
                c = nchoosek(mat(:,1),2);
                [r,~]=size(c);
                for i=1:r
                    % start and end of pair 1
                    s1=mat(c(i,1),2);
                    e1=mat(c(i,1),3);
                    % start and end of pair 2
                    s2=mat(c(i,2),2);
                    e2=mat(c(i,2),3);
                    min_distance_between_i_j_classes(j,i) =min(min(dist_mat{j}(s1:e1,s2:e2)));        
                end    
                 external(1,j) = min(min_distance_between_i_j_classes(j,:));
        end
       
%            optimization_functions.optimize_fmincon(dist_mat , internal, external);
%          optimization_functions.optimize_linprog(dist_mat , internal, external);
         optimization_functions.optimize_python_gekko(dist_mat , internal, external);
%         
    end 
    
      % optimization using fmincon
    function optimize_fmincon(dist_mat , internal, external)
        
    f = @(x) -1*( x(1)*external(1) + x(2)*external(2) + x(3)*external(3) + x(4)*external(4) +  x(5)*external(5) +  x(6)*external(6) + x(7)*external(7) + x(8)*external(8) + x(9)*external(9) + x(10)*external(10) + x(11)*external(11)+x(12)*external(12)+x(13)*external(13)) ...
    / ( x(1)*internal(1) + x(2)*internal(2) + x(3)*internal(3) + x(4)*internal(4) +  x(5)*internal(5) +  x(6)*internal(6) + x(7)*internal(7) + x(8)*internal(8) + x(9)*internal(9) + x(10)*internal(10) + x(11)*internal(11)+x(12)*internal(12)+x(13)*internal(13));
              
    % f =-1.*[m1/n1 m2/n2 m3/n3 m4/n4 m5/n5 m6/n6 m7/n7 m8/n8 m9/n9 m10/n10 m11/n11];
    mat=[];
    for i=1:100
        x0 = randi([0 9], [1, 11]);
        % s=0;
        % s1=3;
        % x0 = [s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11];
        % lb = [low,low,low,low,low,low,low,low,low,low,low];
        % ub = [high,high,high,high,high,high,high,high,high,high,high];

        g =[-1];
            while any(g<0) 
                lb= randi([0 9], [1, 11]);
                ub= randi([0 9], [1, 11]);
                g= ub-lb;
            end 
        % b = zeros(1,11);
        b= [ 5 5 5 5 5 5 5 5 5 5 5 5 5];
        nonlcon = [];
        options = optimoptions('fmincon','MaxFunctionEvaluations', 100000,'MaxIterations', 10000, 'MaxFunctionEvaluations',10000000 );

        % x= fmincon(f,x0,[],[],[],[],lb,ub,nonlcon, options);
        x= fmincon(f,x0,[],[],[],[],lb,ub);
        disp(x);
        disp(['Final Objective: ' num2str(f(x))])
        mat(i,:) = [f(x) x];        
    end
    a= sortrows(mat, 1);
    x= a(1,2:end);
%     a = int64(mat);
%     [r,c] =size(a);
%     x=[];
%     for k=2:c    
%         x(1,k-1)=mode(a(:,k));   
%     end 
    combination_of_distances.combine_matrices_with_weights(x,dist_mat);      
    end 
    
    % optimization using linprog 
    function optimize_linprog(dist_mat , internal, external)
        f =-1.*[external(1)/internal(1) external(2)/internal(2) external(3)/internal(3)...
                external(4)/internal(4) external(5)/internal(5) external(6)/internal(6)...
                external(7)/internal(7) external(8)/internal(8) external(9)/internal(9)...
                external(10)/internal(10) external(11)/internal(11) external(12)/internal(12) external(13)/internal(13) ];
        % b = zeros(1,11);
        b= [ 5 5 5 5 5 5 5 5 5 5 5];
        x= linprog(f,eye(11),b);
        combination_of_distances.combine_matrices_with_weights(x,dist_mat);
        disp(x);
        disp(['Final Objective: ' num2str(f(x))])
    end 
    
    function optimize_python_gekko(dist_mat , internal, external)
        % Initialize model
        m = py.gekko.GEKKO();
        % Initialize Variables
        x1 = m.Var(pyargs('value',1,'lb',1,'ub',9));
        x2 = m.Var(pyargs('value',5,'lb',1,'ub',9));
        x3 = m.Var(pyargs('value',5,'lb',1,'ub',9));
        x4 = m.Var(pyargs('value',1,'lb',1,'ub',9));
        x5 = m.Var(pyargs('value',1,'lb',1,'ub',9));
        x6 = m.Var(pyargs('value',5,'lb',1,'ub',9));
        x7 = m.Var(pyargs('value',5,'lb',1,'ub',9));
        x8 = m.Var(pyargs('value',1,'lb',1,'ub',9));
        x9 = m.Var(pyargs('value',1,'lb',1,'ub',9));
        x10 = m.Var(pyargs('value',5,'lb',1,'ub',9));
        x11 = m.Var(pyargs('value',5,'lb',1,'ub',9));
%         x12 = m.Var(pyargs('value',5,'lb',1,'ub',9));
%         x13 = m.Var(pyargs('value',5,'lb',1,'ub',9));
        
        m.Obj ( -1*( x1*external(1) + x2*external(2) + x3*external(3) + x4*external(4) +  x5*external(5) +  x6*external(6) + x7*external(7) + x8*external(8) + x9*external(9) + x10*external(10) + x11*external(11)) ...
        / ( x1*internal(1) + x2*internal(2) + x3*internal(3) + x4*internal(4) +  x5*internal(5) +  x6*internal(6) + x7*internal(7) + x8*internal(8) + x9*internal(9) + x10*internal(10) + x11*internal(11)) );

        
%         m.Obj ( -1*( x1*external(1) + x2*external(2) + x3*external(3) + x4*external(4) +  x5*external(5) +  x6*external(6) + x7*external(7) + x8*external(8) + x9*external(9) + x10*external(10) + x11*external(11)+ x12*external(12) + x13*external(13)) ...
%         / ( x1*internal(1) + x2*internal(2) + x3*internal(3) + x4*internal(4) +  x5*internal(5) +  x6*internal(6) + x7*internal(7) + x8*internal(8) + x9*internal(9) + x10*internal(10) + x11*internal(11)+ x12*internal(12) + x13*internal(13)) );

         m.solve();
        % Extract values from Python lists using curly brackets
        disp(['x1: ' num2str(x1.VALUE{1})]);
        disp(['x2: ' num2str(x2.VALUE{1})]);
        disp(['x3: ' num2str(x3.VALUE{1})]);
        disp(['x4: ' num2str(x4.VALUE{1})]);
        disp(['x5: ' num2str(x5.VALUE{1})]);
        disp(['x6: ' num2str(x6.VALUE{1})]);
        disp(['x7: ' num2str(x7.VALUE{1})]);
        disp(['x8: ' num2str(x8.VALUE{1})]);
        disp(['x9: ' num2str(x9.VALUE{1})]);
        disp(['x10: ' num2str(x10.VALUE{1})]);
        disp(['x11: ' num2str(x11.VALUE{1})]);
%         disp(['x12: ' num2str(x12.VALUE{1})]);
%         disp(['x13: ' num2str(x13.VALUE{1})]);

%         x= [ x1.VALUE{1} x2.VALUE{1} x3.VALUE{1} x4.VALUE{1} x5.VALUE{1} ...
%             x6.VALUE{1} x7.VALUE{1} x8.VALUE{1} x9.VALUE{1} x10.VALUE{1} x11.VALUE{1} x12.VALUE{1} x13.VALUE{1}];
           x= [ x1.VALUE{1} x2.VALUE{1} x3.VALUE{1} x4.VALUE{1} x5.VALUE{1} ...
            x6.VALUE{1} x7.VALUE{1} x8.VALUE{1} x9.VALUE{1} x10.VALUE{1} x11.VALUE{1}];
   
            combination_of_distances.combine_matrices_with_weights(x,dist_mat);
    end 
    
    
    
end
end 