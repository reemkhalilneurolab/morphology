function get_dendrogram(dist_matrix,name)
% load(fullfile(pwd, CreateTree.save_directory('tortuosity_dist_matrix.mat')));
 load(fullfile(pwd, CreateTree.save_directory('swc_vector.mat')));

%  myfunc = @(x,k) clusterdata(x,'linkage','average','maxclust',k);
%   eva = evalclusters(x,myfunc,'CalinskiHarabasz','klist',[1:6]);
% [rows columns] = size(dist_matrix);
% v = [];
% [rows columns] = size(dist_matrix);
%     for i = 1:rows-1
%                 v = [v dist_matrix(i,i+1:columns)];
%     end
% dist_matrix =[];
% dist_matrix = v;
dist_matrix = squareform(dist_matrix);


% eva = evalclusters(dist_matrix,'linkage','CalinskiHarabasz','KList',[1:3]);
% 
% numClusters=eva.OptimalK;
numClusters =3;
%  numClusters =4;
% DM.get_clusters(dist_matrix);
Z=linkage(dist_matrix,'average');
T = cluster(Z,'Maxclust',numClusters);
% cutoff = median([Z(end-2,3) Z(end-1,3)]);
color = Z(end-numClusters+2,3)-eps;

% scrsz = get(0,'ScreenSize');
% figH1=figure('Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
% gscatter(dist_matrix(:,1),dist_matrix(:,2),T);
% set(gcf,'name',name,'numbertitle','off')
% fname = 'C:\Workspace\Repositories\PARS-TMD\save'; %laptop
% fname = 'D:\GitHub\Morphology_Descriptors\save'; %HPC
% fname = 'C:\Users\rkhalil\Documents\GitHub\PARS-TMD\save';
% fname = fullfile(pwd,'\save');
% filename=sprintf('Cluster_%s.jpg',name);
% saveas(gcf, fullfile(fname, filename));
labels = cellstr(swc_vector(4,:));
% [NumData, labels, AllData]  = xlsread('D:\GitHub\Morphology\save\index.xlsx', 'Sheet1','A1:A29');
% labels = transpose(labels);
 for i=1:numel(labels)
  labels(i)= replace( labels(i) , '_' , '-' )  ; 
 end   
  labels = cellfun(@(x) x(1:end-1), labels, 'Uniform', 0);
  
% for i=1:numel(labels)
%     labels(i) =strcat(labels(i),  num2str(i));
% end
 % strip the last two characters _1
 
%  make the figure full screen
 scrsz = get(0,'ScreenSize');
 figH=figure('Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
 
 [H,T,origpos]  = dendrogram(Z,0, 'Labels', labels, 'Orientation','left','ColorThreshold',color);
  set(H,'LineWidth',2)
ax = gca; % get the axes handle
X = get(ax.Children,'XData'); % Get x values of all lines
Y = get(ax.Children,'YData'); % Get y values of all lines
ax = gca; % get the axes handle
lab = ax.YAxis.TickLabels; % get all the labels
loc = ax.YAxis.TickValues; % get all labels location
[ulab,~,lab_ind] = unique(lab); % find all unique labels
clr = lines(numel(ulab)); % make a color map
templab=lab;

c=0.1;
for k = 1:numel(ulab) % for every type of lable     
    ind = strcmp(ulab{k},lab); % find all instances in lab
    x = repelem(ax.XAxis.Limits(1)-0.01,sum(ind)); % make an x position vector
    % place this lable at the same locations with a distinct color:  
%     origpos=flip(origpos);
    myloc=loc(ind);
    for j=1:numel(myloc) 
        if loc(mod(myloc(j),2)==0)
            lab{myloc(j)}=strcat('-------------',num2str(origpos(myloc(j))));    
        else
            lab{myloc(j)}=strcat('-------',num2str(origpos(myloc(j))));  
        end 
%            lab{myloc(j)}=strcat('----',lab{myloc(j)},'----',num2str(origpos(myloc(j))));  
    end

    c=c+0.02;
    NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*c;
    text(NE(1), NE(2), ulab{k},'Color',clr(k,:), 'VerticalAlignment','top', 'HorizontalAlignment','right')
    text(x,loc(ind),lab(ind),'Color',clr(k,:));
end

ax.YAxis.TickLabels = []; % remove the original labels
% replace the original labels with white space, to keep the axes position:
ax.YAxis.TickLabels = repelem('  ',max(cellfun(@numel,lab)));


set(gcf,'name',name,'numbertitle','off')
% % fname = 'C:\Workspace\Repositories\PARS-TMD\save'; %laptop
% fname = 'D:\GitHub\Morphology_Descriptors\save'; %HPC
% % fname = 'C:\Users\rkhalil\Documents\GitHub\PARS-TMD\save';
fname = fullfile(pwd,'\save');
filename=sprintf('Dendrogram_%s.jpg',name);
saveas(gcf, fullfile(fname, filename))

% NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
% SW = [min(xlim) min(ylim)]+[diff(xlim) diff(ylim)]*0.05;
% NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
% SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
% text(NE(1), NE(2), 'NORTHEAST!', 'VerticalAlignment','top', 'HorizontalAlignment','right')
% text(SW(1), SW(2), 'SOUTHWEST!', 'VerticalAlignment','bottom', 'HorizontalAlignment','left')
% text(NW(1), NW(2), 'NORTHWEST!', 'VerticalAlignment','top', 'HorizontalAlignment','left')
% text(SE(1), SE(2), 'SOUTHEAST!', 'VerticalAlignment','bottom', 'HorizontalAlignment','right')

% for k = 1:numel(Y)
%     if Y{k}(1)==fix(Y{k}(1))
%         line(ax,X{k}(1:2),Y{k}(1:2),'Color',clr(lab_ind(Y{k}(1)),:),...
%             'LineWidth',2);
%     end
%     if Y{k}(3)==fix(Y{k}(3))
%         line(ax,X{k}(3:4),Y{k}(3:4),'Color',clr(lab_ind(Y{k}(3)),:),...
%             'LineWidth',2);
%     end
% end

end