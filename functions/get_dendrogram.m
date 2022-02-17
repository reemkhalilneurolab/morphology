function get_dendrogram(dist_matrix,name)

load(fullfile(pwd, CreateTree.save_directory('swc_vector.mat')));
dist_matrix = squareform(dist_matrix);
[classes,~,idx]=unique(swc_vector(4,:)');
numClusters = size(classes,1);
Z=linkage(dist_matrix,'ward');
T = cluster(Z,'Maxclust',numClusters);
color = Z(end-numClusters+2,3)-eps;

labels = cellstr(swc_vector(4,:));
 for i=1:numel(labels)
  labels(i)= replace( labels(i) , '_' , '-' )  ; 
 end   
  labels = cellfun(@(x) x(1:end-1), labels, 'Uniform', 0);
 
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
fname = fullfile(pwd,'\save');
filename=sprintf('Dendrogram_%s.jpg',name);
saveas(gcf, fullfile(fname, filename))


end