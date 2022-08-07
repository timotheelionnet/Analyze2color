function plot_colocalization_stats2(coloc_trk1,coloc_trk2,matrix,trk1,trk2)
traj_idx_col_num = 5;
matches = match_trajectories4_dia(trk1,trk2,1,20);

%plot distances between spots and between trajectories
%'All' trajectories
dist = [];
avgdist = [];
for i=1:size(matches,1) 
    avgdist = [avgdist,matches(i).avgdist];
    curdist = reshape(matches(i).distances,numel(matches(i).distances),1);
    if ~isempty(curdist)
        curdist = reshape(curdist,numel(curdist),1);
        curdist = curdist(curdist~=0);
        dist = [dist;curdist];
    end    
end

%colocalized trajectories only
sdist = [];
savgdist = [];
traj_idx = unique(coloc_trk1(:,7));
coloc_dist = sqrt(sum((coloc_trk1(:,1:2) - coloc_trk2(:,1:2)).^2,2));
for i=1:numel(traj_idx)
    savgdist = [savgdist,nanmean(coloc_dist(coloc_trk1(:,7)==traj_idx(i)))];
    cursdist = coloc_dist(coloc_trk1(:,7)==traj_idx(i));
    if ~isempty(cursdist)
        cursdist = reshape(cursdist,numel(cursdist),1);
        cursdist = cursdist(cursdist~=0);
        cursdist = cursdist(~isnan(cursdist));
        sdist = [sdist;cursdist];
    end    
end

[n,x] = hist(dist,0:0.4:max(dist)+1);
[ns,xs] = hist(sdist,0:0.4:max(sdist)+1);
figure;
set(gcf,'Name','Spot Distances Distributions');
axes('Position',[0.1,0.2,0.35,0.6]);
bar(x,n,'BarWidth',1,'EdgeColor',[0.85,0.85,0.85],'FaceColor',[0.85,0.85,0.85]);
hold;
plot(xs,ns,'Color',[1,0,0],'LineWidth',2);
set(get(gca,'XLabel'),'String','Distance Between Spots (pix)');
set(get(gca,'YLabel'),'String','Frequency');
set(gca,'YLim',[0,max(n)*1.3]);
set(gca,'XLim',[0,max(x)]);
legend('All Trajectories Within 20 pix','Colocalized');

[n,x] = hist(avgdist,0:0.4:max(avgdist)+1);
[ns,xs] = hist(savgdist,0:0.4:max(savgdist)+1);
axes('Position',[0.6,0.2,0.35,0.6]);
bar(x,n,'BarWidth',1,'EdgeColor',[0.85,0.85,0.85],'FaceColor',[0.85,0.85,0.85]);
hold;
plot(xs,ns,'Color',[1,0,0],'LineWidth',2);
set(get(gca,'XLabel'),'String','Avg Distance Between Trajectories (pix)');
set(get(gca,'YLabel'),'String','Frequency');
set(gca,'YLim',[0,max(n)*1.3]);
set(gca,'XLim',[0,max(x)]);
legend('All Trajectories Within 20 pix','Colocalized');


%plot trajectory length distributions
traj1_idx = unique(trk1(:,traj_idx_col_num));
ntraj1 = numel(traj1_idx);
length1 = zeros(ntraj1,1);
for i=1:ntraj1
    curtraj_time = trk1(trk1(:,traj_idx_col_num)==traj1_idx(i),4);
    length1(i) = max(curtraj_time) - min(curtraj_time) +1;
end

lengthco = [];
cotraj1_idx = unique(coloc_trk1(:,traj_idx_col_num));
ncotraj1 = numel(cotraj1_idx);
lengthco = zeros(ncotraj1,1);
for i=1:ncotraj1
    curtraj_time = coloc_trk1(coloc_trk1(:,traj_idx_col_num)==cotraj1_idx(i),4);
    lengthco(i) = max(curtraj_time) - min(curtraj_time) +1;
end

traj2_idx = unique(trk2(:,traj_idx_col_num));
ntraj2 = numel(traj2_idx);
length2 = zeros(ntraj2,1);
for i=1:ntraj2
    curtraj_time = trk2(trk2(:,traj_idx_col_num)==traj2_idx(i),4);
    length2(i) = max(curtraj_time) - min(curtraj_time) +1;
end

[nlength1,xlength1] = hist(length1,0:1:max(length1)+1);
[nlength2,xlength2] = hist(length2,0:1:max(length2)+1);
[nco,xco] = hist(lengthco,0:1:max(lengthco)+1);

figure; hold;
set(gcf,'Name','Trajectory Duration Distributions');
plot(xlength1,nlength1,'Color',[0,0,0],'LineWidth',2);
plot(xlength2,nlength2,'Color',[0,0,1],'LineWidth',2);
plot(xco,nco,'Color',[1,0,0],'LineWidth',2);
set(get(gca,'XLabel'),'String','Trajectory Duration (frames)');
set(get(gca,'YLabel'),'String','Frequency');
legend('All Channel 1','All Channel 2','Colocalization events');

%plot number of colocalized trajectories distribution
nco1 = [];
for i=1:ntraj1
    nco1 = [nco1;numel(find(matrix(:,1)==traj1_idx(i)))];
end
nco2 = [];
for i=1:ntraj2
    nco2 = [nco2;numel(find(matrix(:,2)==traj2_idx(i)))];
end
[nnco1,xnco1] = hist(nco1,0:1:max([nco1;nco2]));
[nnco2,xnco2] = hist(nco2,0:1:max([nco1;nco2]));
figure; 
set(gcf,'Name','Distributions of Colocalization Events Per Trajectory');
axes('Position',[0.1,0.2,0.35,0.6]);
bar(xnco1,nnco1,'BarWidth',1,'EdgeColor',[0,0,0],'FaceColor',[0,0,0]);
set(get(gca,'XLabel'),'String','Number of colocalization events per trajectory');
set(get(gca,'YLabel'),'String','Frequency');
set(gca,'YScale','log');
legend('Channel 1');
axes('Position',[0.6,0.2,0.35,0.6]);
bar(xnco2,nnco2,'BarWidth',1,'EdgeColor',[0,0,1],'FaceColor',[0,0,1]);
set(get(gca,'XLabel'),'String','Number of colocalization events per trajectory');
set(get(gca,'YLabel'),'String','Frequency');
set(gca,'YScale','log');
legend('Channel 2');

end
