function plot_colocalization_visual_output(trk1,coloc_trk2,neg_trk2,wintitle)

trajfig = figure;
set(trajfig,'Name',wintitle);
hold;
plot_trajectories_in_figure(trk1,trajfig,'k');

%remove redundant datapoints in colocalized trajectory array
[~,unique_idx,~] = unique(coloc_trk2(:,1:5),'rows');
coloc_trk2 = coloc_trk2(sort(unique_idx),:);
plot_trajectories_in_figure(coloc_trk2,trajfig,'r');
plot_trajectories_in_figure(neg_trk2(:,1:5),trajfig,'b');
h = gca;
set(gca,'DataAspectRatio',[1,1,1]);

%legend box
legend_str(1) = {'Channel 1'};
legend_str(2) = {'\color[rgb]{1 0 0} Colocalized Channel 2'};
legend_str(3) = {'\color[rgb]{0 0 1} Orphans Channel 2'};

xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
xmax = xlim(2);
ymax = ylim(2);
text(xmax,ymax,legend_str,'HorizontalAlignment','right','VerticalAlignment','top',...
    'BackgroundColor','w','EdgeColor','k');
end