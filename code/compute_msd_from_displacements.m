function msd = compute_msd_from_displacements(disp)

%input: spots positions in colocalized trajectories formatted as follows:
%disp is formatted as follows:
% col1: traj idx
% col2: time interval in frames
% col3: dx in um
% col4: dy in um
% col5: dr in um

%outputs:
% msd is an array that contains the displacements between co-detection events 
%each entry corresponds to a time interval,
%e.g. dispco{10} contains the displacements at intervals of 10 frames
%each dispco entry is an array formatted as follows:
% col1: coloc traj #, 
% col2: value of the time interval (in frames) at which the cdf was computed
% col3: msd in um^2
% col4: std of the square displacement in um^2
% col5: # of datapoints used to compute the msd


msd = [];
traj_idx = unique(disp(:,1));

for i = 1:numel(traj_idx)
    curdisp = disp(disp(:,1)==traj_idx(i),:);
    tvalues = unique(curdisp(:,2));
    for j=1:numel(tvalues)
        curdisp2 = curdisp(curdisp(:,2) == tvalues(j),:);
        msd = [msd;traj_idx(i),tvalues(j),nanmean(curdisp2(:,5).^2),nanstd(curdisp2(:,5).^2),size(curdisp2,1)];
    end
end



end