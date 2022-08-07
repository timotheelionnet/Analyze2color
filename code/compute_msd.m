function [dt,msd,msd_std,ndata] = compute_msd(traj)
%traj = x,y,I,t

traj = traj(~isnan(traj(:,1)),:);

msd = zeros((size(traj,1)-1),1);
msd_std = zeros((size(traj,1)-1),1);
ndata = zeros((size(traj,1)-1),1);

dt = repmat(traj(:,4),1,size(traj,1)) - repmat(traj(:,4)',size(traj,1),1);
tvalues = unique(dt(dt>0));

dx = repmat(traj(:,1),1,size(traj,1)) - repmat(traj(:,1)',size(traj,1),1);
dy = repmat(traj(:,2),1,size(traj,1)) - repmat(traj(:,2)',size(traj,1),1);
dr2 = dx.^2 + dy.^2;

for i=1:numel(tvalues) 
    ndata(i) = sum(dt(:)==tvalues(i));
    dr_tmp = dr2(dt==tvalues(i));
    msd(i) = nanmean(dr_tmp);
    msd_std(i) = nanstd(dr_tmp);
    
end
dt = tvalues;
end