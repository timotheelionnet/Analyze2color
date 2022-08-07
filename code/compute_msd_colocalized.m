function [msdco,msdall,msdsingle] = compute_msd_colocalized(traj,pix_size)

%input: spots positions in colocalized trajectories formatted as follows:
%1 and 2 refer to the two channels
% col1: x1,
% col2: y1,
% col3: I1,
% col4: x2,
% col5: y2,
% col6: I2,
% col7: t,
% col8: traj1 #,
% col9: traj 2#,
% col10: colocalized trajectory #

%outputs:
% dispco is a cell array that contains the displacements between co-detection events 
%each entry corresponds to a time interval,
%e.g. dispco{10} contains the displacements at intervals of 10 frames
%each dispco entry is an array formatted as follows:
% col1: coloc traj #, 
% col2: dx1, 
% col3: dy1, 
% col4: dr1, 
% col5: dx2,
% col6: dy2,
% col7: dr2

%dispall is a cell array that contains the displacements between all events
%regardless of colocalization.
%format identical as dispco

%dispsingle is a cell array that contains the displacements between all events
%that belong to a bonafide trajectory in a single channel only.
%for example for the channel 1, I include all spots within a trajectory from channel 1
%regardless of colocalization events.
%However, I exclude the isolated spots from channel 1 (i.e. not part of an identified
%trajectory) that colocalize with channel 2 trajectories
%format identical as dispco

traj_idx = unique(traj(:,10));

msdco = cell('');
msdall = cell('');
msdsingle = cell('');
traj(:,1:2) = traj(:,1:2)*pix_size;
traj(:,4:6) = traj(:,4:6)*pix_size;

for i = 1:numel(traj_idx)
    %select current trajectory
    tmptraj = traj(traj(:,10)==traj_idx(i),:);
    
    %crop trajectories for the various data sets
    
    %data set for msdall: keep all spots separately in each channel
    %regardless of codetections (including the isolated spots)
    tmptraj1a = [tmptraj(:,1:3),tmptraj(:,7)];
    tmptraj2a = tmptraj(:,4:7); 
    
    %dataset for msdsingles: in each channel, keep only the spots that
    %belong to a bonafide trajectory (exclude the spots that were added 
    %only because they colocalize with a trajectory in the other channel).
    %Include all these spots, regardless of colocalization events
    tmptraj1s = tmptraj1a;
    tmptraj1s( tmptraj(:,8)==-1 ,1:3) = NaN;
    
    tmptraj2s = tmptraj2a;
    tmptraj2s( tmptraj(:,9)==-1,1:3) = NaN; 
    
    %dataset for msdco: include all co-detected spots.
    tmptraj1c = tmptraj1a;
    tmptraj2c = tmptraj2a;
    tmptraj1c( ~logical( (~isnan(tmptraj(:,1))).*(~isnan(tmptraj(:,4)))),1:3) = NaN;
    tmptraj2c( ~logical( (~isnan(tmptraj(:,1))).*(~isnan(tmptraj(:,4)))),1:3) = NaN;
    
    msdco{traj_idx(i),1}(:,1) = compute_msd(tmptraj1c);
    msdco{traj_idx(i),1}(:,2) = compute_msd(tmptraj2c);
    msdall{traj_idx(i),1}(:,1) = compute_msd(tmptraj1a);
    msdall{traj_idx(i),1}(:,2) = compute_msd(tmptraj2a);
    msdsingle{traj_idx(i),1}(:,1) = compute_msd(tmptraj1s);
    msdsingle{traj_idx(i),1}(:,2) = compute_msd(tmptraj2s);
    
    
end



end