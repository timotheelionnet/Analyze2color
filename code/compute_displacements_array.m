function disp = compute_displacements_array(traj,pix_size)

%input: spots positions in trajectory array formatted as follows:
% col1: x1 (pix),
% col2: y1 (pix),
% col3: I1,
% col4: t (frames)
% col5: traj index,


%pix size: size of the pixels in um


%output:
% disp is an array that contains the displacements for all time intervals
% col1: traj #,
% col2: time interval in frames
% col2: dx1, 
% col3: dy1, 
% col4: dr1

%redundancy check: remove any datapoints that belong to the same trajectory
%and have the same time coordinate, (plus any NaN values)
traj = traj(~isnan(traj(:,1)),:);
[~,idx,~] = unique([traj(:,5),traj(:,4)],'rows'); %each time point from a given trajectory should be entered only once
traj = traj(idx,:);

traj_idx = unique(traj(:,5));
disp = [];
for i = 1:numel(traj_idx)
    %select current trajectory
    tmptraj = traj(traj(:,5)==traj_idx(i),:);
    
    %compute displacements
    disp = [disp;get_displacement_from_current_trajectory_array(tmptraj,pix_size)];
    
end

end

