function trk = consolidate_trajectory_IDs(trk)

%changes trajectory IDs so that there is no negative values (or zero). This
%is important because when trajectories are matched to isolated spots, the
%matched sets of spots are given IDs which are the opposite of their
%matched trajectory ID. 
%The present changes ensure there is a strict one-to-one relationship
%between IDs and trajectories

%If indices negatives or equal to zero are found, all IDs are shifted so
%that the lowest ID is 1. Trajectories are then sorted based on index and
%time.

traj_idx_column = 5;

minidx = min(trk(:,traj_idx_column));

%make sure indices are all striclty positive
if minidx <=0
    trk(:,traj_idx_column) = trk(:,traj_idx_column) - minidx +1;
end

%sort indices
trk = sortrows(trk,[traj_idx_column,4]);


end