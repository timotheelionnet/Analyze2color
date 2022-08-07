function displ = get_displacement_from_current_trajectory_array(traj,pix_size)
    
%traj is formatted as follows:
% col1: x1,
% col2: y1,
% col3: I1,
% col4: t,
% col5: traj index,
%ASSUMES A SINGLE TRAJECTORY!
    
%displ is formatted as follows:
% col1: traj idx
% col2: time interval in frames
% col3: dx
% col4: dy
% col5: dr

if isempty(traj)
    displ = [];
    return
end

npts = size(traj,1);
traj_idx = traj(1,5);

%time values
dt = repmat(traj(:,4),1,npts) - repmat(traj(:,4)',npts,1);
tvalues = unique(dt(dt>0));

%extract the dx, dy and dr for all time points
dx = repmat(traj(:,1),1,npts) - repmat(traj(:,1)',npts,1);
dy = repmat(traj(:,2),1,npts) - repmat(traj(:,2)',npts,1);
dr = sqrt(dx.^2 + dy.^2);

nt = numel(tvalues);
displ = [];
for j=1:nt
    disptmp = [];
    ndata = sum(dt(:)== tvalues(j));

    if ndata == 0 
        disptmp = [];
    else
        disptmp(1:ndata,1) = traj_idx;
        disptmp(1:ndata,2) = tvalues(j);
        disptmp(1:ndata,3) = dx(dt(:)== tvalues(j))*pix_size;
        disptmp(1:ndata,4) = dy(dt(:)== tvalues(j))*pix_size;
        disptmp(1:ndata,5) = dr(dt(:)== tvalues(j))*pix_size;
    end

    displ = [displ;disptmp];
end
end
