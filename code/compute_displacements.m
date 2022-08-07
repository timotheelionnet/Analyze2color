function [dispco,dispall,dispsingle,Intensity] = compute_displacements(traj,pix_size)

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

dispco = cell('');
dispall = cell('');
dispsingle = cell('');
Intensity = zeros(numel(traj_idx,2));
for i = 1:numel(traj_idx)
    %select current trajectory
    tmptraj = traj(traj(:,10)==traj_idx(i),:);
    Intensity(i,1) = nanmean(tmptraj(:,3));
    Intensity(i,2) = nanmean(tmptraj(:,6));
    
    %crop trajectories for the various data sets
    
    %data set for dispall: keep all spots separately in each channel
    %regardless of codetections (including the isolated spots)
    tmptraj1 = tmptraj( logical( ~isnan(tmptraj(:,1))),:);
    npts1 = size(tmptraj1,1);
    tmptraj2 = tmptraj( logical( ~isnan(tmptraj(:,4))),:); 
    npts2 = size(tmptraj2,1);
    
    %dataset for dispsingles: in each channel, keep only the spots that
    %belong to a bonafide trajectory (exclude the spots that were added 
    %only because they colocalize with a trajectory in the other channel).
    %Include all these spots, regardless of colocalization events
    tmptraj1s = tmptraj( logical( (tmptraj(:,8)~=-1).*(~isnan(tmptraj(:,1))) ),:);
    npts1s = size(tmptraj1s,1);
    tmptraj2s = tmptraj( logical( (tmptraj(:,9)~=-1).*(~isnan(tmptraj(:,4))) ),:); 
    npts2s = size(tmptraj2s,1);
    
    %dataset for dispco: include all co-detected spots.
    tmptraj = tmptraj( logical( (~isnan(tmptraj(:,1))).*(~isnan(tmptraj(:,4)))),:);
    npts = size(tmptraj,1);
    
    
    %% compute displacements for only the codetected spots %%%%%%%%%%%%
    %find max number of data points for the current trajectory
    dt = repmat(tmptraj(:,7),1,npts) - repmat(tmptraj(:,7)',npts,1);
    tvalues = unique(dt(dt>0));
    nt = numel(tvalues);
    
    %extract the dx, dy and dr for all time points
    dx1 = repmat(tmptraj(:,1),1,npts) - repmat(tmptraj(:,1)',npts,1);
    dy1 = repmat(tmptraj(:,2),1,npts) - repmat(tmptraj(:,2)',npts,1);
    dr1 = sqrt(dx1.^2 + dy1.^2);
    dx2 = repmat(tmptraj(:,4),1,npts) - repmat(tmptraj(:,4)',npts,1);
    dy2 = repmat(tmptraj(:,5),1,npts) - repmat(tmptraj(:,5)',npts,1);
    dr2 = sqrt(dx2.^2 + dy2.^2);
    
    for j=1:nt
        disptmp = [];
        ndata = sum(dt(:)== tvalues(j));
        
        disptmp(1:ndata,1) = traj_idx(i);
        disptmp(1:ndata,2 ) = dx1(dt(:)== tvalues(j))*pix_size;
        disptmp(1:ndata,3 ) = dy1(dt(:)== tvalues(j))*pix_size;
        disptmp(1:ndata,4 ) = dr1(dt(:)== tvalues(j))*pix_size;
        disptmp(1:ndata,5 ) = dx2(dt(:)== tvalues(j))*pix_size;
        disptmp(1:ndata,6 ) = dy2(dt(:)== tvalues(j))*pix_size;
        disptmp(1:ndata,7 ) = dr2(dt(:)== tvalues(j))*pix_size;
        if size(dispco,1) < tvalues(j)
            dispco{tvalues(j),1} = disptmp;
        else
            dispco{tvalues(j),1} = [dispco{tvalues(j),1};disptmp];
        end
    end
    
    %% compute displacements for all spots %%%%%%%%%%%%
    dt1 = repmat(tmptraj1(:,7),1,npts1) - repmat(tmptraj1(:,7)',npts1,1);
    tvalues1 = unique(dt1(dt1>0));
    
    dt2 = repmat(tmptraj2(:,7),1,npts2) - repmat(tmptraj2(:,7)',npts2,1);
    tvalues2 = unique(dt2(dt2>0));
    
    %extract the dx, dy and dr for all time points
    dx1 = repmat(tmptraj1(:,1),1,npts1) - repmat(tmptraj1(:,1)',npts1,1);
    dy1 = repmat(tmptraj1(:,2),1,npts1) - repmat(tmptraj1(:,2)',npts1,1);
    dr1 = sqrt(dx1.^2 + dy1.^2);
    dx2 = repmat(tmptraj2(:,4),1,npts2) - repmat(tmptraj2(:,4)',npts2,1);
    dy2 = repmat(tmptraj2(:,5),1,npts2) - repmat(tmptraj2(:,5)',npts2,1);
    dr2 = sqrt(dx2.^2 + dy2.^2);
    
    tvalues = unique(union(tvalues1,tvalues2));
    nt = numel(tvalues);
    for j=1:nt
        disptmp = [];
        ndata1 = sum(dt1(:)== tvalues(j));
        ndata2 = sum(dt2(:)== tvalues(j));
        
        if ndata1 == 0 && ndata2 == 0
            disptmp = [];
        elseif ndata1 == 0 && ndata2 ~= 0
            disptmp(1:ndata2,1) = traj_idx(i);
            disptmp(1:ndata2,2 ) = NaN;
            disptmp(1:ndata2,3 ) = NaN;
            disptmp(1:ndata2,4 ) = NaN;
            disptmp(1:ndata2,5 ) = dx2(dt2(:)== tvalues(j))*pix_size;
            disptmp(1:ndata2,6 ) = dy2(dt2(:)== tvalues(j))*pix_size;
            disptmp(1:ndata2,7 ) = dr2(dt2(:)== tvalues(j))*pix_size;
            
        elseif ndata2 == 0 && ndata1 ~= 0
            disptmp(1:ndata1,1) = traj_idx(i);
            disptmp(1:ndata1,2 ) = dx1(dt1(:)== tvalues(j))*pix_size;
            disptmp(1:ndata1,3 ) = dy1(dt1(:)== tvalues(j))*pix_size;
            disptmp(1:ndata1,4 ) = dr1(dt1(:)== tvalues(j))*pix_size;
            disptmp(1:ndata1,5 ) = NaN;
            disptmp(1:ndata1,6 ) = NaN;
            disptmp(1:ndata1,7 ) = NaN;
            
        elseif ndata2 ~= 0 && ndata1 ~= 0
            ndata = max(ndata1,ndata2);
            disptmp = NaN(ndata,7);
            disptmp(1:ndata,1) = traj_idx(i);
            disptmp(1:ndata1,2 ) = dx1(dt1(:)== tvalues(j))*pix_size;
            disptmp(1:ndata1,3 ) = dy1(dt1(:)== tvalues(j))*pix_size;
            disptmp(1:ndata1,4 ) = dr1(dt1(:)== tvalues(j))*pix_size;
            disptmp(1:ndata2,5 ) = dx2(dt2(:)== tvalues(j))*pix_size;
            disptmp(1:ndata2,6 ) = dy2(dt2(:)== tvalues(j))*pix_size;
            disptmp(1:ndata2,7 ) = dr2(dt2(:)== tvalues(j))*pix_size;
        end
        
        if size(dispall,1) < tvalues(j)
            dispall{tvalues(j),1} = disptmp;
        else
            dispall{tvalues(j),1} = [dispall{tvalues(j),1};disptmp];
        end
    end
    
    %% compute displacements for trajectory-only spots %%%%%%%%%%%%
    dt1 = repmat(tmptraj1s(:,7),1,npts1s) - repmat(tmptraj1s(:,7)',npts1s,1);
    tvalues1 = unique(dt1(dt1>0));
    
    dt2 = repmat(tmptraj2s(:,7),1,npts2s) - repmat(tmptraj2s(:,7)',npts2s,1);
    tvalues2 = unique(dt2(dt2>0));
    
    %extract the dx, dy and dr for all time points
    dx1 = repmat(tmptraj1s(:,1),1,npts1s) - repmat(tmptraj1s(:,1)',npts1s,1);
    dy1 = repmat(tmptraj1s(:,2),1,npts1s) - repmat(tmptraj1s(:,2)',npts1s,1);
    dr1 = sqrt(dx1.^2 + dy1.^2);
    dx2 = repmat(tmptraj2s(:,4),1,npts2s) - repmat(tmptraj2s(:,4)',npts2s,1);
    dy2 = repmat(tmptraj2s(:,5),1,npts2s) - repmat(tmptraj2s(:,5)',npts2s,1);
    dr2 = sqrt(dx2.^2 + dy2.^2);
    
    tvalues = unique(union(tvalues1,tvalues2));
    nt = numel(tvalues);
    for j=1:nt
        disptmp = [];
        ndata1 = sum(dt1(:)== tvalues(j));
        ndata2 = sum(dt2(:)== tvalues(j));
        
        if ndata1 == 0 && ndata2 == 0
            disptmp = [];
        elseif ndata1 == 0 && ndata2 ~= 0
            disptmp(1:ndata2,1) = traj_idx(i);
            disptmp(1:ndata2,2 ) = NaN;
            disptmp(1:ndata2,3 ) = NaN;
            disptmp(1:ndata2,4 ) = NaN;
            disptmp(1:ndata2,5 ) = dx2(dt2(:)== tvalues(j))*pix_size;
            disptmp(1:ndata2,6 ) = dy2(dt2(:)== tvalues(j))*pix_size;
            disptmp(1:ndata2,7 ) = dr2(dt2(:)== tvalues(j))*pix_size;
            
        elseif ndata2 == 0 && ndata1 ~= 0
            disptmp(1:ndata1,1) = traj_idx(i);
            disptmp(1:ndata1,2 ) = dx1(dt1(:)== tvalues(j))*pix_size;
            disptmp(1:ndata1,3 ) = dy1(dt1(:)== tvalues(j))*pix_size;
            disptmp(1:ndata1,4 ) = dr1(dt1(:)== tvalues(j))*pix_size;
            disptmp(1:ndata1,5 ) = NaN;
            disptmp(1:ndata1,6 ) = NaN;
            disptmp(1:ndata1,7 ) = NaN;
            
        elseif ndata2 ~= 0 && ndata1 ~= 0
            ndata = max(ndata1,ndata2);
            disptmp = NaN(ndata,7);
            disptmp(1:ndata,1) = traj_idx(i);
            disptmp(1:ndata1,2 ) = dx1(dt1(:)== tvalues(j))*pix_size;
            disptmp(1:ndata1,3 ) = dy1(dt1(:)== tvalues(j))*pix_size;
            disptmp(1:ndata1,4 ) = dr1(dt1(:)== tvalues(j))*pix_size;
            disptmp(1:ndata2,5 ) = dx2(dt2(:)== tvalues(j))*pix_size;
            disptmp(1:ndata2,6 ) = dy2(dt2(:)== tvalues(j))*pix_size;
            disptmp(1:ndata2,7 ) = dr2(dt2(:)== tvalues(j))*pix_size;
        end
        
        if size(dispsingle,1) < tvalues(j)
            dispsingle{tvalues(j),1} = disptmp;
        else
            dispsingle{tvalues(j),1} = [dispsingle{tvalues(j),1};disptmp];
        end
    end
    
    
end