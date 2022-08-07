function [coloc_sum,coloc_traj,undecided_sum,undecided_traj,neg_sum,neg_traj,consolidated_trk1,trk1sum,consolidated_trk2,trk2sum] = ...
    match_trajectories3_dia(trk1,trk2,minpts,max_center_dist,max_spot_dist)

%[robust_sum,robust_traj,match_sum,match_traj,consolidated_trk1,trk1sum,consolidated_trk2,trk2sum] = ...
%    match_trajectories2(trk1,trk2,loc1,loc2)

%% inputs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%trk1, trk2 are lists of trajectories for each channel: 
%it includes spots positions, frame number and trajectory assigned, formated as follows:
%col 1: x
%col 2: y
%col 3: I
%col 4: frame number
%col 5: trajectory number

%loc1, loc2 are the full list of spots for each channel (including those not part of any
%trajectory), formated as follows:
%col 1: x
%col 2: y
%col 3: I
%col 4: frame number

%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%robust_sum: the summary of all the robustly matched trajectories (stringent criteria)
%these trajectories can be used for diffusion studies 
%1 row per trajectory
%col1 = trajectory 1 index, 
%col2 = trajectory 2 index, 
%col 3 = start frame 1,
%col 4 = end frame 1, 
%col 5 = start frame 2, 
%col 6 = end frame 2, 
%col 7 = number of overlapping frames, 
%col 8 = number of simultaneous detections, 
%col 9 = dx of trajectories centers, 
%col 10 = dy of trajectories centers

%robust_traj: the list of all the robustly matched trajectories (stringent criteria)
%these trajectories can be used for diffusion studies 
%list of spots that belong to matched trajctories
%1 row per spot
%col 1 = x1
%col 2 = y1
%col 3 = I1
%col 4 = x2
%col 5 = y2
%col 6 = I2
%col 7 = frame num
%col 8 = index of trajectory 1
%col 9 = index of trajectory 2
%col 10 = index of matched trajectory

%match_sum: the summary of all the matched trajectories (inclusive criteria)
%these trajectories are there for ID only, in order to obtain the fraction
%of colocalized trajectories

%match_traj: the list of all the matched trajectories (inclusive criteria)
%these trajectories are there for ID only, in order to obtain the fraction
%of colocalized trajectories

%consolidated_trk 1: same list and format as the input trk1, except that
%missing points within trajectories are included, with NaNs in lieu of
%position/intensity

%trk1sum: summary of the trajectories in channel 1, formatted as follows:
%col 1 = traj idx,
%col 2 = average x, 
%col 3 = average y, 
%col 4 = std x, 
%col 5 = std y, 
%col 6 = start frame, 
%col 7 = end frame,
%col 8 = number of actual spot detections


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters to match trajectories together
%minimum number of frames over which trajectories need to overlap in order to call a
%robust match. Note that the tajectories might have less codetection events
%than their overlap if the tracking is spotty.
%minpts = 3;

%maximum distance between the trajectory centers (in pix)
%max_center_dist = 5;


%maximum distance between codetected spots (in pix)
%max_spot_dist = 2;

%z-score threshold to determine statistically significant co-detections
zscore_thresh = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add to each trajectory the missing points with NaNs in lieu of position/Intensity
%makes comparisons between trajectories easier

%trajectory 1
traj1 = unique(trk1(:,end));
ntraj1 = numel(traj1);

%the consolidated list of spots indexed by trajectories
% I add to the exisiting list place holder spots with NaNs
%col 1 = x1
%col 2 = y1
%col 3 = I1
%col 4 = frame num
%col 5 = traj num
consolidated_trk1 = zeros(size(trk1));

%the array that stores for each trajectory:
%col 1 = traj1 idx,
%col 2 = average x1, 
%col 3 = average y1, 
%col 4 = std x1, 
%col 5 = std y1, 
%col 6 = start frame 1, 
%col 7 = end frame 1,
%col 8 = number of spot detections 1
trk1sum = zeros(ntraj1,8);

for i=1:ntraj1
    
    %extract ith trajectory
    trktmp = trk1(trk1(:,end)== traj1(i),:);
    
    %store average position, std and start/end frames
    trk1sum(i,1) = traj1(i);
    trk1sum(i,2) = mean(trktmp(:,1));
    trk1sum(i,3) = mean(trktmp(:,2));
    trk1sum(i,4) = std(trktmp(:,1));
    trk1sum(i,5) = std(trktmp(:,2));
    trk1sum(i,6) = trktmp(1,4);
    trk1sum(i,7) = trktmp(end,4);
    trk1sum(i,8) = size(trktmp,1);
    
    %build array consolidated_trk1 that will old the full trajectory
    trktmpo = NaN(trktmp(size(trktmp,1),4) - trktmp(1,4)+1,size(trktmp,2));
    trktmpo(:,4) = trktmp(1,4):trktmp(size(trktmp,1),4);
    trktmpo(:,5) = traj1(i);
    for j=1:size(trktmp,1)
        trktmpo(trktmpo(:,4)==trktmp(j,4),:) = trktmp(j,:);
    end

    if i==1
        consolidated_trk1 = trktmpo;
    else
        consolidated_trk1 = [consolidated_trk1;trktmpo];
    end
end

traj1 = unique(trk1sum(:,1));
ntraj1 = numel(traj1);

%repeat for trajectory 2
traj2 = unique(trk2(:,end));
ntraj2 = numel(traj2);

%the consolidated list of spots indexed by trajectories
% I add to the exisiting list place holder spots with NaNs
%col 1 = x1
%col 2 = y1
%col 3 = I1
%col 4 = frame num
%col 5 = traj num
consolidated_trk2 = zeros(size(trk2));

%the array that stores for each trajectory:
%col 1 = traj2 idx,
%col 2 = average x2, 
%col 3 = average y2, 
%col 4 = std x2, 
%col 5 = std y2, 
%col 6 = start frame 2, 
%col 7 = end frame 2,
%col 8 = number of spot detections 2
trk2sum = zeros(ntraj2,8);

for i=1:ntraj2
    
    %extract ith trajectory
    trktmp = trk2(trk2(:,end)== traj2(i),:);
    
    %store average position, std and start/end frames
    trk2sum(i,1) = traj2(i);
    trk2sum(i,2) = mean(trktmp(:,1));
    trk2sum(i,3) = mean(trktmp(:,2));
    trk2sum(i,4) = std(trktmp(:,1));
    trk2sum(i,5) = std(trktmp(:,2));
    trk2sum(i,6) = trktmp(1,4);
    trk2sum(i,7) = trktmp(end,4);
    trk2sum(i,8) = size(trktmp,1);
    
    %build array consolidated_trk2 that will old the full trajectory
    trktmpo = NaN(trktmp(size(trktmp,1),4) - trktmp(1,4)+1,size(trktmp,2));
    trktmpo(:,4) = trktmp(1,4):trktmp(size(trktmp,1),4);
    trktmpo(:,5) = traj2(i);
    for j=1:size(trktmp,1)
        trktmpo(trktmpo(:,4)==trktmp(j,4),:) = trktmp(j,:);
    end

    if i==1
        consolidated_trk2 = trktmpo;
    else
        consolidated_trk2 = [consolidated_trk2;trktmpo];
    end
end

traj2 = unique(trk2sum(:,1));
ntraj2 = numel(traj2);

%% Match Trajectories

%summary of matched trajectories:
%1 row per trajectory
%col1 = traj1 idx, 
%col2 = traj2 idx, 
%col 3 = start frame 1,
%col 4 = end frame 1, 
%col 5 = start frame 2, 
%col 6 = end frame 2, 
%col 7 = number of overlapping frames, 
%col 8 = number of simultaneous detections, 
%col 9 = dx of trajectories centers, 
%col 10 = dy of trajectories centers
match_sum = zeros(1,10);

%list of spots that belong to matched trajctories
%1 row per spot
%col 1 = x1
%col 2 = y1
%col 3 = I1
%col 4 = x2
%col 5 = y2
%col 6 = I2
%col 7 = frame num
%col 8 = index of trajectory 1
%col 9 = index of trajectory 2
%col 10 = index of matched trajectory
matched_trajectories = zeros(1,10);

nmatchtraj = 0;
nmatchspots = 0;
for i=1:ntraj1
    for j=1:ntraj2
        %check for sufficient overlap in time
        traj1idx = find(trk1sum(:,1) == traj1(i));
        traj2idx = find(trk2sum(:,1) == traj2(j));
        overlap = min(trk1sum(traj1idx,7),trk2sum(traj2idx,7)) ...
            - max(trk1sum(traj1idx,6),trk2sum(traj2idx,6)) +1;
        
        %check that positions match
        if overlap >= minpts
            %maxdist_x = min(maxdist*trk1sum(traj1idx,4),maxdist*trk2sum(traj2idx,4));
            %maxdist_y = min(maxdist*trk1sum(traj1idx,5),maxdist*trk2sum(traj2idx,5));
            dx = abs(trk1sum(traj1idx,2) - trk2sum(traj2idx,2));
            dy = abs(trk1sum(traj1idx,3) - trk2sum(traj2idx,3));
            
            % if centers are close enough, look for enough codetected spots
            if sqrt(dx^2 + dy^2) < max_center_dist;
                do_traj_match = 1;
            else
                do_traj_match = 0;
            end
        else
            do_traj_match = 0;
        end
        
        %store trajectories in case of matching event
        if do_traj_match
            nmatchtraj = nmatchtraj +1;
            
            %make sure the start and end frame are identical by filling in
            %NaNs
            trk1tmp = consolidated_trk1( consolidated_trk1(:,5) == traj1(i),: );
            trk2tmp = consolidated_trk2( consolidated_trk2(:,5) == traj2(j),: );
            
            ts = min(trk1sum(traj1idx,6),trk2sum(traj2idx,6));
            te = max(trk1sum(traj1idx,7),trk2sum(traj2idx,7));
            
            if trk1tmp(1,4) > ts
               trk1tmp = [ NaN( trk1tmp(1,4) - ts ,5) ; trk1tmp];
               trk1tmp(1 : (trk1sum(traj1idx,6) - ts),4) = ts : (trk1sum(traj1idx,6) - 1);
               trk1tmp(1 : (trk1sum(traj1idx,6) - ts),5) = traj1(i);
            end
            
            if trk2tmp(1,4) > ts
               trk2tmp = [ NaN( trk2tmp(1,4) - ts ,5) ; trk2tmp];
               trk2tmp(1 : (trk2sum(traj2idx,6) - ts),4) = ts : (trk2sum(traj2idx,6) - 1);
               trk2tmp(1 : (trk2sum(traj2idx,6) - ts),5) = traj2(j);
            end
            
            if trk1tmp(end,4) < te
                trk1tmp = [trk1tmp; NaN( te - trk1tmp(end,4) ,5)];
                n1 = size(trk1tmp,1);
                trk1tmp((n1 - (te - trk1sum(traj1idx,7)) +1) :n1 , 4) = (trk1sum(traj1idx,7)+1):te;
                trk1tmp((n1 - (te - trk1sum(traj1idx,7)) +1) :n1 , 5) = traj1(i);
            end
            
            if trk2tmp(end,4) < te
                trk2tmp = [trk2tmp; NaN( te - trk2tmp(end,4) ,5)];
                n2 = size(trk2tmp,1);
                trk2tmp((n2 - (te - trk2sum(traj2idx,7)) +1) :n2 , 4) = (trk2sum(traj2idx,7)+1):te;
                trk2tmp((n2 - (te - trk2sum(traj2idx,7)) +1) :n2 , 5) = traj2(j);
            end
            
            imatch = nmatchspots + 1 : nmatchspots + te - ts +1;
            
            matched_trajectories(imatch,1:3) = trk1tmp(:,1:3);
            matched_trajectories(imatch,4:6) = trk2tmp(:,1:3);
            matched_trajectories(imatch,7) = trk2tmp(:,4);
            matched_trajectories(imatch,8) = traj1(i);
            matched_trajectories(imatch,9) = traj2(j);
            matched_trajectories(imatch,10) = nmatchtraj;
            
            match_sum(nmatchtraj,1) = traj1(i);
            match_sum(nmatchtraj,2) = traj2(j);
            match_sum(nmatchtraj,3) = trk1sum(traj1idx,6);
            match_sum(nmatchtraj,4) = trk1sum(traj1idx,7);
            match_sum(nmatchtraj,5) = trk2sum(traj2idx,6);
            match_sum(nmatchtraj,6) = trk2sum(traj2idx,7);
            match_sum(nmatchtraj,7) = overlap;
            match_sum(nmatchtraj,8) = sum(~isnan(matched_trajectories(imatch,1).*matched_trajectories(imatch,4)));
            match_sum(nmatchtraj,9) = dx;
            match_sum(nmatchtraj,10) = dy;
            match_sum(nmatchtraj,11) = nmatchtraj;
            nmatchspots = nmatchspots +  trk1sum(traj1idx,7) - trk1sum(traj1idx,6);    
        end
    end
end

%% clean up trajectories of false codetected spots based on a distance threshold
robust_traj = zeros(1,10);
robust_sum = zeros(1,10);
%Spot to spot distance
dr = sqrt((matched_trajectories(:,1) - matched_trajectories(:,4)).^2 +...
    (matched_trajectories(:,2) - matched_trajectories(:,5)).^2);
%spots that are too far away
dr = dr > max_spot_dist;

%go through each trajectory and clean up individual spots that do not colocalize
cleantraj = 0;
cleanspots = 0;
false_traj = zeros(nmatchtraj,1);

for i = 1:nmatchtraj
    idx = matched_trajectories(:,10) == match_sum(i,11);   
    tmptrk = matched_trajectories(idx,:);
    tmpdr = dr(idx);
    tmptrk(tmpdr,1:6) = NaN;
    
    %crop beginning and end of trajectory if they consist of unmatched
    %detections; ignore trajectories that are completely out of range;
    %ignore trajectories with unique meaningful codetection event
        
    %find index of simultaneous detection events
    ipair = find((~isnan(tmptrk(:,1))).*(~isnan(tmptrk(:,4))));
    
    if numel(ipair)>=minpts
        %crop trajectory
        cleantraj = cleantraj+1;
        tmptrk = tmptrk(ipair(1):ipair(end),:);
        tmp_cleanspots = size(tmptrk,1);
        robust_traj(cleanspots+1:cleanspots+tmp_cleanspots,:) = tmptrk;
        robust_traj(cleanspots+1:cleanspots+tmp_cleanspots,10) = cleantraj;
        cleanspots = cleanspots+tmp_cleanspots;
        
        %update trajectory summary info
        robust_sum(cleantraj,1:2) = match_sum(i,1:2);
        robust_sum(cleantraj,3) = tmptrk(1,7);
        robust_sum(cleantraj,4) = tmptrk(end,7);
        robust_sum(cleantraj,5) = tmptrk(1,7);
        robust_sum(cleantraj,6) = tmptrk(end,7);
        robust_sum(cleantraj,7) = size(tmptrk,1);
        robust_sum(cleantraj,8) = sum(~isnan(tmptrk(:,1).*tmptrk(:,4)));
        ctrpos = nanmean(tmptrk(:,1:5));
        robust_sum(cleantraj,9) = abs(ctrpos(1) - ctrpos(4));
        robust_sum(cleantraj,10) = abs(ctrpos(2) - ctrpos(5));
        robust_sum(cleantraj,11) = cleantraj;
    else
        %store index of trajectory to eliminate it later
        false_traj(i) = 1;
    end    
end

%eliminate trajectories with low codetection events
for i = 1:nmatchtraj
    if false_traj(i) == 1
        matched_trajectories = matched_trajectories(matched_trajectories(:,10)~=match_sum(i,11),:);
    end
end
false_traj = logical(false_traj);
match_sum = match_sum(match_sum(~false_traj,11),:);


%% make sure each trajectory portion has only one correspondent
%so far I display warnings only, no steps taken to fix the issue, since I
%haven't encountered it yet
for i = 1:ntraj1
    idx1 = traj1(i);
    tmptraj = robust_sum(robust_sum(:,1)==idx1,:);
    if size(tmptraj,1)>1
        for k=1:size(tmptraj,1)-1
            for l=k+1:size(tmptraj,1)
                overlap = min(tmptraj(k,4),tmptraj(l,4)) ...
                    - max(tmptraj(k,3),tmptraj(l,3)) +1;
                if overlap >=1
                    disp(['WARNING: found overlapped trajectories ch. 1: ',num2str(tmptraj(k,1)),' with ch. 2: ',num2str(tmptraj(:,2)')]);
                end
            end
        end
    end    
end

for i = 1:ntraj2
    idx2 = traj2(i);
    tmptraj = robust_sum(robust_sum(:,2)==idx2,:);
    if size(tmptraj,1)>1
        for k=1:size(tmptraj,1)-1
            for l=k+1:size(tmptraj,1)
                overlap = min(tmptraj(k,6),tmptraj(l,6)) ...
                    - max(tmptraj(k,5),tmptraj(l,5)) +1;
                if overlap >=1
                    disp(['WARNING: found overlapped trajectories ch. 2: ',num2str(tmptraj(k,2)),' with ch. 1: ',num2str(tmptraj(:,1)')]);
                end
            end
        end
    end    
end

%% Detect non-robust trajectory matches
%Browse through non-robustly matched trajectories to find potential 
%matches with individual spots from the complementary channel

%Match orphan trajectories from channel 1 with spots in channel 2
%select orphan trajectories from channel 1, 
%i.e. trajectories without a match

ntraj1 = unique(trk1sum(:,1));
nrobust1 = unique(robust_sum(:,1));
orphans1 = ntraj1(~ismember(ntraj1,nrobust1));

match_traj = [];
match_sum = [];

nmatchtraj = 0;
nrobusttraj = size(robust_sum,1);
for i=1:numel(orphans1)
    %select current orphan trajectory
    tmptrk = consolidated_trk1(consolidated_trk1(:,5)==orphans1(i),:);
    
    %crop out from the spots array spots that do not overlap temporally with any spot in the
    %trajectory
    tmpspots = loc2(logical( (loc2(:,4)>=tmptrk(1,4)).*(loc2(:,4)<=tmptrk(end,4))) ,:);
    
    %try to find matches
    match_found = 0;
    nmatchspots = 0;
    match_traj_tmp = [];
    for j=1:size(tmptrk,1)
        %extract spots detected at the current frame
        tmpspots2 = tmpspots(tmpspots(:,4) == tmptrk(j,4),:);

        %compute distances
        dist = sqrt( (tmpspots2(:,1) - tmptrk(j,1)).^2 + (tmpspots2(:,2) - tmptrk(j,2)).^2 );

        %find spots within match radius
        if sum(dist <= max_spot_dist) > 1
            [~,idx] = min(dist);
        elseif sum(dist <= max_spot_dist) == 1
            idx = find(dist <= max_spot_dist);
        else
            idx = [];
        end

        %fill in trajectory array if a spot was found
        if ~isempty(idx)
            match_found = match_found +1;
            if match_found == 1
                nmatchtraj = nmatchtraj + 1;
                ts = tmptrk(j,4);
                te = tmptrk(j,4);
            else
                te = tmptrk(j,4);
            end
            nmatchspots = nmatchspots + 1;
            match_traj_tmp(j,1:3) = tmptrk(j,1:3);
            match_traj_tmp(j,4:6) = tmpspots2(idx,1:3);
            match_traj_tmp(j,7) = tmptrk(j,4);
            match_traj_tmp(j,8) = tmptrk(j,5);
            match_traj_tmp(j,9) = -1;
            match_traj_tmp(j,10) = nrobusttraj+nmatchtraj;
        else
            match_traj_tmp(j,1:3) = tmptrk(j,1:3);
            match_traj_tmp(j,4:6) = NaN;
            match_traj_tmp(j,7) = tmptrk(j,4);
            match_traj_tmp(j,8) = tmptrk(j,5);
            match_traj_tmp(j,9) = -1;
            match_traj_tmp(j,10) = NaN;
        end
    end
    
    if match_found >= 1
        match_traj_tmp(:,10) = max(match_traj_tmp(:,10));
        match_traj = [match_traj;match_traj_tmp];
        match_sum(nmatchtraj,1) = tmptrk(j,5);
        match_sum(nmatchtraj,2) = -1;
        match_sum(nmatchtraj,3) = tmptrk(1,4);
        match_sum(nmatchtraj,4) = tmptrk(end,4);
        match_sum(nmatchtraj,5) = -1;
        match_sum(nmatchtraj,6) = -1;
        match_sum(nmatchtraj,7) = te - ts +1;
        match_sum(nmatchtraj,8) = match_found;
        idx = match_traj(:,10) == nrobusttraj+nmatchtraj;
        dx = nanmean(  (match_traj(idx,1) - match_traj(idx,4))   );
        dy = nanmean(  (match_traj(idx,2) - match_traj(idx,5))   );
        match_sum(nmatchtraj,9) = abs(dx);
        match_sum(nmatchtraj,10) = abs(dy);
        match_sum(nmatchtraj,11) = nrobusttraj+nmatchtraj;
        
    end
end


%Match orphan trajectories from channel 2 with spots in channel 1
%select orphan trajectories from channel 2, 
%i.e. trajectories without a match
ntraj2 = unique(trk2sum(:,1));
nrobust2 = unique(robust_sum(:,2));
orphans2 = ntraj2(~ismember(ntraj2,nrobust2));

for i=1:numel(orphans2)
    %select current orphan trajectory
    tmptrk = consolidated_trk2(consolidated_trk2(:,5)==orphans2(i),:);
    
    %crop out from the spots array spots that do not overlap temporally with any spot in the
    %trajectory
    tmpspots = loc1(logical( (loc1(:,4)>=tmptrk(1,4)).*(loc1(:,4)<=tmptrk(end,4))) ,:);
    
    %try to find matches
    match_found = 0;
    nmatchspots = 0;
    match_traj_tmp = [];
    for j=1:size(tmptrk,1)
        %extract spots detected at the current frame
        tmpspots2 = tmpspots(tmpspots(:,4) == tmptrk(j,4),:);
        
        %compute distances
        dist = sqrt( (tmpspots2(:,1) - tmptrk(j,1)).^2 + (tmpspots2(:,2) - tmptrk(j,2)).^2 );
        
        %find spots within match radius
        if sum(dist <= max_spot_dist) > 1
            [~,idx] = min(dist);
        elseif sum(dist <= max_spot_dist) == 1
            idx = find(dist <= max_spot_dist);
        else
            idx = [];
        end
        
        %fill in trajectory arra if a spot was found
        if ~isempty(idx)
            match_found = match_found +1;
            if match_found == 1
                nmatchtraj = nmatchtraj + 1;
                ts = tmptrk(j,4);
                te = tmptrk(j,4);
            else
                te = tmptrk(j,4);
            end
            nmatchspots = nmatchspots + 1;
            match_traj_tmp(j,1:3) = tmpspots2(idx,1:3);
            match_traj_tmp(j,4:6) =tmptrk(j,1:3);
            match_traj_tmp(j,7) = tmptrk(j,4);
            match_traj_tmp(j,8) = -1;
            match_traj_tmp(j,9) = tmptrk(j,5);
            match_traj_tmp(j,10) = nrobusttraj+nmatchtraj;
        else
            match_traj_tmp(j,1:3) = NaN;
            match_traj_tmp(j,4:6) =tmptrk(j,1:3);
            match_traj_tmp(j,7) = tmptrk(j,4);
            match_traj_tmp(j,8) = -1;
            match_traj_tmp(j,9) = tmptrk(j,5);
            match_traj_tmp(j,10) = NaN;
        end        
    end
    
    if match_found >= 1
        %make sure the trajectory index is correct
        match_traj_tmp(:,10) = max(match_traj_tmp(:,10));
        match_traj = [match_traj;match_traj_tmp];
        match_sum(nmatchtraj,1) = -1;
        match_sum(nmatchtraj,2) = tmptrk(j,5);
        match_sum(nmatchtraj,3) = -1;
        match_sum(nmatchtraj,4) = -1;
        match_sum(nmatchtraj,5) = tmptrk(1,4);
        match_sum(nmatchtraj,6) = tmptrk(end,4);
        match_sum(nmatchtraj,7) = te - ts +1;
        match_sum(nmatchtraj,8) = match_found;
        idx = match_traj(:,10) == nrobusttraj+nmatchtraj;
        dx = nanmean(  (match_traj(idx,1) - match_traj(idx,4))   );
        dy = nanmean(  (match_traj(idx,2) - match_traj(idx,5))   );
        match_sum(nmatchtraj,9) = abs(dx);
        match_sum(nmatchtraj,10) = abs(dy);
        match_sum(nmatchtraj,11) = nrobusttraj+nmatchtraj;
    end
end

%% Statistical Significance Check For Non-Robust Trajectories
%I only keep matching events when they are significantly more
%frequent than random co-detections
if ~isempty(match_sum)
    ImgArea = (maxx - minx+1)*(maxy - miny+1);
    DetectionCrossSection = pi*max_spot_dist^2;
    FalsePosProb = DetectionCrossSection / ImgArea;

    %compute # of spots detected / frame in each channel
    frames = unique([loc1(:,4);loc2(:,4)]);
    for i=1:numel(frames)
        Nspots(i,1) = sum(loc1(:,4) == frames(i));
        Nspots(i,2) = sum(loc2(:,4) == frames(i));
    end

    %compare the number of co-detected spots to what is expected by chance
    %quantify it by a "z-score" 
    FalsePos = [];
    nfalse = 0;

    for i=1:size(match_sum,1)
        %number of codetections
        tmptrk = match_traj(match_traj(:,10)==match_sum(i,11),:);
        NumDetections = sum(~isnan(tmptrk(:,1).*tmptrk(:,4)));

        %remove frames where the trajectory is lost
        if match_sum(i,2) == -1
            tmptrk = tmptrk(~isnan(tmptrk(:,1)),:);
        else
            tmptrk = tmptrk(~isnan(tmptrk(:,4)),:);
        end

        %compute the cumulated probability of spot detection over the
        %trajectory
        DetProb = 0;
        if match_sum(i,2) == -1
            for j=1:size(tmptrk,1)
                DetProb = DetProb + FalsePosProb*Nspots(frames == tmptrk(j,7),2);
            end
        else
            for j=1:size(tmptrk,1)
                DetProb = DetProb + FalsePosProb*Nspots(frames == tmptrk(j,7),1);
            end
        end

        %assuming Poisson, compute a "z-score" relative to the null hypothesis
        %(null = spurious detection)
        if DetProb == 0
            zscore = Inf;
        else
            zscore = (NumDetections - DetProb)/sqrt(DetProb);
        end
        match_sum(i,12) = zscore;

        if zscore < zscore_thresh || NumDetections < minpts
            nfalse = nfalse + 1;
            FalsePos(nfalse) = match_sum(i,11);
        end    
    end

    %exclude trajectories with a low z-score 
    [~,real_traj_idx] = setdiff(match_sum(:,11),FalsePos);
    match_sum_clean = match_sum(real_traj_idx,:);
    match_traj_clean = [];
    for i=1:numel(real_traj_idx)
        match_traj_clean = [match_traj_clean;...
            match_traj(match_traj(:,10)== match_sum(real_traj_idx(i),11),:)];    
    end
    
else
    match_sum_clean = [];
    match_traj_clean = [];
    FalsePos = [];
    nfalse = 0;
end

%% Format results

%colocalized trajectories (robust and non-robust)
%combine robust w/ non-robust
coloc_sum = [[robust_sum,zeros(size(robust_sum,1),1)];match_sum_clean];
coloc_traj = [robust_traj;match_traj_clean];

%undecided trajectories (matched spots, but with low z-scores)
undecided_sum = [];
undecided_traj = [];
for i=1:numel(FalsePos)
    undecided_sum = [undecided_sum ; match_sum(match_sum(:,11) == FalsePos(i),:)];
    undecided_traj = [undecided_traj ; match_traj(match_traj(:,10)==FalsePos(i),:)];
end

%noncolocalized trajectories
matchidx1 = union( unique(coloc_sum(:,1)) , unique(undecided_sum(:,1)) );
matchidx1 = matchidx1(matchidx1~=-1);
orphanidx1 = setdiff(trk1sum(:,1) ,matchidx1);

matchidx2 = union( unique(coloc_sum(:,2)) , unique(undecided_sum(:,2)) );
matchidx2 = matchidx2(matchidx2~=-1);
orphanidx2 = setdiff(trk2sum(:,1) ,matchidx2);

neg_sum = [];
neg_traj= [];
nspots = 0;
ntraj = 0;
for i=1:numel(orphanidx1)
    ntraj = ntraj+1;
    tmptrk = consolidated_trk1(consolidated_trk1(:,5)==orphanidx1(i),:);
    neg_sum(ntraj,1) = tmptrk(1,5);
    neg_sum(ntraj,2) = -1;
    neg_sum(ntraj,3) = tmptrk(1,4);
    neg_sum(ntraj,4) = tmptrk(end,4);
    neg_sum(ntraj,5) = -1;
    neg_sum(ntraj,6) = -1;
    neg_sum(ntraj,7) = tmptrk(end,4) - tmptrk(1,4) +1;
    neg_sum(ntraj,8) = 0;
    neg_sum(ntraj,9) = NaN;
    neg_sum(ntraj,10) = NaN;
    neg_sum(ntraj,11) = ntraj;
    
    neg_traj(nspots+1:nspots+size(tmptrk,1),1:3) = tmptrk(:,1:3);
    neg_traj(nspots+1:nspots+size(tmptrk,1),4:6) = NaN;
    neg_traj(nspots+1:nspots+size(tmptrk,1),7) = tmptrk(:,4);
    neg_traj(nspots+1:nspots+size(tmptrk,1),8) = tmptrk(:,5);
    neg_traj(nspots+1:nspots+size(tmptrk,1),9) = -1;
    neg_traj(nspots+1:nspots+size(tmptrk,1),10) = ntraj;        
    nspots = nspots+size(tmptrk,1);
end

for i=1:numel(orphanidx2)
    ntraj = ntraj+1;
    tmptrk = consolidated_trk2(consolidated_trk2(:,5)==orphanidx2(i),:);
    neg_sum(ntraj,1) = -1;
    neg_sum(ntraj,2) = tmptrk(1,5);
    neg_sum(ntraj,3) = -1;
    neg_sum(ntraj,4) = -1;
    neg_sum(ntraj,5) = tmptrk(1,4);
    neg_sum(ntraj,6) = tmptrk(end,4);
    neg_sum(ntraj,7) = 0;
    neg_sum(ntraj,8) = tmptrk(end,4) - tmptrk(1,4) +1;
    neg_sum(ntraj,9) = NaN;
    neg_sum(ntraj,10) = NaN;
    neg_sum(ntraj,11) = ntraj;
    
    neg_traj(nspots+1:nspots+size(tmptrk,1),1:3) = NaN;
    neg_traj(nspots+1:nspots+size(tmptrk,1),4:6) = tmptrk(:,1:3);
    neg_traj(nspots+1:nspots+size(tmptrk,1),7) = tmptrk(:,4);
    neg_traj(nspots+1:nspots+size(tmptrk,1),8) = -1;
    neg_traj(nspots+1:nspots+size(tmptrk,1),9) = tmptrk(:,5);
    neg_traj(nspots+1:nspots+size(tmptrk,1),10) = ntraj;        
    nspots = nspots+size(tmptrk,1);
end

%% Display Results Summary in Pop-up window
ntraj1 = size(trk1sum,1);
ntraj1_coloc = numel(unique(coloc_sum(coloc_sum(:,1)~=-1,1)));
ntraj1_undecided  = numel(unique(undecided_sum(undecided_sum(:,1)~=-1,1)));
ntraj1_neg  = numel(unique(neg_sum(neg_sum(:,1)~=-1,1)));

ntraj2 = size(trk2sum,1);
ntraj2_coloc = numel(unique(coloc_sum(coloc_sum(:,2)~=-1,2)));
ntraj2_undecided  = numel(unique(undecided_sum(undecided_sum(:,2)~=-1,2)));
ntraj2_neg  = numel(unique(neg_sum(neg_sum(:,2)~=-1,2)));

x = sprintf(['Channel 1: ',num2str(ntraj1_coloc),' colocalized trajectories / ',...
    num2str(ntraj1),' total\n',num2str(ntraj1_undecided),' undecided trajectories,\n',...
    num2str(ntraj1_neg),' trajectories without a match\n\n',...
    'Channel 2: ',num2str(ntraj2_coloc),' colocalized trajectories / ',...
    num2str(ntraj2),' total\n',num2str(ntraj2_undecided),' undecided trajectories,\n',...
    num2str(ntraj2_neg),' trajectories without a match']);
dispwin('results',x);

end

