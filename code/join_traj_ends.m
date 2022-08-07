function [jtrk,junctions] = join_traj_ends(trk,dist_threshold,max_time_gap)
%joins together trajectories whose ends are less than max_time_gap frames
%apart and within a radius of dist_threshold

%trk should be formatted as a 5 column array [x,y,t,I,traj_idx]
%where each row is a datapoint; spots are ordered in increasing trajectory
%index and increasing time wihtin each trajectory

%jtrk output connects trajectories that satisfy the criteria
%replaces the column 5 by a new trajectory index
%each connected trajectory is given a new index
%the original traj_idx column is duplicated in column 6 for reference


%junction is an array that stores the list of trajectory junctions made by
%the program

%find all end - start connections
[junctions_end,dr_end,dt_end] = find_all_junction_candidates_endref(trk,dist_threshold,max_time_gap);
%if nothing is found get out of here
if isempty(junctions_end)
    jtrk = trk;
    junctions = unique(trk(:,5));
    return
elseif size(junctions_end,2)<2
    jtrk = trk;
    junctions = unique(trk(:,5));
    return
end
%select optimal start candidate for each trajectory end
match_junctions_end = connect_ends_with_best_start(junctions_end,dr_end,dt_end);

%find all start - end connections
[junctions_start,dr_start,dt_start] = find_all_junction_candidates_startref(trk,dist_threshold,max_time_gap);
%select optimal end candidate for each trajectory start
match_junctions_start = connect_ends_with_best_start(junctions_start,dr_start,dt_start);

%join trajectories only if they are mutual best candidates
junctions = [];
numjunctions = 0;
for i=1:size(match_junctions_end,1)
    if match_junctions_end(i,2) ~=0
        if match_junctions_start(find(match_junctions_start(:,1) == match_junctions_end(i,2)),2) == match_junctions_end(i,1)
                numjunctions = numjunctions+1;
                junctions(numjunctions,1:2) = match_junctions_end(i,1:2);
        end 
    end
end

%patch trajectories together
junctions = patch_junctions(junctions);
traj_idx = unique(trk(:,5));
ntraj = numel(traj_idx);
jtrk = [];
for i=1:ntraj
    k = find(traj_idx(i)==junctions(:,1));
    if ~isempty(k)
        for j=1:size(junctions,2)
            if (j==1) || (junctions(k,j) ~= 0)
                curtrk = trk(trk(:,5)== junctions(k,j),:);
                curtrk = [curtrk,curtrk(:,5)];
                curtrk(:,5) = max(traj_idx) + k;
                jtrk = [jtrk;curtrk];
            end
        end
    else
        if isempty(find(traj_idx(i)==junctions,1))
            curtrk = trk(trk(:,5)== traj_idx(i),:);
            jtrk = [jtrk;curtrk,curtrk(:,5)];
        end
    end
    
end

end

function junctions = patch_junctions(junctions)

i = 1;
while i <= size(junctions,1)
    
    tmin = junctions(:,1);
    tmax = max(junctions(:,2:end),[],2);
    curtmin = tmin(i);
    curtmax = tmax(i);
    
    if sum(tmax == curtmin)
        [imax,jmax] = find(junctions(:,2:end) == curtmin);
        junctions(imax,jmax+2) = junctions(i,2);
        junctions = junctions(~ismember(1:size(junctions,1),i),:);
        i = 1;
    elseif sum(tmin == curtmax)
        imin = find(junctions(:,1) == curtmax);
        [imax,jmax] = find(junctions(:,2:end) == curtmax);
        junctions(imax,jmax+2) = junctions(imin,2);
        junctions = junctions(~ismember(1:size(junctions,1),imin),:);
        i = 1;
    else
        i = i+1;
    end
end

end

function match_junctions = connect_ends_with_best_start(junctions,dr,dt)

match_junctions = junctions(:,1:2);
if size(junctions,2)>=3
    %find trajectory ends that have multiple compatible starts
    double_starts = find(sum(junctions~=0,2)>=3);    
    for i=1:numel(double_starts)       
        %find trajectory starts that are closest in time
        curdt = dt(double_starts(i),:);
        dtmin = min(curdt(curdt>0));
        ndoubles = find(curdt == dtmin);
        
        %if more than one compatible match, find closest trajectory
        if numel(ndoubles)>1
            curdr = dr(double_starts(i),:);
            curdr(curdt~=dtmin) = Inf;
            drmin = min(curdr);
            curmatch = find(curdr == drmin);
            match_junctions(double_starts(i),2) = junctions(double_starts(i),curmatch);
        else
            match_junctions(double_starts(i),2) = junctions(double_starts(i),ndoubles); 
        end    
    end    
end

end

function [matches,dr,dt] = find_all_junction_candidates_endref(trk,dist_threshold,max_time_gap)

%finds all trajectories whose ends are less than max_time_gap frames
%apart and within a radius of dist_threshold

%input:
%trk is trajectory array as a column array [x,y,I,t,traj index] where each
%row is a datapoint, order by increasing time

%outputs:
%matches: first column is traj idx of junction end, 
%other columns are indices of matched trajectory starts

%dr: first column is 0, 
%other columns are the [endpoint - startpoint] distances corresponding to the indices in "matches"

%dt: first column is 0, 
%other columns are the [endpoint - startpoint] time gaps corresponding to the indices in "matches"


%generate array of starts and ends
traj_idx = unique(trk(:,5));
ntraj = numel(traj_idx);

trk_start = zeros(ntraj,size(trk,2));
trk_end = zeros(ntraj,size(trk,2));
for i=1:ntraj
    tmptraj = trk(trk(:,5)==traj_idx(i),:);
    trk_start(i,:) = tmptraj(1,:);
    trk_end(i,:) = tmptraj(end,:);    
end

%find all possible matches between ends and starts
matches = [];%first column is traj idx of junction end, other columns are indices of matched trajectory starts
dr = [];%first column is 0, other columns are the [endpoint - startpoint] distances corresponding to the indices in "matches"
dt = [];%first column is 0, other columns are the [endpoint - startpoint] time gaps corresponding to the indices in "matches"
for i=1:ntraj
    
    %select time points that are compatible in time
    curdt = trk_start(:,4) - trk_end(i,4);
    poss_matches_idx = logical( ( curdt>0 ) .* ( curdt<= max_time_gap ) );
    cur_trk_start = trk_start(poss_matches_idx,:);
    curdt = curdt(poss_matches_idx);
    
    %select time points that are compatible in space
    curdr = sum((cur_trk_start(:,1:2) - repmat(trk_end(i,1:2),size(cur_trk_start,1),1) ).^2,2);
    poss_matches_idx = logical(curdr<=dist_threshold^2);
    poss_matches = cur_trk_start(poss_matches_idx,:);
    curdr = curdr(poss_matches_idx);
    curdt = curdt(poss_matches_idx);
    if ~isempty(poss_matches)
        matches(i,1) = traj_idx(i);
        dr(i,1) = 0;
        dt(i,1) = 0;
        matches(i,2:(size(poss_matches,1)+1)) = poss_matches(:,5)';
        dr(i,2:(size(poss_matches,1)+1)) = curdr;    
        dt(i,2:(size(poss_matches,1)+1)) = curdt; 
    else
        matches(i,1) = traj_idx(i);
        dr(i,1) = 0;
        dt(i,1) = 0;
    end    
end

end

function [matches,dr,dt] = find_all_junction_candidates_startref(trk,dist_threshold,max_time_gap)

%finds all trajectories whose ends are less than max_time_gap frames
%apart and within a radius of dist_threshold

%input:
%trk is trajectory array as a column array [x,y,I,t,traj index] where each
%row is a datapoint, order by increasing time

%outputs:
%matches: first column is traj idx of junction start, 
%other columns are indices of matched trajectory ends

%dr: first column is 0, 
%other columns are the [endpoint - startpoint] distances corresponding to the indices in "matches"

%dt: first column is 0, 
%other columns are the [endpoint - startpoint] time gaps corresponding to the indices in "matches"


%generate array of starts and ends
traj_idx = unique(trk(:,5));
ntraj = numel(traj_idx);

trk_start = zeros(ntraj,size(trk,2));
trk_end = zeros(ntraj,size(trk,2));
for i=1:ntraj
    tmptraj = trk(trk(:,5)==traj_idx(i),:);
    trk_start(i,:) = tmptraj(1,:);
    trk_end(i,:) = tmptraj(end,:);    
end

%find all possible matches between ends and starts
matches = [];%first column is traj idx of junction end, other columns are indices of matched trajectory starts
dr = [];%first column is 0, other columns are the [endpoint - startpoint] distances corresponding to the indices in "matches"
dt = [];%first column is 0, other columns are the [endpoint - startpoint] time gaps corresponding to the indices in "matches"
for i=1:ntraj
    
    %select time points that are compatible in time
    curdt = trk_start(i,4) - trk_end(:,4);
    poss_matches_idx = logical( ( curdt>0 ) .* ( curdt<= max_time_gap ) );
    cur_trk_end = trk_end(poss_matches_idx,:);
    curdt = curdt(poss_matches_idx);
    
    %select time points that are compatible in space
    curdr = sum((repmat(trk_start(i,1:2),size(cur_trk_end,1),1) - cur_trk_end(:,1:2) ).^2,2);
    poss_matches_idx = logical(curdr<=dist_threshold^2);
    poss_matches = cur_trk_end(poss_matches_idx,:);
    curdr = curdr(poss_matches_idx);
    curdt = curdt(poss_matches_idx);
    if ~isempty(poss_matches)
        matches(i,1) = traj_idx(i);
        dr(i,1) = 0;
        dt(i,1) = 0;
        matches(i,2:(size(poss_matches,1)+1)) = poss_matches(:,5)';
        dr(i,2:(size(poss_matches,1)+1)) = curdr;    
        dt(i,2:(size(poss_matches,1)+1)) = curdt; 
    else
        matches(i,1) = traj_idx(i);
        dr(i,1) = 0;
        dt(i,1) = 0;
    end    
end

end