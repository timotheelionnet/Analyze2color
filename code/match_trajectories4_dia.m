function matches = match_trajectories4_dia(trk1,trk2,minpts,max_spot_dist)

%% inputs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%trk1, trk2 are lists of trajectories for each channel: 
%it includes spots positions, frame number and trajectory assigned, formated as follows:
%col 1: x
%col 2: y
%col 3: I
%col 4: frame number
%col 5: trajectory number

traj1 = unique(trk1(:,end));
ntraj1 = numel(traj1);

%compile in a cell array the list of channel 2 trajectories that match each
%channel 1 trajectory.
%each entry of the cell contains an array of values which are the indices of the
%matched trajectories
%match criterion: more than minpts positions are within max_spot_dist

for i=1:ntraj1
    matches(i,1) = match_traj_idx(trk1,trk2,traj1(i),minpts,max_spot_dist); 
end

end

function matches = match_traj_idx(trk1,trk2,idx,minpts,max_spot_dist)
    %build structure
    matches = struct('idx1',{},'idx2',{},'npts',{},'avgdist',{},'avgcolocdist',{},'distances',{},'colocpts',{});
    
    %select the ith trajectory from channel 1
    curtraj1 = trk1(trk1(:,end)==idx,:);
    tmin = nanmin(curtraj1(:,4));
    tmax = nanmax(curtraj1(:,4));
    
    %filter out the trajectories in channel 2 that do not overlap in time
    curtraj2 = trk2( logical((trk2(:,4)>=tmin).*(trk2(:,4)<=tmax)),:);
    trajidx2 = unique(curtraj2(:,end));
    %save('C:\junk\test.txt','curtraj2','-ascii');
    
    %consolidate each channel 2 trajectory so that it has datapoints from tmin to
    %tmax (fill missing datapoints w/ NaN)
    curtraj1 = consolidate_traj(curtraj1,tmin,tmax);
    %save('C:\junk\test1.txt','curtraj1','-ascii');
    curtraj2 = consolidate_traj(curtraj2,tmin,tmax);
    %save('C:\junk\test2.txt','curtraj2','-ascii');
    
    %replicate array 1 for easier calculations
    trajidx2 = unique(curtraj2(:,end));
    ntraj2 = numel(trajidx2);   
    curtraj1 = repmat(curtraj1,ntraj2,1);
    
    %compute distances
    dist = sqrt(sum((curtraj1(:,1:2) - curtraj2(:,1:2)).^2,2));
    
    %extract score
    coloc = dist <= max_spot_dist;
    matches(1).idx1 = idx;
    matches.idx2 = [];
    matches.npts = [];
    matches.avgdist = [];
    matches.avgcolocdist = [];
    matches.distances = [];
    matches.colocpts = [];
    for i = 1:ntraj2
        if sum(coloc(curtraj2(:,end) == trajidx2(i))) >= minpts
            matches.idx2 = [matches.idx2,trajidx2(i)];
            matches.npts = [matches.npts,sum(coloc(curtraj2(:,end) == trajidx2(i)))];
            matches.avgdist = [matches.avgdist,nanmean(dist(curtraj2(:,end) == trajidx2(i)))];
            matches.avgcolocdist = [matches.avgcolocdist,...
                nanmean(dist(logical(coloc.*(curtraj2(:,end) == trajidx2(i)))))];    
            matches.colocpts = [matches.colocpts,coloc(curtraj2(:,end) == trajidx2(i))];
            matches.distances = [matches.distances,dist(curtraj2(:,end) == trajidx2(i))];
        end
    end
end

function trko = consolidate_traj(trki,tmin,tmax)
    %pad each trajectory in trajectory list trki so that
    %there are datapoints for each frame between tmin and tmax (inclusive)
    %missing datapoints are replaced by NaNs
    trajidx = unique(trki(:,end));
    ntraj = numel(trajidx);
    trko = [];
    for i=1:ntraj
        curtrki = trki(trki(:,end)==trajidx(i),:);
        if size(curtrki,1) <= tmax - tmin +1
            curtraj = NaN(tmax - tmin +1,size(trki,2));
            curtraj(:,4) = tmin:tmax;
            curtraj(:,5) = trajidx(i);
            for j=1:size(curtrki,1)
                curtraj(curtraj(:,4) == curtrki(j,4),:) = curtrki(j,:);
            end
        else
            curtraj = curtrki;
        end
        trko = [trko;curtraj];
    end
end