function screened_matches = screen_trajectories(matches,traj_threshold_method, thresh_value)
    
    %build useful arrays    
    idx1 = [matches.idx1];   %array that stores the list of matched trajectories from channel 1     
    matrix = [];    %array that stores the full list of trajectory pairs    
    avgdist = [];   %array that stores the list of avg distances between trajectory pairs
    for i=1:size(matches,1)
        avgdist = [avgdist,matches(i).avgdist];
        
        if ~isempty(matches(i).idx2)
            curmatrix = [matches(i).idx1*ones(numel(matches(i).idx2),1),(matches(i).idx2)'];
            matrix = [matrix;curmatrix];
        end    
    end

    %screen average distances
    switch traj_threshold_method
        case 'manual'
            positives = avgdist < thresh_value;
        case 'auto'
            positives = otsu(avgdist)==1;
    end
    
    %cleanup matches structure to retain positives only
    positives_matrix = matrix(positives,1:2);
    
    screened_matches = struct('idx1',{},'idx2',{},'npts',{},'avgdist',{},'avgcolocdist',{},'colocpts',{});
    screened_matches(1).idx1 = [];
    screened_matches.idx2 = [];
    screened_matches.npts = [];
    screened_matches.avgdist = [];
    screened_matches.avgcolocdist = [];
    screened_matches.colocpts = [];
    screened_matches.distances = [];
    for i=1:size(matches,1)
        screened_matches(i,1).idx1 = matches(i,1).idx1;
    end     
    
    for i =1:size(positives_matrix,1)
        curidx1 = find(idx1 == positives_matrix(i,1));
        curidx2 = find(matches(curidx1).idx2 == positives_matrix(i,2));
        
        screened_matches(curidx1,1).idx1 = positives_matrix(i,1);
        screened_matches(curidx1,1).idx2 = [screened_matches(curidx1).idx2,positives_matrix(i,2)];
        
        screened_matches(curidx1,1).npts = [screened_matches(curidx1).npts,matches(curidx1).npts(curidx2)];
        screened_matches(curidx1,1).avgdist = [screened_matches(curidx1).avgdist,matches(curidx1).avgdist(curidx2)];
        screened_matches(curidx1,1).avgcolocdist = [screened_matches(curidx1).avgcolocdist,matches(curidx1).avgcolocdist(curidx2)];
        screened_matches(curidx1,1).colocpts = [screened_matches(curidx1).colocpts,matches(curidx1).colocpts(:,curidx2)];
        screened_matches(curidx1,1).distances = [screened_matches(curidx1).distances,matches(curidx1).distances(:,curidx2)];
        
    end
    
    %remove empty entries
    
    %screened_matches = screened_matches(~cellfun(@isempty,{screened_matches.idx1}));
    
end