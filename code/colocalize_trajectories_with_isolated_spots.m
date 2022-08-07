function [matrix,neg_list1,neg_list2,...
        coloc_trk1,coloc_trk2,neg_trk1,neg_trk2,trk1,trk2] = ...
        colocalize_trajectories_with_isolated_spots(...
        loc1,loc2,trk1,trk2,...
        neg_list1,neg_list2,coloc_trk1,coloc_trk2,neg_trk1,neg_trk2,...
        spot_coloc_zscore_thresh,mininum_number_of_isolated_colocalizations,max_spot_dist,parfname)

%detects for each orphan trajectory in neg_trk1 array, whether lone spots in loc2 array colocalize with it.
%if the number of detected spots is statistically significant, generate a
%channel 2 trajectory based on these spots and updates the trajectory arrays.
%Does the same thing with neg_trk2 array vs loc1 spots.

%Note that the recreated trajectories are given ID indices equal to the
%opposite of the trajectory with which they colocalize.

%Parameters:
%spot_coloc_zscore_thresh: z score threshold for a statistically
%significant colocalization frequency, based on the spatiotemporal
%frequency of spot detection
%mininum_number_of_isolated_colocalizations: minimum number of spot
%detections in one trajectory to call a colocalization event
%max_spot_dist: distance between channels sufficient to call a
%colocalization event
%parfname: name of the file holding the parameters of the detected spots
%including the size of the image, used to compute the detection frequency.

%the column in which the trajectory indices are stored
traj_idx_col_num = 5;     

%compute the area of each frame in pix^2
%it is used to compuote the probability of random matches
if exist(parfname,'file')
    [xmax,ymax] = set_exclusion_limits(parfname,0);
else
    xmax = max([max(loc1(:,1)),max(loc2(:,1))]);
    ymax = max([max(loc1(:,2)),max(loc2(:,2))]);
end

%colocalize orphan trajectories in channel 1 with isolated spots in
%channel 2
[coloc_trk1,coloc_trk2,neg_list1,neg_trk1,trk2] = match_trajectory_with_orphan_spots(...
    neg_trk1,loc2,neg_list1,coloc_trk1,coloc_trk2,trk2,...
    traj_idx_col_num,max_spot_dist,spot_coloc_zscore_thresh,...
    mininum_number_of_isolated_colocalizations,xmax,ymax);    
    
%colocalize orphan trajectories in channel 2 with isolated spots in
%channel 1
[coloc_trk2,coloc_trk1,neg_list2,neg_trk2,trk1] = match_trajectory_with_orphan_spots(...
    neg_trk2,loc1,neg_list2,coloc_trk2,coloc_trk1,trk1,...
    traj_idx_col_num,max_spot_dist,spot_coloc_zscore_thresh,...
    mininum_number_of_isolated_colocalizations,xmax,ymax);    

%update colocalization matrix
coloc_list1 = coloc_trk1(:,traj_idx_col_num);
coloc_list1(isnan(coloc_list1)) = Inf;
coloc_list2 = coloc_trk2(:,traj_idx_col_num);
coloc_list2(isnan(coloc_list2)) = Inf;
matrix = unique([coloc_list1, coloc_list2],'rows');
matrix(matrix == Inf) = NaN;
neg_list1 = neg_list1(~ismember(neg_list1,matrix(:,1)));
neg_list2 = neg_list2(~ismember(neg_list2,matrix(:,2)));

end



function [coloc_trk1,coloc_trk2,neg_list1,neg_trk1,trk2] = match_trajectory_with_orphan_spots(...
    neg_trk1,loc2,neg_list1,coloc_trk1,coloc_trk2,trk2,...
    traj_idx_col_num,max_spot_dist,zscore_thresh,minpts,xmax,ymax)
  
orphans1_idx = unique(neg_trk1(:,traj_idx_col_num));

%the maximum colocalization event index (all new colocalization events will
%be given sequential indices starting from that value +1)
nmatchtraj = max(coloc_trk1(:,7));

for i=1:numel(orphans1_idx) %loop over orphan trajectories and find potential matching spots
    
    %select current orphan trajectory
    tmptrk = neg_trk1(neg_trk1(:,traj_idx_col_num)==orphans1_idx(i),:);
    
    %crop out from the spots array the spots that do not overlap temporally
    %with the current trajectory
    tmpspots = loc2(logical( (loc2(:,4)>=tmptrk(1,4)).*(loc2(:,4)<=tmptrk(end,4))) ,:);
    
    %loop over each datapoint from the trajectory try to find matches
    match_found = 0;
    nmatchspots = 0;
    coloc_trk1_tmp = [];
    coloc_trk2_tmp = [];
    for j=1:size(tmptrk,1)
        %extract spots detected at the current frame
        tmpspots2 = tmpspots(tmpspots(:,4) == tmptrk(j,4),:);

        %compute distances with the current datapoint form the orphan
        %trajectory
        dist = sqrt( (tmpspots2(:,1) - tmptrk(j,1)).^2 + (tmpspots2(:,2) - tmptrk(j,2)).^2 );

        %find closest spot within accepted radius
        if sum(dist <= max_spot_dist) >= 1
            [~,idx] = min(dist);
        else
            idx = [];
        end

        %fill in trajectory array if a spot was found
        if ~isempty(idx)
            match_found = match_found +1;
            if match_found == 1
                nmatchtraj = nmatchtraj + 1;
                coloc_trk1_tmp = tmptrk;
                coloc_trk1_tmp(:,6) = 0;
                coloc_trk1_tmp(j,6) = 1;
                coloc_trk1_tmp(:,7) = nmatchtraj;

                coloc_trk2_tmp = NaN*zeros(size(coloc_trk1_tmp));
                coloc_trk2_tmp(:,4) = coloc_trk1_tmp(:,4);
                %give the colocalized spot a trajectory ID that is equal to the opposite the ID of the corresponding
                %trajectory
                coloc_trk2_tmp(:,5) = -coloc_trk1_tmp(:,5);
                coloc_trk2_tmp(:,6:7) = coloc_trk1_tmp(:,6:7);
                coloc_trk2_tmp(j,1:3) = tmpspots2(idx,1:3);
            else
                coloc_trk1_tmp(j,6) = 1;
                coloc_trk1_tmp(j,7) = nmatchtraj;
                coloc_trk2_tmp(j,6:7) = coloc_trk1_tmp(j,6:7);
                coloc_trk2_tmp(j,1:4) = tmpspots2(idx,1:4);
            end
            nmatchspots = nmatchspots + 1;
        end
    end
    
    if match_found>=1
        %if matches were found, test whether the colocalization frequency is statistically significant
        is_significant = compute_match_significance(...
            coloc_trk1_tmp,loc2,zscore_thresh,minpts,max_spot_dist,xmax,ymax);
        
        %if significant, add current matches to the colocalized arrays
        %add newly created trajectory to the channel 2 list
        if is_significant
            coloc_trk1 = [coloc_trk1;coloc_trk1_tmp];
            coloc_trk2 = [coloc_trk2;coloc_trk2_tmp];
            trk2 = [trk2;coloc_trk2_tmp(:,1:5)];
        end
    end
end    
    
%remove positive hits from the orphan trajectory list  
neg_trk1 = neg_trk1(~ismember(neg_trk1(:,traj_idx_col_num),coloc_trk1(:,traj_idx_col_num)),:);    
end



