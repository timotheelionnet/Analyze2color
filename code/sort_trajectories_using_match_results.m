function [matrix,neg_list1,neg_list2,coloc_trk1,coloc_trk2,neg_trk1,neg_trk2] = ...
    sort_trajectories_using_match_results(trk1,trk2,matches,max_spot_dist)

%coloc_trk1/2: list of datapoints that belong to colocalized trajectories
%each row is one datapoint from one time point
%columns:
%1:x
%2:y
%3:Intensity (or place holder)
%4:time
%5:trajectory 1 index
%6:colocalization flag (1 if the datapoint satisfies the colocalication
%criterion; 0 otherwise)
%7: colocalized trajectory index
%IMPOTANT NOTE: trajectories are always stored entirely, i.e. one trajectory from channel 1 might
%be stored twice if it is matched by 2 trajectories from channel 2

%neg_trk1/2: list of datapoints that belong to orphan trajectories
%each row is one datapoint from one time point
%columns:
%1:x
%2:y
%3:Intensity (or place holder)
%4:time
%5: trajectory 1 index


%compute the list of paired trajectories
%each row lists a pair of indices (the IDs of the channel 1 and channel 2
%trajectories that colocalize)
matrix = [];
for i=1:size(matches,1)
    if ~isempty(matches(i).idx2)
        curmatrix = [matches(i).idx1*ones(numel(matches(i).idx2),1),(matches(i).idx2)'];
        matrix = [matrix;curmatrix];
    end    
end

%compile colocalized trajectories
coloc_trk1 = [];
coloc_trk2 = [];
ncoloc = 1;
for i=1:size(matches,1)
    curtrk1 = trk1(trk1(:,end) == matches(i).idx1,:);
    for j=1:numel(matches(i).idx2)
        curtrk2 = trk2(trk2(:,end) == matches(i).idx2(j),:);
        
        tmin = min([min(curtrk1(:,4)),min(curtrk2(:,4))]);
        tmax = max([max(curtrk1(:,4)),max(curtrk2(:,4))]);
        
        trk1o = consolidate_traj(curtrk1,tmin,tmax);
        trk2o = consolidate_traj(curtrk2,tmin,tmax);
        
        colocalization_flag = sqrt( sum( (trk1o(:,1:2) - trk2o(:,1:2)).^2 ,2)  );
        colocalization_flag = colocalization_flag < max_spot_dist;
        trk1o(:,6) = colocalization_flag;
        trk1o(:,7) = ncoloc;
        trk2o(:,6) = colocalization_flag;
        trk2o(:,7) = ncoloc;
        coloc_trk1 = [coloc_trk1;trk1o];
        coloc_trk2 = [coloc_trk2;trk2o]; 
        
        ncoloc = ncoloc+1;
        
    end
end

%compile orphan trajectories
coloc_idx1 = unique(matrix(:,1));
all_idx1 = unique(trk1(:,end));
neg_trk1 = [];
neg_list1 = [];
for i=1:numel(all_idx1)
    k = find(all_idx1(i)==coloc_idx1,1);
    if isempty(k)
        neg_list1 = [neg_list1;all_idx1(i)];
        neg_trk1 = [neg_trk1;trk1(trk1(:,end)==all_idx1(i),:)];
    end    
end

coloc_idx2 = unique(matrix(:,2));
all_idx2 = unique(trk2(:,end));
neg_trk2 = [];
neg_list2 = [];
for i=1:numel(all_idx2)
    k = find(all_idx2(i)==coloc_idx2,1);
    if isempty(k)
        neg_list2 = [neg_list2;all_idx2(i)];
        neg_trk2 = [neg_trk2;trk2(trk2(:,end)==all_idx2(i),:)];
    end    
end

end