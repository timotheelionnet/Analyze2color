function is_significant = compute_match_significance(trk1,loc2,zscore_thresh,minpts,max_spot_dist,xmax,ymax)
%compute the statistical significance of the matches found:
%could they have happened by random occurence of isolated spots?

%compute the probability 
ImgArea = (xmax - 1)*(ymax - 1);
DetectionCrossSection = pi*max_spot_dist^2;
FalsePosProb = DetectionCrossSection / ImgArea;

%compute # of spots detected / frame in each channel
frames = unique(loc2(:,4));
for i=1:numel(frames)
    nspots2(i) = sum(loc2(:,4) == frames(i));
end

%compare the number of co-detected spots to what is expected by chance
%quantify it by a "z-score" 
NumDetections = sum(trk1(:,6)==1);

%compute the cumulated probability of spot detection over the
%trajectory
frames2 = ismember(frames,trk1(:,4));
DetProb = sum( FalsePosProb * nspots2(frames2) );

%compute z-score equivalent for Poisson distributed events
%zscore_thresh = 2 corresponds to events that are 5% or less likely to
%happen by chance
%zscore_thresh = 3: 1.5% or less 
%zscore_thresh = 4: 0.5% or less 
if DetProb == 0
    zscore = Inf;
else
    zscore = (NumDetections - DetProb)/sqrt(DetProb);
end

is_significant = (zscore >= zscore_thresh).*(NumDetections >= minpts);
end