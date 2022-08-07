function cdf = compute_cdf_single_traj_array(disp,npts)

% disp is an array that contains the displacements for all time intervals
% col1: traj #,
% col2: time interval in frames
% col2: dx1, 
% col3: dy1, 
% col4: dr1

%npts is the number of bins each histogram is set to have

%output: cdf is an array that contains all the cumulative displacement distributions
%cdf is formatted as follows:
%col 1: trajectory index 
%col 2: value of the time interval (in frames) of the histogram to which the cdf value belongs
%col 3: cdf displacement value (histogram x, um)
%col 4: cdf displacement frequency (histogram y, normalized to 1)
%col 5: number of datapoints that were used to build the cdf

cdf = [];
tvalues = unique(disp(:,2));

for i=1:numel(tvalues)    
    spotstmp = disp(disp(:,2)==tvalues(i),:);
    traj = unique(spotstmp(:,1));
    for j=1:numel(traj)  
        spots = spotstmp(spotstmp(:,1) == traj(j),:);
        spots = spots(~isnan(spots(:,2)) ,:  );

        n = []; x = [];
        ndata = numel(spots(:,4));
        if ndata >= 1
            dx = ( max(spots(:,4)) - min(spots(:,4)) ) /npts;
            if dx == 0
                if max(spots(:,4)) > 0
                    [n,x] = hist(spots(:,4),0.4999*max(spots(:,4)):max(spots(:,4)):1.4999*max(spots(:,4)));
                end
            else
                [n,x] = hist(spots(:,4),min(spots(:,4))-dx:dx:max(spots(:,4)) + dx);
            end
            n = cumsum(n)/sum(n);
        end
        cdf_tmp = [repmat(traj(j),size(x')),repmat(i,size(x')),x',n',repmat(ndata,size(x'))];  
        cdf = [cdf;cdf_tmp];
    end
    
end


end