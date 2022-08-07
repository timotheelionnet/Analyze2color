function cdf = compute_cdf_single_traj_array2(disp,npts)

% disp is an array that contains the displacements for all time intervals
% col1: traj idx
% col2: time interval in frames
% col3: dx
% col4: dy
% col5: dr

%npts is the number of bins each histogram is set to have

%output: cdf is an array that contains all the cumulative displacement distributions
%cdf is formatted as follows:
%col 1: trajectory index 
%col 2: value of the time interval (in frames) of the histogram to which the cdf value belongs
%col 3: cdf displacement value (histogram x, um)
%col 4: cdf displacement frequency (histogram y, normalized to 1)
%col 5: number of datapoints that were used to build the cdf

%update from v1: translated the x coordinate so that the data points fall
%on the right end of each bin rather than at the bin center; this avoids
%bin size dependent artefacts.


cdf = [];
tvalues = unique(disp(:,2));

for i=1:numel(tvalues)    
    spotstmp = disp(disp(:,2)==tvalues(i),:);
    traj = unique(spotstmp(:,1));
    for j=1:numel(traj)  
        spots = spotstmp(spotstmp(:,1) == traj(j),:);
        spots = spots(~isnan(spots(:,2)) ,:  );

        n = []; x = [];
        ndata = numel(spots(:,5));
        if ndata >= 1
            dx = ( max(spots(:,5)) - min(spots(:,5)) ) /npts;
            if dx == 0
                if max(spots(:,5)) > 0
                    [n,x] = hist(spots(:,5),0.4999*max(spots(:,5)):max(spots(:,5)):1.4999*max(spots(:,5)));
                    %displace the datapoint to the right end of the bin;
                    %otherwise using a datapoint centered w/ respect to the
                    %bin with a cumulative distribution is inaccurate and
                    %generate bin size artefacts
                    x = x + max(spots(:,5))/2;
                end
            else
                [n,x] = hist(spots(:,5),min(spots(:,5))-dx:dx:max(spots(:,5)) + dx);
                x = x+dx/2;
            end
            n = cumsum(n)/sum(n);
        end
        cdf_tmp = [repmat(traj(j),size(x')),repmat(i,size(x')),x',n',repmat(ndata,size(x'))];  
        cdf = [cdf;cdf_tmp];
    end
    
end


end