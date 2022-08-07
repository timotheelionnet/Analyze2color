function D = fit_cdf_array(cdf,disp,time_interval)

%cdf: cumulative displacement histogram cdf is formatted as follows:
%col 1: trajectory index 
%col 2: value of the time interval (in frames) of the histogram to which the cdf value belongs
%col 3: cdf displacement value (histogram x), in um
%col 4: cdf displacement frequency (histogram x), in um
%col 5: number of datapoints that were used to build the cdf

% disp is an array that contains the displacements for all time intervals
% col1: traj #,
% col2: time interval in frames
% col3: dx1, in um 
% col4: dy1, in um
% col5: dr1

%time interval: value of the time interval between frames in s

% D is an array that contains the displacements for all time intervals
% col1: traj #,
% col2: time interval in frames
% col2: computed diffusion coefficient, in um^2/s
% col3: dy1
% col4: dr1

D = [];

traj_idx = unique(cdf(:,1));

for j=1:numel(traj_idx)
    curcdf = cdf(cdf(:,1)== traj_idx(j),:);
    tvalues = unique(curcdf(:,2));
    for i=1:numel(tvalues)
        curcdf2 = curcdf(curcdf(:,2) == tvalues(i),:);
        t = tvalues(i)*time_interval;
        
        if curcdf2(1,5)>3 
            %if large number of datapoints were used for the histogram,
            %proceed to NLLS fit
            
            %use sqrt(N) as statistical weight
            weights = 1./sqrt(curcdf2(1,5)*curcdf2(:,4));
            weights(curcdf2(:,4)==0) = 1;
            
            %NLLS fit
            [D_tmp,gof_tmp] = fit_cdf(curcdf2(:,3),curcdf2(:,4),weights,t);
            rsq_tmp = gof_tmp.rsquare;
        else
            
            %if only 3 datapoints or less, fit will be poorly constrained;
            %compute D using the avg displacement instead
            disp_idx = logical((disp(:,1) == traj_idx(j)).*(disp(:,2) ==tvalues(i)));
            spots = disp(disp_idx,:);
            D_tmp = (nanmean(spots(:,5)))^2/(4*t);
            rsq_tmp = NaN;
        end
        
        D = [D;traj_idx(j),tvalues(i),D_tmp,curcdf2(1,5),rsq_tmp];
    end
end


end