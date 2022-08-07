function [cdf1,cdf2,ndata1,ndata2] = compute_cdf(disp,npts)
%input: disp is a cell array that contains the displacements between points along trajectories 
%each entry corresponds to a time interval,
%e.g. dispco{10} contains the displacements at intervals of 10 frames
%each dispco entry is an array formatted as follows:
% col1: coloc traj #, 
% col2: dx1, 
% col3: dy1, 
% col4: dr1, 
% col5: dx2,
% col6: dy2,
% col7: dr2

%output: cdf is a  cell array that contains all the cumulative displacement distributions
%cdf{i,j} store the cumulative distribution of the trajectory with index j,
%sampled with an interval i times the unit interval.
%cdf{i,j} is formatted as follows:
%col 1: trajectory index (i.e. j)
%col 2: cdf channel 1 x
%col 3: cdf channel 1 y
%col 4: cdf channel 2 x
%col 5: cdf channel 2 y

%ndata1 is the number of data points used to construct the coresponding cdf
%ndata1(i,j) is the number of intervals of length i within traj j

cdf1 = cell('');
cdf2 = cell('');
ndata1 = [];
ndata2 = [];
for i=1:numel(disp)
    
    spotstmp = disp{i,1};
    if ~isempty(spotstmp)
        traj = unique(spotstmp(:,1));
        for j=1:numel(traj)  
            spots = spotstmp(spotstmp(:,1) == traj(j),:);
            spots1 = spots(~isnan(spots(:,2)) ,:  );

            n1 = []; x1 = [];
            ndata1(i,traj(j)) = numel(spots1(:,4));
            if ndata1(i,traj(j))>= 1
                dx = ( max(spots1(:,4)) - min(spots1(:,4)) ) /npts;
                if dx == 0
                    if max(spots1(:,4)) > 0
                        [n1,x1] = hist(spots1(:,4),0.4999*max(spots1(:,4)):max(spots1(:,4)):1.4999*max(spots1(:,4)));
                    end
                else
                    [n1,x1] = hist(spots1(:,4),min(spots1(:,4))-dx:dx:max(spots1(:,4)) + dx);
                end
                n1 = cumsum(n1)/sum(n1);
            end
            cdf1{i,traj(j)} = [repmat(traj(j),size(x1')),x1',n1'];

            spots2 = spots(~isnan(spots(:,5)),:   );
            n2 = []; x2 = [];
            ndata2(i,traj(j)) = numel(spots2(:,7));
            if ndata2(i,traj(j))>= 1
                dx = ( max(spots2(:,7)) - min(spots2(:,7)) ) /npts;
                if dx == 0
                    if max(spots2(:,7)) > 0
                        [n2,x2] = hist(spots2(:,7),0.4999*max(spots2(:,7)):max(spots2(:,7)):1.4999*max(spots2(:,7)));
                    end
                else
                    [n2,x2] = hist(spots2(:,7),min(spots2(:,7))-dx:dx:max(spots2(:,7)) + dx);
                end
                n2 = cumsum(n2)/sum(n2);
            end

            cdf2{i,traj(j)} = [repmat(traj(j),size(x2')),x2',n2'];
        end
    end
end

end