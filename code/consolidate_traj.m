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