function plot_trajectories_in_figure(loc,fh,varargin)

%loc should be a file containing trajectories formatted as an array:
%[x, y, ..., traj index]

%optional argument is a character coresponding to the color of the curves
%default is random assignment of colors, a separate color per trajectory

traj = unique(loc(:,5));
ntraj = numel(traj);


figure(fh);
if ~ishold
    hold;
end

R = 0:0.1:1;

for i=1:ntraj
    
    curR = R(mod(i,numel(R))+1);
    curG = R(mod(ceil(i/2),numel(R))+1);
    curB = R(mod(ceil(i/3),numel(R))+1);
    
    spots = loc(loc(:,end) == traj(i),:);
    spots = spots(~isnan(spots(:,1)),:);
    if size(spots,2)>=5
        tmin = num2str(min(spots(:,4)));
        tmax = num2str(max(spots(:,4)));
    end
    if numel(varargin) == 0
        plot(spots(:,1),spots(:,2),'Color',[curR,curG,curB],'DisplayName',['trajectory ',num2str(traj(i)),...
            '; tmin = ',tmin,'; tmax = ',tmax]);
    else
        if ischar(varargin{1})
            trajcolor = varargin{1};
            plot(spots(:,1),spots(:,2),'Color',trajcolor,'DisplayName',['trajectory ',num2str(traj(i)),...
            '; tmin = ',tmin,'; tmax = ',tmax]);
        else
            plot(spots(:,1),spots(:,2),'Color',[curR,curG,curB],'DisplayName',['trajectory ',num2str(traj(i)),...
            '; tmin = ',tmin,'; tmax = ',tmax]);
        end
    end
end



end