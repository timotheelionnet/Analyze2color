function Analyze2color()

%analysis parameters

params = set_parameters();
if ~isstruct(params)
    return
else
    pix_size = params.pix_size;
    time_interval = params.time_interval;
    exclusion_radius = params.exclusion_radius;
    minpts = params.min_overlap;
    max_center_dist = params.max_center_dist;
    max_spot_dist = params.max_spot_dist;
end
% pix_size = 0.107;
% time_interval = 0.01;
% exclusion_radius = 5;
%min_overlap = 3;


%% Load data
user_entry = 1;
if user_entry
    
%% load files Channel 1
[loc1name,loc1dirname,fidx] = uigetfile('*.loc','Load First Channel Loc Result','\');
if fidx == 0
    return
end

[~,loc1name_no_ext,~] = fileparts(loc1name);
trk1name = [loc1name_no_ext,'.trk'];
listing1 = dir(loc1dirname);

trkfound = 0;
for i=1:numel(listing1)
    if strcmp(listing1(i).name,trk1name)
        trkfound = 1;
        trk1dirname = loc1dirname;
    end
end

if ~trkfound
    [trk1name,trk1dirname,fidx] = uigetfile('*.trk','Load First Channel Trk Result',loc1dirname);
    if fidx == 0
        return
    end    
end


par1name = [loc1name_no_ext,'.par'];
parfound = 0;
for i=1:numel(listing1)
    if strcmp(listing1(i).name,par1name)
        parfound = 1;
        par1dirname = loc1dirname;
    end
end

if ~parfound
    [par1name,par1dirname,fidx] = uigetfile('*.par','Load Localization Parameter File Channel 1',loc1dirname);
    if fidx == 0
        return
    end    
end

%% load files Channel 2
[loc2name,loc2dirname,fidx] = uigetfile('*.loc','Load Second Channel Loc Result',loc1dirname);
if fidx == 0
    return
end

[~,loc2name_no_ext,~] = fileparts(loc2name);
trk2name = [loc2name_no_ext,'.trk'];
listing2 = dir(loc2dirname);

trkfound = 0;
for i=1:numel(listing2)
    if strcmp(listing2(i).name,trk2name)
        trkfound = 1;
        trk2dirname = loc2dirname;
    end
end

if ~trkfound
    [trk2name,trk2dirname,fidx] = uigetfile('*.trk','Load Second Channel Trk Result',loc2dirname);
    if fidx == 0
        return
    end    
end


save_dirname = uigetdir(loc2dirname,'Saving Directory') ;
else
    loc1name = 'batch_spot_detection.loc';
    loc1dirname = 'Z:\Timothee Lionnet\data\Jeff Trick Data\120618\MS2\registered';
    loc1name_no_ext = 'batch_spot_detection';
    trk1name = 'MS2_batch_spot_detection_2_1_5.trk';
    trk1dirname = 'Z:\Timothee Lionnet\data\Jeff Trick Data\120618\MS2\registered';
    
    loc2name = 'batch_spot_detection.loc';
    loc2dirname = 'Z:\Timothee Lionnet\data\Jeff Trick Data\120618\PP7\registered';
    loc2name_no_ext = 'batch_spot_detection';
    trk2name = 'PP7_batch_spot_detection_2_1_5.trk';
    trk2dirname = 'Z:\Timothee Lionnet\data\Jeff Trick Data\120618\PP7\registered';
    
    save_dirname = 'D:\junk\';
end

loc1 = load(fullfile(loc1dirname,loc1name));
trk1 = load(fullfile(trk1dirname,trk1name));
loc2 = load(fullfile(loc2dirname,loc2name));
trk2 = load(fullfile(trk2dirname,trk2name));


%% match trajectories
boundaries = exclusion_radius*ones(1,4);
[boundaries(2),boundaries(4)] = set_exclusion_limits(fullfile(par1dirname,par1name),exclusion_radius);


fh = dispwin('Status','Matching Trajectories...');
[coloc_sum,coloc_traj,undecided_sum,~,neg_sum,neg_traj,~,trk1sum,~,trk2sum] = ...
    match_trajectories3(trk1,trk2,loc1,loc2,boundaries,minpts,max_center_dist,max_spot_dist);


%% compute displacement histograms and diffusion properties
fh = dispwin('Status','Computing Diffusion Kinetics...',fh);

[coloc_res,coloc_dispco,coloc_dispall,coloc_disp] = compute_displacements_and_fit_figure2(coloc_traj,pix_size,time_interval,'Colocalized Trajectories Diffusion Coefficient','co');

[neg_res,neg_dispco,neg_dispall,neg_disp] = compute_displacements_and_fit_figure2(neg_traj,pix_size,time_interval,'Isolated Trajectories Diffusion Coefficient','single');
%save('D:\junk\logJeff.mat','coloc_sum','coloc_res','coloc_traj','coloc_dispco','coloc_dispall');

%% Launch Trajectory Browsers
mine_traj3(coloc_res,coloc_sum,coloc_traj,time_interval,'co','Colocalized Trajectories Browser');

neg_sum1 = neg_sum(neg_sum(:,1)~=-1,:);
neg_traj1 = neg_traj(neg_traj(:,8)~=-1,:);
mine_traj3(neg_res,neg_sum1,neg_traj1,time_interval,'single','Isolated Trajectories Browser: Channel 1');

neg_sum2 = neg_sum(neg_sum(:,2)~=-1,:);
neg_traj2 = neg_traj(neg_traj(:,9)~=-1,:);
mine_traj3(neg_res,neg_sum2,neg_traj2,time_interval,'single','Isolated Trajectories Browser: Channel 2');


%% save results
% data
if isempty(coloc_sum) || isempty(coloc_res.D1c) || isempty(coloc_res.D2c)
    coloc_list = 'no trajectories found';
else
    coloc_list = coloc_sum(:,11);
    
    coloc_list(:,2) = coloc_res.D1c(1,coloc_sum(:,11))';
    coloc_list(:,3) = coloc_res.D2c(1,coloc_sum(:,11))';
end
save(fullfile(save_dirname,[loc1name_no_ext,'_',loc2name_no_ext,'.co']),'coloc_list','-ascii');

if isempty(neg_sum) || isempty(neg_res.D1c) || isempty(neg_res.D2c)
    neg_list = 'no trajectories found';
else
    neg_list = neg_sum(:,11);
    neg_list(:,2) = neg_res.D1s(1,neg_sum(:,11))';
    neg_list(:,3) = neg_res.D2s(1,neg_sum(:,11))';
end
save(fullfile(save_dirname,[loc1name_no_ext,'_',loc2name_no_ext,'.neg']),'neg_list','-ascii');

save(fullfile(save_dirname,[loc1name_no_ext,'_',loc2name_no_ext,'.res']),'coloc_res','neg_res','coloc_sum','neg_sum','coloc_traj','neg_traj','trk1sum','trk2sum','coloc_disp','neg_disp');

% summary file
summary_filename = fullfile(save_dirname,[loc1name_no_ext,'_match_results.txt']);
fid1 = fopen(summary_filename,'w');

NL = sprintf('\r\n');

file_str = ['Loc File Channel 1: ',fullfile(loc1dirname,loc1name),NL,...
    'Trk File Channel 1: ',fullfile(trk1dirname,trk1name),NL,...
    'Loc File Channel 2: ',fullfile(loc2dirname,loc2name),NL,...
    'Trk File Channel 2: ',fullfile(trk2dirname,trk2name),NL,NL];


params_str = ['Pixel size: ',num2str(pix_size),' um',NL,...
    'Time between frames: ',num2str(time_interval),' s',NL,...
    'Exclusion radius: ',num2str(exclusion_radius),' pix',NL,...
    'Colocalized trajectories minimum overlap :',num2str(minpts),' frames' ,NL,...
    'Maximum distance allowed between colocalized trajectories centers: ',num2str(max_center_dist),' pix',NL,...
    'Maximum distance allowed between colocalized spots: ',num2str(max_spot_dist), ' pix',NL,NL];

if isempty(coloc_sum)
    ntraj1_coloc = 0;
    ntraj2_coloc = 0;
    total_coloc_traj = 0;
else
    ntraj1_coloc = sum( unique(coloc_sum(:,1)) ~= -1);
    ntraj2_coloc = sum( unique(coloc_sum(:,2)) ~= -1);
    total_coloc_traj = size(coloc_sum,1);
end

if isempty(undecided_sum)
    ntraj1_undecided = 0;
    ntraj2_undecided = 0;
else
    ntraj1_undecided = sum( unique(undecided_sum(:,1)) ~= -1);
    ntraj2_undecided = sum( unique(undecided_sum(:,2)) ~= -1);
end

if isempty(neg_sum)
    ntraj1_neg = 0;
    ntraj2_neg = 0;    
else
    ntraj1_neg = sum( unique(neg_sum(:,1)) ~= -1);
    ntraj2_neg = sum( unique(neg_sum(:,2)) ~= -1);
end

if isempty(trk1sum)
    ntraj1 = 0;
else
    ntraj1 = size(trk1sum,1);
end

if isempty(trk2sum)
    ntraj2 = 0;
else
    ntraj2 = size(trk2sum,1);
end

ftraj1_coloc = 100*ntraj1_coloc/ntraj1;
ftraj1_undecided = 100*ntraj1_undecided/ntraj1;
ftraj1_neg = 100*ntraj1_neg/ntraj1;

ftraj2_coloc = 100*ntraj2_coloc/ntraj2;
ftraj2_undecided = 100*ntraj2_undecided/ntraj2;
ftraj2_neg = 100*ntraj2_neg/ntraj2;
        
match_str = [num2str(total_coloc_traj),' colocalized trajectories',NL,NL,...
    'Channel 1: ',...
    num2str(ntraj1),' total trajectories',NL,...
    num2str(ntraj1_coloc),' colocalized trajectories (',num2str(ftraj1_coloc),'%)',NL,...
    num2str(ntraj1_undecided),' undecided trajectories (',num2str(ftraj1_undecided),'%)',NL,...
    num2str(ntraj1_neg),' trajectories without a match (',num2str(ftraj1_neg),'%)',NL,NL,...
    'Channel 2: ',...
    num2str(ntraj2),' total trajectories',NL,...
    num2str(ntraj2_coloc),' colocalized trajectories (',num2str(ftraj2_coloc),'%)',NL,...
    num2str(ntraj2_undecided),' undecided trajectories (',num2str(ftraj2_undecided),'%)',NL,...
    num2str(ntraj2_neg),' trajectories without a match (',num2str(ftraj2_neg),'%)',NL,NL,NL];    

if isempty(coloc_res.D1c)
    ColocAvgD1 = NaN;
    ColocStdD1 = NaN;
    ColocMinD1 = NaN;
    ColocMaxD1 = NaN;
else
    ColocAvgD1 = nanmean(coloc_res.D1c(1,:));
    ColocStdD1 = nanstd(coloc_res.D1c(1,:));
    ColocMinD1 = min(coloc_res.D1c(1,:));
    ColocMaxD1 = max(coloc_res.D1c(1,:));
end

if isempty(coloc_res.D2c)
    ColocAvgD2 = NaN;
    ColocStdD2 = NaN;
    ColocMinD2 = NaN;
    ColocMaxD2 = NaN;
else
    ColocAvgD2 = nanmean(coloc_res.D2c(1,:));
    ColocStdD2 = nanstd(coloc_res.D2c(1,:));
    ColocMinD2 = min(coloc_res.D2c(1,:));
    ColocMaxD2 = max(coloc_res.D2c(1,:));
end

if isempty(neg_res.D1s)
    NegAvgD1 = NaN;
    NegStdD1 = NaN;
    NegMinD1 = NaN;
    NegMaxD1 = NaN;
else
    NegAvgD1 = nanmean(neg_res.D1s(1,:));
    NegStdD1 = nanstd(neg_res.D1s(1,:));
    NegMinD1 = min(neg_res.D1s(1,:));
    NegMaxD1 = max(neg_res.D1s(1,:));
end

if isempty(neg_res.D2s)
    NegAvgD2 = NaN;
    NegStdD2 = NaN;
    NegMinD2 = NaN;
    NegMaxD2 = NaN;
else
    NegAvgD2 = nanmean(neg_res.D2s(1,:));
    NegStdD2 = nanstd(neg_res.D2s(1,:));
    NegMinD2 = min(neg_res.D2s(1,:));
    NegMaxD2 = max(neg_res.D2s(1,:));
end

D_str = ['Diffusion Coefficients (calculated using 1 frame displacement histograms):',NL,NL,...
    'Colocalized Trajectories',NL,...
    'Channel 1 Average D: ',num2str(ColocAvgD1),' um2s-1',NL,...
    'Channel 1 std D: ',num2str(ColocStdD1),' um2s-1',NL,...
    'Channel 1 min D: ',num2str(ColocMinD1),' um2s-1',NL,...
    'Channel 1 max D: ',num2str(ColocMaxD1),' um2s-1',NL,NL,...
    'Channel 2 Average D: ',num2str(ColocAvgD2),' um2s-1',NL,...
    'Channel 2 std D: ',num2str(ColocStdD2),' um2s-1',NL,...
    'Channel 2 min D: ',num2str(ColocMinD2),' um2s-1',NL,...
    'Channel 2 max D: ',num2str(ColocMaxD2),' um2s-1',NL,NL,NL,...
    'Isolated Trajectories',NL,...
    'Channel 1 Average D: ',num2str(NegAvgD1),' um2s-1',NL,...
    'Channel 1 std D: ',num2str(NegStdD1),' um2s-1',NL,...
    'Channel 1 min D: ',num2str(NegMinD1),' um2s-1',NL,...
    'Channel 1 max D: ',num2str(NegMaxD1),' um2s-1',NL,NL,...
    'Channel 2 Average D: ',num2str(NegAvgD2),' um2s-1',NL,...
    'Channel 2 std D: ',num2str(NegStdD2),' um2s-1',NL,...
    'Channel 2 min D: ',num2str(NegMinD2),' um2s-1',NL,...
    'Channel 2 max D: ',num2str(NegMaxD2),' um2s-1',NL,NL];

fprintf(fid1,'%s%s%s%s \r\n',file_str,params_str,match_str,D_str);
fclose(fid1);


close(fh);

end



function params = set_parameters()
params = 0;
%setting the figure up
fh = figure('Units','characters',...
                  'NumberTitle','off',...
                  'Name','Set Parameters',...
                  'Position',[10 10 80 20],...
                  'Visible','off'); 
              

%Pixel Size
uicontrol('Parent',fh,...
    'Units','characters',...
    'Style','text','String','Pixel Size (um): ',...
    'Position',[2 18.1 20 1.4]);

uicontrol('Parent',fh,...
    'Units','characters',...
    'Tag','hpix_size',...
    'Style','edit','String',num2str(0.107),...
    'Position',[22 18.1 10 1.4]);

%Time Interval
uicontrol('Parent',fh,...
    'Units','characters',...
    'Style','text','String','Time between frames (s): ',...
    'Position',[2 16.5 30 1.4]);

uicontrol('Parent',fh,...
    'Units','characters',...
    'Tag','htime_interval',...
    'Style','edit','String',num2str(0.010),...
    'Position',[32 16.5 10 1.4]);

%Exclusion Radius
uicontrol('Parent',fh,...
    'Units','characters',...
    'Style','text','String','Size of Exclusion Region around image edges (pix): ',...
    'Position',[2 14 40 2.2]);

uicontrol('Parent',fh,...
    'Units','characters',...
    'Tag','hexclusion_radius',...
    'Style','edit','String',num2str(5),...
    'Position',[42 14 10 2.2]);

%Minimum Overlap
uicontrol('Parent',fh,...
    'Units','characters',...
    'Style','text','String','Matched Trajectories Minimum Number of Codetections : ',...
    'Position',[2 11.5 40 2.2]);

uicontrol('Parent',fh,...
    'Units','characters',...
    'Tag','hmin_overlap',...
    'Style','edit','String',num2str(3),...
    'Position',[42 11.5 10 2.2]);

%Distance Between Trajectory centers
uicontrol('Parent',fh,...
    'Units','characters',...
    'Style','text','String','Maximum Distance allowed between trajectory centers (pix): ',...
    'Position',[2 9 40 2.2]);

uicontrol('Parent',fh,...
    'Units','characters',...
    'Tag','hmax_center_dist',...
    'Style','edit','String',num2str(5),...
    'Position',[42 9 10 2.2]);

%Distance between spots
uicontrol('Parent',fh,...
    'Units','characters',...
    'Style','text','String','Maximum Distance allowed between codetected spots (pix): ',...
    'Position',[2 6.5 40 2.2]);

uicontrol('Parent',fh,...
    'Units','characters',...
    'Tag','hmax_spot_dist',...
    'Style','edit','String',num2str(2),...
    'Position',[42 6.5 10 2.2]);

%Action buttons
uicontrol('Parent',fh,...
    'Units','characters',...
    'Tag','hnext',...
    'Style','pushbutton',...
    'String','Next',...
    'Position',[60 0.5 10 2]);

%Action buttons
uicontrol('Parent',fh,...
    'Units','characters',...
    'Tag','hcancel',...
    'Style','pushbutton',...
    'String','Cancel',...
    'Position',[10 0.5 10 2]);

h = guihandles(fh);

set(h.hnext,'Callback',{@nextcancel_action,fh});
set(h.hcancel,'Callback',{@nextcancel_action,fh});

    function nextcancel_action(src,eventdata,fh)
        if strcmp(get(src,'Tag'),'hcancel')
            setappdata(fh,'params',0);
            uiresume
        elseif strcmp(get(src,'Tag'),'hnext')
            h = guihandles(fh);
            params = struct();
            params.pix_size = str2double(get(h.hpix_size,'String'));
            params.time_interval = str2double(get(h.htime_interval,'String'));
            params.exclusion_radius = str2double(get(h.hexclusion_radius,'String'));
            params.min_overlap = str2double(get(h.hmin_overlap,'String'));
            params.max_center_dist = str2double(get(h.hmax_center_dist,'String'));
            params.max_spot_dist = str2double(get(h.hmax_spot_dist,'String'));
            setappdata(fh,'params',params);
            uiresume
        end
    end

set(fh,'Visible','on');
uiwait

params = getappdata(fh,'params');
close(fh);
end
