function Analyze2color_diatrack2()
tic
%use only first 2 columns of diatrack trajectories
use_2_columns_only = 1;
user_entry = 1;

%% Enter Parameters & Load data

%analysis parameters
params = set_parameters();
if ~isstruct(params)
    return
else
    pix_size = params.pix_size;
    time_interval = params.time_interval;
    exclusion_radius = params.exclusion_radius;
    minpts = params.min_overlap;
    max_spot_dist = params.max_spot_dist;
    if params.max_traj_dist >=0
        traj_screen_method = 'manual';
        traj_screen_tresh_value = params.max_traj_dist;
    else
        traj_screen_method = 'auto';
        traj_screen_tresh_value = 0;
    end
    perform_diffusion_analysis = params.perform_diffusion_analysis;
    x1shift = params.x_shift;
    y1shift = params.y_shift;
end
% pix_size = 0.107;
% time_interval = 0.01;
% exclusion_radius = 5;
%min_overlap = 3;

if user_entry  
    % load files Channel 1
    [dia1fname,dia1dirname,fidx] = uigetfile('*.txt','Load First Channel Diatrack Result','\');
    if fidx == 0
        return
    end
    % load files Channel 2
    [dia2fname,dia2dirname,fidx] = uigetfile('*.txt','Load Second Channel Diatrack Result',dia1dirname);
    if fidx == 0
        return
    end
    save_dirname = uigetdir(dia1dirname,'Saving Directory') ;
else
%     dia1fname = 'batch_spot_detection.loc';
%     dia1dirname = 'Z:\Timothee Lionnet\data\Jeff Trick Data\120618\MS2\registered';
%     dia1name_no_ext = 'batch_spot_detection';
%     trk1name = 'MS2_batch_spot_detection_2_1_5.trk';
%     trk1dirname = 'Z:\Timothee Lionnet\data\Jeff Trick Data\120618\MS2\registered';    
%     dia2fname = 'batch_spot_detection.loc';
%     dia2dirname = 'Z:\Timothee Lionnet\data\Jeff Trick Data\120618\PP7\registered';
%     dia2name_no_ext = 'batch_spot_detection';
%     trk2name = 'PP7_batch_spot_detection_2_1_5.trk';
%     trk2dirname = 'Z:\Timothee Lionnet\data\Jeff Trick Data\120618\PP7\registered';     
%     save_dirname = 'D:\junk\';
end

%load and parse diatrack output
trk1 = import_diatrack_data_XXL(fullfile(dia1dirname,dia1fname),use_2_columns_only);
[~,dia1name_no_ext,~] = fileparts(dia1fname);
trk2 = import_diatrack_data_XXL(fullfile(dia2dirname,dia2fname),use_2_columns_only);
[~,dia2name_no_ext,~] = fileparts(dia2fname);


%add shift in the coordinates
trk1(:,1) = trk1(:,1) + x1shift;
trk1(:,2) = trk1(:,2) + y1shift;


%% match trajectories
fh = dispwin('Status','Matching Trajectories, first screen...');
%first screen 
    % for each trajectory from channel one, finds all trajectories in channel 2
    % that satisfy the distance criterion (max_spot_dist) for at least minpts datapoints.
    matches = match_trajectories4_dia(trk1,trk2,minpts,max_spot_dist);

fh = dispwin('Status','Matching Trajectories, second screen...',fh);
%second screen 
    %measures the average distance between datapoints from each pair of matched trajectories
    %excludes trajectories based on either
    % 1) an automatic threshold on the avg distance (otsu's method)
    % 2) a manual threshold on the avg distance
    screened_matches = screen_trajectories(matches,traj_screen_method,traj_screen_tresh_value);

%TO DO: INSERT HERE ONE_TO_ONE SCREENING TO EXCLUDE DOUBLE MATCHES (IF
%NEEDED)

%create list of matches and negative trajectories
    [matrix,neg_list1,neg_list2,coloc_trk1,coloc_trk2,neg_trk1,neg_trk2] = ...
        sort_trajectories_using_match_results(trk1,trk2,screened_matches,max_spot_dist);

%plot msg and figures for visual control
    %plot spot distances histograms: screened versus first inclusive data
    plot_colocalization_stats(screened_matches,matrix,trk1,trk2);

    %plot all trajectories; channel 1 is in black, colocalized channel 2 is
    %red, orphan channel 2 is blue
    plot_colocalization_visual_output(trk1,coloc_trk2,neg_trk2,'Trajectories Map');

    %display results message
    x = sprintf(['Channel 1: ',num2str(numel(unique(coloc_trk1(:,5)))),' colocalized trajectories / ',...
    num2str(num2str(numel(unique(trk1(:,5))))),' total\n',...
    num2str(num2str(numel(unique(neg_trk1(:,5))))),' trajectories without a match\n\n',...
    'Channel 2: ',num2str(numel(unique(coloc_trk2(:,5)))),' colocalized trajectories / ',...
    num2str(num2str(numel(unique(trk2(:,5))))),' total\n',...
    num2str(num2str(numel(unique(neg_trk2(:,5))))),' trajectories without a match']);
    dispwin('results',x);

%save results
    save_colocalization_results(save_dirname,dia1dirname,dia1fname,dia2dirname,dia2fname,...
        params,matrix,neg_list1,neg_list2,trk1,trk2,coloc_trk1,coloc_trk2,neg_trk1,neg_trk2);

%TO DO: CLEAR VARIABLES THAT WON'T BE USED LATER

%% compute diffusion properties
%TO DO: REORDER LATER TO AVOID MEMORY SATURATION BY CLEARING VARIABLES
if perform_diffusion_analysis
    %compute displacements within each trajectory datapoints (units of frames, um)
    fh = dispwin('Status','Diffusion Characterization: Computing Displacements...',fh);
    [disp_co1,disp_co2,disp_all1,disp_all2,disp_neg1,disp_neg2] = compute_all_displacements(...
        coloc_trk1,coloc_trk2,neg_trk1,neg_trk2,pix_size);    

    %compute CDF (units of frames, um)
    fh = dispwin('Status','Diffusion Characterization: Computing CDFs...',fh);
    npts = 10;
    [cdf_co1,cdf_co2,cdf_all1,cdf_all2,cdf_neg1,cdf_neg2] = compute_all_CDFs(...
        disp_co1,disp_co2,disp_all1,disp_all2,disp_neg1,disp_neg2,npts);        

    %compute MSD (units of um2, frames)
    fh = dispwin('Status','Diffusion Characterization: Computing MSDs...',fh);
    [msd_co1,msd_co2,msd_all1,msd_all2,msd_neg1,msd_neg2] = compute_all_MSDs(...
        disp_co1,disp_co2,disp_all1,disp_all2,disp_neg1,disp_neg2);

    %Fit CDFs to obtain D coefficients (units of um2.s-1; 
    %each time interval used to compute D is stored in units of frames)
    fh = dispwin('Status','Diffusion Characterization: fitting CDFs...',fh);
    [D_co1,D_co2,D_all1,D_all2,D_neg1,D_neg2] = fit_all_Ds(...
        cdf_co1,cdf_co2,cdf_all1,cdf_all2,cdf_neg1,cdf_neg2,...
        disp_co1,disp_co2,disp_all1,disp_all2,disp_neg1,disp_neg2,...
        time_interval);

    %% plot diffusion analysis results

    %Plot D vs. time delay
        plot_D_vs_time_delay_array(D_co1,D_co2,time_interval,'Colocalized Trajectories');
        plot_D_vs_time_delay_array(D_neg1,D_neg2,time_interval,'Isolated Trajectories');

    %Plot Histogram of D values at t=5 frames
        nframes = 5;
        nbins = 100;
        plot_D_histograms(D_co1,D_co2,D_neg1,D_neg2,nbins,nframes,'Diffusion Coefficients');


    %% save kinetic analysis results   
    save_diffusion_analysis_results(save_dirname,dia1dirname,dia1fname,dia2dirname,dia2fname,...
        disp_co1,disp_co2,disp_all1,disp_all2,disp_neg1,disp_neg2,...
        cdf_co1,cdf_co2,cdf_all1,cdf_all2,cdf_neg1,cdf_neg2,...
        msd_co1,msd_co2,msd_all1,msd_all2,msd_neg1,msd_neg2,...
        D_co1,D_co2,D_all1,D_all2,D_neg1,D_neg2);
else
    msd_co1 = [];
    msd_co2 = [];
    cdf_co1 = [];
    cdf_co2 = [];
    D_co1 = [];
    D_co2 = [];
    msd_neg1 = [];
    msd_neg2 = [];
    cdf_neg1 = [];
    cdf_neg2 = [];
    D_neg1 = [];
    D_neg2 = [];
end
close(fh);


%% Launch Trajectory Browsers
mine_traj3_array(coloc_trk1,coloc_trk2,msd_co1,msd_co2,cdf_co1,cdf_co2,D_co1,D_co2,...
    time_interval,'Colocalized Trajectories Browser');
mine_traj3_array(neg_trk1,[],msd_neg1,[],cdf_neg1,[],D_neg1,[],...
    time_interval,'Isolated Trajectories Browser: Channel 1');
mine_traj3_array([],neg_trk2,[],msd_neg2,[],cdf_neg2,[],D_neg2,...
    time_interval,'Isolated Trajectories Browser: Channel 2');




t = toc


end


function save_diffusion_analysis_results(save_dirname,dia1dirname,dia1fname,dia2dirname,dia2fname,...
    disp_co1,disp_co2,disp_all1,disp_all2,disp_neg1,disp_neg2,...
    cdf_co1,cdf_co2,cdf_all1,cdf_all2,cdf_neg1,cdf_neg2,...
    msd_co1,msd_co2,msd_all1,msd_all2,msd_neg1,msd_neg2,...
    D_co1,D_co2,D_all1,D_all2,D_neg1,D_neg2)

    [~,dia1name_no_ext,~] = fileparts(dia1fname);
    [~,dia2name_no_ext,~] = fileparts(dia2fname);
    
    %save displacements data
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_disp_co1.txt']),'disp_co1','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_disp_co2.txt']),'disp_co2','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_disp_all1.txt']),'disp_all1','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_disp_all2.txt']),'disp_all2','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_disp_neg1.txt']),'disp_neg1','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_disp_neg2.txt']),'disp_neg2','-ascii');

    %save cdf data
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_cdf_co1.txt']),'cdf_co1','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_cdf_co2.txt']),'cdf_co2','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_cdf_all1.txt']),'cdf_all1','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_cdf_all2.txt']),'cdf_all2','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_cdf_neg1.txt']),'cdf_neg1','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_cdf_neg2.txt']),'cdf_neg2','-ascii');
    
    %save msd data
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_msd_co1.txt']),'msd_co1','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_msd_co2.txt']),'msd_co2','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_msd_all1.txt']),'msd_all1','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_msd_all2.txt']),'msd_all2','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_msd_neg1.txt']),'msd_neg1','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_msd_neg2.txt']),'msd_neg2','-ascii');

    %save D data
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_D_co1.txt']),'D_co1','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_D_co2.txt']),'D_co2','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_D_all1.txt']),'D_all1','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_D_all2.txt']),'D_all2','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_D_neg1.txt']),'D_neg1','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_D_neg2.txt']),'D_neg2','-ascii');

    %append diffusion data to match file
    summary_filename = fullfile(save_dirname,[dia1name_no_ext,'_coloc_results.txt']);
    fid1 = fopen(summary_filename,'a');
    NL = sprintf('\r\n');
    idx_t1 = D_co1(:,2) == 1;
    ColocAvgD1 = nanmean(D_co1(idx_t1,3));
    ColocStdD1 = nanstd(D_co1(idx_t1,3));
    ColocMinD1 = min(D_co1(idx_t1,3));
    ColocMaxD1 = max(D_co1(idx_t1,3));
    idx_t2 = D_co2(:,2) == 1;
    ColocAvgD2 = nanmean(D_co2(idx_t2,3));
    ColocStdD2 = nanstd(D_co2(idx_t2,3));
    ColocMinD2 = min(D_co2(idx_t2,3));
    ColocMaxD2 = max(D_co2(idx_t2,3));
    idx_t1 = D_neg1(:,2) == 1;
    NegAvgD1 = nanmean(D_neg1(idx_t1,3));
    NegStdD1 = nanstd(D_neg1(idx_t1,3));
    NegMinD1 = min(D_neg1(idx_t1,3));
    NegMaxD1 = max(D_neg1(idx_t1,3));
    idx_t2 = D_neg2(:,2) == 1;
    NegAvgD2 = nanmean(D_neg2(idx_t2,3));
    NegStdD2 = nanstd(D_neg2(idx_t2,3));
    NegMinD2 = min(D_neg2(idx_t2,3));
    NegMaxD2 = max(D_neg2(idx_t2,3));
    
    D_str = ['Diffusion Coefficients (calculated using 1 frame displacement data):',NL,NL,...
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

    fprintf(fid1,'%s \r\n',D_str);
    fclose(fid1);
end


function save_colocalization_results(save_dirname,dia1dirname,dia1fname,dia2dirname,dia2fname,...
    params,matrix,neg_list1,neg_list2,trk1,trk2,coloc_trk1,coloc_trk2,neg_trk1,neg_trk2)
    
    [~,dia1name_no_ext,~] = fileparts(dia1fname);
    [~,dia2name_no_ext,~] = fileparts(dia2fname);
    
    %save colocalization matrix
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_coloc_list.txt']),'matrix','-ascii');

    %save lists of non colocalized trajectories
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_neg1_list.txt']),'neg_list1','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_neg2_list.txt']),'neg_list2','-ascii');
    
    %save colocalized trajectories
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_coloc1.trk']),'coloc_trk1','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_coloc2.trk']),'coloc_trk2','-ascii');
    
    %save non colocalized trajectories
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_neg1.trk']),'neg_trk1','-ascii');
    save(fullfile(save_dirname,[dia1name_no_ext,'_',dia2name_no_ext,'_neg2.trk']),'neg_trk2','-ascii');
    
    %summary file
    summary_filename = fullfile(save_dirname,[dia1name_no_ext,'_coloc_results.txt']);
    fid1 = fopen(summary_filename,'w');

    NL = sprintf('\r\n');

    file_str = ['Diatrack File Channel 1: ',fullfile(dia1dirname,dia1fname),NL,...
        'Diatrack File Channel 2: ',fullfile(dia2dirname,dia2fname),NL,NL];


    params_str = ['Pixel size: ',num2str(params.pix_size),' um',NL,...
        'Time between frames: ',num2str(params.time_interval),' s',NL,...
        'Exclusion radius: ',num2str(params.exclusion_radius),' pix',NL,...
        'Colocalized trajectories minimum overlap :',num2str(params.min_overlap),' frames' ,NL,...
        'Maximum distance allowed between colocalized traj',num2str(params.max_traj_dist),' frames' ,NL,...
        'Maximum distance allowed between colocalized spots: ',num2str(params.max_spot_dist), ' pix',NL,NL];
    
    total_coloc_traj = size(matrix,1);
    ntraj1 = numel(unique(trk1(:,5)));
    ntraj1_coloc = numel(unique(coloc_trk1(:,5)));
    ntraj1_neg = numel(unique(neg_trk1(:,5)));
    ntraj2 = numel(unique(trk2(:,5)));
    ntraj2_coloc = numel(unique(coloc_trk2(:,5)));
    ntraj2_neg = numel(unique(neg_trk2(:,5)));

    ftraj1_coloc = 100*ntraj1_coloc/ntraj1;
    ftraj1_neg = 100*ntraj1_neg/ntraj1;

    ftraj2_coloc = 100*ntraj2_coloc/ntraj2;
    ftraj2_neg = 100*ntraj2_neg/ntraj2;

    match_str = [num2str(total_coloc_traj),' colocalized trajectories',NL,NL,...
        'Channel 1: ',...
        num2str(ntraj1),' total trajectories',NL,...
        num2str(ntraj1_coloc),' colocalized trajectories (',num2str(ftraj1_coloc),'%)',NL,...
        num2str(ntraj1_neg),' trajectories without a match (',num2str(ftraj1_neg),'%)',NL,NL,...
        'Channel 2: ',...
        num2str(ntraj2),' total trajectories',NL,...
        num2str(ntraj2_coloc),' colocalized trajectories (',num2str(ftraj2_coloc),'%)',NL,...
        num2str(ntraj2_neg),' trajectories without a match (',num2str(ftraj2_neg),'%)',NL,NL,NL];    

    fprintf(fid1,'%s%s%s \r\n',file_str,params_str,match_str);
    fclose(fid1);
end



function params = set_parameters()
params = 0;
%setting the figure up
fh = figure('Units','characters',...
                  'NumberTitle','off',...
                  'Name','Set Parameters',...
                  'Position',[10 10 80 20],...
                  'Visible','off'); 

%Perform Diffusion Analysis Switch
% uicontrol('Parent',fh,...
%     'Units','characters',...
%     'Style','text','String','Perform Diffusion Analysis: ',...
%     'Position',[2 55 20 1.4]);

uicontrol('Parent',fh,...
    'Units','characters',...
    'Tag','hperform_diffusion_analysis',...
    'Style','radiobutton','String','Diffusion Analysis?',...
    'Position',[55 15 24 4]);              
              
              
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
    'Style','edit','String',num2str(0),...
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

%Distance between spots
uicontrol('Parent',fh,...
    'Units','characters',...
    'Style','text','String','Maximum Distance allowed between codetected spots (pix): ',...
    'Position',[2 9 40 2.2]);

uicontrol('Parent',fh,...
    'Units','characters',...
    'Tag','hmax_spot_dist',...
    'Style','edit','String',num2str(3),...
    'Position',[42 9 10 2.2]);

%Distance between trajectories for stringent screen
uicontrol('Parent',fh,...
    'Units','characters',...
    'Style','text','String','Maximum Distance allowed between trajectories (pix; -1 for auto): ',...
    'Position',[2 6.5 40 2.2]);

uicontrol('Parent',fh,...
    'Units','characters',...
    'Tag','hmax_traj_dist',...
    'Style','edit','String',num2str(4),...
    'Position',[42 6.5 10 2.2]);

%Shift Channel 1 Coordinates
uicontrol('Parent',fh,...
    'Units','characters',...
    'Style','text','String','X shift (new X1 = old X1 + X shift): ',...
    'Position',[2 4 20 2.2]);

uicontrol('Parent',fh,...
    'Units','characters',...
    'Tag','hx_shift',...
    'Style','edit','String',num2str(-0.4),...
    'Position',[22 4 10 2.2]);

uicontrol('Parent',fh,...
    'Units','characters',...
    'Style','text','String','Y shift (new Y1 = old Y1 + Y shift): ',...
    'Position',[34 4 20 2.2]);

uicontrol('Parent',fh,...
    'Units','characters',...
    'Tag','hy_shift',...
    'Style','edit','String',num2str(-2.2),...
    'Position',[54 4 10 2.2]);


%Action buttons
uicontrol('Parent',fh,...
    'Units','characters',...
    'Tag','hnext',...
    'Style','pushbutton',...
    'String','Next',...
    'Position',[60 0.5 10 2]);

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
            params.max_spot_dist = str2double(get(h.hmax_spot_dist,'String'));
            params.max_traj_dist = str2double(get(h.hmax_traj_dist,'String'));
            params.x_shift = str2double(get(h.hx_shift,'String'));
            params.y_shift = str2double(get(h.hy_shift,'String'));
            params.perform_diffusion_analysis = get(h.hperform_diffusion_analysis,'Value');
            setappdata(fh,'params',params);
            uiresume
        end
    end

set(fh,'Visible','on');
uiwait

params = getappdata(fh,'params');
close(fh);
end
