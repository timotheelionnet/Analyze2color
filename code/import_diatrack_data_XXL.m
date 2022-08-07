function trk = import_diatrack_data_XXL(fname,use_2_columns_only)

%trk is the list of trajectories: 
%it includes spots positions, frame number and trajectory assigned, formated as follows:
%col 1: x
%col 2: y
%col 3: I (placeholder column of zeros for diatrack files; useful to keep
%format compatibility with Localize data)
%col 4: frame number
%col 5: trajectory number

%fname is the filename (including path) of the diatrack file
%use_2_columns_only is an option (set to 0 or 1) that when set to 1, uses
%only the first 2 columns of each three columns set corresponding to a
%unique trajectory in diatrack

%% load data as a cell array

%load file, extract lines from the files
fid = fopen(fname);
diadata = textscan(fid,'%s','delimiter','\t');
fclose(fid);
disp('file loaded');

diadata = vertcat(diadata{:});

%remove first line (header)
diadata = diadata(2:end,:);

%find the number of columns using the delimiters (multiple adjacent white
%spaces)
x = strfind(diadata{1}, ' ');
y = diff([0,x])~=1;

%parse each line to extract values in between spaces,
%store in fdiadata cell array
fdiadata = cellfun(@textscan,diadata,repmat({'%s'},size(diadata)),repmat({sum(y)},size(diadata)),...
    repmat({'delimiter'},size(diadata)),repmat({' '},size(diadata)),...
    repmat({'MultipleDelimsAsOne'},size(diadata)),repmat({1},size(diadata)));
clear('diadata');
fdiadata = cellfun(@transpose,fdiadata,'UniformOutput',0);
%fdiadata = cellfun(@squeeze,fdiadata,'UniformOutput',0);

%consolidate fdiadata (combine nested 1D cells into one 2D cell)
x= cell(size(fdiadata,1), size(fdiadata{1,1},2));
for i=1:size(x,1)
  x(i,:) = fdiadata{i,1};      
end
fdiadata = x;
clear('x');

%% stack trajectories vertically, and remove zeros that pad the end of trajectories

if use_2_columns_only
    %remove every third column (if there is no information there)
    idx = 1:size(fdiadata,2);
    idx = rem(idx,3)~=0;
    fdiadata = fdiadata(:,idx);
    
    %stack trajectories vertically
    ntraj = size(fdiadata,2)/2;
    tmptraj = mat2cell(fdiadata,size(fdiadata,1),2*ones(ntraj,1));
    clear('fdiadata');
    tmptraj = reshape(tmptraj,ntraj,1);
    
    x= cell(size(tmptraj,1)*size(tmptraj{1,1},1), size(tmptraj{1,1},2));
    for i=1:size(tmptraj,1)
      x((i-1)*size(tmptraj{1,1},1) +1 : i*size(tmptraj{1,1},1),:) = tmptraj{i,1};      
    end
    tmptraj = x;
    clear('x');
    
    %find and remove vacant entries that pad the end of trajectories
    idx = cellfun(@strcmp,tmptraj,repmat({'0.00'},size(tmptraj)));
    idx = ~logical(idx(:,1).*idx(:,2));
    tmptraj = tmptraj(idx,:);
    clear('idx');
    
    %convert strings to double
    tmptraj = cellfun(@str2double,tmptraj);
else
    %keep all columns 
    ntraj = size(fdiadata,2)/3;
    tmptraj = mat2cell(fdiadata,size(fdiadata,1),3*ones(ntraj,1));
    clear('fdiadata');
    tmptraj = reshape(tmptraj,ntraj,1);
    
    %stack trajectories vertically
    x= cell(size(tmptraj,1)*size(tmptraj{1,1},1), size(tmptraj{1,1},2));
    for i=1:size(tmptraj,1)
      x((i-1)*size(tmptraj{1,1},1) +1 : i*size(tmptraj{1,1},1),:) = tmptraj{i,1};      
    end
    tmptraj = x;
    clear('x');
    
    %find and remove vacant entries that pad the end of trajectories
    idx = cellfun(@strcmp,tmptraj,repmat({'0.00'},size(tmptraj)));
    idx = ~logical(idx(:,1).*idx(:,2));
    tmptraj = tmptraj(idx,:);
    clear('idx');
    
    %convert strings to double
    tmptraj = cellfun(@str2double,tmptraj);
end

fdiadata = tmptraj;
clear('tmptraj');

%% add Intensity, time column and trajectory index column

%isolate the rows that store the starting frame of each trajectory
idx = (fdiadata(:,1)~=0) .* (fdiadata(:,2)==0);

time_increments = (1:size(fdiadata,1))'; %a timer column

% reset the timer to the correct frame at the start of each trajectory
reset_times = fdiadata(:,1).*idx - time_increments.*idx;
for k=2:length(reset_times)
  if reset_times(k)==0
    reset_times(k)=reset_times(k-1);
  end
end
tcolumn = time_increments + reset_times;

%use first timeframe = time zero convention
tcolumn = tcolumn -2;
if use_2_columns_only
    fdiadata = [fdiadata,zeros(size(tcolumn)),tcolumn,cumsum(idx~=0)];
else
    fdiadata = [fdiadata,tcolumn,cumsum(idx~=0)];
end
%remove the rows that store the starting frame of each trajectory
trkraw = fdiadata(~idx,:);

clear('idx','time_increments','reset_times','tcolumn','fdiadata');
disp('file parsed');



%% clean up bogus trajectories
%identified by at least a pair of integer coordinates

%remove all integer coordinates pairs 
trk = trkraw(~logical((mod(trkraw(:,1),1)==0).*(mod(trkraw(:,2),1)==0)),:);

%REMOVED, COULD BE USEFUL TO KEEP INDICES CONSISTENCY WITH THE DIATRACK
%FILE FOR IGOR USERS
%(rewrite trajectory index in case entire trajectories were deleted)
% old_traj_idx = unique(trk(:,end));
% idx = 1:numel(old_traj_idx)-1;
% 
% new_traj_index = 1:max(old_traj_idx);
% new_traj_index(old_traj_idx(1:numel(old_traj_idx))) = idx(1:numel(old_traj_idx));
% 
% trk(:,end) = new_traj_index(trk(:,end));

end