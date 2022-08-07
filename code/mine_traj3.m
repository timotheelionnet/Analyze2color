function mine_traj3(datares,trajsum,traj,time_interval,traj_plot_mode,wintitle)

if isempty(trajsum) || isempty(traj) || isempty(datares)
   dispwin('Warning','Empty data set, Trajectory Browser not loaded') ;
   return
end

%create figure and axes
fh = create_figure_for_traj_mining(datares,trajsum,traj,time_interval,wintitle,traj_plot_mode);

set(fh,'Visible','on');
h = guihandles(fh);

set(h.hprev,'Callback', {@change_traj,fh,datares,trajsum,traj,time_interval,traj_plot_mode});
set(h.hnext,'Callback', {@change_traj,fh,datares,trajsum,traj,time_interval,traj_plot_mode});
set(h.hTrajIdx,'Callback', {@change_traj,fh,datares,trajsum,traj,time_interval,traj_plot_mode});

end

function change_traj(src,eventdata,fh,datares,trajsum,traj,time_interval,traj_plot_mode)
    srctag= get(src,'Tag');
    if strcmp(srctag,'hnext')
        info = getappdata(fh,'info');
        idx = info.idx;
        idx = idx +1;
    elseif strcmp(srctag,'hprev')
        info = getappdata(fh,'info');
        idx = info.idx;
        idx = idx -1;
    elseif strcmp(srctag,'hTrajIdx')
        h = guihandles(fh);
        TrajIdx = str2double(get(h.hTrajIdx,'String'));
        idx = find(trajsum(:,11) == TrajIdx);
        if isempty(idx)
            if TrajIdx > max(trajsum(:,11))
                idx = size(trajsum,1);
            elseif TrajIdx < min(trajsum(:,11))
                idx = 1;
            else
                while isempty(idx) && TrajIdx <= max(trajsum(:,11))
                    TrajIdx = TrajIdx +1;
                    idx = find(trajsum(:,11) == TrajIdx);
                end
            end    
        end
    else
       disp('wrong source');
    end
    
    update_info(fh,datares,trajsum,idx,traj_plot_mode);
    update_plots(fh,datares,trajsum,time_interval,traj,idx,traj_plot_mode);
end

function fh = create_figure_for_traj_mining(datares,trajsum,traj,time_interval,wintitle,traj_plot_mode)

% fh = figure(...
%               'Units','characters',...
%               'MenuBar','none',...
%               'Toolbar','none',...
%               'NumberTitle','off',...
%               'Name',wintitle,...
%               'Position',[20 3 200 55],...
%               'Visible','off'); 
fh = figure(...
              'Units','characters',...
              'NumberTitle','off',...
              'Name',wintitle,...
              'Position',[20 3 200 60],...
              'Visible','off'); 


%main panels          
msdpanel = uipanel(fh,'Title','MSD','Units','characters',...
             'Position',[1 37 100 18],'TitlePosition','centertop');
         
Ipanel = uipanel(fh,'Title','Intensity','Units','characters',...
             'Position',[1 23 100 14],'TitlePosition','centertop');
         
Dpanel = uipanel(fh,'Title','apparent D','Units','characters',...
             'Position',[111 5 100 50],'TitlePosition','centertop');
         
Trajpanel = uipanel(fh,'Title','Traj ','Units','characters',...
             'Position',[1 5 100 18],'TitlePosition','centertop');

%main axes         
axes('Parent',msdpanel,'Units','characters','Position',[12,3,80,12],'Tag','hmsd'); 
xlabel('Time (s)');
ylabel('M.S.D. (um2)');
set(gca,'XGrid','on');
set(gca,'YGrid','on');

axes('Parent',Dpanel,'Units','characters','Position',[12,30,72,18],'Tag','hCDF'); 
xlabel('Displacement (um)');
ylabel('CDF');
set(gca,'XGrid','on');
set(gca,'YGrid','on');

axes('Parent',Dpanel,'Units','characters','Position',[12,10,80,17],'Tag','hD'); 
xlabel('Delay (s)');
ylabel('Apparent D (um2s-1)');
set(gca,'XGrid','on');
set(gca,'YGrid','on');

axes('Parent',Ipanel,'Units','characters','Position',[12,3,80,8],'Tag','hI'); 
xlabel('time (s)');
ylabel('Intensity (AU)');
set(gca,'XGrid','on');
set(gca,'YGrid','on');

axes('Parent',Trajpanel,'Units','characters','Position',[45,4,50,12],'Tag','htraj');     
xlabel('X Position (pix)');
ylabel('Y Position (pix)');
set(gca,'XGrid','on');
set(gca,'YGrid','on');
axis ij

%text displays
%colocalized traj info
uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Style','text','String','Coloc Traj Idx: ',...
    'Position',[1,15,16,1.2]); 

uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Tag','hcoloctrajnum',...
    'Style','text','String','0',...
    'Position',[16,15,4,1.2]);

uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Style','text','String','Ctr X:',...
    'Position',[1,13.5,4,1.2]); 

uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Tag','hctrX',...
    'Style','text','String','0',...
    'Position',[7,13.5,8,1.2]);

uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Style','text','String','Ctr Y: ',...
    'Position',[1,12,6,1.2]); 

uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Tag','hctrY',...
    'Style','text','String','0',...
    'Position',[7,12,8,1.2]);

%traj 1 info
uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Style','text','String','Ch.1 Traj Idx:',...
    'Position',[1,10,14,1.2]); 
uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Tag','htraj1num',...
    'Style','text','String','0',...
    'Position',[15,10,4,1.2]); 

uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Style','text','String','Ch.1 Traj tstart: ',...
    'Position',[1,8.5,16,1.2]); 
uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Tag','htraj1ts',...
    'Style','text','String','0',...
    'Position',[17,8.5,4,1.2]); 

uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Style','text','String','Ch.1 Traj tend: ',...
    'Position',[1,7,15,1.2]); 
uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Tag','htraj1te',...
    'Style','text','String','0',...
    'Position',[16,7,4,1.2]);


%traj 2 info
uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Style','text','String','Ch.2: ',...
    'Position',[20,10,6,1.2]); 
uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Tag','htraj2num',...
    'Style','text','String','0',...
    'Position',[26,10,4,1.2]); 

uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Style','text','String','Ch.2: ',...
    'Position',[24,8.5,6,1.2]); 
uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Tag','htraj2ts',...
    'Style','text','String','0',...
    'Position',[30,8.5,4,1.2]); 

uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Style','text','String','Ch.2: ',...
    'Position',[22,7,6,1.2]); 
uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Tag','htraj2te',...
    'Style','text','String','0',...
    'Position',[28,7,4,1.2]);

%overlap info
uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Style','text','String','Time Overlap: ',...
    'Position',[1,5,16,1.2]); 
uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Tag','hoverlap',...
    'Style','text','String','0',...
    'Position',[17,5,4,1.2]);

uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Style','text','String','Codetections: ',...
    'Position',[1,3.5,16,1.2]); 
uicontrol('Parent',Trajpanel,...
    'Units','characters',...
    'Tag','hcodetect',...
    'Style','text','String','0',...
    'Position',[17,3.5,4,1.2]);


%D apparent info
uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Style','text','String','Channel 1 D (single; 10 ms): ',...
    'Position',[1,5,40,1.2]); 
uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Tag','hCDF1s_10',...
    'Style','text','String','0',...
    'Position',[42,5,8,1.2]);

uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Style','text','String','Channel 1 D (single; 50 ms): ',...
    'Position',[1,3.5,40,1.2]); 
uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Tag','hCDF1s_50',...
    'Style','text','String','0',...
    'Position',[42,3.5,8,1.2]);

uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Style','text','String','Channel 1 D (single; 100 ms): ',...
    'Position',[1,2,40,1.2]); 
uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Tag','hCDF1s_100',...
    'Style','text','String','0',...
    'Position',[42,2,8,1.2]);

uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Style','text','String','Channel 2: ',...
    'Position',[51,5,15,1.2]); 
uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Tag','hCDF2s_10',...
    'Style','text','String','0',...
    'Position',[66,5,8,1.2]);

uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Style','text','String','Channel 2: ',...
    'Position',[51,3.5,15,1.2]); 
uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Tag','hCDF2s_50',...
    'Style','text','String','0',...
    'Position',[66,3.5,8,1.2]);

uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Style','text','String','Channel 2: ',...
    'Position',[51,2,15,1.2]); 
uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Tag','hCDF2s_100',...
    'Style','text','String','0',...
    'Position',[66,2,8,1.2]);


% Next/Prev pushbuttons
uicontrol(fh,'Style','pushbutton',...
                    'Tag','hprev',...
                    'Units','characters',...
                    'String','prev',...
                    'Position',[60 0.4 22 3]);
                
uicontrol(fh,'Style','pushbutton',...
                    'Tag','hnext',...
                    'Units','characters',...
                    'String','next',...
                    'Position',[116 0.4 22 3]);

% Trajectory Index Prompt to access directly arbitrary indices
uicontrol(fh,'Style','edit',...
                    'Tag','hTrajIdx',...
                    'Units','characters',...
                    'String','1',...
                    'Position',[99 1 6 1.5]);

idx = 1;
info = struct();
setappdata(fh,'info',info);
update_info(fh,datares,trajsum,idx,traj_plot_mode);
update_plots(fh,datares,trajsum,time_interval,traj,idx,traj_plot_mode);

end

function update_plots(fh,datares,trajsum,time_interval,traj,idx,traj_plot_mode)

%make sure index is within boundaries
if idx <= 1
    idx = 1;
end
if idx > size(trajsum,1)
    idx = size(trajsum,1);
end

h = guihandles(fh);

%plot trajectories
axes(h.htraj);
cla
if ~ishold(h.htraj)
    hold(h.htraj)
end
trajtmp = traj(traj(:,10) == trajsum(idx,11),:);
plot(trajtmp(:,2),trajtmp(:,1),'Color','r');
scatter(trajtmp(:,2),trajtmp(:,1),6,trajtmp(:,7),'filled');
plot(trajtmp(:,5),trajtmp(:,4),'Color','g');
scatter(trajtmp(:,5),trajtmp(:,4),6,trajtmp(:,7),'filled');
colorbar;

%plot Intensities
axes(h.hI);
cla
if ~ishold(h.hI)
    hold(h.hI)
end
t_axis = traj(traj(:,10) == trajsum(idx,11),7);
Itmp1 = traj(traj(:,10) == trajsum(idx,11),3);
Itmp2 = traj(traj(:,10) == trajsum(idx,11),6);
plot(t_axis,Itmp1,'Color','r');
scatter(t_axis,Itmp1,6,t_axis,'filled');
plot(t_axis,Itmp2,'Color','g');
scatter(t_axis,Itmp2,6,t_axis,'filled');
colorbar;

%plot MSD
axes(h.hmsd);
cla
if ~ishold(h.hmsd)
    hold(h.hmsd)
end

if strcmp(traj_plot_mode,'single')
    msd1 = datares.msds{trajsum(idx,11),1}(:,1);
    msd2 = datares.msds{trajsum(idx,11),1}(:,2);
elseif strcmp(traj_plot_mode,'co')
    msd1 = datares.msdc{trajsum(idx,11),1}(:,1);
    msd2 = datares.msdc{trajsum(idx,11),1}(:,2);
end
ts = time_interval*(1:size(msd1,1));
plot(ts,msd1,':ro');
plot(ts,msd2,':go');
xlim([0 max(ts)]);
ylim([0 max([max(msd1) max(msd2)])]);

%plot cdfs
axes(h.hCDF);
cla
if ~ishold(h.hCDF)
    hold(h.hCDF)
end

xmax = 0;
%10 ms data
if strcmp(traj_plot_mode,'single')
    cdf1 = datares.cdf1s{1,trajsum(idx,11)};
    cdf2 = datares.cdf2s{1,trajsum(idx,11)};
elseif strcmp(traj_plot_mode,'co')
    cdf1 = datares.cdf1c{1,trajsum(idx,11)};
    cdf2 = datares.cdf2c{1,trajsum(idx,11)};
end
if ~isempty(cdf1)
    plot(cdf1(:,2),cdf1(:,3),':rs','LineWidth',2);
    xmax = max([xmax;cdf1(:,2)]);
end
if ~isempty(cdf2)
    plot(cdf2(:,2),cdf2(:,3),':gs','LineWidth',2);
    xmax = max([xmax;cdf2(:,2)]);
end

%50 ms data
if strcmp(traj_plot_mode,'single')
    cdf1 = datares.cdf1s{5,trajsum(idx,11)};
    cdf2 = datares.cdf2s{5,trajsum(idx,11)};
elseif strcmp(traj_plot_mode,'co')
    cdf1 = datares.cdf1c{5,trajsum(idx,11)};
    cdf2 = datares.cdf2c{5,trajsum(idx,11)};
end
if ~isempty(cdf1)
    plot(cdf1(:,2),cdf1(:,3),':ro','LineWidth',2);
    xmax = max([xmax;cdf1(:,2)]);
end
if ~isempty(cdf2)
    plot(cdf2(:,2),cdf2(:,3),':go','LineWidth',2);
    xmax = max([xmax;cdf2(:,2)]);
end

%100 ms data
if strcmp(traj_plot_mode,'single')
    cdf1 = datares.cdf1s{10,trajsum(idx,11)};
    cdf2 = datares.cdf2s{10,trajsum(idx,11)};
elseif strcmp(traj_plot_mode,'co')
    cdf1 = datares.cdf1c{10,trajsum(idx,11)};
    cdf2 = datares.cdf2c{10,trajsum(idx,11)};
end
if ~isempty(cdf1)
    plot(cdf1(:,2),cdf1(:,3),':rd','LineWidth',2);
    xmax = max([xmax;cdf1(:,2)]);
end
if ~isempty(cdf2)
    plot(cdf2(:,2),cdf2(:,3),':gd','LineWidth',2);
    xmax = max([xmax;cdf2(:,2)]);
end

%fit and plot D values
if strcmp(traj_plot_mode,'single')
    D1plot = datares.D1s;
    D2plot = datares.D2s;
elseif strcmp(traj_plot_mode,'co')
    D1plot = datares.D1c;
    D2plot = datares.D2c;
end

if isnan(xmax)
    xmax = 1;
end

if xmax ~=0
    xfit = 0:(xmax/100):xmax;
    if ~isempty(D1plot)
        if D1plot(1,trajsum(idx,11))~=0
            yfit = 1 - exp( - xfit.^2/(4*D1plot(1,trajsum(idx,11))*1*time_interval) );
            plot(xfit,yfit,'-rs','MarkerSize',3);
        end
    end
    if ~isempty(D2plot)
        if D2plot(1,trajsum(idx,11))~=0
            yfit = 1 - exp( - xfit.^2/(4*D2plot(1,trajsum(idx,11))*1*time_interval) );
            plot(xfit,yfit,'-gs','MarkerSize',3);
        end
    end
    if size(D1plot,1) >= 5
        if D1plot(5,trajsum(idx,11))~=0
            yfit = 1 - exp( - xfit.^2/(4*D1plot(5,trajsum(idx,11))*5*time_interval) );
            plot(xfit,yfit,'-ro','MarkerSize',3);
        end
    end
    if size(D2plot,1) >= 5
        if D2plot(5,trajsum(idx,11))~=0
            yfit = 1 - exp( - xfit.^2/(4*D2plot(5,trajsum(idx,11))*5*time_interval) );
            plot(xfit,yfit,'-go','MarkerSize',3);
        end
    end
    if size(D1plot,1) >= 10
        if D1plot(10,trajsum(idx,11))~=0
            yfit = 1 - exp( - xfit.^2/(4*D1plot(10,trajsum(idx,11))*10*time_interval) );
            plot(xfit,yfit,'-rd','MarkerSize',3);
        end
    end

    if size(D2plot,1) >= 10
        if D2plot(10,trajsum(idx,11))~=0
            yfit = 1 - exp( - xfit.^2/(4*D2plot(10,trajsum(idx,11))*10*time_interval) );
            plot(xfit,yfit,'-gd','MarkerSize',3);
        end
    end
else
    xmax = 1;
end
xlim([0 xmax]);


%plot D
axes(h.hD);
cla
if ~ishold(h.hD)
    hold(h.hD)
end

ymax = 0;
if sum(~isnan(D1plot(:,trajsum(idx,11))))~=0
    k1 = find(~isnan(D1plot(:,trajsum(idx,11)) ));
    D1tmp = D1plot(k1,trajsum(idx,11));
    x1tmp = time_interval*k1;
    plot(x1tmp,D1tmp,'-ro');
    ymax = max([max(D1tmp) ymax]);
end

if sum(~isnan(D2plot(:,trajsum(idx,11))))~=0
    k2 = find(~isnan(D2plot(:,trajsum(idx,11)) ));
    D2tmp = D2plot(k2,trajsum(idx,11));
    x2tmp = time_interval*(k2);
    plot(x2tmp,D2tmp,'-go');
    ymax = max([ymax max(D2tmp)]);
end
if ymax == 0, ymax = 1; end
ylim([0,ymax]);
end

function update_info(fh,datares,trajsum,idx,traj_plot_mode)

if strcmp(traj_plot_mode,'single')
    D1info = datares.D1s;
    D2info = datares.D2s;
elseif strcmp(traj_plot_mode,'co')
    D1info = datares.D1c;
    D2info = datares.D2c;
end
%make sure index is within boundaries
if idx <= 1
    idx = 1;
end
if idx > size(trajsum,1)
    idx = size(trajsum,1);
end

info = getappdata(fh,'info');
info.idx = idx;
info.coloctraj = trajsum(idx,11);
info.traj1 = trajsum(idx,1);
info.traj2 = trajsum(idx,2);
info.s1 = trajsum(idx,3);
info.e1 = trajsum(idx,4);
info.s2 = trajsum(idx,5);
info.e2 = trajsum(idx,6);
info.overlap = trajsum(idx,7);
info.simult = trajsum(idx,8);
info.Ctrx = trajsum(idx,9);
info.Ctry = trajsum(idx,10);


if size(D1info,1) >= 1
    info.D1s1 = D1info(1,trajsum(idx,11));
else
    info.D1s1 = NaN;
end

if size(D2info,1) >= 1
    info.D2s1 = D2info(1,trajsum(idx,11));
else
    info.D2s1 = NaN;
end

if size(D1info,1) >= 5
    info.D1s5 = D1info(5,trajsum(idx,11));
else
    info.D1s5 = NaN;
end

if size(D2info,1) >= 5
    info.D2s5 = D2info(5,trajsum(idx,11));
else
    info.D2s5 = NaN;
end

if size(D1info,1) >= 10
    info.D1s10 = D1info(10,trajsum(idx,11));
else
    info.D1s10 = NaN;
end

if size(D2info,1) >= 10
    info.D2s10 = D2info(10,trajsum(idx,11));
else
    info.D2s10 = NaN;
end

setappdata(fh,'info',info);

h = guihandles(fh);
set(h.hTrajIdx,'String',num2str(trajsum(idx,11)));
set(h.hcoloctrajnum,'String',num2str(info.coloctraj));
set(h.hctrX,'String',num2str(info.Ctrx));
set(h.hctrY,'String',num2str(info.Ctry));
set(h.htraj1num,'String',num2str(info.traj1));
set(h.htraj1ts,'String',num2str(info.s1));
set(h.htraj1te,'String',num2str(info.e1));
set(h.htraj2num,'String',num2str(info.traj2));
set(h.htraj2ts,'String',num2str(info.s2));
set(h.htraj2te,'String',num2str(info.e2));
set(h.hoverlap,'String',num2str(info.overlap));
set(h.hcodetect,'String',num2str(info.simult));

set(h.hCDF1s_10,'String',num2str(info.D1s1));
set(h.hCDF1s_50,'String',num2str(info.D1s5));
set(h.hCDF1s_100,'String',num2str(info.D1s10));
set(h.hCDF2s_10,'String',num2str(info.D2s1));
set(h.hCDF2s_50,'String',num2str(info.D2s5));
set(h.hCDF2s_100,'String',num2str(info.D2s10));


end
