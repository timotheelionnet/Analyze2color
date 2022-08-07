function mine_traj3_array(trk1,trk2,MSD1,MSD2,cdf1,cdf2,D1,D2,time_interval,wintitle)

%create figure and axes
fh = create_figure_for_traj_mining(trk1,trk2,MSD1,MSD2,cdf1,cdf2,D1,D2,time_interval,wintitle);

set(fh,'Visible','on');
h = guihandles(fh);

set(h.hprev,'Callback', {@change_traj,fh,trk1,trk2,MSD1,MSD2,cdf1,cdf2,D1,D2,time_interval});
set(h.hnext,'Callback', {@change_traj,fh,trk1,trk2,MSD1,MSD2,cdf1,cdf2,D1,D2,time_interval});
set(h.hTrajIdx,'Callback', {@change_traj,fh,trk1,trk2,MSD1,MSD2,cdf1,cdf2,D1,D2,time_interval});

end

function change_traj(src,eventdata,fh,trk1,trk2,MSD1,MSD2,cdf1,cdf2,D1,D2,time_interval)
    srctag= get(src,'Tag');
    info = getappdata(fh,'info');
    if strcmp(srctag,'hnext')
        info = getappdata(fh,'info');
        idx = info.idx;
        k = find(info.idx_list == idx);
        if k < numel(info.idx_list)
            idx = info.idx_list(k+1);
        end
    elseif strcmp(srctag,'hprev')
        info = getappdata(fh,'info');
        idx = info.idx;
        k = find(info.idx_list == idx);
        if k > 1
            idx = info.idx_list(k-1);
        end     
    elseif strcmp(srctag,'hTrajIdx')
        h = guihandles(fh);
        TrajIdx = str2double(get(h.hTrajIdx,'String'));
        k = find(info.idx_list == TrajIdx);
        if isempty(k)
            if TrajIdx > max(info.idx_list)
                TrajIdx = max(info.idx_list);
            elseif TrajIdx < min(info.idx_list)
                TrajIdx = min(info.idx_list);
            else
                while isempty(k) && TrajIdx <= max(info.idx_list)
                    TrajIdx = TrajIdx +1;
                    k = find(info.idx_list == TrajIdx);
                end
            end  
        end
        idx = TrajIdx;
    else
       disp('wrong source');
    end
    
    update_info(fh,trk1,trk2,D1,D2,idx);
    update_plots(fh,trk1,trk2,MSD1,MSD2,cdf1,cdf2,D1,D2,time_interval,idx);
end

function fh = create_figure_for_traj_mining(trk1,trk2,MSD1,MSD2,cdf1,cdf2,D1,D2,time_interval,wintitle)

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
    'Tag','hCDF1s_1',...
    'Style','text','String','0',...
    'Position',[42,5,8,1.2]);

uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Style','text','String','Channel 1 D (single; 50 ms): ',...
    'Position',[1,3.5,40,1.2]); 
uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Tag','hCDF1s_5',...
    'Style','text','String','0',...
    'Position',[42,3.5,8,1.2]);

uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Style','text','String','Channel 1 D (single; 100 ms): ',...
    'Position',[1,2,40,1.2]); 
uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Tag','hCDF1s_10',...
    'Style','text','String','0',...
    'Position',[42,2,8,1.2]);

uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Style','text','String','Channel 2: ',...
    'Position',[51,5,15,1.2]); 
uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Tag','hCDF2s_1',...
    'Style','text','String','0',...
    'Position',[66,5,8,1.2]);

uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Style','text','String','Channel 2: ',...
    'Position',[51,3.5,15,1.2]); 
uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Tag','hCDF2s_5',...
    'Style','text','String','0',...
    'Position',[66,3.5,8,1.2]);

uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Style','text','String','Channel 2: ',...
    'Position',[51,2,15,1.2]); 
uicontrol('Parent',Dpanel,...
    'Units','characters',...
    'Tag','hCDF2s_10',...
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

info = struct();
if isempty(trk2) 
   info.traj_plot_mode = 'single traj 1';
   info.idx_list = unique(trk1(:,5));
else
    if isempty(trk1)
        info.traj_plot_mode = 'single traj 2';
        info.idx_list = unique(trk2(:,5));
    else
        info.traj_plot_mode = 'colocalized';
        info.idx_list = unique(trk1(:,7));
    end
end
idx = info.idx_list(1);
setappdata(fh,'info',info);
update_info(fh,trk1,trk2,D1,D2,idx);
update_plots(fh,trk1,trk2,MSD1,MSD2,cdf1,cdf2,D1,D2,time_interval,idx);

end

function update_plots(fh,trk1,trk2,MSD1,MSD2,cdf1,cdf2,D1,D2,time_interval,idx)

info = getappdata(fh,'info');

%make sure index is within boundaries
if idx <= min(info.idx_list)
    idx = min(info.idx_list);
end
if idx > max(info.idx_list)
    idx = max(info.idx_list);
end

h = guihandles(fh);

%plot trajectories
axes(h.htraj);
cla
if ~ishold(h.htraj)
    hold(h.htraj)
end
if ~isempty(trk1)
    traj1tmp = trk1(trk1(:,5) == info.traj1,:);
    plot(traj1tmp(:,2),traj1tmp(:,1),'Color','k','LineWidth',2);
    scatter(traj1tmp(:,2),traj1tmp(:,1),6,traj1tmp(:,4),'filled');
end
if ~isempty(trk2)
    traj2tmp = trk2(trk2(:,5) == info.traj2,:);
    plot(traj2tmp(:,2),traj2tmp(:,1),'Color',[0.5,0.8,1],'LineWidth',2);
    scatter(traj2tmp(:,2),traj2tmp(:,1),6,traj2tmp(:,4),'filled');
end
colorbar;

%plot Intensities
axes(h.hI);
cla
if ~ishold(h.hI)
    hold(h.hI)
end
if ~isempty(trk1)
    t_axis = trk1(trk1(:,5) == info.traj1,4);
    Itmp1 = trk1(trk1(:,5) == info.traj1,3);
    plot(t_axis,Itmp1,'Color','k');
    scatter(t_axis,Itmp1,6,t_axis,'filled');
else
    t_axis = trk2(trk2(:,5) == info.traj2,4);
    Itmp2 = trk2(trk2(:,5) == info.traj2,3);
    plot(t_axis,Itmp2,'Color',[0.5,0.8,1]);
    scatter(t_axis,Itmp2,6,t_axis,'filled');
end
colorbar;

%plot MSD
axes(h.hmsd);
cla
if ~ishold(h.hmsd)
    hold(h.hmsd)
end

ts = [];
msd_tmp = [];
if ~isempty(MSD1)
    msd1_tmp = MSD1(MSD1(:,1) == info.traj1,3);
    ts1 = time_interval*(1:size(msd1_tmp,1));
    plot(ts1,msd1_tmp,':ko');
    ts = [ts,ts1];
    msd_tmp = [msd_tmp;msd1_tmp];
end
if ~isempty(MSD2) && ~isempty(trk2)
    msd2_tmp = MSD2(MSD2(:,1) == info.traj2,3);
    ts2 = time_interval*(1:size(msd2_tmp,1));
    plot(ts2,msd2_tmp,':o','Color',[0.5,0.8,1]);
    ts = [ts,ts2];
    msd_tmp = [msd_tmp;msd2_tmp];
end

if isempty(ts)
    ts = 1;
end
if isempty(msd_tmp)
    msd_tmp = 1;
end
xlim([0 max(ts)]);
ylim([0 max(msd_tmp)]);

%plot cdfs
axes(h.hCDF);
cla
if ~ishold(h.hCDF)
    hold(h.hCDF)
end

xmax = 0;
%1 frame data
if ~isempty(cdf1)
    idxtmp = logical((cdf1(:,1) == info.traj1).*((cdf1(:,2) == 1)));
    cdf1_tmp = cdf1(idxtmp,3:4);
    if ~isempty(cdf1_tmp)
        plot(cdf1_tmp(:,1),cdf1_tmp(:,2),':ks','LineWidth',2);
        xmax = max([xmax;cdf1_tmp(:,1)]);
    end
end
if ~isempty(cdf2)
    idxtmp = logical((cdf2(:,1) == info.traj2).*((cdf2(:,2) == 1)));
    cdf2_tmp = cdf2(idxtmp,3:4);
    if ~isempty(cdf2_tmp)
        plot(cdf2_tmp(:,1),cdf2_tmp(:,2),':s','Color',[0.5,0.8,1],'LineWidth',2);
        xmax = max([xmax;cdf2_tmp(:,1)]);
    end
end

% %5 frames data
% if ~isempty(cdf1)
%     idxtmp = logical((cdf1(:,1) == info.traj1).*((cdf1(:,2) == 5)));
%     cdf1_tmp = cdf1(idxtmp,3:4);
%     if ~isempty(cdf1_tmp)
%         plot(cdf1_tmp(:,1),cdf1_tmp(:,2),':ko','LineWidth',2);
%         xmax = max([xmax;cdf1_tmp(:,1)]);
%     end
% end
% if ~isempty(cdf2)
%     idxtmp = logical((cdf2(:,1) == info.traj2).*((cdf2(:,2) == 5)));
%     cdf2_tmp = cdf2(idxtmp,3:4);
%     if ~isempty(cdf2_tmp)
%         plot(cdf2_tmp(:,1),cdf2_tmp(:,2),':bo','LineWidth',2);
%         xmax = max([xmax;cdf2_tmp(:,1)]);
%     end
% end

%10 frames data
if ~isempty(cdf1)
    idxtmp = logical((cdf1(:,1) == info.traj1).*((cdf1(:,2) == 10)));
    cdf1_tmp = cdf1(idxtmp,3:4);
    if ~isempty(cdf1_tmp)
        plot(cdf1_tmp(:,1),cdf1_tmp(:,2),':kd','LineWidth',2);
        xmax = max([xmax;cdf1_tmp(:,1)]);
    end
end
if ~isempty(cdf2)
    idxtmp = logical((cdf2(:,1) == info.traj2).*((cdf2(:,2) == 10)));
    cdf2_tmp = cdf2(idxtmp,3:4);
    if ~isempty(cdf2_tmp)
        plot(cdf2_tmp(:,1),cdf2_tmp(:,2),':d','Color',[0.5,0.8,1],'LineWidth',2);
        xmax = max([xmax;cdf2_tmp(:,1)]);
    end
end

%fit and plot D values
if ~isempty(D1)
    D1plot = D1(D1(:,1)== info.traj1,:);   
else
    D1plot = [];
end
if ~isempty(D2)
    D2plot = D2(D2(:,1)== info.traj2,:);  
else
    D2plot = [];
end

if isnan(xmax)
    xmax = 1;
end

if xmax ~=0
    xfit = 0:(xmax/100):xmax;
    if ~isempty(D1plot)
        if D1plot(D1plot(:,2)==1,3)~=0
            yfit = 1 - exp( - xfit.^2/(4*D1plot(D1plot(:,2)==1,3)*1*time_interval) );
            plot(xfit,yfit,'-ks','MarkerSize',3);
        end
        %     if ~isempty(D1plot(D1plot(:,2)==5,3))
        %         if D1plot(D1plot(:,2)==5,3)~=0
        %             yfit = 1 - exp( - xfit.^2/(4*D1plot(D1plot(:,2)==5,3)*5*time_interval) );
        %             plot(xfit,yfit,'-ko','MarkerSize',3);
        %         end
        %     end
        if ~isempty(D1plot(D1plot(:,2)==10,3))
            if D1plot(D1plot(:,2)==10,3)~=0
                yfit = 1 - exp( - xfit.^2/(4*D1plot(D1plot(:,2)==10,3)*10*time_interval) );
                plot(xfit,yfit,'-kd','MarkerSize',3);
            end
        end
        
    end
    if ~isempty(D2plot)
        if D2plot(D2plot(:,2)==1,3)~=0
            yfit = 1 - exp( - xfit.^2/(4*D2plot(D2plot(:,2)==1,3)*1*time_interval) );
            plot(xfit,yfit,'-s','Color',[0.5,0.8,1],'MarkerSize',3);
        end
        %     if ~isempty(D2plot(D2plot(:,2)==5,3))
        %         if D2plot(D2plot(:,2)==5,3)~=0
        %             yfit = 1 - exp( - xfit.^2/(4*D2plot(D2plot(:,2)==5,3)*5*time_interval) );
        %             plot(xfit,yfit,'-bo','MarkerSize',3);
        %         end
        %     end
        if ~isempty(D2plot(D2plot(:,2)==10,3))
            if D2plot(D2plot(:,2)==10,3)~=0
                yfit = 1 - exp( - xfit.^2/(4*D2plot(D2plot(:,2)==10,3)*10*time_interval) );
                plot(xfit,yfit,'-d','Color',[0.5,0.8,1],'MarkerSize',3);
            end
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
if ~isempty(D1)
    if sum(~isnan(D1plot(:,3)))~=0
        k1 = find(~isnan(D1plot(:,3)));
        D1tmp = D1plot(k1,3);
        x1tmp = D1plot(k1,2);
        plot(x1tmp,D1tmp,'-ko');
        ymax = max([max(D1tmp) ymax]);
    end
end
if ~isempty(D2)
    if sum(~isnan(D2plot(:,3)))~=0
        k2 = find(~isnan(D2plot(:,3)));
        D2tmp = D2plot(k2,3);
        x2tmp = D2plot(k2,2);
        plot(x2tmp,D2tmp,'-o','Color',[0.5,0.8,1]);
        ymax = max([ymax max(D2tmp)]);
    end
end
if ymax == 0, ymax = 1; end
ylim([0,ymax]);
end

function update_info(fh,trk1,trk2,D1,D2,idx)
info = getappdata(fh,'info');

%make sure index is within boundaries
if idx <= min(info.idx_list)
    idx = min(info.idx_list);
end
if idx > max(info.idx_list)
    idx = max(info.idx_list);
end


info.idx = idx;
switch info.traj_plot_mode 
    case 'colocalized'
        info.coloctraj = idx;
        k1 = find(trk1(:,7)==idx);
        if ~isempty(k1)
            info.traj1 = trk1(k1(1),5);
            info.s1 = trk1(k1(1),4);
            info.e1 = trk1(k1(end),4);
            info.overlap = numel(k1);
            info.simult = sum(trk1(k1,6));
            info.Ctrx = nanmean(trk1(k1,1));
            info.Ctry = nanmean(trk1(k1,2));
        else
            info.traj1 = 0;
            info.s1 = 0;
            info.e1 = 0;
            info.overlap = 0;
            info.simult = 0;
            info.Ctrx = [];
            info.Ctry = [];
        end
        
        k2 = find(trk2(:,7)==idx);
        if ~isempty(k2)
            info.traj2 = trk2(k2(1),5);
            info.s2 = trk2(k2(1),4);
            info.e2 = trk2(k2(end),4);
            info.simult = sum(trk1(k1,6));
            info.Ctrx = nanmean( [ mean(trk1(k1,1)), mean(trk2(k2,1))] );
            info.Ctry = nanmean( [ mean(trk1(k1,2)), mean(trk2(k2,2))] );
        else
            info.traj2 = 0;
            info.s2 = 0;
            info.e2 = 0;
            info.simult = 0;
            info.Ctrx = nanmean( [ info.Ctrx, mean(trk2(k2,1))] );
            info.Ctry = nanmean( [ info.Ctrx, mean(trk2(k2,2))] );
        end
        
        if isempty(info.Ctrx)
            info.Ctrx = 0;
            info.Ctry = 0;
        end
    case 'single traj 1'
        info.coloctraj = 0;
        info.traj1 = idx;
        info.traj2 = 0;
        k1 = find(trk1(:,5)==idx);
        if ~isempty(k1)
            info.s1 = trk1(k1(1),4);
            info.e1 = trk1(k1(end),4);
            info.Ctrx = nanmean(trk1(k1,1));
            info.Ctry = nanmean(trk1(k1,2));
        else
            info.s1 = 0;
            info.e1 = 0;
            info.Ctrx = 0;
            info.Ctry = 0;
        end
        
        info.s2 = 0;
        info.e2 = 0;
        info.overlap = 0;
        info.simult = 0;
        
    case 'single traj 2'
        info.coloctraj = 0;
        info.traj1 = 0;
        info.traj2 = idx;
        k2 = find(trk2(:,5)==idx);
        if ~isempty(k2)
            info.s2 = trk2(k2(1),4);
            info.e2 = trk2(k2(end),4);
            info.Ctrx = mean(trk2(k2,1));
            info.Ctry = mean(trk2(k2,2));
        else
            info.s2 = 0;
            info.e2 = 0;
            info.Ctrx = 0;
            info.Ctry = 0;
        end
        info.s1 = 0;
        info.e1 = 0;
        info.overlap = 0;
        info.simult = 0;
        
end

if ~isempty(D1)
    curD1 = D1(D1(:,1)==info.traj1,:);
    k1 = find(curD1(:,2)==1);
    if ~isempty(k1)
        info.D1s1 = curD1(k1,3);
    else
        info.D1s1 = NaN;
    end
    k5 = find(curD1(:,2)==5);
    if ~isempty(k5)
        info.D1s5 = curD1(k5,3);
    else
        info.D1s5 = NaN;
    end
    k10 = find(curD1(:,2)==10);
    if ~isempty(k10)
        info.D1s10 = curD1(k10,3);
    else
        info.D1s10 = NaN;
    end   
else
    info.D1s1 = NaN;
    info.D1s5 = NaN;
    info.D1s10 = NaN; 
end
if ~isempty(D2)
    curD2 = D2(D2(:,1)==info.traj2,:);
    k1 = find(curD2(:,2)==1);
    if ~isempty(k1)
        info.D2s1 = curD2(k1,3);
    else
        info.D2s1 = NaN;
    end
    k5 = find(curD2(:,2)==5);
    if ~isempty(k5)
        info.D2s5 = curD2(k5,3);
    else
        info.D2s5 = NaN;
    end
    k10 = find(curD2(:,2)==10);
    if ~isempty(k10)
        info.D2s10 = curD2(k10,3);
    else
        info.D2s10 = NaN;
    end  
else
    info.D2s1 = NaN;
    info.D2s5 = NaN;
    info.D2s10 = NaN;
end

setappdata(fh,'info',info);

h = guihandles(fh);

switch info.traj_plot_mode 
    case 'colocalized'
        set(h.hTrajIdx,'String',num2str(info.coloctraj));
    case 'single traj 1'
        set(h.hTrajIdx,'String',num2str(info.traj1));
    case 'single traj 2'
        set(h.hTrajIdx,'String',num2str(info.traj2));
end
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
set(h.hCDF1s_1,'String',num2str(info.D1s1));
set(h.hCDF1s_5,'String',num2str(info.D1s5));
set(h.hCDF1s_10,'String',num2str(info.D1s10));
set(h.hCDF2s_1,'String',num2str(info.D2s1));
set(h.hCDF2s_5,'String',num2str(info.D2s5));
set(h.hCDF2s_10,'String',num2str(info.D2s10));



end
