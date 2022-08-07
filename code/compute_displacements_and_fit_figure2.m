function [res,dispco,dispall,dispsingle] = compute_displacements_and_fit_figure2(traj,pix_size,interval,wintitle,traj_plot_mode)
npts = 10;

%compute displacements
[dispco,dispall,dispsingle,Intensity] = compute_displacements(traj,pix_size);
% dispco is a cell array that contains the displacements between co-detection events 
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

%compute cdfs
[cdf1c,cdf2c,n1c,n2c] = compute_cdf(dispco,npts);
[nx,ny] = size(cdf1c);
if  ny < max(traj(:,10))
    cdf1c(:,ny+1:max(traj(:,10))) = repmat({[]},nx,max(traj(:,10))-ny);
end
[nx,ny] = size(cdf2c);
if  ny < max(traj(:,10))
    cdf2c(:,ny+1:max(traj(:,10))) = repmat({[]},nx,max(traj(:,10))-ny);
end
[nx,ny] = size(n1c);
if  ny < max(traj(:,10))
    n1c(:,ny+1:max(traj(:,10))) = zeros(nx,max(traj(:,10))-ny);
end
[nx,ny] = size(n2c);
if  ny < max(traj(:,10))
    n2c(:,ny+1:max(traj(:,10))) = zeros(nx,max(traj(:,10))-ny);
end

[cdf1a,cdf2a,n1a,n2a] = compute_cdf(dispall,npts);
[nx,ny] = size(cdf1a);
if  ny < max(traj(:,10))
    cdf1a(:,ny+1:max(traj(:,10))) = repmat({[]},nx,max(traj(:,10))-ny);
end
[nx,ny] = size(cdf2a);
if  ny < max(traj(:,10))
    cdf2a(:,ny+1:max(traj(:,10))) = repmat({[]},nx,max(traj(:,10))-ny);
end
[nx,ny] = size(n1a);
if  ny < max(traj(:,10))
    n1a(:,ny+1:max(traj(:,10))) = zeros(nx,max(traj(:,10))-ny);
end
[nx,ny] = size(n2a);
if  ny < max(traj(:,10))
    n2a(:,ny+1:max(traj(:,10))) = zeros(nx,max(traj(:,10))-ny);
end

[cdf1s,cdf2s,n1s,n2s] = compute_cdf(dispsingle,npts);
[nx,ny] = size(cdf1s);
if  ny < max(traj(:,10))
    cdf1s(:,ny+1:max(traj(:,10))) = repmat({[]},nx,max(traj(:,10))-ny);
end
[nx,ny] = size(cdf2s);
if  ny < max(traj(:,10))
    cdf2s(:,ny+1:max(traj(:,10))) = repmat({[]},nx,max(traj(:,10))-ny);
end
[nx,ny] = size(n1s);
if  ny < max(traj(:,10))
    n1s(:,ny+1:max(traj(:,10))) = zeros(nx,max(traj(:,10))-ny);
end
[nx,ny] = size(n2s);
if  ny < max(traj(:,10))
    n2s(:,ny+1:max(traj(:,10))) = zeros(nx,max(traj(:,10))-ny);
end

%compute msds
[msdco,msdall,msdsingle] = compute_msd_colocalized(traj,pix_size);

disp('start fitting...');

%perform the fit on colocalized trajectories (stringent definition: only
%considering colocalized detections)
D1c = [];
gof1c = cell('');
rsq1c = [];
for i=1:size(cdf1c,1)
    t = interval*i;
    for j=1:size(cdf1c,2)
        if ~isempty(cdf1c{i,j})
            if n1c(i,j)>3
                weights = 1./sqrt(n1c(i,j)*cdf1c{i,j}(:,3));
                weights(cdf1c{i,j}(:,3)==0) = 1;
                [D1c(i,j),gof1c{i,j}] = fit_cdf(cdf1c{i,j}(:,2),cdf1c{i,j}(:,3),weights,t);
                rsq1c(i,j) = gof1c{i,j}.rsquare;
            else
                spots = dispco{i,1}(dispco{i,1}(:,1)==j,:);
                D1c(i,j) = (nanmedian(spots(:,4)))^2/(4*t);
                gof1c{i,j} = [];
                rsq1c(i,j) = NaN;
            end
        else
            D1c(i,j) = NaN;
            gof1c{i,j} = [];
            rsq1c(i,j) = NaN;
        end
    end
end
D2c = [];
gof2c = cell('');
rsq2c = [];
for i=1:size(cdf2c,1)
    t = interval*i;
    for j=1:size(cdf2c,2)
        if ~isempty(cdf2c{i,j})
            if n2c(i,j)>npts
                weights = 1./sqrt(n2c(i,j)*cdf2c{i,j}(:,3));
                weights(cdf2c{i,j}(:,3)==0) = 1;
                [D2c(i,j),gof2c{i,j}] = fit_cdf(cdf2c{i,j}(:,2),cdf2c{i,j}(:,3),weights,t);
                rsq2c(i,j) = gof2c{i,j}.rsquare;
            else
                spots = dispco{i,1}(dispco{i,1}(:,1)==j,:);
                D2c(i,j) = (nanmedian(spots(:,7)))^2/(4*t);
                gof2c{i,j} = [];
                rsq2c(i,j) = NaN;
            end
        else
            D2c(i,j) = NaN;
            gof2c{i,j} = [];
            rsq2c(i,j) = NaN;
        end
    end
end

%perform the fit on colocalized trajectories (inclusive definition)
D1a = [];
gof1a = cell('');
rsq1a = [];
for i=1:size(cdf1a,1)
    t = interval*i;
    for j=1:size(cdf1a,2)
        if ~isempty(cdf1a{i,j})
            if n1a(i,j)>npts
                weights = 1./sqrt(n1a(i,j)*cdf1a{i,j}(:,3));
                weights(cdf1a{i,j}(:,3)==0) = 1;
                [D1a(i,j),gof1a{i,j}] = fit_cdf(cdf1a{i,j}(:,2),cdf1a{i,j}(:,3),weights,t);
                rsq1a(i,j) = gof1a{i,j}.rsquare;
            else
                spots = dispall{i,1}(dispall{i,1}(:,1)==j,:);
                D1a(i,j) = (nanmedian(spots(:,4)))^2/(4*t);
                gof1a{i,j} = [];
                rsq1a(i,j) = NaN;
            end
        else
            D1a(i,j) = NaN;
            gof1a{i,j} = [];
            rsq1a(i,j) = NaN;
        end
    end
end
D2a = [];
gof2a = cell('');
rsq2a = [];
for i=1:size(cdf2a,1)
    t = interval*i;
    for j=1:size(cdf2a,2)
        if ~isempty(cdf2a{i,j})
            if n2a(i,j)>npts
                weights = 1./sqrt(n2a(i,j)*cdf2a{i,j}(:,3));
                weights(cdf2a{i,j}(:,3)==0) = 1;
                [D2a(i,j),gof2a{i,j}] = fit_cdf(cdf2a{i,j}(:,2),cdf2a{i,j}(:,3),weights,t);
                rsq2a(i,j) = gof2a{i,j}.rsquare;
            else
                spots = dispall{i,1}(dispall{i,1}(:,1)==j,:);
                D2a(i,j) = (nanmedian(spots(:,7)))^2/(4*t);
                gof2a{i,j} = [];
                rsq2a(i,j) = NaN;
            end
        else
            D2a(i,j) = NaN;
            gof2a{i,j} = [];
            rsq2a(i,j) = NaN;
        end
    end
end

%perform the fit on single trajectories
D1s = [];
gof1s = cell('');
rsq1s = [];
for i=1:size(cdf1s,1)
    t = interval*i;
    for j=1:size(cdf1s,2)
        if ~isempty(cdf1s{i,j})
            if n1s(i,j)>npts
                weights = 1./sqrt(n1s(i,j)*cdf1s{i,j}(:,3));
                weights(cdf1s{i,j}(:,3)==0) = 1;
                [D1s(i,j),gof1s{i,j}] = fit_cdf(cdf1s{i,j}(:,2),cdf1s{i,j}(:,3),weights,t);
                rsq1s(i,j) = gof1s{i,j}.rsquare;
            else
                spots = dispsingle{i,1}(dispsingle{i,1}(:,1)==j,:);
                D1s(i,j) = (nanmedian(spots(:,4)))^2/(4*t);
                gof1s{i,j} = [];
                rsq1s(i,j) = NaN;
            end
        else
            D1s(i,j) = NaN;
            gof1s{i,j} = [];
            rsq1s(i,j) = NaN;
        end
    end
end
D2s = [];
gof2s = cell('');
rsq2s = [];
for i=1:size(cdf2s,1)
    t = interval*i;
    for j=1:size(cdf2s,2)
        if ~isempty(cdf2s{i,j})
            if n2s(i,j)>npts
                weights = 1./sqrt(n2s(i,j)*cdf2s{i,j}(:,3));
                weights(cdf2s{i,j}(:,3)==0) = 1;
                [D2s(i,j),gof2s{i,j}] = fit_cdf(cdf2s{i,j}(:,2),cdf2s{i,j}(:,3),weights,t);
                rsq2s(i,j) = gof2s{i,j}.rsquare;
            else
                spots = dispsingle{i,1}(dispsingle{i,1}(:,1)==j,:);
                D2s(i,j) = (nanmedian(spots(:,7)))^2/(4*t);
                gof2s{i,j} = [];
                rsq2s(i,j) = NaN;
            end
            
        else
            D2s(i,j) = NaN;
            gof2s{i,j} = [];
            rsq2s(i,j) = NaN;
        end
    end
end


%compute median values
D1c_med = median_nonzeros(D1c);
D2c_med = median_nonzeros(D2c);
D1a_med = median_nonzeros(D1a);
D2a_med = median_nonzeros(D2a);
D1s_med = median_nonzeros(D1s);
D2s_med = median_nonzeros(D2s);

D1c_std = std_nonzeros(D1c);
D2c_std = std_nonzeros(D2c);
D1a_std = std_nonzeros(D1a);
D2a_std = std_nonzeros(D2a);
D1s_std = std_nonzeros(D1s);
D2s_std = std_nonzeros(D2s);
 

%store results in a structure
res = struct();
res.I1 = Intensity(:,1)';
res.I2 = Intensity(:,2)';

res.msdc = msdco;
clear('msdco');
res.msda = msdall;
clear('msdall');
res.msds = msdsingle;
clear('msdsingle');

res.cdf1c = cdf1c;
clear('cdf1c');
res.cdf1a = cdf1a;
clear('cdf1a');
res.cdf1s = cdf1s;
clear('cdf1s');
res.cdf2c = cdf2c;
clear('cdf2c');
res.cdf2a = cdf2a;
clear('cdf2a');
res.cdf2s = cdf2s;
clear('cdf2s');

res.n1c = n1c;
clear('n1c');
res.n1a = n1a;
clear('n1a');
res.n1s = n1s;
clear('n1s');
res.n2c = n2c;
clear('n2c');
res.n2a = n2a;
clear('n2a');
res.n2s = n2s;
clear('n2s');

res.D1c = D1c;
clear('D1c');
res.D1a = D1a;
clear('D1a');
res.D1s = D1s;
clear('D1s');
res.D2c = D2c;
clear('D2c');
res.D2a = D2a;
clear('D2a');
res.D2s = D2s;
clear('D2s');

res.gof1c = gof1c;
clear('gof1c');
res.gof1a = gof1a;
clear('gof1a');
res.gof1s = gof1s;
clear('gof1s');
res.rsq1c = rsq1c;
clear('rsq1c');
res.rsq1a = rsq1a;
clear('rsq1a');
res.rsq1s = rsq1s;
clear('rsq1s');
res.gof2c = gof2c;
clear('gof2c');
res.gof2a = gof2a;
clear('gof2a');
res.gof2s = gof2s;
clear('gof2s');
res.rsq2c = rsq2c;
clear('rsq2c');
res.rsq2a = rsq2a;
clear('rsq2a');
res.rsq2s = rsq2s;
clear('rsq2s');

res.D1c_med = D1c_med;
clear('D1c_med');
res.D1a_med = D1a_med;
clear('D1a_med');
res.D1s_med = D1s_med;
clear('D1s_med');
res.D2c_med = D2c_med;
clear('D2c_med');
res.D2a_med = D2a_med;
clear('D2a_med');
res.D2s_med = D2s_med;
clear('D2s_med')

res.D1c_std = D1c_std;
clear('D1c_std');
res.D1a_std = D1a_std;
clear('D1a_std');
res.D1s_std = D1s_std;
clear('D1s_std');
res.D2c_std = D2c_std;
clear('D2c_std');
res.D2a_std = D2a_std;
clear('D2a_std');
res.D2s_std = D2s_std;
clear('D2s_std')

%plot the apparent D value vs time delay

figure; hold;
set(gcf,'Name',wintitle);
ymax = 0;
ymin = 0;
if strcmp(traj_plot_mode,'single')
    errorbar(interval*(1:numel(res.D1s_med)),res.D1s_med,res.D1s_std,'-or');
    errorbar(interval*(1:numel(res.D2s_med)),res.D2s_med,res.D2s_std,'-og');
    ymax = max(  [max(res.D1s_med), max(res.D2s_med)]  );
    ymin = min(  [min(res.D1s_med), max(res.D2s_med)]  );
elseif strcmp(traj_plot_mode,'co')
    errorbar(interval*(1:numel(res.D1c_med)),res.D1c_med,res.D1c_std,'-or');
    errorbar(interval*(1:numel(res.D2c_med)),res.D2c_med,res.D2c_std,'-og');
    ymax = max(  [max(res.D1c_med), max(res.D2c_med)]  );
    ymin = min(  [min(res.D1c_med), max(res.D2c_med)]  );    
end

ymax = ymax + 0.25*(ymax-ymin);
ymin = ymin - 0.25*(ymax-ymin);
if ymax <= ymin
    ymin = min(0,ymin);
    ymax = max(1,ymax);
end
%ylim([ymin,ymax]);
ylim([-0.1,1]);
set(gca,'XGrid','on');
set(gca,'YGrid','on');
xlabel('Delay (s)');
ylabel('Median Apparent D (um2s-1)');



end