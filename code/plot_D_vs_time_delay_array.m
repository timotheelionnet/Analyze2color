function plot_D_vs_time_delay_array(D1,D2,time_interval,wintitle)

%compute median and std of D for each time interval
D1med = [];
if ~isempty(D1)
    tvalues = unique(D1(:,2));
    D1med = zeros(numel(tvalues),3);
    for i=1:numel(tvalues)
        D1med(i,1) = tvalues(i)*time_interval;
        curD1 = D1(D1(:,2) == tvalues(i),3);
        curD1 = curD1(logical((~isnan(curD1)).*(curD1~=0)));
        D1med(i,2) = median(curD1);
        D1med(i,3) = std(curD1);
    end
end

D2med = [];
if ~isempty(D2)
    tvalues = unique(D2(:,2));
    D2med = zeros(numel(tvalues),3);
    for i=1:numel(tvalues)
        D2med(i,1) = tvalues(i)*time_interval;
        curD2 = D2(D2(:,2) == tvalues(i),3);
        curD2 = curD2(logical((~isnan(curD2)).*(curD2~=0)));
        D2med(i,2) = median(curD2);
        D2med(i,3) = std(curD2);
    end
end

%plot the apparent D value vs time delay
figure; hold;
set(gcf,'Name',wintitle);
ymax = 0;
ymin = 0;
if ~isempty(D1med)
    errorbar(D1med(:,1),D1med(:,2),D1med(:,3),'-ok');
    ymax = max(  [max(D1med(:,2)), ymax]  );
    ymin = min(  [min(D1med(:,2)), ymin]  );
end
if ~isempty(D2med)
    errorbar(D2med(:,1),D2med(:,2),D2med(:,3),'-ob');
    ymax = max(  [max(D2med(:,2)), ymax]  );
    ymin = min(  [min(D2med(:,2)), ymin]  );
end 
ymax = ymax + 0.25*(ymax-ymin);
ymin = ymin - 0.25*(ymax-ymin);

if ymax <= ymin
    ymin = min(0,ymin);
    ymax = max(1,ymax);
end

%ylim([-0.1,1]);
set(gca,'XGrid','on');
set(gca,'YGrid','on');
xlabel('Delay (s)');
ylabel('Median Apparent D (um2s-1)');
legend('Channel 1','Channel 2');

end