function plot_D_histograms(D_co1,D_co2,D_neg1,D_neg2,nbins,nframes,wintitle)

figure;
set(gcf,'Name',wintitle);
hold;

%select the D coefs that have been calculated using nframes as the time
%interval
D_co1 = D_co1(D_co1(:,2)==nframes,3);
D_co2 = D_co2(D_co2(:,2)==nframes,3);
D_neg1 = D_neg1(D_neg1(:,2)==nframes,3);
D_neg2 = D_neg2(D_neg2(:,2)==nframes,3);

%generate corresponding histograms
Dmax = max([max(D_co1(:)),max(D_co2(:)),max(D_neg1(:)),max(D_neg2(:))]);
bins = 0:Dmax/nbins:Dmax;
[nco1,xco1] = hist(D_co1,bins);
[nco2,xco2] = hist(D_co2,bins);
[nneg1,xneg1] = hist(D_neg1,bins);
[nneg2,xneg2] = hist(D_neg2,bins);

%plot histograms
plot(xco1,nco1/sum(nco1),'Color',[0,0,0],'LineWidth',2);
plot(xneg1,nneg1/sum(nneg1),'Color',[0,0,0],'LineStyle',':','LineWidth',2);
plot(xco2,nco2/sum(nco2),'Color',[1,0,0],'LineWidth',2);
plot(xneg2,nneg2/sum(nneg2),'Color',[0,0,1],'LineWidth',2);

set(get(gca,'XLabel'),'String','Diffusion Coefficient (\mum^2.s^{-1})');
set(get(gca,'YLabel'),'String','Frequency');

legend('Channel 1 Colocalized','Channel 1 Isolated',...
    'Channel 2 Colocalized','Channel 2 Isolated');

end
    