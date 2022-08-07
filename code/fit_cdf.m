function [D,gof] = fit_cdf(xdata,ydata,weights,t)

k = find(ydata<=0.5);
D0 = xdata(k(end))^2/(4*t);
xdata = xdata/sqrt(4*t);

s = fitoptions('Method','NonlinearLeastSquares',...
               'Display','off',...
               'Lower',1e-8,...
               'Upper',15,...
               'Startpoint',D0,...
               'Weights',weights,...
               'TolX',0.0001,...
               'MaxIter',100);
         
f = fittype('1 - exp ( -x.^2/a)','options',s);
[c,gof] = fit(xdata,ydata,f);  
D = coeffvalues(c);


end