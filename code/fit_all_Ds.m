function [D_co1,D_co2,D_all1,D_all2,D_neg1,D_neg2] = fit_all_Ds(...
    cdf_co1,cdf_co2,cdf_all1,cdf_all2,cdf_neg1,cdf_neg2,...
    disp_co1,disp_co2,disp_all1,disp_all2,disp_neg1,disp_neg2,...
    time_interval)

D_co1 = fit_cdf_array(cdf_co1,disp_co1,time_interval);
D_co2 = fit_cdf_array(cdf_co2,disp_co2,time_interval);
D_all1 = fit_cdf_array(cdf_all1,disp_all1,time_interval);
D_all2 = fit_cdf_array(cdf_all2,disp_all2,time_interval);
D_neg1 = fit_cdf_array(cdf_neg1,disp_neg1,time_interval);
D_neg2 = fit_cdf_array(cdf_neg2,disp_neg2,time_interval);

    
end