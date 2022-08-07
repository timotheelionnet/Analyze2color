function [cdf_co1,cdf_co2,cdf_all1,cdf_all2,cdf_neg1,cdf_neg2] = compute_all_CDFs(...
    disp_co1,disp_co2,disp_all1,disp_all2,disp_neg1,disp_neg2,npts) 

    cdf_co1 = compute_cdf_single_traj_array(disp_co1,npts);
    cdf_co2 = compute_cdf_single_traj_array(disp_co2,npts);
    cdf_all1 = compute_cdf_single_traj_array(disp_all1,npts);
    cdf_all2 = compute_cdf_single_traj_array(disp_all2,npts);
    cdf_neg1 = compute_cdf_single_traj_array(disp_neg1,npts);
    cdf_neg2 = compute_cdf_single_traj_array(disp_neg2,npts);
    
end