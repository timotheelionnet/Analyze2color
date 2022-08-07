function [msd_co1,msd_co2,msd_all1,msd_all2,msd_neg1,msd_neg2] = compute_all_MSDs(...
    disp_co1,disp_co2,disp_all1,disp_all2,disp_neg1,disp_neg2)

    msd_co1 = compute_msd_from_displacements(disp_co1);
    msd_co2 = compute_msd_from_displacements(disp_co2);
    msd_all1 = compute_msd_from_displacements(disp_all1);
    msd_all2 = compute_msd_from_displacements(disp_all2);
    msd_neg1 = compute_msd_from_displacements(disp_neg1);
    msd_neg2 = compute_msd_from_displacements(disp_neg2);
    
end