function [disp_co1,disp_co2,disp_all1,disp_all2,disp_neg1,disp_neg2] = compute_all_displacements(...
    coloc_trk1,coloc_trk2,neg_trk1,neg_trk2,pix_size)

disp_neg1 = compute_displacements_array(neg_trk1,pix_size);
    disp_neg2 = compute_displacements_array(neg_trk2,pix_size);

    %these arrays compile only the colocalized datapoints from each trajectory
    idx_co = logical( (~isnan(coloc_trk1(:,1))).*(~isnan(coloc_trk2(:,1))));
    coloc_trk1_co = coloc_trk1( idx_co,:);
    coloc_trk2_co = coloc_trk2( idx_co,:);

    disp_co1 = compute_displacements_array(coloc_trk1_co,pix_size);
    disp_co2 = compute_displacements_array(coloc_trk2_co,pix_size);

    %these arrays compile all datapoints from trajectories that feature
    %colocalization events
    coloc_trk1_all = coloc_trk1(~isnan(coloc_trk1(:,1)),:);
    [~,idx1,~] = unique([coloc_trk1_all(:,5),coloc_trk1_all(:,4)],'rows'); %each time point from a given trajectory should be entered only once
    coloc_trk1_all = coloc_trk1_all(idx1,:);

    coloc_trk2_all = coloc_trk2(~isnan(coloc_trk2(:,1)),:);
    [~,idx2,~] = unique([coloc_trk2_all(:,5),coloc_trk2_all(:,4)],'rows'); %each time point from a given trajectory should be entered only once
    coloc_trk2_all = coloc_trk2_all(idx2,:);

    disp_all1 = compute_displacements_array(coloc_trk1_all,pix_size);
    disp_all2 = compute_displacements_array(coloc_trk2_all,pix_size);

end