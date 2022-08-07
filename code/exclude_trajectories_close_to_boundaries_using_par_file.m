function trk = exclude_trajectories_close_to_boundaries_using_par_file(trk,exclusion_radius,filename)

[dirname,fname_noext,extension] = fileparts(filename);

if strcmp(extension,'par3') || strcmp(extension,'par')
    [xmax,ymax] = set_exclusion_limits(filename,exclusion_radius);
else
    %guess name of .par file
    parfname = fullfile(dirname,[fname_noext,'.par3']);
    if exist(parfname,'file')
        [xmax,ymax] = set_exclusion_limits(parfname,exclusion_radius);
    else
        parfname = fullfile(dirname,[fname_noext,'.par']);
        if exist(parfname,'file')
            [xmax,ymax] = set_exclusion_limits(parfname,exclusion_radius);
        else
            disp('no .par file found: could not exclude trajectories close from border');
            return
        end
    end
end

%exclude out of boundaries spots
included_spots = (trk(:,1)>exclusion_radius).*(trk(:,1)<xmax);
included_spots = included_spots.*(trk(:,2)>exclusion_radius).*(trk(:,2)<ymax);

trk = trk(logical(included_spots),:);
end