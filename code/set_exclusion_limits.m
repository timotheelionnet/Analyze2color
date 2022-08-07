function [xmax,ymax] = set_exclusion_limits(parfname,exclusion_radius)

par = textread(parfname,'%c')';

kx = strfind(par,'nx');
ky = strfind(par,'ny');
kz = strfind(par,'nz');

xmax = str2double(par(kx+3:ky-1))-exclusion_radius;

if isempty(kz)
    kTRACK = strfind(par,'TRACK.PRO');
    if isempty(kTRACK)
        ymax = str2double(par(ky+3:end))-exclusion_radius;
    else
        ymax = str2double(par(ky+3:kTRACK-1))-exclusion_radius;
    end
else
    ymax = str2double(par(ky+3:kz-1))-exclusion_radius;
end


end
