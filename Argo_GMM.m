% Argo GMM
addpath("C:\Users\ajs82292\Desktop\Research\Matlab\Source\seawater")
addpath("C:\Users\ajs82292\Desktop\Research\Matlab\Source\gsw_matlab_v3_06_16")
load Argo_profiles.mat

% remove profiles that don't have values between 300-900 m
profiles = profiles(arrayfun(@(p) numel(p.depth) == 601, profiles));

% extract variables
lat = [profiles.lat];
lon = [profiles.lon];

pres = gsw_p_from_z()

clear profiles

%convert into  absolute sal and potential temp
for i = 1:size(lat,1)

end

% standardize