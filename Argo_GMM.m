% Argo GMM
cd("C:\Users\ajs82292\Desktop\Research\Matlab\Source\Southern_Ocean")
addpath("C:\Users\ajs82292\Desktop\Research\Matlab\Source\seawater")
addpath(genpath('C:\Users\ajs82292\Desktop\Research\Matlab\Source\gsw_matlab_v3_06_16'));
load Argo_profiles.mat



% remove profiles that don't have all values between 300-900 m

keep = true(1, length(profiles));

for i = 1:length(profiles)
    if any(isnan(profiles(i).temp)) || any(isnan(profiles(i).salt))
        keep(i) = false;
    end
end

profiles = profiles(keep);  % Keep only the valid profiles

profiles = profiles(arrayfun(@(p) numel(p.depth) == 601, profiles));

% Convert into SA and potential temp
for i = 1:length(profiles)
    presData = gsw_p_from_z(-1*profiles(i).depth,profiles(i).lat) ; % gives pressure from depth and lattitude
    profiles(i).pressure = presData;
    absSalData(:,i) = gsw_SA_from_SP(profiles(i).salt,profiles(i).pressure,profiles(i).lon,profiles(i).lat) ;
    profiles(i).abs_sal = absSalData(:,i) ;
    potTempData(:,i) = gsw_pt_from_t(profiles(i).abs_sal,profiles(i).temp,profiles(i).pressure) ; % p_reff = 0 since not specified
    profiles(i).pot_temp = potTempData(:,i) ;
end

% standardize across all ARGO floats (z-score per m across all ARGO data)

mean_sal = mean(absSalData, 2);
std_sal = std(absSalData,0,2) ;
mean_temp = mean(potTempData,2) ;
std_temp = std(potTempData,0,2) ;

profiles = rmfield(profiles, {'temp', 'salt'});
clear presData i

