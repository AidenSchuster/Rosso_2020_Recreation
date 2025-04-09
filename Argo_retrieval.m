% retrieve ARGO Data

% base URL
baseURL = 'https://www.ncei.noaa.gov/data/oceans/argo/gadr/data/indian/';

% start and end year/month
startYear = 2010; startMonth = 12;
endYear = 2018; endMonth = 9;

saveDir = 'C:\Users\ajs82292\Desktop\Research\Matlab\Source\Southern_Ocean';

% Define web options (increase timeout)
options = weboptions('Timeout', 60);

% Loop through years and months
for year = startYear:endYear
    for month = 1:12
        % Skip months before Dec 2010 in the first year
        if year == 2010 && month < startMonth
            continue;
        end
        % Stop when reaching Sept 2018
        if year == 2018 && month > endMonth
            break;
        end
        
        % Format the URL for the specific year/month
        folderURL = sprintf('%s%04d/%02d/', baseURL, year, month);
        
        fprintf('Accessing: %s\n', folderURL);
        
        % Try to read the directory page (some months may be missing)
        try
            html = webread(folderURL);
        catch
            fprintf('Skipping missing directory: %s\n', folderURL);
            continue;
        end
        
        % Extract file names (NetCDF files)
        fileNames = regexp(html, '(?<=href=")[^"]+\.nc', 'match');
        
        % Create a subdirectory for this month/year
        saveSubDir = fullfile(saveDir, sprintf('%04d_%02d', year, month));
        if ~exist(saveSubDir, 'dir')
            mkdir(saveSubDir);
        end
        
        % Download each file
        for i = 1:length(fileNames)
            fileUrl = [folderURL, fileNames{i}];  % Full URL of the file
            savePath = fullfile(saveSubDir, fileNames{i});  % Save path
            
            % Check if file already exists to avoid re-downloading
            if exist(savePath, 'file')
                fprintf('Skipping (already exists): %s\n', fileNames{i});
                continue;
            end
            
            % Retry mechanism for downloading
        maxRetries = 10 ;
       for attempt = 1:maxRetries
    try
        websave(savePath, fileUrl, options);
        fprintf('Successfully downloaded: %s\n', fileNames{i});
        break; % Exit loop if successful
    catch ME
        fprintf('Error on attempt %d: %s\n', attempt, ME.message);
        pause(5); % Wait before retrying
    end
        end
            fprintf('Downloading: %s\n', fileNames{i});
            websave(savePath, fileUrl,options);
        end
    end
end

disp('All downloads completed!');

%% Load all information
% Setup
% Initialize the structure for storing profiles
addpath('C:\Users\ajs82292\Desktop\Research\Matlab\Source\seawater');
profiles = struct('lat', {}, 'lon', {}, 'temp', {}, 'salt', {}, 'depth', {});

% Loop through files
fileList = dir('C:\Users\ajs82292\Desktop\Research\Matlab\Source\Southern_Ocean/**/*.nc');  % Recursively find all files

for k = 1:length(fileList)
    try
        fpath = fullfile(fileList(k).folder, fileList(k).name);  % Full path to the .nc file
        clear pres lat lon salt temp temp_qc sal_qc
        % Read necessary data from the NetCDF file
        pres = ncread(fpath, 'pres');
        temp = ncread(fpath, 'temp');
        salt = ncread(fpath, 'psal');
        temp_qc = ncread(fpath, 'temp_qc');  % Quality control flags for temperature
        sal_qc = ncread(fpath, 'psal_qc');  % Quality control flags for salinity
        lat = ncread(fpath, 'latitude');  % Latitude of the profile
        lon = ncread(fpath, 'longitude');  % Longitude of the profile

        if numel(lat) > 1
        % Option 1: Use the first element
        lat = lat(1);
        lon = lon(1) ;
        end
        
        temp_qc = temp_qc(:);
        sal_qc = sal_qc(:);
        pres = pres(:) ;
        salt = salt(:);
        temp = temp(:) ;
        % Clean NaNs and invalid QC values (QC should be 1, 5, or 8)
        valid = ismember(double(temp_qc) - 48, [1, 5, 8]) & ismember(double(sal_qc) - 48, [1, 5, 8]);
        
        if sum(valid) < 2
            fprintf('Skipping %s: not enough valid data for interpolation\n', fileList(k).name);
            continue  % Or 'continue' if inside a loop
        end

        % Remove NaNs and invalid QC data
        pres = pres(valid);
        temp = temp(valid);
        salt = salt(valid);
        temp_qc = temp_qc(valid);
        sal_qc = sal_qc(valid);
        %convert pressure to depth

        depth = sw_dpth(pres,lat) ;

        [depth, idx_unique] = unique(depth);  % Get unique values and their indices
        temp = temp(idx_unique);  % Apply the same indexing to the temperature
        salt = salt(idx_unique) ;
        if max(depth) < 400, continue; end  % Skip if not deep enough
        
        try
        % Interpolate the profile data at 1-meter intervals
        interpZ = (1:max(depth))';  % Interpolated depths from 1m to the maximum depth
        tempInterp = interp1(depth, temp, interpZ, 'linear', NaN);
        saltInterp = interp1(depth, salt, interpZ, 'linear', NaN);
        
       catch ME
    if contains(ME.message, 'Interpolation requires at least two sample points') || ...
       contains(ME.message, 'Sample points must be unique')
        fprintf('⚠️ Pausing on file: %s\nError: %s\n',fileList(k).name , ME.message);
        keyboard
    else
        rethrow(ME)
    end
        end

        % Extract data between 300m and 900m after interpolation
        depthRange = interpZ >= 300 & interpZ <= 900;
        tempInterp = tempInterp(depthRange);
        saltInterp = saltInterp(depthRange);
        
        % Check if there is valid data within the 300-900m range
        if any(~isnan(tempInterp)) && any(~isnan(saltInterp))  % Only save if both are valid in this range
            % Add this profile to the structure
            profiles(end+1).lat = lat;
            profiles(end).lon = lon;
            profiles(end).temp = tempInterp;
            profiles(end).salt = saltInterp;
            profiles(end).depth = interpZ(depthRange);  % Store the depth in this range
        end
        
    catch ME
        fprintf('Skipping %s: %s\n', fileList(k).name, ME.message);
        continue;
    end
end

% The 'profiles' struct now contains all the valid profiles

%% 
% index for the proper boundaries
test_profiles = profiles ;
lat_vals = [test_profiles.lat];
lon_vals = [test_profiles.lon];

% Logical index for filtering
keep_idx = lat_vals <= -30 & lon_vals > 0 & lon_vals < 180;

test_profiles_filtered = test_profiles(keep_idx);

%% Save
profiles = test_profiles_filtered ;

save("Argo_profiles.mat", 'profiles', '-v7.3');

