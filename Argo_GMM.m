% Argo GMM
cd("C:\Users\ajs82292\Desktop\Research\Matlab\Source\Southern_Ocean")
addpath("C:\Users\ajs82292\Desktop\Research\Matlab\Script\Southern_Ocean")
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

run = 2 ;
if run == 1
% Convert into SA and potential temp
for i = 1:length(profiles)
    presData = gsw_p_from_z(-1*profiles(i).depth,profiles(i).lat) ; % gives pressure from depth and lattitude
    profiles(i).pressure = presData;
    absSalData(:,i) = gsw_SA_from_SP(profiles(i).salt,profiles(i).pressure,profiles(i).lon,profiles(i).lat) ;
    profiles(i).abs_sal = absSalData(:,i) ;
    potTempData(:,i) = gsw_pt_from_t(profiles(i).abs_sal,profiles(i).temp,profiles(i).pressure) ; % p_reff = 0 since not specified
    profiles(i).pot_temp = potTempData(:,i) ;
end

profiles = rmfield(profiles, {'temp', 'salt','depth'});

clear keep presData
% standardize across all ARGO floats (z-score per m across all ARGO data) (m intervals not decibar like in the paper but shouldn't matter much)

mean_sal = mean(absSalData, 2);
std_sal = std(absSalData,0,2) ;
mean_temp = mean(potTempData,2) ;
std_temp = std(potTempData,0,2) ;

%z_sal = NaN(size(absSalData));
%z_temp = NaN(size(potTempData));

z_sal = (absSalData - mean_sal) ./ std_sal;  
z_sal = z_sal'; %invert
z_temp = (potTempData - mean_temp) ./ std_temp;
z_temp = z_temp';
clear i mean_sal mean_temp std_temp std_sal absSalData potTempData
save z_sal.mat z_sal
save z_temp.mat z_temp
end
load z_sal.mat
load z_temp.mat
clear i keep run
%% Run PCA

[coeff_sal, score_sal, ~ , ~] = pca(z_sal,'Centered','off');
[coeff_temp, score_temp, ~ , ~] = pca(z_temp,'Centered','off');

% Select Training and Testing indices

lat = [profiles.lat] ;
lon  = [profiles.lon] ;
grid_res = 0.1 ; % 0.1x0.1 boxes
lat_edges = -90:grid_res:-30;
lon_edges = 0:grid_res:180;
[~, lat_bin] = histc(lat, lat_edges);
[~, lon_bin] = histc(lon, lon_edges);
grid_id = lat_bin + (lon_bin - 1) * length(lat_edges);
unique_grid_ids = unique(grid_id);

% Initialize train_idx as a logical array of all zeros
train_idx = false(size(lat));  % Logical array of the same size as the number of profiles

% For each grid cell, randomly select one index for training
for i = 1:length(unique_grid_ids)
    idx_in_cell = find(grid_id == unique_grid_ids(i));  % Find profiles in the current grid cell
    if ~isempty(idx_in_cell)
        % Set the corresponding index in train_idx to true (selected for training)
        train_idx(idx_in_cell(randi(length(idx_in_cell)))) = true;
    end
end

% Count how many profiles were selected for training
fprintf('Using %.2f%% of total profiles for training (%d out of %d profiles).\n', ...
    100 * sum(train_idx) / length(lat), sum(train_idx), length(lat));
train_idx = train_idx';
%% Train GMM
k = 9 ; % number of clusters (same as used in Rosso paper)
first_eof_cons = 1 ;
num_eofs = 1+first_eof_cons ; % only used first two eofs
feature_matrix = [score_temp(train_idx,first_eof_cons:num_eofs),score_sal(train_idx,first_eof_cons:num_eofs)] ;

run = 2 ; % to save a given iteration of the model
if run == 1
options = statset('MaxIter', 2000, 'Display', 'final');  % Increase iterations to 500
indepen_model = fitgmdist(feature_matrix, k, 'Options', options,'Replicates', 1, 'RegularizationValue', 0.004);
save Rosso_model.mat indepen_model
end
load Rosso_model.mat
cluster_labels = cluster(indepen_model, feature_matrix);
cluster_probs = posterior(indepen_model, feature_matrix);

% apply to test data (already calculated PCA's all togehter)
test_feature_matrix = [score_temp(~train_idx,first_eof_cons:num_eofs),score_sal(~train_idx,first_eof_cons:num_eofs)] ;
test_labels = cluster(indepen_model, test_feature_matrix);
test_probs = posterior(indepen_model, test_feature_matrix);

%% Plotting
aspect_ratio = cosd(50) ; % Aspect Ratio at 65 N
unique_clusters = unique(cluster_labels);
cmap = [
    0.1, 0.2, 0.8;   % Blue 
    0.8, 0.1, 0.2;   % Red
    0.1, 0.8, 0.2;   % Green
    0.8, 0.8, 0.1;   % Yellow
    0.1, 0.8, 0.8;   % Cyan
    0.8, 0.1, 0.8;   % Magenta
    0.5, 0.5, 0.5;   % Gray
    0.9, 0.5, 0.2;   % Orange
    0.4, 0.3, 0.6;   % Navy
    0.6, 0.3, 0.1    % Brown
];
%cluster_color_map = containers.Map(unique_clusters, num2cell(cmap, 2));

load coastlines  % loads 'coastlat' and 'coastlon'
figure
hold on
plot(coastlon, coastlat, 'k')
gscatter(lon(train_idx), lat(train_idx), cluster_labels, cmap, 'o', 2.8, 'filled');
gscatter(lon(~train_idx),lat(~train_idx),test_labels,cmap,'o',2.8,'filled')
axis equal
xlabel('Longitude')
ylabel('Latitude')
title('World Coastline')
xlim([0,180])
ylim([-72,-28])
daspect([1 aspect_ratio 1])

%% Plot Sal and Temp avgs and std devs

%Plot Jones type plot (need actual salinity and temp for shape)
sal_cell = {profiles.salt};  
sal_mat = cell2mat(sal_cell)' ;

temp_cell = {profiles.temp};  
temp_mat = cell2mat(temp_cell)' ;

train_sal_mat = sal_mat(train_idx,:) ;
train_temp_mat = temp_mat(train_idx,:) ;

clear sal_cell temp_cell sal_mat temp_mat

n_clusters = length(unique_clusters);
n_rows = ceil(sqrt(n_clusters));
n_cols = ceil(n_clusters / n_rows);

%Sal
figure;
DepInterval_custom = 300:900 ;

for i = 1:length(unique_clusters)
    nexttile;

    cluster_id = unique_clusters(i);
    cluster_idx = cluster_labels == cluster_id;
    %sal_cluster = with_sal_test(cluster_idx,starting_depth:ending_depth) ;
    sal_cluster = train_sal_mat(cluster_idx,:) ;
end


for i = 1:n_clusters
    cluster_id = unique_clusters(i);
    cluster_color = cmap(i, :);  
    % Get high-probability profiles for this cluster
    cluster_mask = (cluster_labels == cluster_id) ;
    cluster_profiles = train_sal_mat(cluster_mask, :);
    cluster_labels_subset = cluster_labels(cluster_idx);
    if size(cluster_profiles, 1) < 5
        continue;
    end

    mean_profile = mean(cluster_profiles, 1,'omitnan');
    std_profile = std(cluster_profiles, 0, 1,'omitnan');
    n_profiles = size(cluster_profiles, 1);
    n_show = max(1, round(0.25 * n_profiles));
    show_idx = randperm(n_profiles, n_show);
    cluster_color = cmap(cluster_id,:);
    subplot(n_rows, n_cols, i);
    hold on;

    for k = 1:n_show
        plot(cluster_profiles(show_idx(k), :), DepInterval_custom, ...
            'Color', [0.6 0.6 0.6 0.3]);
    end

    plot(mean_profile, DepInterval_custom, 'Color',cluster_color, 'LineWidth', 3);
    plot(mean_profile + std_profile, DepInterval_custom,':', 'Color',cluster_color);
    plot(mean_profile - std_profile, DepInterval_custom, ':', 'Color',cluster_color);
    plot(mean_profile + 2*std_profile, DepInterval_custom, ':', 'Color',cluster_color);
    plot(mean_profile - 2*std_profile, DepInterval_custom, ':', 'Color',cluster_color);

    title(sprintf('Cluster %d', cluster_id));
    xlabel('Salinity');
    ylabel('Depth');
    axis ij;
end
sgtitle('Salinity by Cluster');


% Temp
figure
for i = 1:length(unique_clusters)
    nexttile;

    cluster_id = unique_clusters(i);
    cluster_idx = cluster_labels == cluster_id;
    %sal_cluster = with_sal_test(cluster_idx,starting_depth:ending_depth) ;
    temp_cluster = train_temp_mat(cluster_idx,:) ;
end


for i = 1:n_clusters
    cluster_id = unique_clusters(i);
    cluster_color = cmap(i, :);  
    % Get high-probability profiles for this cluster
    cluster_mask = (cluster_labels == cluster_id) ;
    cluster_profiles = train_temp_mat(cluster_mask, :);
    cluster_labels_subset = cluster_labels(cluster_idx);
    if size(cluster_profiles, 1) < 5
        continue;
    end

    mean_profile = mean(cluster_profiles, 1,'omitnan');
    std_profile = std(cluster_profiles, 0, 1,'omitnan');
    n_profiles = size(cluster_profiles, 1);
    n_show = max(1, round(0.25 * n_profiles));
    show_idx = randperm(n_profiles, n_show);
    cluster_color = cmap(cluster_id,:);
    subplot(n_rows, n_cols, i);
    hold on;

    for k = 1:n_show
        plot(cluster_profiles(show_idx(k), :), DepInterval_custom, ...
            'Color', [0.6 0.6 0.6 0.3]);
    end

    plot(mean_profile, DepInterval_custom, 'Color',cluster_color, 'LineWidth', 3);
    plot(mean_profile + std_profile, DepInterval_custom,':', 'Color',cluster_color);
    plot(mean_profile - std_profile, DepInterval_custom, ':', 'Color',cluster_color);
    plot(mean_profile + 2*std_profile, DepInterval_custom, ':', 'Color',cluster_color);
    plot(mean_profile - 2*std_profile, DepInterval_custom, ':', 'Color',cluster_color);

    title(sprintf('Cluster %d', cluster_id));
    xlabel('Temperature');
    ylabel('Depth');
    axis ij;
end
sgtitle('Temperature by Cluster');