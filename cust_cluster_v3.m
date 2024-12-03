clear; close all;


%% Input data
data = load('D:\matlab\mmWaveRadar\time_distance.mat');
tdmat = data.time_distance;
non_zero_rows = any(tdmat ~= 0, 2);non_zero_cols = any(tdmat ~= 0, 1);
tdmat = tdmat(non_zero_rows, non_zero_cols);tdmat(:,1) = []; % Ignore first column
[r,c] = size(tdmat);
indices = repmat(1:r,1,c);readings = reshape(tdmat,1,[]);

figure(1);
plot(indices, readings, 'x', 'DisplayName', 'Original Data');
title('Unclassified');legend('show');grid on;


%% Finding Optimal K Value (Elbow Method)
k_values = 1:5;
wcss = zeros(size(k_values));

for k = k_values
    [~, ~, sumd] = kmeans(readings', k, 'Replicates', 5);
    wcss(k) = sum(sumd); % Sum of distances within clusters
end

wcss_norm = (wcss - min(wcss)) / (max(wcss) - min(wcss));
start_point = [1, wcss_norm(1)];
end_point = [length(k_values), wcss_norm(end)];

% Distances to straight line
distances = zeros(1, length(wcss_norm));
for i = 1:length(wcss_norm)
    point = [i, wcss_norm(i)];
    distances(i) = abs(det([end_point - start_point; point - start_point])) / norm(end_point - start_point);
end

[~, elbow_idx] = max(distances);
optimal_k = k_values(elbow_idx);


%% Performing Clustering
[idx_kmeans, C_kmeans] = kmeans(readings', optimal_k); % optimal_k
type_1 = zeros(1,r*c);type_2 = zeros(1,r*c);
type_1(idx_kmeans==1) = readings(idx_kmeans==1);type_2(idx_kmeans==2) = readings(idx_kmeans==2);
type_1(type_1==0) = nan;type_2(type_2==0) = nan;

figure(2);
gscatter(indices, readings, idx_kmeans,  'br', 'xx');
xlabel('Index');ylabel('Value');title('K-Means Clustered');
legend('Cluster 1', 'Cluster 2');legend('show');grid on;


%% Fitting Curves to Clusters (Experimental)
unique_indices = unique(indices);
agg_type_1 = arrayfun(@(x) mean(type_1(indices == x), 'omitnan'), unique_indices);
agg_type_2 = arrayfun(@(x) mean(type_2(indices == x), 'omitnan'), unique_indices);
interp_type_1 = interp1(unique_indices, agg_type_1, unique_indices, 'linear', 'extrap');
interp_type_2 = interp1(unique_indices, agg_type_2, unique_indices, 'linear', 'extrap');
filled_1 = fillmissing(interp_type_1, 'makima'); filled_2 = fillmissing(interp_type_2, 'makima'); % Fill method: spline, pchip, makima
sm1 = smoothdata(filled_1, "movmean", 15,'includenan'); sm2 = smoothdata(filled_2, "movmean", 15,'includenan'); % Inbuilt smooth func
f1 = find(~isnan(agg_type_1), 1, 'first');l1 = find(~isnan(agg_type_1), 1, 'last');
f2 = find(~isnan(agg_type_2), 1, 'first');l2 = find(~isnan(agg_type_2), 1, 'last');

figure(3);
plot(1:r, interp_type_1, 'x', 'DisplayName', 'Interpolated T1');hold on;
plot(1:r, interp_type_2, 'x', 'DisplayName', 'Interpolated T2');
plot(f1:l1, sm1(f1:l1), '-', 'DisplayName', 'Smoothened Line (T1)');
plot(f2:l2, sm2(f2:l2), '-', 'DisplayName', 'Smoothened Line (T2)');hold off;
legend('show');grid on;