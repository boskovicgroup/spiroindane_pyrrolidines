%% load the data and convert to arrays
ctrl1 = readtable('../zarkolab/publications/bmc_spiroindane_pyrrolidines/Ctrl1.csv');
ctrl1_array = table2array(ctrl1(:,4:end));
ctrl2 = readtable('../zarkolab/publications/bmc_spiroindane_pyrrolidines/Ctrl2.csv');
ctrl2_array = table2array(ctrl2(:,4:end));
ctrl3 = readtable('../zarkolab/publications/bmc_spiroindane_pyrrolidines/Ctrl3.csv');
ctrl3_array = table2array(ctrl3(:,4:end));


%% Calculate magnitudes of fingerprints and average per replicate of controls
ctrl1_magn = sqrt(sum(ctrl1_array.^2,2));
ctrl2_magn = sqrt(sum(ctrl2_array.^2,2));
ctrl3_magn = sqrt(sum(ctrl3_array.^2,2));
ctrl_magns = [ctrl1_magn; ctrl2_magn; ctrl3_magn];
ctrl_magns_ave = mean(ctrl_magns,1);
%% Remove outlier treatments from all replicates
ctrl1_array(ctrl1_magn > 1e3 | ctrl2_magn > 1e4 | ctrl3_magn > 1e4,:) = [];
ctrl2_array(ctrl1_magn > 1e3 | ctrl2_magn > 1e4 | ctrl3_magn > 1e4,:) = [];
ctrl3_array(ctrl1_magn > 1e3 | ctrl2_magn > 1e4 | ctrl3_magn > 1e4,:) = [];
%% Recalculate magnitudes
ctrl1_magn = sqrt(sum(ctrl1_array.^2,2));
ctrl2_magn = sqrt(sum(ctrl2_array.^2,2));
ctrl3_magn = sqrt(sum(ctrl3_array.^2,2));
ctrl_magns = [ctrl1_magn; ctrl2_magn; ctrl3_magn];
ctrl_magns_ave = mean(ctrl_magns,1);
%% Normalize by subtracting the mean per row
ctrl1_norm = [];
for i = 1:size(ctrl1_array,1)
    ctrl1_norm(i,:)= ctrl1_array(i,:)-mean(ctrl1_array(i,:));
end
ctrl2_norm = [];
for i = 1:size(ctrl2_array,1)
    ctrl2_norm(i,:)= ctrl2_array(i,:)-mean(ctrl2_array(i,:));
end
ctrl3_norm = [];
for i = 1:size(ctrl3_array,1)
    ctrl3_norm(i,:)= ctrl3_array(i,:)-mean(ctrl3_array(i,:));
end
%% How correlated are the control compounds between pairs of replicates?
% ctrl1 vs ctrl2
ctrl1_ctrl2_cor = [];
for i = 1:size(ctrl1_norm)
    ctrl1_ctrl2_cor(i) = ctrl1_norm(i,1:end)*ctrl2_norm(i,1:end)'/(norm(ctrl1_norm(i,1:end))*norm(ctrl2_norm(i,1:end)));
end;
% ctrl1 vs ctr3
ctrl1_ctrl3_cor = [];
for i = 1:size(ctrl1_norm)
    ctrl1_ctrl3_cor(i) = ctrl1_norm(i,1:end)*ctrl3_norm(i,1:end)'/(norm(ctrl1_norm(i,1:end))*norm(ctrl3_norm(i,1:end)));
end;
% ctrl 2 vs ctr3
ctrl2_ctrl3_cor = [];
for i = 1:size(ctrl1_norm)
    ctrl2_ctrl3_cor(i) = ctrl2_norm(i,1:end)*ctrl3_norm(i,1:end)'/(norm(ctrl2_norm(i,1:end))*norm(ctrl3_norm(i,1:end)));
end;

p1 = plot(ctrl1_magn, 'LineWidth',3);hold on;
p2 = plot(ctrl2_magn, 'LineWidth',3);
p3 = plot(ctrl3_magn, 'LineWidth',3); hold off;
axis([0 124 0 600]);
h = [p1;p2;p3];
legend(h, 'Controls, replicate 1', 'Controls, replicate 2', 'Controls, replicate 3');
xlabel('Treatment')
ylabel('Fingerprint magnitude')
title('Control compounds produce reproducible fingerprints')

mean(ctrl1_ctrl2_cor)
mean(ctrl2_ctrl3_cor)
mean(ctrl1_ctrl3_cor)


%% Perform SVD on control compounds
for i = 1:size(ctrl1_array,1)
    ctrl1_norm(i,:)= ctrl1_array(i,:)-mean(ctrl1_array(i,:));
end
[ctrl1_U, ctrl1_S, ctrl1_V] = svd(ctrl1_norm);


for i = 1:size(ctrl2_array,1)
    ctrl2_norm(i,:)= ctrl2_array(i,:)-mean(ctrl2_array(i,:));
end
[ctrl2_U, ctrl2_S, ctrl2_V] = svd(ctrl2_norm);


for i = 1:size(ctrl3_array,1)
    ctrl3_norm(i,:)= ctrl3_array(i,:)-mean(ctrl3_array(i,:));
end
[ctrl3_U, ctrl3_S, ctrl3_V] = svd(ctrl3_norm);

% Normalize the sigma values for controls
ctrl1_s_sum = sum(sum(ctrl1_S.^2));
ctrl1_s_norm = ctrl1_S.^2/ctrl1_s_sum;
ctrl2_s_norm(1,1)

ctrl2_s_sum = sum(sum(ctrl2_S.^2));
ctrl2_s_norm = ctrl2_S.^2/ctrl2_s_sum;
ctrl2_s_norm(1,1)

ctrl3_s_sum = sum(sum(ctrl3_S.^2));
ctrl3_s_norm = ctrl3_S.^2/ctrl3_s_sum;
ctrl3_s_norm(1,1)

%% Correlation between magnitudes and reproducibility
% get magnitude vectors and correlation vectors for controls
ctrl1_magn = sqrt(sum(ctrl1_norm.^2,2))
ctrl2_magn = sqrt(sum(ctrl2_norm.^2,2))
ctrl3_magn = sqrt(sum(ctrl3_norm.^2,2))
ctrl_magns = [ctrl1_magn ctrl2_magn ctrl3_magn]
ctrl_magns_ave = mean(ctrl_magns,2)
ctrl_magns_ave(ctrl1_magn > 1e3 | ctrl2_magn > 1e3 | ctrl3_magn > 1e3) = []
ctrl_cors = [ctrl1_ctrl2_cor; ctrl1_ctrl3_cor; ctrl2_ctrl3_cor];
ctrl_cors_ave = mean(ctrl_cors,1)
scatter(ctrl_cors_ave, ctrl_magns_ave)
xlabel('Reproducibility')
ylabel('Magnitude')
title('Control fingerprints')
%% Load the test compounds
rep1 = readtable('../zarkolab/publications/bmc_spiroindane_pyrrolidines/Rep1.csv');
rep1_array = table2array(rep1(:,4:end));
rep2 = readtable('../zarkolab/publications/bmc_spiroindane_pyrrolidines/Rep2.csv');
rep2_array = table2array(rep2(:,4:end));
% rep3 = readtable('../zarkolab/publications/bmc_spiroindane_pyrrolidines/Rep3.csv');
% rep3_array = table2array(rep3(:,4:end));

%% Calculate magnitudes for the test compounds
rep1_magn = sqrt(sum(rep1_array.^2,2));
rep2_magn = sqrt(sum(rep2_array.^2,2));
%rep3_magn = sqrt(sum(rep3_array.^2,2));
rep_magns = [rep1_magn rep2_magn rep3_magn];
rep_magns_ave = mean(rep_magns,2);
%% Check for outliers and remove from test compounds
rep1_array(rep1_magn > 1e3 | rep2_magn > 1e3,:) = [];
rep2_array(rep1_magn > 1e3 | rep2_magn > 1e3,:) = [];
%rep3_array(rep1_magn > 1e3 | rep2_magn > 1e3,:) = [];

%% Recalculate magnitudes
rep1_magn = sqrt(sum(rep1_array.^2,2));
rep2_magn = sqrt(sum(rep2_array.^2,2));
%rep3_magn = sqrt(sum(rep3_array.^2,2));
rep_magns = [rep1_magn rep2_magn];
rep_magns_ave = mean(rep_magns,2);
%% Normalization and SVD

for i = 1:size(rep1_array,1)
    rep1_norm(i,:)= rep1_array(i,:)-mean(rep1_array(i,:));
end
[rep1_U, rep1_S, repl1_V] = svd(rep1_norm);


for i = 1:size(rep2_array,1)
    rep2_norm(i,:)= rep2_array(i,:)-mean(rep2_array(i,:));
end
[rep2_U, rep2_S, rep2_V] = svd(rep2_norm);

% 
% for i = 1:size(rep3_array,1)
%     rep3_norm(i,:)= rep3_array(i,:)-mean(rep3_array(i,:));
% end
% [rep3_U, rep3_S, rep3_V] = svd(rep3_norm);

% Normalize the sigma values for controls
rep1_s_sum = sum(sum(rep1_S.^2));
rep1_s_norm = rep1_S.^2/rep1_s_sum;
rep1_s_norm(1,1)

rep2_s_sum = sum(sum(rep2_S.^2));
rep2_s_norm = rep2_S.^2/rep2_s_sum;
rep2_s_norm(1,1)

% rep3_s_sum = sum(sum(rep3_S.^2));
% rep3_s_norm = rep3_S.^2/rep3_s_sum;
% rep3_s_norm(1,1)

%% %% How are replicates correlated within test compounds?

% rep1 vs rep2
for i = 1:size(rep1_norm)
    rep1_rep2_cor(i) = rep1_norm(i,1:end)*rep2_norm(i,1:end)'/(norm(rep1_norm(i,1:end))*norm(rep2_norm(i,1:end)));
end;

% % rep1 vs rep3
% for i = 1:size(rep1_norm)
%     rep1_rep3_cor(i) = rep1_norm(i,1:end)*rep3_norm(i,1:end)'/(norm(rep1_norm(i,1:end))*norm(rep3_norm(i,1:end)));
% end;
% 
% % rep2 vs rep3
% for i = 1:size(rep1_norm)
%     rep2_rep3_cor(i) = rep2_norm(i,1:end)*rep3_norm(i,1:end)'/(norm(rep2_norm(i,1:end))*norm(rep3_norm(i,1:end)));
% end;

mean(rep1_rep2_cor)
%mean(rep2_rep3_cor)
%mean(rep1_rep3_cor)
rep_cors = [rep1_rep2_cor];
rep_cors_ave = mean(rep_cors, 1)';
%%
figure
hold on
plot(rep_cors_ave, rep_magns_ave, 'or', 'MarkerFaceColor','r');
plot(ctrl_cors_ave, ctrl_magns_ave, 'ob');
axis([0 1 0 400]);
legend('Test compounds', 'Controls');
text(rep_cors_ave(rep_magns_ave>50 & rep_cors_ave > 0.5)+0.01, rep_magns_ave(rep_magns_ave>50 & rep_cors_ave >0.5)+5, table2array(rep1(rep_magns_ave>50 & rep_cors_ave>0.5,2)));
xlabel('Fingerprint reproducibility');
ylabel('Fingerprint magnitude');
title('Compound activity')
%% Plot magnitudes
p1 = plot(rep1_magn, 'LineWidth', 3);hold on;
p2 = plot(rep2_magn, 'LineWidth', 3);
%p3 = plot(rep3_magn); hold off;
axis([0 181 0 300]);
h = [p1;p2];
legend(h, 'This collection, replicate 1', 'This collection, replicate 2');
xlabel('Treatment')
ylabel('Fingerprint magnitude')
%% What cut-off for correlation in all replicates to choose for activity?
% reproducible set of controls at 0.6 correlation at 3 replicates
reproducible_set_ctrl = ctrl1(ctrl1_ctrl2_cor > 0.6 & ctrl1_ctrl3_cor>0.6 & ctrl2_ctrl3_cor>0.6, 2:3)
reproducible_set_test = rep1(rep1_rep2_cor > 0.5 & rep1_rep3_cor > 0.5 & rep2_rep3_cor>0.5, 2:3)
