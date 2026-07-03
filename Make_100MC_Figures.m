%% Make_100MC_Figures.m
%  Regenerates fig_method_comparison.png from the 100-MC results
%  (plos_revision_results.mat). Plotting code copied verbatim from
%  Tdoa_full_run.m; only the data source differs. No EKF runs.
%  Runtime: ~5 seconds.
%
%  NOTE: this is the ONLY figure that depends on the Monte Carlo count.
%  fig_cdf_errors uses a single fixed-seed representative run (rng 2042)
%  and does NOT need regeneration. fig_R1Q6_tag_geometry.png was already
%  produced by PLOS_Reviewer_Sims.m Part P5.

clear; clc; close all;
load('plos_revision_results.mat');   % res2_pos, res2_hdg (8 x 3 x 100)
res_pos = res2_pos;  res_hdg = res2_hdg;
fprintf('Loaded 100-MC results: %d runs\n', size(res_pos,3));

mnames  = {'GPS pos only','GPS+NavIC pos only',...
           '4-ant GPS att','4-ant GPS+NavIC att',...
           '1-tag UWB+GPS','Multi-tag UWB only',...
           'Multi-tag+GPS','Multi-tag+GPS+NavIC (proposed)'};
snames  = {'S1:Nominal','S2:Shadow','S3:Outage'};

% === BAR CHART — All methods x All scenarios (identical to Tdoa_full_run.m) ===
figure('Position',[100 100 1000 500]);
subplot(1,2,1);
pos_data = zeros(8,3);
for mi=1:8; for si=1:3; pos_data(mi,si)=mean(res_pos(mi,si,:)); end; end
sel = [1 3 4 5 6 7 8]; % GPS, GPS4ant, GPS+NavIC4ant, 1tag, MTonly, MT+GPS, Proposed
b = bar(pos_data(sel,:), 'grouped');
b(1).FaceColor = [0.3 0.6 0.9]; b(2).FaceColor = [0.9 0.6 0.2]; b(3).FaceColor = [0.8 0.2 0.2];
set(gca,'XTickLabel',mnames(sel),'XTickLabelRotation',35,'FontSize',8);
ylabel('Position RMSE (m)'); legend(snames,'Location','northwest','FontSize',7);
grid on; set(gca,'YScale','log'); ylim([0.5 300]);
title('(a) Position RMSE Comparison');

subplot(1,2,2);
hdg_data = zeros(8,3);
for mi=1:8; for si=1:3; hdg_data(mi,si)=mean(res_hdg(mi,si,:)); end; end
b = bar(hdg_data(sel,:), 'grouped');
b(1).FaceColor = [0.3 0.6 0.9]; b(2).FaceColor = [0.9 0.6 0.2]; b(3).FaceColor = [0.8 0.2 0.2];
set(gca,'XTickLabel',mnames(sel),'XTickLabelRotation',35,'FontSize',8);
ylabel('Heading RMSE (\circ)'); legend(snames,'Location','northwest','FontSize',7);
grid on; set(gca,'YScale','log'); ylim([0.05 200]);
title('(b) Heading RMSE Comparison');
saveas(gcf,'fig_method_comparison.png');
fprintf('Saved fig_method_comparison.png (100-MC data)\n');

% Sanity check against the manuscript's Table 5 values:
fprintf('\nSanity check (must match Table 5):\n');
fprintf('  Proposed S1: pos %.2f m, hdg %.3f deg  (expect 1.53 / 0.257)\n', ...
    mean(res_pos(8,1,:)), mean(res_hdg(8,1,:)));
fprintf('  1-tag   S1: pos %.2f m, hdg %.3f deg  (expect 38.34 / 61.741)\n', ...
    mean(res_pos(5,1,:)), mean(res_hdg(5,1,:)));
