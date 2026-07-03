%% TDOA_Attitude_NavIC_Fusion.m
%  ═══════════════════════════════════════════════════════════════════════════
%  Attitude Observability from Hyperbolic Range-Difference Geometry:
%  TDOA-Based Multi-Tag Rigid Body Pose with GPS+NavIC Augmentation
%
%  Parts:
%    A: Vessel trajectory + platform 6-DOF motion
%    B: GPS+NavIC satellite constellation + visibility model
%    C: Multi-antenna GNSS attitude model (4 antennas)
%    D: Multi-tag UWB TDOA measurement model
%    E: 12-state EKF (GNSS attitude + UWB attitude fusion)
%    F: 8-method head-to-head (3 scenarios × 30 MC)
%    G: Attitude observability verification (Theorems 1-3)
%    H: GDAP validation + tag placement analysis
%    I: GEO attitude anchoring demonstration
%    J: Complementary fusion bound validation
%    K: Platform compensation (4 strategies × 3 sea states)
%    L: Maritime RAIM (auto-calibrated)
%    M: UTIL validation with multi-tag injection
%
%  Real Datasets Used:
%    - NavIC: Kaggle NavIC_Dataset (86,400 epochs, 7 PRNs, 24h)
%    - GPS:   4-antenna RINEX 2.x (Master + 3 Slaves, Germany)
%    - UWB:   UTIL dataset parameters (σ=0.333m, 8 anchors)
%
%  Author: Navya B.R., REVA University
%  Date:   March 2026
%  ═══════════════════════════════════════════════════════════════════════════
clear; clc; close all;
rng(42);
tic_total = tic;  % Start timing
fprintf('╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║  TDOA Attitude + GPS+NavIC Fusion — FULL PAPER RUN              ║\n');
fprintf('║  (30 MC × 300s — expect ~2-4 hours)                             ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n\n');

%% ═══ LOAD REAL DATASET PARAMETERS ═══
use_real = false;
if exist('real_dataset_params.mat','file')
    load('real_dataset_params.mat');
    use_real = true;
    fprintf('  ✓ Loaded REAL parameters from NavIC + GPS datasets\n');
    fprintf('    NavIC: %d sats, σ_PR=%.3f m, C/N0=%.1f dB-Hz\n', ...
        real_par.navic_n_sats, real_par.navic_sig_pr, real_par.navic_mean_cn0);
    fprintf('    GPS:   σ_CP=%.4f m, baselines %.1f-%.1fm, %d sats\n', ...
        real_par.gps_sig_carrier, min(real_par.gps_baselines), ...
        max(real_par.gps_baselines), real_par.gps_n_sats);
    fprintf('\n');
else
    fprintf('  ⚠ real_dataset_params.mat not found\n');
    fprintf('    Run Real_Dataset_Integration.m first to use real data\n');
    fprintf('    Continuing with assumed parameters...\n\n');
end

%% ═══ PARAMETERS ═══
par.dt = 0.05; par.T = 300; par.N = round(par.T/par.dt);
par.t = (0:par.N-1)*par.dt; par.NMC = 30;

% Process noise
par.sig_a = 0.30; par.sig_vz = 0.10;
par.sig_om = 0.01; par.sig_rp = 0.005;

% GNSS position
par.sig_gnss = 1.5; par.gnss_rate = 1.0;
par.gnss_mp_infl = 6.0; par.gnss_mp_rng = 200; par.gnss_shadow_rng = 80;

% GNSS attitude (multi-antenna)
par.gnss_ant_body = [40 0 15; -40 0 8; 0 -10 12; 0 10 12]; % 4 antennas on ship
par.lambda_L1 = 0.1903;  % GPS L1 wavelength (m)

% === REAL DATA OVERRIDE: GPS carrier phase noise ===
if use_real
    par.sig_carrier = real_par.gps_sig_carrier;  % From real RINEX data
else
    par.sig_carrier = 0.003;  % Assumed value
end

% UWB TDOA
par.sig_tdoa = 0.50; par.pdrop = 0.10;
par.pnlos = 0.15; par.sig_nlos = 1.5;
par.uwb_max_rng = 600;

% === REAL DATA: NavIC pseudorange noise ===
if use_real
    par.navic_sig_pr = real_par.navic_sig_pr;  % From Kaggle NavIC dataset
else
    par.navic_sig_pr = 0.72;  % Typical NavIC L5
end

% EKF
par.kappa_max = 30; par.beta_geo = 0.5;
par.tau_nis = 200; par.gamma_nis = 5;

% Ship tags (body frame)
par.tags_body = [40 0 8; -40 0 4; 0 -15 12];
par.Ntags = 3;
par.gnss_body = [10; 0; 15];

% Platform anchors (dock body frame) — FPSO 160×50 m
par.anch_body = [0 0 5; 80 0 5; -80 0 5; 0 25 5; 0 -25 5; 40 12.5 35];
par.M = 6; par.ref = 1; par.plat_pos = [0;0;0];

% ═══ WAVE & VESSEL PARAMETERS (JONSWAP-based) ═══
% Sea state (default: moderate, SS4)
par.Hs = 1.5;          % Significant wave height (m)
par.Tp = 10;           % Peak wave period (s)
par.gamma_js = 3.3;    % JONSWAP peakedness (3.3 = North Sea)
par.n_waves = 30;      % Frequency components for spectral decomposition
par.wave_dir = 0;      % Dominant wave direction (rad, 0=head seas)
par.wave_spread = 15;  % Directional spreading (deg)

% Dock vessel parameters (FPSO-type, 160m × 50m × 15m draft)
par.dock_L = 160;      % Length (m)
par.dock_B = 50;       % Beam (m)
par.dock_T = 15;       % Draft (m)
par.dock_disp = 80000; % Displacement (tonnes)
par.dock_GM   = 4.0;   % Transverse metacentric height (m) — typical FPSO
par.dock_GML  = 200;   % Longitudinal metacentric height (m)
par.dock_Cb   = par.dock_disp / (1.025 * par.dock_L * par.dock_B * par.dock_T); % Block coefficient
par.dock_drift = 0.3;  % Mooring drift amplitude (m)

% Dock natural periods — COMPUTED from vessel parameters
rho_sw = 1025;  % seawater density (kg/m³)
g = 9.81;
dock_m = par.dock_disp * 1000;  % kg
dock_Aw = par.dock_Cb * par.dock_L * par.dock_B;  % waterplane area (m²)
dock_a33 = 0.6 * dock_m;  % heave added mass ≈ 60% of displacement
par.dock_Tn_heave = 2*pi*sqrt((dock_m + dock_a33) / (rho_sw * g * dock_Aw));
dock_k_roll = 0.4 * par.dock_B;  % roll radius of gyration ≈ 0.4·B
par.dock_Tn_roll  = 2*pi * dock_k_roll / sqrt(g * par.dock_GM);
dock_k_pitch = 0.25 * par.dock_L;  % pitch radius of gyration ≈ 0.25·L
par.dock_Tn_pitch = 2*pi * dock_k_pitch / sqrt(g * par.dock_GML);
par.dock_zeta_heave = 0.10; % Damping ratio: heave (viscous + radiation)
par.dock_zeta_roll  = 0.05; % Damping ratio: roll (bilge keels etc.)
par.dock_zeta_pitch = 0.08; % Damping ratio: pitch

% Ship vessel parameters (supply vessel ~80m × 18m × 5m draft)
par.ship_L = 80;       % Length (m)
par.ship_B = 18;       % Beam (m)
par.ship_T = 5;        % Draft (m)
par.ship_disp = 4000;  % Displacement (tonnes)
par.ship_GM   = 1.5;   % Transverse metacentric height (m) — supply vessel
par.ship_GML  = 80;    % Longitudinal metacentric height (m)
par.ship_Cb   = par.ship_disp / (1.025 * par.ship_L * par.ship_B * par.ship_T);

% Ship natural periods — COMPUTED from vessel parameters
ship_m = par.ship_disp * 1000;  % kg
ship_Aw = par.ship_Cb * par.ship_L * par.ship_B;
ship_a33 = 0.5 * ship_m;  % heave added mass ≈ 50% for smaller vessel
par.ship_Tn_heave = 2*pi*sqrt((ship_m + ship_a33) / (rho_sw * g * ship_Aw));
ship_k_roll = 0.4 * par.ship_B;
par.ship_Tn_roll  = 2*pi * ship_k_roll / sqrt(g * par.ship_GM);
ship_k_pitch = 0.25 * par.ship_L;
par.ship_Tn_pitch = 2*pi * ship_k_pitch / sqrt(g * par.ship_GML);
par.ship_zeta_heave = 0.12; % Damping ratio: heave
par.ship_zeta_roll  = 0.06; % Damping ratio: roll
par.ship_zeta_pitch = 0.10; % Damping ratio: pitch

fprintf('  Vessel natural periods (computed from mass/dimensions):\n');
fprintf('    Dock (%dt, %dm×%dm): Tn_heave=%.1fs Tn_roll=%.1fs Tn_pitch=%.1fs\n',...
    par.dock_disp, par.dock_L, par.dock_B, ...
    par.dock_Tn_heave, par.dock_Tn_roll, par.dock_Tn_pitch);
fprintf('    Ship (%dt, %dm×%dm):  Tn_heave=%.1fs Tn_roll=%.1fs Tn_pitch=%.1fs\n',...
    par.ship_disp, par.ship_L, par.ship_B, ...
    par.ship_Tn_heave, par.ship_Tn_roll, par.ship_Tn_pitch);

% Location: Indian Ocean (19°N, 72°E — off Mumbai)
par.lat = 19; par.lon = 72;

fprintf('Tags:%d Anchors:%d MC:%d Location:%.0f°N,%.0f°E\n',...
    par.Ntags,par.M,par.NMC,par.lat,par.lon);
fprintf('σ_carrier=%.4fm (source:%s) | σ_NavIC_PR=%.3fm (source:%s)\n\n',...
    par.sig_carrier, tern(use_real,'REAL RINEX','assumed'), ...
    par.navic_sig_pr, tern(use_real,'REAL Kaggle','assumed'));

%% ═══ PART A: TRAJECTORY + PLATFORM ═══
fprintf('--- Part A: Trajectory + Platform (JONSWAP wave model) ---\n');
% Generate wave field ONCE — shared between ship and dock (same ocean!)
waves = gen_waves(par);
traj = gen_traj(par, waves);
plat = gen_plat(par, waves);
fprintf('  Ship: (%.0f,%.0f)→(%.1f,%.1f) m  Heading: %.1f°→%.1f°\n',...
    traj.pos(1,1),traj.pos(2,1),traj.pos(1,end),traj.pos(2,end),...
    rad2deg(traj.att(1,1)),rad2deg(traj.att(1,end)));
fprintf('  Sea state: Hs=%.1fm Tp=%.0fs (JONSWAP γ=%.1f)\n',...
    par.Hs, par.Tp, par.gamma_js);
fprintf('  Ship  motion — heave: ±%.2fm  roll: ±%.1f°  pitch: ±%.1f°\n',...
    max(abs(traj.pos(3,:))), rad2deg(max(abs(traj.att(2,:)))), ...
    rad2deg(max(abs(traj.att(3,:)))));
fprintf('  Dock  motion — heave: ±%.2fm  roll: ±%.1f°  pitch: ±%.1f°\n',...
    max(abs(plat.hv)), rad2deg(max(abs(plat.rl))), ...
    rad2deg(max(abs(plat.pt))));

%% ═══ PART B: GPS+NavIC SATELLITE MODEL ═══
fprintf('\n--- Part B: GPS+NavIC Satellite Constellation ---\n');

% === REAL DATA OVERRIDE: Use real NavIC az/el if available ===
if use_real && exist('sat_navic_real','var') && ~isempty(sat_navic_real)
    [sat_gps, ~] = gen_satellites(par);  % GPS still from model
    sat_navic = sat_navic_real;          % NavIC from Kaggle dataset
    fprintf('  GPS satellites: %d (modelled)\n', size(sat_gps,1));
    fprintf('  NavIC satellites: %d (FROM REAL KAGGLE DATA)\n', size(sat_navic,1));
    for si = 1:size(sat_navic,1)
        fprintf('    NavIC PRN: az=%.1f° el=%.1f°\n', sat_navic(si,1), sat_navic(si,2));
    end
else
    [sat_gps, sat_navic] = gen_satellites(par);
    fprintf('  GPS satellites: %d | NavIC satellites: %d (modelled)\n', ...
        size(sat_gps,1), size(sat_navic,1));
end

% Compute visibility vs distance to dock
dists_vis = [500 400 300 200 150 100 80 50 30 10];
nvis_gps = zeros(size(dists_vis));
nvis_navic = zeros(size(dists_vis));
nvis_total = zeros(size(dists_vis));
pdop_gps = zeros(size(dists_vis));
pdop_total = zeros(size(dists_vis));

for di = 1:length(dists_vis)
    d = dists_vis(di);
    mask = compute_dock_mask(d, par);
    [ng,pg] = count_visible(sat_gps, mask, par);
    [nn,~]  = count_visible(sat_navic, mask, par);
    [nt,pt] = count_visible([sat_gps; sat_navic], mask, par);
    nvis_gps(di) = ng; nvis_navic(di) = nn;
    nvis_total(di) = nt; pdop_gps(di) = pg; pdop_total(di) = pt;
end

fprintf('  At 500m: GPS=%d NavIC=%d Total=%d PDOP=%.1f\n',...
    nvis_gps(1),nvis_navic(1),nvis_total(1),pdop_total(1));
fprintf('  At  80m: GPS=%d NavIC=%d Total=%d PDOP=%.1f\n',...
    nvis_gps(7),nvis_navic(7),nvis_total(7),pdop_total(7));
fprintf('  At  30m: GPS=%d NavIC=%d Total=%d PDOP=%.1f\n',...
    nvis_gps(9),nvis_navic(9),nvis_total(9),pdop_total(9));

% Fig 1: Satellite visibility
figure('Position',[100 100 900 400]);
subplot(1,2,1);
plot(dists_vis, nvis_gps, 'b-o', 'LineWidth',2,'MarkerFaceColor','b'); hold on;
plot(dists_vis, nvis_total, 'r-s', 'LineWidth',2,'MarkerFaceColor','r');
plot(dists_vis, nvis_navic, 'g--^', 'LineWidth',1.5,'MarkerFaceColor','g');
yline(5,'k:','Min for attitude','LineWidth',1.5);
set(gca,'XDir','reverse','FontSize',11);
xlabel('Distance to dock (m)'); ylabel('Visible satellites');
legend('GPS only','GPS+NavIC','NavIC only','Location','southwest');
title('Satellite Visibility vs Distance'); grid on;
subplot(1,2,2);
plot(dists_vis, pdop_gps, 'b-o', 'LineWidth',2,'MarkerFaceColor','b'); hold on;
plot(dists_vis, pdop_total, 'r-s', 'LineWidth',2,'MarkerFaceColor','r');
set(gca,'XDir','reverse','FontSize',11);
xlabel('Distance to dock (m)'); ylabel('PDOP');
legend('GPS only','GPS+NavIC','Location','northeast');
title('PDOP vs Distance'); grid on;
saveas(gcf,'fig_satellite_visibility.png');

%% ═══ PART C: MULTI-ANTENNA GNSS ATTITUDE MODEL ═══
fprintf('\n--- Part C: Multi-antenna GNSS attitude ---\n');

dists_att = linspace(500, 10, 50);
sig_hdg_gps = zeros(size(dists_att));
sig_hdg_navic = zeros(size(dists_att));
amb_success_gps = zeros(size(dists_att));
amb_success_navic = zeros(size(dists_att));

B_surge = norm(par.gnss_ant_body(1,:)-par.gnss_ant_body(2,:));
for di = 1:length(dists_att)
    d = dists_att(di);
    mask = compute_dock_mask(d, par);
    [ng,~] = count_visible(sat_gps, mask, par);
    [nt,~] = count_visible([sat_gps;sat_navic], mask, par);
    mp = 1 + 5*exp(-d/80);
    sig_cp = par.sig_carrier * mp;
    sig_hdg_gps(di) = rad2deg(sig_cp / (B_surge * sqrt(max(ng,1))));
    sig_hdg_navic(di) = rad2deg(sig_cp / (B_surge * sqrt(max(nt,1))));
    amb_success_gps(di) = max(0, 1 - exp(-0.8*(ng-5)));
    amb_success_navic(di) = max(0, 1 - exp(-0.8*(nt-5)));
end

fprintf('  Baseline: %.0f m | σ_carrier=%.4f m (source:%s)\n',...
    B_surge, par.sig_carrier, tern(use_real,'REAL','assumed'));
fprintf('  At 500m: σ_ψ=%.4f° | At 80m: σ_ψ=%.4f°(GPS) %.4f°(+NavIC)\n',...
    sig_hdg_gps(1), sig_hdg_gps(find(dists_att<=80,1)), sig_hdg_navic(find(dists_att<=80,1)));

% Fig 2: GNSS attitude accuracy
figure('Position',[100 100 900 400]);
subplot(1,2,1);
semilogy(dists_att, sig_hdg_gps, 'b-', 'LineWidth',2); hold on;
semilogy(dists_att, sig_hdg_navic, 'r-', 'LineWidth',2);
yline(1,'k:','1° target','LineWidth',1.5);
set(gca,'XDir','reverse','FontSize',11);
xlabel('Distance to dock (m)'); ylabel('Heading σ (°)');
legend('GPS 4-ant','GPS+NavIC 4-ant','Location','northwest'); grid on;
title('GNSS Attitude Accuracy');
subplot(1,2,2);
plot(dists_att, amb_success_gps*100, 'b-', 'LineWidth',2); hold on;
plot(dists_att, amb_success_navic*100, 'r-', 'LineWidth',2);
yline(99,'k:','99% target','LineWidth',1.5);
set(gca,'XDir','reverse','FontSize',11);
xlabel('Distance to dock (m)'); ylabel('Ambiguity resolution success (%)');
legend('GPS','GPS+NavIC','Location','southwest'); grid on;
title('Integer Ambiguity Resolution');
saveas(gcf,'fig_gnss_attitude.png');

%% ═══ PART D/E: MEASUREMENT + EKF READY ═══
fprintf('\n--- Parts D/E: Measurement + EKF engines ready ---\n');
mtest = gen_meas(traj, plat, par, 1, sat_gps, sat_navic, 'gps_navic');
ng = sum([mtest.gnss.valid]); nu = sum([mtest.uwb.valid]);
na = sum([mtest.gnss_att.valid]);
fprintf('  S1: %d GNSS pos, %d GNSS att, %d UWB epochs\n', ng, na, nu);

%% ═══ PART F: 8-METHOD HEAD-TO-HEAD ═══
fprintf('\n--- Part F: 8-method Head-to-Head Monte Carlo ---\n');
burn = round(2/par.dt); NMC = par.NMC;
fprintf('  (8 methods × 3 scenarios × %d MC × %.0fs = %d EKF runs)\n', ...
    NMC, par.T, 8*3*NMC);
tic_F = tic;

methods = {'gps_pos','gpsnav_pos','gps_4ant','gpsnav_4ant',...
           'uwb1tag_gps','mtag_only','mtag_gps','mtag_gpsnav'};
mnames  = {'GPS pos only','GPS+NavIC pos only',...
           '4-ant GPS att','4-ant GPS+NavIC att',...
           '1-tag UWB+GPS','Multi-tag UWB only',...
           'Multi-tag+GPS','Multi-tag+GPS+NavIC (proposed)'};
snames  = {'S1:Nominal','S2:Shadow','S3:Outage'};

res_pos = zeros(8,3,NMC);
res_hdg = zeros(8,3,NMC);
res_roll = zeros(8,3,NMC);
res_pitch = zeros(8,3,NMC);

for si = 1:3
    fprintf('  Scenario %d:\n', si);
    for mc = 1:NMC
        if mod(mc,10)==0; fprintf('    MC %d/%d\n', mc, NMC); end
        rng(mc*1000+si);
        for mi = 1:8
            if any(mi==[2 4 8]); gnss_mode='gps_navic'; else; gnss_mode='gps_only'; end
            mm = gen_meas(traj, plat, par, si, sat_gps, sat_navic, gnss_mode);
            [est,~] = run_ekf(traj, plat, mm, par, methods{mi}, 'perfect');
            ep = traj.pos(:,burn+1:end) - est.pos(:,burn+1:end);
            res_pos(mi,si,mc) = sqrt(mean(sum(ep.^2,1)));
            eh = wrapToPi(traj.att(1,burn+1:end) - est.att(1,burn+1:end));
            res_hdg(mi,si,mc) = rad2deg(sqrt(mean(eh.^2)));
            er = wrapToPi(traj.att(2,burn+1:end) - est.att(2,burn+1:end));
            res_roll(mi,si,mc) = rad2deg(sqrt(mean(er.^2)));
            ept = wrapToPi(traj.att(3,burn+1:end) - est.att(3,burn+1:end));
            res_pitch(mi,si,mc) = rad2deg(sqrt(mean(ept.^2)));
        end
    end
end

% Print tables
fprintf('\n╔══════════════════════════════════════════════════════════════════════╗\n');
fprintf('║          POSITION RMSE (metres, %d MC)                              ║\n', NMC);
fprintf('╠══════════════════════════════════════════════════════════════════════╣\n');
for mi=1:8
    fprintf('║ %-28s│',mnames{mi});
    for si=1:3; fprintf(' %6.2f±%.2f │',mean(res_pos(mi,si,:)),std(res_pos(mi,si,:))); end
    fprintf('\n');
end
fprintf('╚══════════════════════════════════════════════════════════════════════╝\n');

fprintf('\n╔══════════════════════════════════════════════════════════════════════╗\n');
fprintf('║          HEADING RMSE (degrees, %d MC)                              ║\n', NMC);
fprintf('╠══════════════════════════════════════════════════════════════════════╣\n');
for mi=1:8
    fprintf('║ %-28s│',mnames{mi});
    for si=1:3; fprintf(' %6.2f±%.2f │',mean(res_hdg(mi,si,:)),std(res_hdg(mi,si,:))); end
    fprintf('\n');
end
fprintf('╚══════════════════════════════════════════════════════════════════════╝\n');

fprintf('\n╔══════════════════════════════════════════════════════════════════════╗\n');
fprintf('║          ROLL/PITCH RMSE (degrees, %d MC)                           ║\n', NMC);
fprintf('╠══════════════════════════════════════════════════════════════════════╣\n');
for mi=[6 7 8]  % only multi-tag methods can estimate roll/pitch
    fprintf('║ %-28s│',mnames{mi});
    for si=1:3
        fprintf(' R:%.2f P:%.2f │', mean(res_roll(mi,si,:)), mean(res_pitch(mi,si,:)));
    end
    fprintf('\n');
end
fprintf('╚══════════════════════════════════════════════════════════════════════╝\n');

% Full 6-DOF summary for proposed method
fprintf('\n  ┌── PROPOSED METHOD: Full 6-DOF Summary ──────────────────────────┐\n');
for si=1:3
    fprintf('  │ %s: pos=%.2fm hdg=%.2f° roll=%.2f° pitch=%.2f° │\n', ...
        snames{si}, mean(res_pos(8,si,:)), mean(res_hdg(8,si,:)), ...
        mean(res_roll(8,si,:)), mean(res_pitch(8,si,:)));
end
fprintf('  └──────────────────────────────────────────────────────────────────┘\n');

for si=2:3
    p_1tag = mean(res_pos(5,si,:)); p_prop = mean(res_pos(8,si,:));
    h_gps4 = mean(res_hdg(3,si,:)); h_nav4 = mean(res_hdg(4,si,:));
    fprintf('\n  %s: Pos vs 1-tag: %.1f%% | Hdg GPS4ant vs +NavIC: %.2f° vs %.2f°\n',...
        snames{si}, (p_1tag-p_prop)/p_1tag*100, h_gps4, h_nav4);
end
fprintf('\n  Part F completed in %.1f minutes\n', toc(tic_F)/60);

% ═══ PUBLICATION FIGURES (IEEE standard) ═══

% === FIG 3: BAR CHART — All methods × All scenarios (THE comparison figure) ===
figure('Position',[100 100 1000 500]);
% Position RMSE bar chart
subplot(1,2,1);
pos_data = zeros(8,3);
for mi=1:8; for si=1:3; pos_data(mi,si)=mean(res_pos(mi,si,:)); end; end
% Only show key methods (skip redundant ones)
sel = [1 3 4 5 6 7 8]; % GPS, GPS4ant, GPS+NavIC4ant, 1tag, MTonly, MT+GPS, Proposed
b = bar(pos_data(sel,:), 'grouped');
b(1).FaceColor = [0.3 0.6 0.9]; b(2).FaceColor = [0.9 0.6 0.2]; b(3).FaceColor = [0.8 0.2 0.2];
set(gca,'XTickLabel',mnames(sel),'XTickLabelRotation',35,'FontSize',8);
ylabel('Position RMSE (m)'); legend(snames,'Location','northwest','FontSize',7);
grid on; set(gca,'YScale','log'); ylim([0.5 300]);
title('(a) Position RMSE Comparison');

% Heading RMSE bar chart
subplot(1,2,2);
hdg_data = zeros(8,3);
for mi=1:8; for si=1:3; hdg_data(mi,si)=mean(res_hdg(mi,si,:)); end; end
b = bar(hdg_data(sel,:), 'grouped');
b(1).FaceColor = [0.3 0.6 0.9]; b(2).FaceColor = [0.9 0.6 0.2]; b(3).FaceColor = [0.8 0.2 0.2];
set(gca,'XTickLabel',mnames(sel),'XTickLabelRotation',35,'FontSize',8);
ylabel('Heading RMSE (°)'); legend(snames,'Location','northwest','FontSize',7);
grid on; set(gca,'YScale','log'); ylim([0.05 200]);
title('(b) Heading RMSE Comparison');
saveas(gcf,'fig_method_comparison.png');

% === FIG 4: CDF — Position and Heading error (standard in IEEE papers) ===
% Run one MC for each key method in S1, store full error vectors
figure('Position',[100 100 1000 400]);
rng(2042);
cdf_methods = [3 4 5 6 7 8]; % GPS4ant, +NavIC4ant, 1tag, MTonly, MT+GPS, Proposed
cdf_colors = [0 0.4 0.8; 0 0.7 0.3; 0.8 0.5 0; 0.6 0 0.6; 0.2 0.2 0.2; 1 0 0];
cdf_styles = {'--','-.',':','--','-','-'};
cdf_widths = [1.2 1.2 1.2 1.2 1.5 2.5];

subplot(1,2,1); hold on;
for ci = 1:length(cdf_methods)
    mi = cdf_methods(ci);
    if any(mi==[2 4 8]); gm='gps_navic'; else; gm='gps_only'; end
    mm_c = gen_meas(traj,plat,par,1,sat_gps,sat_navic,gm);
    [est_c,~] = run_ekf(traj,plat,mm_c,par,methods{mi},'perfect');
    pos_err = sqrt(sum((traj.pos(:,burn+1:end)-est_c.pos(:,burn+1:end)).^2,1));
    [f,x] = ecdf(pos_err);
    plot(x, f*100, cdf_styles{ci}, 'Color', cdf_colors(ci,:), 'LineWidth', cdf_widths(ci));
end
xlabel('Position Error (m)'); ylabel('CDF (%)');
xline(2,'k:','2m target','FontSize',8);
legend(mnames(cdf_methods),'Location','southeast','FontSize',7);
grid on; xlim([0 30]); set(gca,'FontSize',10);
title('(a) Position Error CDF — S1:Nominal');

subplot(1,2,2); hold on;
for ci = 1:length(cdf_methods)
    mi = cdf_methods(ci);
    if any(mi==[2 4 8]); gm='gps_navic'; else; gm='gps_only'; end
    mm_c = gen_meas(traj,plat,par,1,sat_gps,sat_navic,gm);
    [est_c,~] = run_ekf(traj,plat,mm_c,par,methods{mi},'perfect');
    hdg_err = rad2deg(abs(wrapToPi(traj.att(1,burn+1:end)-est_c.att(1,burn+1:end))));
    [f,x] = ecdf(hdg_err);
    plot(x, f*100, cdf_styles{ci}, 'Color', cdf_colors(ci,:), 'LineWidth', cdf_widths(ci));
end
xlabel('Heading Error (°)'); ylabel('CDF (%)');
xline(1,'k:','1° target','FontSize',8);
legend(mnames(cdf_methods),'Location','southeast','FontSize',7);
grid on; xlim([0 10]); set(gca,'FontSize',10);
title('(b) Heading Error CDF — S1:Nominal');
saveas(gcf,'fig_cdf_errors.png');

% === FIG 5: 6-DOF TIME-SERIES — Scenario 3 with all attitude angles ===
% Shows position + heading + roll + pitch errors over time
% KEY: demonstrates that TDOA gives ALL THREE attitude angles
rng(3042);
mm_ts = gen_meas(traj,plat,par,3,sat_gps,sat_navic,'gps_navic');
% Run 4 key methods for comparison
ts_mi = [3 4 6 8]; % GPS4ant, +NavIC, MTonly, Proposed
ts_names = mnames(ts_mi);
ts_colors = [0 0.4 0.8; 0 0.7 0.3; 0.6 0 0.6; 1 0 0];
ts_widths = [1.2 1.2 1.2 2.0];

figure('Position',[100 100 1000 800]);
ts_results = cell(1,4);
for ci = 1:4
    mi = ts_mi(ci);
    if any(mi==[2 4 8]); gm='gps_navic'; else; gm='gps_only'; end
    mm_t = gen_meas(traj,plat,par,3,sat_gps,sat_navic,gm);
    [ts_results{ci},~] = run_ekf(traj,plat,mm_t,par,methods{mi},'perfect');
end

% Position error
subplot(4,1,1); hold on;
for ci=1:4
    plot(par.t, sqrt(sum((traj.pos-ts_results{ci}.pos).^2,1)), '-', ...
        'Color', ts_colors(ci,:), 'LineWidth', ts_widths(ci));
end
ylabel('Pos Error (m)'); grid on;
lg = legend(ts_names,'Location','northeast','FontSize',7); set(lg,'AutoUpdate','off');
yl=ylim;
patch([120 180 180 120],[0 0 yl(2) yl(2)],[.5 .5 .5],'FaceAlpha',.1,'EdgeColor','none','HandleVisibility','off');
patch([180 300 300 180],[0 0 yl(2) yl(2)],[1 .3 .3],'FaceAlpha',.1,'EdgeColor','none','HandleVisibility','off');
text(150,yl(2)*0.85,'Shadow','FontSize',9,'HorizontalAlignment','center','FontWeight','bold');
text(240,yl(2)*0.85,'GNSS Outage','FontSize',9,'HorizontalAlignment','center','Color','r','FontWeight','bold');
title('Scenario 3: Full 6-DOF Error Time-Series'); set(gca,'FontSize',10);

% Heading error
subplot(4,1,2); hold on;
for ci=1:4
    plot(par.t, rad2deg(abs(wrapToPi(traj.att(1,:)-ts_results{ci}.att(1,:)))), '-', ...
        'Color', ts_colors(ci,:), 'LineWidth', ts_widths(ci));
end
ylabel('Heading Error (°)'); grid on; ylim([0 min(120, max(ylim))]);
yl=ylim;
patch([120 180 180 120],[0 0 yl(2) yl(2)],[.5 .5 .5],'FaceAlpha',.1,'EdgeColor','none','HandleVisibility','off');
patch([180 300 300 180],[0 0 yl(2) yl(2)],[1 .3 .3],'FaceAlpha',.1,'EdgeColor','none','HandleVisibility','off');
set(gca,'FontSize',10);

% Roll error
subplot(4,1,3); hold on;
for ci=1:4
    plot(par.t, rad2deg(abs(wrapToPi(traj.att(2,:)-ts_results{ci}.att(2,:)))), '-', ...
        'Color', ts_colors(ci,:), 'LineWidth', ts_widths(ci));
end
ylabel('Roll Error (°)'); grid on; ylim([0 min(40, max(ylim))]);
yl=ylim;
patch([120 180 180 120],[0 0 yl(2) yl(2)],[.5 .5 .5],'FaceAlpha',.1,'EdgeColor','none','HandleVisibility','off');
patch([180 300 300 180],[0 0 yl(2) yl(2)],[1 .3 .3],'FaceAlpha',.1,'EdgeColor','none','HandleVisibility','off');
set(gca,'FontSize',10);

% Pitch error
subplot(4,1,4); hold on;
for ci=1:4
    plot(par.t, rad2deg(abs(wrapToPi(traj.att(3,:)-ts_results{ci}.att(3,:)))), '-', ...
        'Color', ts_colors(ci,:), 'LineWidth', ts_widths(ci));
end
ylabel('Pitch Error (°)'); xlabel('Time (s)'); grid on; ylim([0 min(30, max(ylim))]);
yl=ylim;
patch([120 180 180 120],[0 0 yl(2) yl(2)],[.5 .5 .5],'FaceAlpha',.1,'EdgeColor','none','HandleVisibility','off');
patch([180 300 300 180],[0 0 yl(2) yl(2)],[1 .3 .3],'FaceAlpha',.1,'EdgeColor','none','HandleVisibility','off');
set(gca,'FontSize',10);
saveas(gcf,'fig_6dof_timeseries.png');

% === FIG 6: ATTITUDE FUSION — The KEY figure proving fusion benefit ===
% Use S3 (outage) to show what happens when GNSS fails:
%   - GNSS attitude: works early, then FAILS completely
%   - UWB multi-tag: works throughout (from TDOA geometry)
%   - Fused: best of both + robust
rng(4042);
mm_af = gen_meas(traj,plat,par,3,sat_gps,sat_navic,'gps_navic');
[est_gnss4,~] = run_ekf(traj,plat,mm_af,par,'gpsnav_4ant','perfect');
mm_af2 = gen_meas(traj,plat,par,3,sat_gps,sat_navic,'gps_navic');
[est_uwb,~]  = run_ekf(traj,plat,mm_af2,par,'mtag_only','perfect');
mm_af3 = gen_meas(traj,plat,par,3,sat_gps,sat_navic,'gps_navic');
[est_fused,~]= run_ekf(traj,plat,mm_af3,par,'mtag_gpsnav','perfect');

figure('Position',[100 100 1000 700]);
% Heading: plot UNWRAPPED angles to avoid ±180° jump
subplot(3,1,1); hold on;
true_hdg_uw = rad2deg(unwrap(traj.att(1,:)));
gnss_hdg_uw = rad2deg(unwrap(est_gnss4.att(1,:)));
uwb_hdg_uw  = rad2deg(unwrap(est_uwb.att(1,:)));
fused_hdg_uw = rad2deg(unwrap(est_fused.att(1,:)));
plot(par.t, true_hdg_uw, 'k-', 'LineWidth',1.2);
plot(par.t, gnss_hdg_uw, 'b-', 'LineWidth',1.0);
plot(par.t, uwb_hdg_uw, '-', 'Color',[0.6 0 0.6], 'LineWidth',1.0);
plot(par.t, fused_hdg_uw, 'r-', 'LineWidth',1.5);
ylabel('Heading ψ (°)'); grid on;
lg = legend('True','GNSS 4-ant (GPS+NavIC)','UWB multi-tag only','Proposed (fused)',...
    'Location','best','FontSize',7); set(lg,'AutoUpdate','off');
yl=ylim;
patch([120 180 180 120],[yl(1) yl(1) yl(2) yl(2)],[.5 .5 .5],'FaceAlpha',.08,'EdgeColor','none','HandleVisibility','off');
patch([180 300 300 180],[yl(1) yl(1) yl(2) yl(2)],[1 .3 .3],'FaceAlpha',.08,'EdgeColor','none','HandleVisibility','off');
title('Attitude Tracking: GNSS vs UWB vs Fused — S3:GNSS Outage'); set(gca,'FontSize',10);

% Roll: plot actual angles (small enough to see)
subplot(3,1,2); hold on;
plot(par.t, rad2deg(traj.att(2,:)), 'k-', 'LineWidth',1.2);
plot(par.t, rad2deg(est_gnss4.att(2,:)), 'b-', 'LineWidth',1.0);
plot(par.t, rad2deg(est_uwb.att(2,:)), '-', 'Color',[0.6 0 0.6], 'LineWidth',1.0);
plot(par.t, rad2deg(est_fused.att(2,:)), 'r-', 'LineWidth',1.5);
ylabel('Roll φ (°)'); grid on; set(gca,'FontSize',10);
yl=ylim;
patch([120 180 180 120],[yl(1) yl(1) yl(2) yl(2)],[.5 .5 .5],'FaceAlpha',.08,'EdgeColor','none','HandleVisibility','off');
patch([180 300 300 180],[yl(1) yl(1) yl(2) yl(2)],[1 .3 .3],'FaceAlpha',.08,'EdgeColor','none','HandleVisibility','off');

% Pitch: plot actual angles
subplot(3,1,3); hold on;
plot(par.t, rad2deg(traj.att(3,:)), 'k-', 'LineWidth',1.2);
plot(par.t, rad2deg(est_gnss4.att(3,:)), 'b-', 'LineWidth',1.0);
plot(par.t, rad2deg(est_uwb.att(3,:)), '-', 'Color',[0.6 0 0.6], 'LineWidth',1.0);
plot(par.t, rad2deg(est_fused.att(3,:)), 'r-', 'LineWidth',1.5);
ylabel('Pitch θ (°)'); xlabel('Time (s)'); grid on; set(gca,'FontSize',10);
yl=ylim;
patch([120 180 180 120],[yl(1) yl(1) yl(2) yl(2)],[.5 .5 .5],'FaceAlpha',.08,'EdgeColor','none','HandleVisibility','off');
patch([180 300 300 180],[yl(1) yl(1) yl(2) yl(2)],[1 .3 .3],'FaceAlpha',.08,'EdgeColor','none','HandleVisibility','off');
saveas(gcf,'fig_attitude_fusion.png');

% === FIG 7: ATTITUDE ERROR COMPARISON — zoomed to show fusion benefit ===
% Uses same S3 results as Fig 6 above
figure('Position',[100 100 1000 500]);
t_eval = burn+1:par.N;  % evaluation period (after burn-in)

subplot(1,3,1); hold on;
eh_g = rad2deg(abs(wrapToPi(traj.att(1,t_eval)-est_gnss4.att(1,t_eval))));
eh_u = rad2deg(abs(wrapToPi(traj.att(1,t_eval)-est_uwb.att(1,t_eval))));
eh_f = rad2deg(abs(wrapToPi(traj.att(1,t_eval)-est_fused.att(1,t_eval))));
plot(par.t(t_eval), eh_g, 'b-', 'LineWidth',1); 
plot(par.t(t_eval), eh_u, '-', 'Color',[0.6 0 0.6], 'LineWidth',1);
plot(par.t(t_eval), eh_f, 'r-', 'LineWidth',1.5);
ylabel('|Heading Error| (°)'); xlabel('Time (s)');
legend('GNSS 4-ant','UWB multi-tag','Fused','FontSize',7,'Location','northwest');
title(sprintf('(a) Heading\nGNSS=%.1f° UWB=%.2f° Fused=%.2f°',...
    sqrt(mean(eh_g.^2)), sqrt(mean(eh_u.^2)), sqrt(mean(eh_f.^2))));
grid on; set(gca,'FontSize',10); ylim([0 min(120,max(ylim))]);

subplot(1,3,2); hold on;
er_g = rad2deg(abs(wrapToPi(traj.att(2,t_eval)-est_gnss4.att(2,t_eval))));
er_u = rad2deg(abs(wrapToPi(traj.att(2,t_eval)-est_uwb.att(2,t_eval))));
er_f = rad2deg(abs(wrapToPi(traj.att(2,t_eval)-est_fused.att(2,t_eval))));
plot(par.t(t_eval), er_g, 'b-', 'LineWidth',1);
plot(par.t(t_eval), er_u, '-', 'Color',[0.6 0 0.6], 'LineWidth',1);
plot(par.t(t_eval), er_f, 'r-', 'LineWidth',1.5);
ylabel('|Roll Error| (°)'); xlabel('Time (s)');
title(sprintf('(b) Roll\nGNSS=%.1f° UWB=%.2f° Fused=%.2f°',...
    sqrt(mean(er_g.^2)), sqrt(mean(er_u.^2)), sqrt(mean(er_f.^2))));
grid on; set(gca,'FontSize',10); ylim([0 min(40,max(ylim))]);

subplot(1,3,3); hold on;
ep_g = rad2deg(abs(wrapToPi(traj.att(3,t_eval)-est_gnss4.att(3,t_eval))));
ep_u = rad2deg(abs(wrapToPi(traj.att(3,t_eval)-est_uwb.att(3,t_eval))));
ep_f = rad2deg(abs(wrapToPi(traj.att(3,t_eval)-est_fused.att(3,t_eval))));
plot(par.t(t_eval), ep_g, 'b-', 'LineWidth',1);
plot(par.t(t_eval), ep_u, '-', 'Color',[0.6 0 0.6], 'LineWidth',1);
plot(par.t(t_eval), ep_f, 'r-', 'LineWidth',1.5);
ylabel('|Pitch Error| (°)'); xlabel('Time (s)');
title(sprintf('(c) Pitch\nGNSS=%.1f° UWB=%.2f° Fused=%.2f°',...
    sqrt(mean(ep_g.^2)), sqrt(mean(ep_u.^2)), sqrt(mean(ep_f.^2))));
grid on; set(gca,'FontSize',10); ylim([0 min(30,max(ylim))]);
sgtitle('Attitude Error Comparison — S3:GNSS Outage');
saveas(gcf,'fig_attitude_error_comparison.png');

% === FIG 8: TRAJECTORY + SYSTEM LAYOUT ===
figure('Position',[100 100 800 600]);
hold on;
% Ship trajectory
plot(traj.pos(1,:), traj.pos(2,:), 'b-', 'LineWidth',2);
plot(traj.pos(1,1), traj.pos(2,1), 'go', 'MarkerSize',12, 'MarkerFaceColor','g');
plot(traj.pos(1,end), traj.pos(2,end), 'rs', 'MarkerSize',12, 'MarkerFaceColor','r');
% Dock anchors
for j=1:par.M
    plot(par.anch_body(j,1), par.anch_body(j,2), 'k^', 'MarkerSize',10, 'MarkerFaceColor','k');
    text(par.anch_body(j,1)+3, par.anch_body(j,2)+3, sprintf('A%d',j-1), 'FontSize',9);
end
% Dock outline
dock_x = [-80 80 80 -80 -80]; dock_y = [-25 -25 25 25 -25];
plot(dock_x, dock_y, 'k-', 'LineWidth',2);
% Shadow/outage zones
theta_circ = linspace(0,2*pi,100);
plot(80*cos(theta_circ), 80*sin(theta_circ), 'r:', 'LineWidth',1.5);
text(60, -60, 'Shadow (80m)', 'Color','r', 'FontSize',9);
% Distance markers
for d_mark = [100 200 300 400 500]
    plot(d_mark*cos(theta_circ), d_mark*sin(theta_circ), ':', 'Color',[0.7 0.7 0.7]);
    text(d_mark*0.7, d_mark*0.7, sprintf('%dm',d_mark), 'Color',[0.5 0.5 0.5], 'FontSize',8);
end
% Phase labels
text(350, 60, 'Phase 1: Approach', 'FontSize',10, 'Color','b');
text(100, -40, 'Phase 2: Shadow', 'FontSize',10, 'Color',[0.8 0.5 0]);
text(20, 30, 'Phase 3: Docking', 'FontSize',10, 'Color','r');
xlabel('X (m)'); ylabel('Y (m)');
legend('Ship trajectory','Start','End','Dock anchors','Location','northeast');
title('System Layout: Ship Approach to Floating Dry Dock');
axis equal; grid on; set(gca,'FontSize',11);
saveas(gcf,'fig_system_layout.png');

%% ═══ PART G: ATTITUDE OBSERVABILITY ═══
fprintf('\n--- Part G: Attitude observability (Theorems 1-3) ---\n');
ntags_test = [1 2 3 4];
tag_cfgs = {[40 0 8]; [40 0 8; -40 0 4]; [40 0 8; -40 0 4; 0 -15 12];
            [40 0 8; -40 0 4; 0 -15 12; 0 15 6]};
p_test = [100;0;0]; att_test = [pi;0;0];
for ti=1:4
    [cp,ch,rk] = compute_crlb_6dof(p_test, att_test, tag_cfgs{ti}, par);
    % Classification: rank alone is insufficient — check CRLB quality
    if rk < 6
        obs_str = 'UNOBSERVABLE';
    elseif rad2deg(ch) > 100 || cp > 1000
        obs_str = 'ILL-CONDITIONED';  % numerically full-rank but practically singular
    else
        obs_str = 'OBSERVABLE';
    end
    fprintf('  %d tags: rank=%d (%s) CRLB_pos=%.2f m CRLB_hdg=%.2f°\n',...
        ntags_test(ti), rk, obs_str, cp, rad2deg(ch));
end

fprintf('  Verifying across 10000 random configs...\n');
rank_counts = zeros(4,1);
for trial = 1:10000
    for ti = 1:4
        tags_r = (rand(ntags_test(ti),3)-0.5).*[80 30 20];
        [cp_r,ch_r,rk] = compute_crlb_6dof(p_test, att_test, tags_r, par);
        % Count as observable only if rank=6 AND CRLB is practically useful
        if rk >= 6 && rad2deg(ch_r) < 100 && cp_r < 1000
            rank_counts(ti) = rank_counts(ti)+1;
        end
    end
end
for ti=1:4
    fprintf('    %d tags: rank=6 in %.1f%% of configs\n', ntags_test(ti), rank_counts(ti)/100);
end

figure('Position',[100 100 700 350]);
bar(ntags_test, rank_counts/100, 0.6, 'FaceColor',[.2 .6 .9]);
ylabel('Configs achieving rank 6 (%)'); xlabel('Number of tags');
title('Attitude Observability vs Number of Tags');
set(gca,'FontSize',12); grid on;
saveas(gcf,'fig_observability.png');

%% ═══ PART H: GDAP VALIDATION ═══
fprintf('\n--- Part H: GDAP vs distance + tag baseline ---\n');
dists_gdap = linspace(500, 10, 50);
baselines = [20 40 80];
gdap_all = zeros(length(baselines), length(dists_gdap));

for bi = 1:length(baselines)
    L = baselines(bi);
    tags_scaled = [L/2 0 8; -L/2 0 4; 0 -L/4 12];
    for di = 1:length(dists_gdap)
        [~,ch,~] = compute_crlb_6dof([dists_gdap(di);0;0], [pi;0;0], tags_scaled, par);
        gdap_all(bi,di) = rad2deg(ch);
    end
end

ratio_test = dists_gdap / 40;
gdap_40 = gdap_all(2,:);
scale_factor = median(gdap_40 ./ ratio_test);
fprintf('  GDAP decomposition: GDAP ≈ %.3f × (d/L_tag)\n', scale_factor);
fprintf('  Correlation: r = %.3f\n', corr(ratio_test', gdap_40'));

figure('Position',[100 100 800 400]);
subplot(1,2,1);
for bi=1:3; plot(dists_gdap, gdap_all(bi,:), '-', 'LineWidth',2); hold on; end
yline(1,'k:','1° target','LineWidth',1.5);
set(gca,'XDir','reverse','FontSize',11);
xlabel('Distance (m)'); ylabel('GDAP — heading bound (°)');
legend(arrayfun(@(b) sprintf('L_{tag}=%dm',b), baselines, 'Uni',0),'Location','northwest');
title('GDAP vs Distance'); grid on;
subplot(1,2,2);
plot(ratio_test, gdap_40, 'bo', 'MarkerSize',4); hold on;
plot(ratio_test, scale_factor*ratio_test, 'r-', 'LineWidth',2);
xlabel('d / L_{tag}'); ylabel('GDAP (°)');
legend('Numerical','Linear fit','Location','northwest');
title('GDAP Linearity'); grid on; set(gca,'FontSize',11);
saveas(gcf,'fig_gdap.png');

%% ═══ PART I: GEO ATTITUDE ANCHORING ═══
fprintf('\n--- Part I: GEO attitude anchoring ---\n');
% Progressive satellite blockage: remove low-elevation GPS sats
% Compare GPS-only vs GPS+NavIC heading when MEO sats are blocked
% This demonstrates NavIC's high-elevation benefit

n_block = 0:2:min(12, size(sat_gps,1));
hdg_rmse_gps = zeros(size(n_block));
hdg_rmse_navic = zeros(size(n_block));
NMC_I = 20;

for bi = 1:length(n_block)
    nb = n_block(bi);
    rmse_g = zeros(1,NMC_I); rmse_n = zeros(1,NMC_I);
    for mc = 1:NMC_I
        rng(mc*7000+bi);
        
        % Create reduced GPS constellation (remove lowest elevation sats)
        sg = sat_gps;
        if nb > 0 && nb <= size(sg,1)
            [~,idx] = sort(sg(:,2));  % sort by elevation ascending
            sg(idx(1:min(nb,size(sg,1))), :) = [];  % remove lowest
        end
        
        % GPS-only: run EKF with reduced GPS constellation
        mm_g = gen_meas(traj, plat, par, 1, sg, sat_navic, 'gps_only');
        [est_g, ~] = run_ekf(traj, plat, mm_g, par, 'gps_4ant', 'perfect');
        eh_g = wrapToPi(traj.att(1,burn+1:end) - est_g.att(1,burn+1:end));
        rmse_g(mc) = rad2deg(sqrt(mean(eh_g.^2)));
        
        % GPS+NavIC: same reduced GPS + all NavIC satellites
        mm_n = gen_meas(traj, plat, par, 1, sg, sat_navic, 'gps_navic');
        [est_n, ~] = run_ekf(traj, plat, mm_n, par, 'gpsnav_4ant', 'perfect');
        eh_n = wrapToPi(traj.att(1,burn+1:end) - est_n.att(1,burn+1:end));
        rmse_n(mc) = rad2deg(sqrt(mean(eh_n.^2)));
    end
    hdg_rmse_gps(bi) = mean(rmse_g);
    hdg_rmse_navic(bi) = mean(rmse_n);
    fprintf('  Block %2d MEO: GPS=%d sats hdg=%.2f° | +NavIC=%d sats hdg=%.2f° (%.0f%% better)\n',...
        nb, size(sg,1), hdg_rmse_gps(bi), size(sg,1)+size(sat_navic,1), ...
        hdg_rmse_navic(bi), (hdg_rmse_gps(bi)-hdg_rmse_navic(bi))/max(hdg_rmse_gps(bi),0.01)*100);
end

fprintf('\n  NavIC high-elevation satellites (from %s):\n', tern(use_real,'REAL DATA','model'));
for gi = 1:min(3,size(sat_navic,1))
    fprintf('    NavIC-%d: az=%.1f° el=%.1f°\n', gi, sat_navic(gi,1), sat_navic(gi,2));
end

figure('Position',[100 100 700 400]);
plot(n_block, hdg_rmse_gps, 'b-o', 'LineWidth',2,'MarkerFaceColor','b'); hold on;
plot(n_block, hdg_rmse_navic, 'r-s', 'LineWidth',2,'MarkerFaceColor','r');
xlabel('MEO satellites blocked'); ylabel('Heading RMSE (°)');
legend('GPS only','GPS+NavIC (high-el anchored)','Location','northwest');
title('GEO/High-Elevation Attitude Anchoring (EKF-based)');
grid on; set(gca,'FontSize',12);
saveas(gcf,'fig_geo_anchoring.png');

%% ═══ PART J: COMPLEMENTARY FUSION BOUND ═══
fprintf('\n--- Part J: Complementary fusion bound ---\n');
dists_comp = linspace(500, 10, 50);
crlb_gnss = zeros(size(dists_comp));
crlb_uwb = zeros(size(dists_comp));
crlb_fused = zeros(size(dists_comp));

for di = 1:length(dists_comp)
    d = dists_comp(di);
    mask = compute_dock_mask(d, par);
    [ns,~] = count_visible([sat_gps;sat_navic], mask, par);
    
    % === REALISTIC near-dock GNSS carrier-phase degradation ===
    % Multipath: steel dock causes severe carrier multipath at close range
    %   - Open sky (>300m):  σ_cp ≈ 0.003 m (clean carrier phase)
    %   - Moderate (100-300m): σ_cp ≈ 0.01-0.03 m (emerging multipath)
    %   - Shadow (<100m):    σ_cp ≈ 0.05-0.15 m (severe multipath from dock)
    %   - Very close (<30m): σ_cp ≈ 0.10-0.20 m (near-field diffraction)
    mp_carrier = 1 + 15*exp(-d/60) + 30*exp(-d/20);  % aggressive multipath model
    sig_cp = par.sig_carrier * mp_carrier;
    
    % Ambiguity resolution: fails in high-multipath near dock
    %   - Success requires ADOP < 0.12 cycles
    %   - ADOP degrades with multipath, fewer sats, poor geometry
    %   - Below ~80m: ambiguity resolution becomes unreliable
    if ns >= 5
        adop_factor = sig_cp / par.sig_carrier;  % how much worse than clean
        p_amb = max(0, 1 - exp(-0.8*(ns-5))) * exp(-0.5*(adop_factor-1));
    else
        p_amb = 0;
    end
    
    % GNSS heading CRLB: σ_ψ = σ_cp / (B × √(N-1))
    B_gnss = 80;  % baseline (m)
    if ns >= 5 && p_amb > 0.01
        sig_hdg_gnss = sig_cp / (B_gnss * sqrt(max(ns-1,1)));  % radians
        % Weight by ambiguity success: effective CRLB = σ/p_amb
        % (when ambiguity fails, heading is random → infinite CRLB)
        effective_sig = sig_hdg_gnss / sqrt(p_amb);
        J_gnss = 1 / effective_sig^2;
    else
        J_gnss = 0;
    end
    
    % UWB attitude CRLB (from multi-tag TDOA geometry)
    [~,ch_rad,~] = compute_crlb_6dof([d;0;0],[pi;0;0],par.tags_body,par);
    ch_deg = rad2deg(ch_rad);
    if ch_deg < 100
        J_uwb = 1 / ch_rad^2;
    else
        J_uwb = 0;
    end
    
    % Store individual CRLBs
    if J_gnss > 1e-10
        crlb_gnss(di) = rad2deg(1/sqrt(J_gnss));
    else
        crlb_gnss(di) = 180;
    end
    crlb_uwb(di) = ch_deg;
    
    % Fused CRLB: J_fused = J_gnss + J_uwb (FIM additivity)
    J_total = J_gnss + J_uwb;
    if J_total > 1e-10
        crlb_fused(di) = rad2deg(1/sqrt(J_total));
    else
        crlb_fused(di) = 180;
    end
end

% Find crossover distance (where GNSS and UWB CRLBs are equal)
% Avoid comparing where one is at 180° (invalid)
valid_both = crlb_gnss < 170 & crlb_uwb < 170;
if any(valid_both)
    diffs = abs(crlb_gnss - crlb_uwb);
    diffs(~valid_both) = Inf;
    [~, cross_idx] = min(diffs);
else
    cross_idx = round(length(dists_comp)/2);
end
fprintf('  Crossover distance: %.0f m\n', dists_comp(cross_idx));
fprintf('  At crossover: GNSS=%.2f° UWB=%.2f° Fused=%.2f°\n',...
    crlb_gnss(cross_idx), crlb_uwb(cross_idx), crlb_fused(cross_idx));
improvement = (1 - crlb_fused(cross_idx)/min(crlb_gnss(cross_idx),crlb_uwb(cross_idx)))*100;
fprintf('  Improvement at crossover: %.0f%% vs best single\n', improvement);

% Also print at key distances
for dd = [500 200 100 50 30]
    [~,idx] = min(abs(dists_comp - dd));
    fprintf('    d=%3dm: GNSS=%.2f° UWB=%.2f° Fused=%.2f° (p_amb≈%.0f%%)\n',...
        dd, crlb_gnss(idx), crlb_uwb(idx), crlb_fused(idx), ...
        max(0,1-exp(-0.8*(count_visible([sat_gps;sat_navic],compute_dock_mask(dd,par),par)-5)))*...
        exp(-0.5*(par.sig_carrier*(1+15*exp(-dd/60)+30*exp(-dd/20))/par.sig_carrier-1))*100);
end

figure('Position',[100 100 700 450]);
semilogy(dists_comp, crlb_gnss, 'b-', 'LineWidth',2); hold on;
semilogy(dists_comp, crlb_uwb, 'g-', 'LineWidth',2);
semilogy(dists_comp, crlb_fused, 'r-', 'LineWidth',2.5);
xline(dists_comp(cross_idx),'k:',sprintf('Crossover %.0fm',dists_comp(cross_idx)),'LineWidth',1.5);
yline(1,'k--','1°','LineWidth',1);
yline(0.5,'k:','0.5° target','LineWidth',1);
set(gca,'XDir','reverse','FontSize',12);
xlabel('Distance to dock (m)'); ylabel('Heading CRLB (°)');
legend('GNSS attitude only','UWB attitude only','Fused (proposed)','Location','north');
title('Complementary Attitude Fusion Bound (Proposition 5)'); grid on;
ylim([0.01 200]);
saveas(gcf,'fig_complementary_bound.png');

%% ═══ PART K: PLATFORM COMPENSATION ═══
fprintf('\n--- Part K: Platform compensation ---\n');
% Sea states defined by Hs and Tp (JONSWAP spectrum)
ss = struct('Hs',{0.5, 1.5, 3.0}, 'Tp',{7, 10, 13}, ...
            'nm',{'Calm (SS3)', 'Moderate (SS4)', 'Rough (SS5-6)'});
comp_m={'perfect','gnss_plat','seastate','none'};
comp_n={'Dock GNSS','Dock pos only','Sea-state','None'};
comp_res=zeros(4,3);
for si=1:3
    ps=par; ps.Hs=ss(si).Hs; ps.Tp=ss(si).Tp;
    waves_k = gen_waves(ps);  % new wave field for this sea state
    pl=gen_plat(ps, waves_k);
    fprintf('  %s (Hs=%.1fm Tp=%.0fs): dock heave=±%.2fm roll=±%.1f°\n', ...
        ss(si).nm, ss(si).Hs, ss(si).Tp, ...
        max(abs(pl.hv)), rad2deg(max(abs(pl.rl))));
    for ci=1:4
        rt=zeros(1,20);
        for mc=1:20
            rng(mc*4000+si*100+ci);
            mm=gen_meas(traj,pl,ps,3,sat_gps,sat_navic,'gps_navic');
            [e,~]=run_ekf(traj,pl,mm,ps,'mtag_gpsnav',comp_m{ci});
            ep=traj.pos(:,burn+1:end)-e.pos(:,burn+1:end);
            rt(mc)=sqrt(mean(sum(ep.^2,1)));
        end
        comp_res(ci,si)=mean(rt);
    end
    fprintf('    dock_gnss=%.2f none=%.2f m\n',comp_res(1,si),comp_res(4,si));
end

figure('Position',[100 100 600 400]);
imagesc(comp_res); colormap(flipud(hot)); cb=colorbar; cb.Label.String='RMSE (m)';
set(gca,'XTick',1:3,'XTickLabel',{sprintf('Calm\nHs=%.1fm',ss(1).Hs), ...
    sprintf('Moderate\nHs=%.1fm',ss(2).Hs), sprintf('Rough\nHs=%.1fm',ss(3).Hs)});
set(gca,'YTick',1:4,'YTickLabel',comp_n);
for ci=1:4; for si=1:3
    text(si,ci,sprintf('%.2f',comp_res(ci,si)),'HorizontalAlignment','center',...
        'FontWeight','bold','FontSize',11);
end; end
title('Platform Compensation Performance'); set(gca,'FontSize',12);
saveas(gcf,'fig_compensation.png');

%% ═══ PART L: RAIM ═══
fprintf('\n--- Part L: Maritime RAIM ---\n');

% Step 1: Calibrate NIS distribution from fault-free runs
nis_cal = [];
nis_dof = [];  % track measurement dimension per epoch
for mc = 1:15
    rng(mc*9999);
    mm = gen_meas(traj, plat, par, 1, sat_gps, sat_navic, 'gps_navic');
    [~, dg] = run_ekf(traj, plat, mm, par, 'mtag_gpsnav', 'perfect');
    valid_nis = dg.nis(dg.nis > 0);
    nis_cal = [nis_cal, valid_nis];
end

% Step 2: Determine threshold
% Method: Use empirical distribution with conservative percentile
% Then add safety margin to control Pfa < 0.05
nis_median = median(nis_cal);
nis_std = std(nis_cal);
nis_p99 = prctile(nis_cal, 99);
nis_p999 = prctile(nis_cal, 99.9);

% Threshold = max of 99.9th percentile and (median + 5σ)
tau_r = max(nis_p999, nis_median + 5*nis_std);
% Additional safety: at least 3× median (prevents Pfa >> 0.05)
tau_r = max(tau_r, 3 * nis_median);
fprintf('  NIS calibration: median=%.1f std=%.1f P99=%.1f P99.9=%.1f\n', ...
    nis_median, nis_std, nis_p99, nis_p999);
fprintf('  Selected threshold: tau=%.1f (%.1f× median)\n', tau_r, tau_r/nis_median);

% Step 3: Fault detection using SLIDING WINDOW
% Detect = window of W consecutive epochs where NIS > tau
% This reduces false alarms from single-epoch spikes
W_detect = 3;  % detection window: 3 consecutive high NIS epochs (150 ms)

fmag = [0 2 5 8 12 20]; NMC_L = 50;
Pd = zeros(size(fmag)); TTA = Pd; Pfa_per_epoch = 0;

for fi = 1:length(fmag)
    fm = fmag(fi); dc = 0; ts = 0; n_fa_epochs = 0; n_total_epochs = 0;
    for mc = 1:NMC_L
        rng(mc*5000+fi);
        mm = gen_meas(traj, plat, par, 1, sat_gps, sat_navic, 'gps_navic');
        kf = round(60/par.dt);
        % Inject fault from kf onward
        for k = kf:par.N
            if mm.uwb(k).valid
                for ti = 1:length(mm.uwb(k).z_tags)
                    if ~isempty(mm.uwb(k).z_tags{ti})
                        mm.uwb(k).z_tags{ti}(1) = mm.uwb(k).z_tags{ti}(1) + fm;
                    end
                end
            end
        end
        [~, dg] = run_ekf(traj, plat, mm, par, 'mtag_gpsnav', 'perfect');
        
        % Sliding window detection
        nis_seq = dg.nis(kf:end);
        detected = false;
        for wi = W_detect:length(nis_seq)
            if all(nis_seq(wi-W_detect+1:wi) > tau_r)
                if ~detected
                    detected = true;
                    ts = ts + (wi - W_detect + 1) * par.dt;
                end
                break;
            end
        end
        if detected; dc = dc + 1; end
        
        % Count false alarms (only for fm=0)
        if fm == 0
            n_total_epochs = n_total_epochs + length(nis_seq);
            % Count windows where all W exceed threshold
            for wi = W_detect:length(nis_seq)
                if all(nis_seq(wi-W_detect+1:wi) > tau_r)
                    n_fa_epochs = n_fa_epochs + 1;
                end
            end
        end
    end
    Pd(fi) = dc / NMC_L;
    TTA(fi) = ts / max(dc, 1);
    if fm == 0
        Pfa_per_epoch = n_fa_epochs / max(n_total_epochs, 1);
        fprintf('  Pfa=%.4f (per window, W=%d epochs)\n', Pd(fi), W_detect);
        fprintf('  Pfa per epoch=%.6f\n', Pfa_per_epoch);
    else
        fprintf('  Fault %2dm: Pd=%.3f TTA=%.1fs\n', fm, Pd(fi), TTA(fi));
    end
end

figure('Position',[100 100 800 350]);
subplot(1,2,1); plot(fmag, Pd*100, 'bo-', 'LineWidth',2, 'MarkerFaceColor','b');
yline(90,'r--','90%'); ylim([0 105]); xlabel('Fault magnitude (m)'); ylabel('Detection rate (%)');
title(sprintf('Fault Detection (W=%d, \\tau=%.0f)', W_detect, tau_r));
grid on; set(gca,'FontSize',11);
subplot(1,2,2); vi = fmag>0 & Pd>0;
if any(vi); plot(fmag(vi), TTA(vi), 'rs-', 'LineWidth',2, 'MarkerFaceColor','r'); end
xlabel('Fault magnitude (m)'); ylabel('Time to Alert (s)');
title('Time to Alert'); grid on; set(gca,'FontSize',11);
saveas(gcf,'fig_raim.png');

%% ═══ PART M: UTIL VALIDATION ═══
fprintf('\n--- Part M: UTIL validation ---\n');
pu=par; pu.T=123; pu.N=round(pu.T/pu.dt); pu.t=(0:pu.N-1)*pu.dt;
pu.sig_tdoa=0.333; pu.pdrop=0.005; pu.pnlos=0.006; pu.sig_nlos=0.5;
pu.anch_body=[-2.5 -2.5 0; 2.5 -2.5 0; -2.5 2.5 0; 2.5 2.5 0;
              -2.5 -2.5 2.5; 2.5 -2.5 2.5; -2.5 2.5 2.5; 2.5 2.5 2.5];
pu.M=8; pu.ref=1; pu.plat_pos=[0;0;0];
pu.tags_body=[1.5 0 0.3; -1.5 0 0.3; 0 -1 0.3]; pu.Ntags=3;
pu.gnss_ant_body=par.gnss_ant_body; pu.gnss_body=[0;0;.5];
pu.sig_carrier=par.sig_carrier; pu.lambda_L1=par.lambda_L1;
% Wave parameters for UTIL scale (indoor lab — minimal motion)
pu.Hs=0.01; pu.Tp=7; pu.gamma_js=3.3; pu.n_waves=20;
pu.dock_Tn_heave=5; pu.dock_Tn_roll=5; pu.dock_Tn_pitch=5;
pu.dock_zeta_heave=0.1; pu.dock_zeta_roll=0.1; pu.dock_zeta_pitch=0.1;
pu.dock_drift=0.01; pu.dock_L=5; pu.dock_B=5;
pu.ship_Tn_heave=3; pu.ship_Tn_roll=3; pu.ship_Tn_pitch=3;
pu.ship_zeta_heave=0.1; pu.ship_zeta_roll=0.1; pu.ship_zeta_pitch=0.1;
pu.ship_L=3; pu.ship_B=1; pu.ship_T=0.5; pu.ship_disp=10;

tu.pos=zeros(3,pu.N); tu.vel=tu.pos; tu.att=zeros(3,pu.N);
for k=1:pu.N; tk=pu.t(k);
    tu.pos(:,k)=[2*sin(2*pi/60*tk); 2*sin(4*pi/60*tk); 1+.3*sin(2*pi/40*tk)];
    tu.att(1,k)=atan2(cos(4*pi/60*tk)*4*pi/60, cos(2*pi/60*tk)*2*pi/60);
    if k>1; tu.vel(:,k)=(tu.pos(:,k)-tu.pos(:,k-1))/pu.dt; end
end
bu=round(2/pu.dt);
mcfg=struct('Hs',{0, 0.5, 1.5},'Tp',{7,8,10},'nm',{'Static','Moderate','Rough'});
L_util = norm(pu.tags_body(1,:)-pu.tags_body(2,:));
d_util = 2;
gdap_pred = rad2deg(pu.sig_tdoa * d_util / (L_util * sqrt(pu.M-1)));
fprintf('  UTIL GDAP prediction: %.1f° (L=%.1fm, d=%.0fm)\n', gdap_pred, L_util, d_util);

for ci=1:3
    pc=pu; pc.Hs=mcfg(ci).Hs; pc.Tp=mcfg(ci).Tp;
    if pc.Hs == 0; pc.Hs = 0.01; end  % avoid zero for spectrum
    waves_m = gen_waves(pc);
    plc=gen_plat(pc, waves_m); rng(42);
    mm=gen_meas(tu,plc,pc,1,sat_gps,sat_navic,'gps_navic');
    [ec,~]=run_ekf(tu,plc,mm,pc,'mtag_only','perfect');
    ep=tu.pos(:,bu+1:end)-ec.pos(:,bu+1:end);
    eh=wrapToPi(tu.att(1,bu+1:end)-ec.att(1,bu+1:end));
    fprintf('  %-12s pos=%.1fcm hdg=%.2f° (pred=%.1f°)\n', mcfg(ci).nm,...
        sqrt(mean(sum(ep.^2,1)))*100, rad2deg(sqrt(mean(eh.^2))), gdap_pred);
end

%% ═══ SUMMARY ═══
fprintf('\n╔════════════════════════════════════════════════════════════════════╗\n');
fprintf('║  FULL PAPER RUN COMPLETE                                         ║\n');
fprintf('╠════════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Total runtime: %.1f minutes (%.0f seconds)                     ║\n', toc(tic_total)/60, toc(tic_total));
if use_real
fprintf('║  ★ REAL DATA PARAMETERS USED:                                    ║\n');
fprintf('║    NavIC: Kaggle NavIC_Dataset (86,400 epochs, σ_PR=%.3fm)     ║\n', par.navic_sig_pr);
fprintf('║    GPS:   4-antenna RINEX (baselines %.1f-%.1fm, σ_CP=%.4fm) ║\n', ...
    min(real_par.gps_baselines), max(real_par.gps_baselines), par.sig_carrier);
fprintf('║    UWB:   UTIL dataset parameters (σ=0.333m)                    ║\n');
end
fprintf('║  Part B: GPS+NavIC constellation (%d+%d sats)                    ║\n', size(sat_gps,1), size(sat_navic,1));
fprintf('║  Part F: 8 methods × 3 scenarios × %d MC                        ║\n', NMC);
fprintf('║  Part G: Observability (Theorems 1-3, 10000 random configs)      ║\n');
fprintf('║  Part H: GDAP (3 baselines × 50 distances)                      ║\n');
fprintf('║  Part I: GEO anchoring (7 blockage levels × %d MC)              ║\n', NMC_I);
fprintf('║  Part J: Complementary bound (50 distances)                      ║\n');
fprintf('║  Part K: Compensation (4 × 3 × 20 MC)                           ║\n');
fprintf('║  Part L: RAIM (%d faults × %d MC)                               ║\n', length(fmag), NMC_L);
fprintf('║  Part M: UTIL validation (3 motion configs)                      ║\n');
fprintf('║                                                                  ║\n');
fprintf('║  Figures saved:                                                   ║\n');
fprintf('║    fig_system_layout        — Dock + trajectory + anchors        ║\n');
fprintf('║    fig_satellite_visibility — Sat count + PDOP vs distance       ║\n');
fprintf('║    fig_gnss_attitude        — GNSS heading σ vs distance         ║\n');
fprintf('║    fig_method_comparison    — Bar chart: 7 methods × 3 scenarios ║\n');
fprintf('║    fig_cdf_errors           — CDF of position + heading errors   ║\n');
fprintf('║    fig_6dof_timeseries      — Pos+Hdg+Roll+Pitch error vs time  ║\n');
fprintf('║    fig_attitude_fusion      — True vs GNSS vs UWB vs Fused att  ║\n');
fprintf('║    fig_attitude_error_comp  — Hdg/Roll/Pitch error comparison   ║\n');
fprintf('║    fig_observability        — Rank-6 achievability vs N_tags     ║\n');
fprintf('║    fig_gdap                 — GDAP vs distance + linearity      ║\n');
fprintf('║    fig_geo_anchoring        — GPS vs NavIC under blockage (EKF) ║\n');
fprintf('║    fig_complementary_bound  — GNSS vs UWB vs Fused CRLB        ║\n');
fprintf('║    fig_compensation         — Platform comp heatmap (3 SS)      ║\n');
fprintf('║    fig_raim                 — Fault detection + time-to-alert   ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════╝\n');

% ═══ SAVE ALL RESULTS ═══
fprintf('\nSaving workspace...\n');
save('full_run_results.mat', 'res_pos', 'res_hdg', 'res_roll', 'res_pitch', ...
    'mnames', 'snames', ...
    'rank_counts', 'ntags_test', 'gdap_all', 'dists_gdap', 'baselines', ...
    'scale_factor', 'hdg_rmse_gps', 'hdg_rmse_navic', 'n_block', ...
    'crlb_gnss', 'crlb_uwb', 'crlb_fused', 'dists_comp', ...
    'comp_res', 'comp_n', 'ss', 'Pd', 'TTA', 'fmag', 'tau_r', 'W_detect', ...
    'nvis_gps', 'nvis_navic', 'nvis_total', 'pdop_gps', 'pdop_total', 'dists_vis', ...
    'sig_hdg_gps', 'sig_hdg_navic', 'dists_att', ...
    'par', 'sat_gps', 'sat_navic');
fprintf('  Saved: full_run_results.mat\n');
fprintf('  (Contains all tables, figures data, and parameters for paper)\n');
fprintf('\nDone! Runtime: %.1f minutes\n', toc(tic_total)/60);

%% ═══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
%  ═══════════════════════════════════════════════════════════════════════════

function s = tern(cond, a, b)
    if cond; s = a; else; s = b; end
end

function traj = gen_traj(par, waves)
    N=par.N; dt=par.dt; t=par.t;
    pos=zeros(3,N); vel=zeros(3,N); att=zeros(3,N);
    pos(:,1)=[500;50;0]; att(1,1)=pi;
    
    % Wave-induced ship motion using JONSWAP + ship RAOs
    % Physical scaling: ship responds MORE than dock (lighter, shorter)
    heave_rms = (par.Hs/4) * 0.9;                     % supply vessel heave ~90% of Hs/4
    roll_rms  = deg2rad((par.Hs / par.ship_B) * 30);   % ~30×(Hs/B) degrees — more than dock
    pitch_rms = deg2rad((par.Hs / par.ship_L) * 30);   % ~30×(Hs/L) degrees — more than dock
    surge_rms = par.Hs * 0.15;                          % surge perturbation
    sway_rms  = par.Hs * 0.10;                          % sway perturbation
    
    ship_hv = compute_rao_response(waves, par.ship_Tn_heave, par.ship_zeta_heave, heave_rms, t);
    ship_rl = compute_rao_response(waves, par.ship_Tn_roll, par.ship_zeta_roll, roll_rms, t);
    ship_pt = compute_rao_response(waves, par.ship_Tn_pitch, par.ship_zeta_pitch, pitch_rms, t);
    ship_sg = compute_rao_response(waves, par.ship_Tn_heave*1.3, 0.15, surge_rms, t);
    ship_sw = compute_rao_response(waves, par.ship_Tn_roll*0.8, 0.10, sway_rms, t);
    
    for k=2:N; tk=t(k);
        % Guidance: approach path toward dock
        if tk<60; th=atan2(-pos(2,k-1),-pos(1,k-1)); sp=4;
        elseif tk<120; th=atan2(-pos(2,k-1),-pos(1,k-1)); sp=max(2,4-(tk-60)/60*2);
        elseif tk<180; th=atan2(-pos(2,k-1),-pos(1,k-1)); sp=max(1,2-(tk-120)/60);
        elseif tk<240; fr=(tk-180)/60; th=pi+fr*(-pi/2); sp=max(.5,1-fr*.5);
        else; th=att(1,k-1); sp=.2; end
        dh=wrapToPi(th-att(1,k-1)); om=max(min(.3*dh,.05),-.05);
        
        % Heading = guidance + wave-induced yaw
        att(1,k)=wrapToPi(att(1,k-1)+om*dt);
        % Roll and pitch from wave spectrum (NOT simple sine anymore)
        att(2,k) = ship_rl(k);   % wave-induced roll
        att(3,k) = ship_pt(k);   % wave-induced pitch
        
        cs=norm(vel(1:2,k-1)); if cs<.01; cs=sp; end
        ns=cs+.1*(sp-cs)*dt;
        vel(1,k)=ns*cos(att(1,k)) + (ship_sg(min(k,N))-ship_sg(max(k-1,1)))/dt*0.1;
        vel(2,k)=ns*sin(att(1,k)) + (ship_sw(min(k,N))-ship_sw(max(k-1,1)))/dt*0.1;
        vel(3,k) = (ship_hv(min(k,N))-ship_hv(max(k-1,1)))/dt; % heave velocity
        pos(1,k) = pos(1,k-1) + vel(1,k)*dt;
        pos(2,k) = pos(2,k-1) + vel(2,k)*dt;
        pos(3,k) = ship_hv(k);  % heave displacement directly
    end
    traj.pos=pos; traj.vel=vel; traj.att=att;
    traj.ship_hv=ship_hv; traj.ship_rl=ship_rl; traj.ship_pt=ship_pt;
end

function plat = gen_plat(par, waves)
    t=par.t; N=par.N;
    
    % Dock 6-DOF response using RAOs
    % amp_scale = desired RMS amplitude (peak ≈ 2.5 × RMS for broadband)
    %
    % Physical scaling with Hs:
    %   Heave RMS ≈ (Hs/4) × heave_transfer  (heave_transfer ~0.5-0.8 for FPSO)
    %   Roll  RMS ≈ (Hs/B) × roll_factor      (rad, empirical)
    %   Pitch RMS ≈ (Hs/L) × pitch_factor     (rad, empirical)
    
    heave_rms = (par.Hs/4) * 0.6;                    % ~0.6 transfer for moored FPSO
    roll_rms  = deg2rad((par.Hs / par.dock_B) * 40);  % empirical: ~40×(Hs/B) degrees
    pitch_rms = deg2rad((par.Hs / par.dock_L) * 30);  % empirical: ~30×(Hs/L) degrees
    yaw_rms   = deg2rad(par.Hs * 0.3);                % slow yaw: ~0.3°/m of Hs
    surge_rms = par.dock_drift * par.Hs;               % mooring drift scales with Hs
    sway_rms  = par.dock_drift * par.Hs * 0.5;
    
    hv = compute_rao_response(waves, par.dock_Tn_heave, par.dock_zeta_heave, heave_rms, t);
    rl = compute_rao_response(waves, par.dock_Tn_roll, par.dock_zeta_roll, roll_rms, t);
    pt = compute_rao_response(waves, par.dock_Tn_pitch, par.dock_zeta_pitch, pitch_rms, t);
    yw = compute_rao_response(waves, 30, 0.2, yaw_rms, t);
    sg = compute_rao_response(waves, 45, 0.15, surge_rms, t);
    sw = compute_rao_response(waves, 55, 0.12, sway_rms, t);
    
    % Compute anchor world positions at each epoch
    aw=zeros(par.M,3,N); da=zeros(par.M,N);
    for k=1:N
        R=eul2r(yw(k),pt(k),rl(k));
        tr=par.plat_pos+[sg(k);sw(k);hv(k)];
        for j=1:par.M
            aw(j,:,k)=(R*par.anch_body(j,:)'+tr)';
            da(j,k)=norm(squeeze(aw(j,:,k))'-par.anch_body(j,:)'-par.plat_pos);
        end
    end
    plat.hv=hv; plat.rl=rl; plat.pt=pt; plat.yw=yw;
    plat.sg=sg; plat.sw=sw;
    plat.aw=aw; plat.dm=mean(da(:)); plat.dx=max(da(:));
end

function waves = gen_waves(par)
    % JONSWAP wave spectrum generation
    % Returns struct with frequencies, amplitudes, and random phases
    % Same phases used by both ship and dock (correlated sea)
    
    nw = par.n_waves;
    wp = 2*pi/par.Tp;  % peak frequency (rad/s)
    
    % Frequency range: 0.3*wp to 3*wp
    w = linspace(0.3*wp, 3*wp, nw);
    dw = w(2)-w(1);
    
    % JONSWAP spectrum: S(w) = α·g²/w⁵ · exp(-5/4·(wp/w)⁴) · γ^r
    % where r = exp(-(w-wp)²/(2·σ²·wp²))
    % σ = 0.07 for w≤wp, 0.09 for w>wp
    g = 9.81;
    alpha = 0.0081;  % Phillips constant (adjusted for Hs)
    
    S = zeros(1, nw);
    for i = 1:nw
        if w(i) <= wp; sigma = 0.07; else; sigma = 0.09; end
        r = exp(-(w(i)-wp)^2 / (2*sigma^2*wp^2));
        S(i) = alpha * g^2 / w(i)^5 * exp(-5/4*(wp/w(i))^4) * par.gamma_js^r;
    end
    
    % Scale spectrum to match desired Hs: Hs = 4·sqrt(m0)
    m0 = trapz(w, S);
    Hs_computed = 4*sqrt(m0);
    S = S * (par.Hs / Hs_computed)^2;
    
    % Wave amplitudes from spectrum
    a = sqrt(2 * S * dw);
    
    % Random phases (SHARED between ship and dock — same ocean!)
    ph = 2*pi*rand(1, nw);
    
    waves.w = w;
    waves.a = a;
    waves.ph = ph;
    waves.nw = nw;
    waves.Hs = par.Hs;
    waves.Tp = par.Tp;
end

function x = compute_rao_response(waves, Tn, zeta, amp_scale, t)
    % Compute vessel response to wave spectrum using simplified RAO
    % RAO(w) = 1 / sqrt((1-(w/wn)²)² + (2·ζ·w/wn)²)
    %
    % amp_scale = desired RMS amplitude of output signal
    % The RAO shapes the frequency content; output is then
    % normalized so that RMS(output) = amp_scale
    
    wn = 2*pi/Tn;  % natural frequency
    N = length(t);
    x = zeros(1, N);
    
    for i = 1:waves.nw
        w = waves.w(i);
        r = w/wn;
        rao_mag = 1 / sqrt((1 - r^2)^2 + (2*zeta*r)^2);
        rao_phase = atan2(-2*zeta*r, 1-r^2);
        x = x + waves.a(i) * rao_mag * sin(w*t + waves.ph(i) + rao_phase);
    end
    
    % Normalize: scale so that RMS = amp_scale
    x_rms = sqrt(mean(x.^2));
    if x_rms > 1e-10
        x = x * (amp_scale / x_rms);
    end
end

function R = eul2r(y,p,r)
    cy=cos(y);sy=sin(y);cp=cos(p);sp=sin(p);cr=cos(r);sr=sin(r);
    R=[cy*cp,cy*sp*sr-sy*cr,cy*sp*cr+sy*sr;sy*cp,sy*sp*sr+cy*cr,sy*sp*cr-cy*sr;-sp,cp*sr,cp*cr];
end

function [sat_gps, sat_navic] = gen_satellites(par)
    sat_gps = [30 45; 85 25; 140 60; 195 35; 230 50; 280 20; 320 40; 355 15];
    % NavIC: 3 GEO (high el) + 4 GSO — these are overridden if real data loaded
    sat_navic = [170 65; 225 70; 310 55; 145 45; 200 50; 260 40; 340 35];
end

function mask = compute_dock_mask(d, par)
    if d > 300; mask.az_block = 0; mask.el_block = 0;
    else
        mask.az_block = min(90, 45 + 45*(1-d/300));
        mask.el_block = min(40, 10 + 30*(1-d/300));
    end
    mask.dock_az = 0;
end

function [nvis, pdop_val] = count_visible(sats, mask, par)
    if isempty(sats); nvis=0; pdop_val=99; return; end
    vis = true(size(sats,1),1);
    for s = 1:size(sats,1)
        az = sats(s,1); el = sats(s,2);
        az_diff = abs(wrapTo180(az - mask.dock_az));
        if az_diff < mask.az_block && el < mask.el_block; vis(s) = false; end
        if el < 10; vis(s) = false; end
    end
    nvis = sum(vis);
    if nvis >= 4
        sv = sats(vis,:);
        H = zeros(nvis, 4);
        for s = 1:nvis
            az_r = deg2rad(sv(s,1)); el_r = deg2rad(sv(s,2));
            H(s,:) = [cos(el_r)*sin(az_r), cos(el_r)*cos(az_r), sin(el_r), 1];
        end
        G = inv(H'*H);
        pdop_val = sqrt(G(1,1)+G(2,2)+G(3,3));
    else; pdop_val = 99; end
end

function meas = gen_meas(traj, plat, par, scen, sat_gps, sat_navic, gnss_mode)
    N=par.N; dt=par.dt; M=par.M; ref=par.ref;
    gp=round(1/par.gnss_rate/dt);
    gnss(N)=struct('valid',false,'z',[],'R',[]);
    gnss_att(N)=struct('valid',false,'z',[],'R',[]);
    uwb(N)=struct('valid',false,'z_tags',{cell(1,par.Ntags)},'aidx_tags',{cell(1,par.Ntags)},'apos',[]);

    for k=1:N
        gnss(k).valid=false; gnss_att(k).valid=false;
        uwb(k).valid=false; uwb(k).z_tags=cell(1,par.Ntags); uwb(k).aidx_tags=cell(1,par.Ntags);

        if mod(k-1,gp)==0
            rp=norm(traj.pos(:,k)-par.plat_pos);
            sg=par.sig_gnss; out=false;
            if rp<par.gnss_shadow_rng
                if scen>=2; out=true; else; sg=sg*par.gnss_mp_infl; end
            elseif rp<par.gnss_mp_rng; sg=sg*par.gnss_mp_infl*(rp/par.gnss_mp_rng); end
            if scen==3 && par.t(k)>180; out=true; end
            if ~out
                R_ship=eul2r(traj.att(1,k),traj.att(2,k),traj.att(3,k));
                gw=R_ship*par.gnss_body+traj.pos(:,k);
                mask_k = compute_dock_mask(rp, par);
                if strcmp(gnss_mode,'gps_navic'); sats_k=[sat_gps;sat_navic];
                else; sats_k=sat_gps; end
                [ns_k, pdop_k] = count_visible(sats_k, mask_k, par);
                
                % Position noise scales with PDOP (more sats → better DOP → less noise)
                pdop_ref = 2.0;  % reference PDOP for sig_gnss
                sig_pos = sg * min(pdop_k/pdop_ref, 5.0);  % cap at 5× inflation
                
                gnss(k).valid=true; gnss(k).z=gw+sig_pos*randn(3,1); gnss(k).R=sig_pos^2*eye(3);
                if ns_k >= 5
                    mp = 1+5*exp(-rp/80);
                    B = norm(par.gnss_ant_body(1,:)-par.gnss_ant_body(2,:));
                    sig_att_h = par.sig_carrier*mp / (B*sqrt(ns_k));
                    sig_att_rp = sig_att_h * 3;
                    p_amb = max(0, 1-exp(-0.8*(ns_k-5)));
                    if rand < p_amb
                        gnss_att(k).valid = true;
                        gnss_att(k).z = traj.att(:,k) + [sig_att_h; sig_att_rp; sig_att_rp].*randn(3,1);
                        gnss_att(k).R = diag([sig_att_h^2, sig_att_rp^2, sig_att_rp^2]);
                    end
                end
            end
        end

        if rand<par.pdrop; continue; end
        a=squeeze(plat.aw(:,:,k));
        R_ship=eul2r(traj.att(1,k),traj.att(2,k),traj.att(3,k));
        any_valid=false;
        for ti=1:par.Ntags
            tag_w=(R_ship*par.tags_body(ti,:)'+traj.pos(:,k))';
            rng_all=zeros(M,1);
            for j=1:M; rng_all(j)=norm(tag_w-a(j,:)); end
            if rng_all(ref)>par.uwb_max_rng; continue; end
            zt=[]; vi=[];
            for j=1:M
                if j==ref||rng_all(j)>par.uwb_max_rng; continue; end
                rd=rng_all(j)-rng_all(ref); bias=0;
                if rand<par.pnlos; bias=max(abs(par.sig_nlos*randn)*.8,0); end
                zt=[zt; rd+bias+par.sig_tdoa*randn]; vi=[vi; j];
            end
            if length(zt)>=3
                uwb(k).z_tags{ti}=zt; uwb(k).aidx_tags{ti}=vi; any_valid=true;
            end
        end
        if any_valid; uwb(k).valid=true; uwb(k).apos=a; end
    end
    meas.gnss=gnss; meas.gnss_att=gnss_att; meas.uwb=uwb;
end

function [est, dg] = run_ekf(traj, plat, meas, par, method, comp)
    N=par.N; dt=par.dt; ref=par.ref; nx=12;
    x=zeros(nx,1); x(1:3)=traj.pos(:,1)+5*randn(3,1);
    x(4:6)=traj.vel(:,1)+randn(3,1); x(7:9)=traj.att(:,1)+.1*randn(3,1);
    P=diag([25 25 25 4 4 4 .1 .05 .05 .01 .01 .01]);
    Qp=par.sig_a^2*[dt^4/4*eye(3),dt^3/2*eye(3);dt^3/2*eye(3),dt^2*eye(3)];
    Qa=diag([par.sig_om^2,par.sig_rp^2,par.sig_rp^2])*dt;
    Qw=diag([par.sig_om^2,par.sig_rp^2,par.sig_rp^2])*dt^2*.1;
    Q=blkdiag(Qp,Qa,Qw);
    est.pos=zeros(3,N); est.vel=zeros(3,N); est.att=zeros(3,N);
    est.pos(:,1)=x(1:3); est.vel(:,1)=x(4:6); est.att(:,1)=x(7:9);
    dg.nis=zeros(1,N);

    for k=2:N
        F=eye(nx); F(1:3,4:6)=dt*eye(3); F(7:9,10:12)=dt*eye(3);
        x=F*x; P=F*P*F'+Q; x(7:9)=wrapToPi(x(7:9));

        skip_gnss = any(strcmp(method,{'mtag_only'}));
        if ~skip_gnss && meas.gnss(k).valid
            R_est=eul2r(x(7),x(8),x(9));
            z_pred=R_est*par.gnss_body+x(1:3);
            Hg=zeros(3,nx); Hg(1:3,1:3)=eye(3);
            for ai=1:3; xp=x;xp(6+ai)=xp(6+ai)+1e-5;
                Rp=eul2r(xp(7),xp(8),xp(9));
                Hg(:,6+ai)=(Rp*par.gnss_body+xp(1:3)-z_pred)/1e-5; end
            nu=meas.gnss(k).z-z_pred; Sg=Hg*P*Hg'+meas.gnss(k).R;
            Kg=P*Hg'/Sg; x=x+Kg*nu;
            P=(eye(nx)-Kg*Hg)*P*(eye(nx)-Kg*Hg)'+Kg*meas.gnss(k).R*Kg';
            x(7:9)=wrapToPi(x(7:9));
        end

        use_gnss_att = any(strcmp(method,{'gps_4ant','gpsnav_4ant','mtag_gps','mtag_gpsnav'}));
        if use_gnss_att && isfield(meas,'gnss_att') && meas.gnss_att(k).valid
            Ha=zeros(3,nx); Ha(1,7)=1; Ha(2,8)=1; Ha(3,9)=1;
            nu_a = wrapToPi(meas.gnss_att(k).z - x(7:9));
            Sa=Ha*P*Ha'+meas.gnss_att(k).R;
            Ka=P*Ha'/Sa; x=x+Ka*nu_a;
            P=(eye(nx)-Ka*Ha)*P*(eye(nx)-Ka*Ha)'+Ka*meas.gnss_att(k).R*Ka';
            x(7:9)=wrapToPi(x(7:9));
        end

        do_tdoa = any(strcmp(method,{'mtag_only','mtag_gps','mtag_gpsnav','uwb1tag_gps'}));
        if do_tdoa && meas.uwb(k).valid
            a_all=meas.uwb(k).apos;
            if strcmp(comp,'none'); a_use=par.anch_body+par.plat_pos';
            elseif strcmp(comp,'seastate'); a_use=par.anch_body+par.plat_pos';
                a_use(:,3)=a_use(:,3)+plat.hv(k);
            elseif strcmp(comp,'gnss_plat'); a_use=a_all;
                for jj=1:size(a_use,1); a_use(jj,:)=a_use(jj,:)+0.1*randn(1,3); end
            else; a_use=a_all; end

            R_est=eul2r(x(7),x(8),x(9));
            ntags_use=par.Ntags;
            if strcmp(method,'uwb1tag_gps'); ntags_use=1; end

            all_nu=[]; all_H=[]; all_R=[];
            for ti=1:ntags_use
                zt=meas.uwb(k).z_tags{ti}; if isempty(zt); continue; end
                aidx=meas.uwb(k).aidx_tags{ti}; if isempty(aidx); continue; end
                m=min(length(zt),length(aidx)); zt=zt(1:m); aidx=aidx(1:m);
                tb=par.tags_body(ti,:)'; tw=R_est*tb+x(1:3);
                h=zeros(m,1); J=zeros(m,nx);
                r0=norm(tw-a_use(ref,:)');
                for ii=1:m; j=aidx(ii);
                    rj=norm(tw-a_use(j,:)'); ej=(tw-a_use(j,:)')/max(rj,.01);
                    e0=(tw-a_use(ref,:)')/max(r0,.01);
                    h(ii)=rj-r0; J(ii,1:3)=(ej-e0)';
                    for ai=1:3; xp=x;xp(6+ai)=xp(6+ai)+1e-5;
                        Rp=eul2r(xp(7),xp(8),xp(9)); twp=Rp*tb+xp(1:3);
                        rjp=norm(twp-a_use(j,:)'); r0p=norm(twp-a_use(ref,:)');
                        J(ii,6+ai)=((rjp-r0p)-h(ii))/1e-5; end
                end
                all_nu=[all_nu;zt-h]; all_H=[all_H;J];
                all_R=blkdiag(all_R,par.sig_tdoa^2*eye(m));
            end

            if ~isempty(all_nu)
                is_ad = any(strcmp(method,{'mtag_gpsnav','mtag_gps'}));
                if is_ad; Jp=all_H(:,1:3);
                    if size(Jp,1)>=3; kap=min(cond(Jp'*Jp),par.kappa_max);
                    else; kap=par.kappa_max; end
                    all_R=all_R*(1+par.beta_geo*log(1+(kap-1))); end
                St=all_H*P*all_H'+all_R;
                nis_v=all_nu'/St*all_nu; dg.nis(k)=nis_v;
                if is_ad && nis_v>par.tau_nis; all_R=all_R*par.gamma_nis;
                    St=all_H*P*all_H'+all_R; end
                Kt=P*all_H'/St; x=x+Kt*all_nu;
                P=(eye(nx)-Kt*all_H)*P*(eye(nx)-Kt*all_H)'+Kt*all_R*Kt';
                x(7:9)=wrapToPi(x(7:9));
            end
        end
        est.pos(:,k)=x(1:3); est.vel(:,k)=x(4:6); est.att(:,k)=x(7:9);
    end
end

function [cp, ch, rk] = compute_crlb_6dof(p, att, tags, par)
    M=par.M; ref=par.ref; R_s=eul2r(att(1),att(2),att(3));
    a=par.anch_body+par.plat_pos'; all_J=[];
    for ti=1:size(tags,1); tw=R_s*tags(ti,:)'+p;
        r0=norm(tw-a(ref,:)');
        for j=1:M; if j==ref; continue; end
            rj=norm(tw-a(j,:)'); ej=(tw-a(j,:)')/max(rj,.01);
            e0=(tw-a(ref,:)')/max(r0,.01); row=zeros(1,6); row(1:3)=(ej-e0)';
            for ai=1:3; ap=att;ap(ai)=ap(ai)+1e-5;
                Rp=eul2r(ap(1),ap(2),ap(3)); twp=Rp*tags(ti,:)'+p;
                rjp=norm(twp-a(j,:)'); r0p=norm(twp-a(ref,:)');
                row(3+ai)=((rjp-r0p)-(rj-r0))/1e-5; end
            all_J=[all_J;row]; end; end
    FIM=all_J'*(1/par.sig_tdoa^2*eye(size(all_J,1)))*all_J;
    rk = rank(FIM);
    if rk>=6; C=inv(FIM); cp=sqrt(trace(C(1:3,1:3))); ch=sqrt(C(4,4));
    else; cp=999; ch=999; end
end

function a = wrapTo180(a)
    a = mod(a+180,360)-180;
end