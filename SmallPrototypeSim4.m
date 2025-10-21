% Lithium Recovery Simulation – Mechanochemical Milling + Water Leach
% FINAL version with physics-aware Monte Carlo and validated kinetics

clear; clc; close all;

%% Fixed Inputs
ore_mass = 100;            % grams of ore
salt = 'NaCl';             % salt type
T_leach_base = 363.15;     % 90°C in Kelvin
t_leach = 3600;            % 1 hour leach

N = 1000;                  % Monte Carlo runs

%% Baseline Constants
x_Li = 0.015;              % 1.5% Li content
t_mill = 2; BPR = 8;       % 2 hr milling, 8:1 BPR
C_Li = 1.2; I = 0.1;

%% === BASELINE SIMULATION ===

conv = milling_conv(t_mill, salt, BPR);
leach = leach_shrinking_core_T(t_leach, T_leach_base);  % validated
ppt = precipitation_yield(C_Li, I);
pen = min(max(1 - 0.15 / (1 + 0.4 * conv), 0.7), 1.0);   % bounded impurity penalty
ppt = ppt * pen;

Li_out = lithium_yield(ore_mass, x_Li, conv, leach, ppt);
fprintf('\n--- BASELINE ---\n');
fprintf('Conv: %.1f%% | Leach: %.1f%% | Ppt: %.1f%% | Li recovered: %.3f g\n', ...
    100*conv, 100*leach, 100*ppt, Li_out);
fprintf('Overall yield vs Li-in-ore: %.1f%%\n\n', 100 * Li_out / (ore_mass * x_Li));

%% === MONTE CARLO SIMULATION ===

% Ore variability (wider CV)
xLi_samples = lognrnd(log(0.015), 0.2, [N 1]);
Li_available = ore_mass .* xLi_samples;

% BPR/time correlation
z = randn(N,1);
BPR_samples = min(max(8 + 1*z, 5), 10);
t_mill_samples = min(max(2 - 0.5*z, 0.5), 3.5);

% Temperature variation (60–100°C)
T_samples = min(max(normrnd(363.15, 5, [N 1]), 333), 373);
C_Li_samples = min(max(normrnd(1.2, 0.2, [N 1]), 0), 2.0);
I_samples = min(max(normrnd(0.1, 0.05, [N 1]), 0), 1.0);

Li_outputs = zeros(N,1);

for i = 1:N
    conv_i = milling_conv(t_mill_samples(i), salt, BPR_samples(i));
    leach_i = leach_shrinking_core_T(t_leach, T_samples(i));
    ppt_i = precipitation_yield(C_Li_samples(i), I_samples(i));

    % Apply impurity penalty (bounded)
    pen = min(max(1 - 0.15 / (1 + 0.4 * conv_i), 0.7), 1.0);
    ppt_i = ppt_i * pen;

    eta_total = conv_i * leach_i * ppt_i;
    Li_outputs(i) = min(Li_available(i), Li_available(i) * eta_total);
end

% Stats
P10 = prctile(Li_outputs, 10);
P50 = prctile(Li_outputs, 50);
P90 = prctile(Li_outputs, 90);

fprintf('--- MONTE CARLO (N = %d) ---\n', N);
fprintf('P10: %.2f g | P50: %.2f g | P90: %.2f g\n', P10, P50, P90);

%% === PLOT ===

figure;
histogram(Li_outputs, 30, 'FaceColor', [0.2 0.5 0.8]); hold on;
xline(Li_out, 'r', 'LineWidth', 2, 'Label', 'Baseline');
xline(ore_mass * x_Li, '--k', 'LineWidth', 1.5, 'Label', 'Li in ore');
xlabel('Recovered Lithium (g)');
ylabel('Frequency');
title('Monte Carlo: Lithium Recovery Distribution');
grid on;

%% === FUNCTIONS ===

function eta_conv = milling_conv(t_mill, salt, BPR)
    BPRsat = 10;
    fBPR = 1 + 0.6 * (BPR ./ BPRsat) ./ (1 + BPR ./ BPRsat);
    s = 1.0 + 0.2 * strcmpi(salt, 'MgCl2');
    k = 0.9 * fBPR * s;   % 1/h
    eta_conv = 1 - exp(-k .* max(t_mill, 0));
    eta_conv = min(max(eta_conv, 0), 1);
end

function eta = leach_shrinking_core_T(t, T)
    % Calibrated constants: 95% yield at 90°C in 1 h
    k0 = 3.3; Ea = 25e3; R = 8.314;
    kL = k0 * exp(-Ea ./ (R * T));
    eta = 1 - exp(-kL .* t);
    eta = min(max(eta, 0), 1);
end

function eta = precipitation_yield(C, I)
    Csat = 0.9 + 0.3 * I;
    S = max(C ./ max(Csat, 1e-6), 0); n = 3;
    eta_kin = S.^n ./ (1 + S.^n);
    eta_imp = exp(-1.2 * I);
    eta = eta_kin .* eta_imp;
    eta = min(max(eta, 0), 1);
end

function Li_out = lithium_yield(ore_mass, x_Li, conv, leach, ppt)
    Li_out = ore_mass .* x_Li .* conv .* leach .* ppt;
end
