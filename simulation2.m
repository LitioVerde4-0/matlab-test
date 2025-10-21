% Lithium Extraction Model with Physics-Based Milling Influence
% Based on LDF Kinetics: dC/dt = k_eff * (C_max - C)
clc; clear; close all

%% USER PARAMETERS
milling_time_hr = 2;          % user-defined milling time [hr]
milling_exponent = 0.5;       % n in R ∝ t_m^-n
R0 = 10e-6;                   % initial particle radius [m]
D_eff = 1e-14;                % effective diffusivity [m^2/s]
mechanism = "diffusion";      % choose "film" or "diffusion"

%% Physical k_eff calculation
R = R0 * milling_time_hr^-milling_exponent;  % updated radius
switch mechanism
    case "film"
        k_eff = 1 / R;  % external mass transfer
    case "diffusion"
        k_eff = D_eff / R^2;  % internal diffusion
end

%% Leaching Setup
C_max = 200;        % max lithium concentration in ppm
C0 = 0;             % initial concentration
tspan = [0 120*60];  % 120 minutes in seconds

%% Solve the ODE
dCdt = @(t, C) k_eff * (C_max - C);
[t, C] = ode45(dCdt, tspan, C0);
fprintf('Final Li⁺ Concentration = %.2f ppm\n', C(end));
fprintf('k_eff = %.4f 1/s\n', k_eff);


%% 3D Animation of Clay Particle (Color = Li in solution)
figure('Color','w')
[Xs, Ys, Zs] = sphere(40);  % higher resolution
r = 1;
Xs = r * Xs; Ys = r * Ys; Zs = r * Zs;

for i = 1:10:length(t)
    green_intensity = C(i)/C_max;  % scale 0 to 1
    cla
    surf(Xs, Ys, Zs, ...
        'FaceColor', [1-green_intensity, 1, 1-green_intensity], ...
        'EdgeColor', 'none');
    axis equal off
    light; lighting gouraud
    camlight headlight
    title(sprintf('Li in Solution: %.1f ppm\nTime = %.0f s', C(i), t(i)), 'FontSize', 14)
    drawnow
end

%% Plot Concentration vs Time
figure;
plot(t/60, C, 'LineWidth', 2); grid on
xlabel('Time (min)'); ylabel('Li⁺ Concentration (ppm)');
title(sprintf('LDF Leaching Model | Milling = %.1f hr (%s-limited)', ...
    milling_time_hr, mechanism))
