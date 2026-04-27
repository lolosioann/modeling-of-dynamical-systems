clear; close all; clc;

%% =========================================================
%  Lab 2 - Part 2: Lyapunov Estimation
% Lolos Ioannis - 10674

%% System parameters
g = 9.81;
l = 1.5;     
c = 0.5;   

%% Observer gains
lambda1 = 20;
lambda2 = 100;

%% Adaptation gains
gamma1 = 5;
gamma2 = 2;

%% Simulation
T     = 40;
dt    = 1e-3;
tspan = 0:dt:T;
u     = @(t) 0.5 * sin(t);

%% ---- a) no disturbance ----
fprintf('Part a)\n');
run_case(g, l, c, gamma1, gamma2, lambda1, lambda2, u, 0, tspan, '(a) no noise');

%% ---- b) errors vs noise amplitude ----
eta0_list = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2];
el_final  = zeros(size(eta0_list));
ec_final  = zeros(size(eta0_list));

for i = 1:numel(eta0_list)
    [el_final(i), ec_final(i)] = run_case_noisy( ...
        g, l, c, gamma1, gamma2, lambda1, lambda2, u, eta0_list(i), tspan);
end

% plot errors vs noise
figure('Name', 'err_vs_eta', 'Position', [100 100 720 380]);
plot(eta0_list, abs(el_final), '-o', 'Color', [0.22 0.45 0.70], 'LineWidth', 1.5); hold on;
plot(eta0_list, abs(ec_final), '-s', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
grid on;
set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--', 'Box', 'on', 'FontSize', 10);
xlabel('$\eta_0$',               'Interpreter', 'latex');
ylabel('Final $|$error$|$',      'Interpreter', 'latex');
title('\textbf{Estimation error vs noise amplitude}', 'Interpreter', 'latex');
legend('$|e_l|$', '$|e_c|$', 'Location', 'northwest', 'Interpreter', 'latex');

%% Estimation for eta_0 = 0.02
mid_idx = ceil(numel(eta0_list) / 2);
fprintf('\nPart b) eta0 = %.3f =====\n', eta0_list(mid_idx));
run_case(g, l, c, gamma1, gamma2, lambda1, lambda2, u, ...
    eta0_list(mid_idx), tspan, ...
    sprintf('(\\beta)\\ \\eta_0 = %.3f', eta0_list(mid_idx)));


%% Helpers

function run_case(g, l, c, g1, g2, lam1, lam2, u, eta0, tspan, tag)
    [tt, X] = simulate(g, l, c, g1, g2, lam1, lam2, u, eta0, tspan);

    theta     = X(:,1);
    theta_hat = X(:,3);
    th1_hat   = X(:,5);
    th2_hat   = X(:,6);

    l_hat = g ./ th1_hat;
    c_hat = th2_hat;

    e_theta = theta - theta_hat;
    e_l     = l_hat - l;
    e_c     = c_hat - c;

    % theta vs theta_hat
    figure('Name', ['theta_' tag], 'Position', [100 100 760 380]);
    plot(tt, theta,     'Color', [0.22 0.45 0.70], 'LineWidth', 1.5); hold on;
    plot(tt, theta_hat, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2, 'LineStyle', '--');
    grid on;
    set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--', 'Box', 'on', 'FontSize', 10);
    xlabel('Time [s]',                        'Interpreter', 'latex');
    ylabel('$\theta,\ \hat{\theta}$ [rad]',   'Interpreter', 'latex');
    title(['Lyapunov estimator'], 'Interpreter', 'latex');
    legend('$\theta(t)$', '$\hat{\theta}(t)$', 'Location', 'best', 'Interpreter', 'latex');

    % Errors
    figure('Name', ['err_' tag], 'Position', [100 100 760 520]);

    subplot(3,1,1);
    plot(tt, e_theta, 'Color', [0.22 0.45 0.70], 'LineWidth', 1.2); grid on;
    ylabel('$e_{\theta}(t)$', 'Interpreter', 'latex');
    title('\textbf{Estimation errors}', 'Interpreter', 'latex');
    set(gca, 'GridAlpha', 0.3, 'Box', 'on');

    subplot(3,1,2);
    plot(tt, e_l, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2); grid on;
    ylabel('$e_l(t)$', 'Interpreter', 'latex');
    set(gca, 'GridAlpha', 0.3, 'Box', 'on');

    subplot(3,1,3);
    plot(tt, e_c, 'Color', [0.10 0.60 0.30], 'LineWidth', 1.2); grid on;
    ylabel('$e_c(t)$',     'Interpreter', 'latex');
    xlabel('Time [s]',     'Interpreter', 'latex');
    set(gca, 'GridAlpha', 0.3, 'Box', 'on');

    fprintf('  l_hat(T) = %.4f  (true = %.3f,  e_l = %+.4f)\n', l_hat(end), l, e_l(end));
    fprintf('  c_hat(T) = %.4f  (true = %.3f,  e_c = %+.4f)\n', c_hat(end), c, e_c(end));
    fprintf('  |e_theta(T)| = %.4e\n\n', abs(e_theta(end)));
end


function [el, ec] = run_case_noisy(g, l, c, g1, g2, lam1, lam2, u, eta0, tspan)
    [~, X] = simulate(g, l, c, g1, g2, lam1, lam2, u, eta0, tspan);

    N       = round(1 / mean(diff(tspan)));
    th1_avg = mean(X(max(1,end-N):end, 5));
    th2_avg = mean(X(max(1,end-N):end, 6));

    el = g / th1_avg - l;
    ec = th2_avg - c;
end


function [tt, X] = simulate(g, l, c, g1, g2, lam1, lam2, u, eta0, tspan)
    
    % [θ; θ'; θ_hat; x2_hat; θ1_hat; θ2_hat]
    theta0 = 0.2;
    x0     = [theta0; 0; theta0; 0; (g/l)*0.9; c*0.9];

    opts    = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    f       = @(t, x) rhs(t, x, g, l, c, g1, g2, lam1, lam2, u, eta0);
    [tt, X] = ode45(f, tspan, x0, opts);
end


function dx = rhs(t, x, g, l, c, g1, g2, lam1, lam2, u, eta0)

    theta   = x(1);
    thetad  = x(2);
    theta_h = x(3);
    x2_h    = x(4);
    th1_h   = x(5);
    th2_h   = x(6);

    eta = eta0 * sin(20*pi*t);

    theta_m  = theta  + eta;
    thetad_m = thetad + eta;
    u_m      = u(t)   + eta;

    theta_ddot = -(g/l)*sin(theta) - c*thetad + u(t);

    e1 = theta_m  - theta_h;
    e2 = thetad_m - x2_h;

    th_h_dot = x2_h  + lam1 * e1;
    x2h_dot  = -th1_h*sin(theta_m) - th2_h*thetad_m + u_m + lam2 * e1;

    dth1 = g1 * e2 * sin(theta_m);
    dth2 = g2 * e2 * thetad_m;

    dx = [thetad;
          theta_ddot;
          th_h_dot;
          x2h_dot;
          dth1;
          dth2];
end