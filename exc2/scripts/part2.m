% Lab02 - Part 2: Lyapunov-based estimation of l, c (pendulum)
% Plant: theta_ddot = -(g/l) sin(theta) - c*theta_dot + u
% Reparametrize: theta1* = g/l, theta2* = c
% Series-parallel identifier with state x2_hat:
%   x2_hat_dot = -theta1_hat*sin(theta) - theta2_hat*theta_dot + u + am*(theta_dot - x2_hat)
% Lyapunov updates:
%   theta1_hat_dot = gamma1 * e * sin(theta)
%   theta2_hat_dot = gamma2 * e * theta_dot
% Recovery: l_hat = g / theta1_hat, c_hat = theta2_hat

clear; close all; clc;

%% Parameters
g   = 9.81;
l   = 1.5;       % true (in [1, 2.5])
c   = 0.5;       % true (in [0.2, 0.8])

gamma1 = 8;
gamma2 = 8;
am     = 5;

T   = 40;
dt  = 1e-3;
tspan = 0:dt:T;

u = @(t) 0.5 * sin(t);

%% ---- (α) noise-free ----
run_case(g, l, c, gamma1, gamma2, am, u, 0, tspan, '(α) no noise');

%% ---- (β) noise sweep ----
eta0_list = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2];
el_final = zeros(size(eta0_list));
ec_final = zeros(size(eta0_list));
for i = 1:numel(eta0_list)
    [el_final(i), ec_final(i)] = run_case_noisy(g, l, c, gamma1, gamma2, am, u, eta0_list(i), tspan);
end

% single example run for mid eta0
mid_idx = ceil(numel(eta0_list)/2);
run_case(g, l, c, gamma1, gamma2, am, u, eta0_list(mid_idx), tspan, ...
    sprintf('(β) eta_0 = %.3f', eta0_list(mid_idx)));

% final error vs noise amplitude
figure('Name','err_vs_eta','Position',[100 100 720 380]);
plot(eta0_list, abs(el_final), '-o','Color',[0.22 0.45 0.70],'LineWidth',1.5); hold on;
plot(eta0_list, abs(ec_final), '-s','Color',[0.85 0.33 0.10],'LineWidth',1.5);
grid on; set(gca,'GridAlpha',0.3,'GridLineStyle','--','Box','on','FontSize',10);
xlabel('\eta_0'); ylabel('Final |error|');
title('\bfParameter estimation error vs noise amplitude','Interpreter','tex');
legend('|e_l|','|e_c|','Location','northwest');

%% ---- helpers ----
function run_case(g, l, c, g1, g2, am, u, eta0, tspan, tag)
    [tt, X] = simulate(g, l, c, g1, g2, am, u, eta0, tspan);
    theta      = X(:,1); thetad  = X(:,2);
    theta_hat  = X(:,3);
    th1_hat    = X(:,5); th2_hat = X(:,6);

    l_hat = g ./ th1_hat;
    c_hat = th2_hat;

    e_theta = theta - theta_hat;
    e_l     = l_hat - l;
    e_c     = c_hat - c;

    % theta vs theta_hat
    figure('Name',['theta_' tag],'Position',[100 100 760 380]);
    plot(tt, theta,     'Color',[0.22 0.45 0.70],'LineWidth',1.5); hold on;
    plot(tt, theta_hat, 'Color',[0.85 0.33 0.10],'LineWidth',1.2,'LineStyle','--');
    grid on; set(gca,'GridAlpha',0.3,'GridLineStyle','--','Box','on','FontSize',10);
    xlabel('Time [s]'); ylabel('\theta, \thetâ [rad]');
    title(['\bfLyapunov estimator — ' tag],'Interpreter','tex');
    legend('\theta(t)','\thetâ(t)','Location','best');

    % errors
    figure('Name',['err_' tag],'Position',[100 100 760 520]);
    subplot(3,1,1);
    plot(tt, e_theta,'Color',[0.22 0.45 0.70],'LineWidth',1.2); grid on;
    ylabel('e_{\theta}(t)'); title(['\bfEstimation errors — ' tag],'Interpreter','tex');
    set(gca,'GridAlpha',0.3,'Box','on');
    subplot(3,1,2);
    plot(tt, e_l,'Color',[0.85 0.33 0.10],'LineWidth',1.2); grid on;
    ylabel('e_l(t)'); set(gca,'GridAlpha',0.3,'Box','on');
    subplot(3,1,3);
    plot(tt, e_c,'Color',[0.10 0.60 0.30],'LineWidth',1.2); grid on;
    ylabel('e_c(t)'); xlabel('Time [s]'); set(gca,'GridAlpha',0.3,'Box','on');

    fprintf('--- %s ---\n', tag);
    fprintf('  l_hat(T)=%.4f (true %.3f, err=%.3e)\n', l_hat(end), l, e_l(end));
    fprintf('  c_hat(T)=%.4f (true %.3f, err=%.3e)\n', c_hat(end), c, e_c(end));
    fprintf('  final|e_theta|=%.3e\n\n', abs(e_theta(end)));
end

function [el, ec] = run_case_noisy(g, l, c, g1, g2, am, u, eta0, tspan)
    [~, X] = simulate(g, l, c, g1, g2, am, u, eta0, tspan);
    th1_hat = X(end,5);  th2_hat = X(end,6);
    % average over last second to reduce HF ripple influence
    N = round(1/mean(diff(tspan)));
    th1_avg = mean(X(max(1,end-N):end,5));
    th2_avg = mean(X(max(1,end-N):end,6));
    el = g/th1_avg - l;
    ec = th2_avg - c;
end

function [tt, X] = simulate(g, l, c, g1, g2, am, u, eta0, tspan)
    % State: [theta; theta_dot; theta_hat; x2_hat; theta1_hat; theta2_hat]
    % theta_hat is integrated as d/dt theta_hat = x2_hat (using estimator state)
    x0 = [0.2; 0; 0.2; 0; g/2.0; 0.1];  % init theta1_hat so l_hat(0)=2
    f  = @(t, x) rhs(t, x, g, l, c, g1, g2, am, u, eta0);
    [tt, X] = ode45(f, tspan, x0);
end

function dx = rhs(t, x, g, l, c, g1, g2, am, u, eta0)
    theta   = x(1); thetad = x(2);
    theta_h = x(3); x2_h   = x(4);
    th1_h   = x(5); th2_h  = x(6);

    eta = eta0 * sin(20*pi*t);

    % noisy measurements used by estimator
    theta_m  = theta  + eta;
    thetad_m = thetad + eta;
    u_m      = u(t)   + eta;

    % plant (true dynamics, no noise)
    theta_ddot = -(g/l)*sin(theta) - c*thetad + u(t);

    % estimator state x2_hat dynamics (uses noisy measurements)
    x2h_dot = -th1_h*sin(theta_m) - th2_h*thetad_m + u_m + am*(thetad_m - x2_h);

    % theta_hat reconstruction
    th_h_dot = x2_h;

    % error for adaptation (uses noisy theta_dot)
    e = thetad_m - x2_h;

    % Lyapunov updates
    dth1 = g1 * e * sin(theta_m);
    dth2 = g2 * e * thetad_m;

    dx = [thetad;
          theta_ddot;
          th_h_dot;
          x2h_dot;
          dth1;
          dth2];
end
