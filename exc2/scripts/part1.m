clear; close all; clc;

%% =========================================================
%  Lab 2 - Part 1: Gradient Descent
% Lolos Ioannis - 10674

%% System parameters
J = 0.025;
k = 0.2;   
b = 1.0;

%% Adaptation gains
gamma1 = 200;
gamma2 = 200;

%% Simulation
T     = 60;
dt    = 1e-3;
tspan = 0:dt:T;
u     = @(t) 0.25 * sin(0.5*pi*t); 

%% ---- a) no disturbance ----
fprintf('Part a)\n');
run_case(J, k, b, gamma1, gamma2, u, @(t) 0, tspan, '(α) d(t) = 0', '_a');

%% ---- b) with disturbance: 0.02sin(2t) ----
fprintf('Part b) d(t) = 0.02 sin(2t)\n');
run_case(J, k, b, gamma1, gamma2, u, @(t) 0.02*sin(2*t), tspan, ...
    '(\beta) d(t) = 0.02 sin(2t)', '_b');


%% HELPER FUNCS

function run_case(J, k, b, g1, g2, u, d, tspan, tag, suffix)
    % The first guesses are 1.5* the true value to
    % clearly see the convergence
    %
    % [φ; φ'; φ_hat; φ'_hat; k_hat; b_hat]
    x0 = [0; 0; 0; 0; k*1.5; b*1.5];

    opts    = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    f       = @(t, x) rhs(t, x, J, k, b, g1, g2, u, d);
    [tt, X] = ode45(f, tspan, x0, opts);

    phi   = X(:,1);  phid  = X(:,2);
    phi_h = X(:,3);
    k_hat = X(:,5);  b_hat = X(:,6);

    e_phi = phi - phi_h;
    e_k   = k_hat - k;
    e_b   = b_hat - b;

    %% phi vs phi_hat plot
    figure('Name', ['phi' suffix], 'Position', [100 100 760 380]);
    plot(tt, phi,   'Color', [0.22 0.45 0.70], 'LineWidth', 1.5); hold on;
    plot(tt, phi_h, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2, 'LineStyle', '--');
    grid on;
    set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--', 'Box', 'on', 'FontSize', 10);
    xlabel('Time [s]');
    ylabel('$\phi,\ \hat{\phi}$ [rad]', 'Interpreter', 'latex');
    title(['Gradient estimator — ' tag]);
    legend({'$\phi(t)$', '$\hat{\phi}(t)$'}, 'Interpreter', 'latex', 'Location', 'best');

    %% errors plot
    figure('Name', ['err' suffix], 'Position', [100 100 760 520]);

    subplot(3,1,1);
    plot(tt, e_phi, 'Color', [0.22 0.45 0.70], 'LineWidth', 1.2); grid on;
    ylabel('e_{\phi}(t)');
    title(['\bfEstimation errors — ' tag], 'Interpreter', 'tex');
    set(gca, 'GridAlpha', 0.3, 'Box', 'on');

    subplot(3,1,2);
    plot(tt, e_k, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2); grid on;
    ylabel('e_k(t)');
    set(gca, 'GridAlpha', 0.3, 'Box', 'on');

    subplot(3,1,3);
    plot(tt, e_b, 'Color', [0.10 0.60 0.30], 'LineWidth', 1.2); grid on;
    ylabel('e_b(t)'); xlabel('Time [s]');
    set(gca, 'GridAlpha', 0.3, 'Box', 'on');

    fprintf('  k_hat(T) = %.4f  (true = %.3f,  e_k = %+.4f)\n', k_hat(end), k, e_k(end));
    fprintf('  b_hat(T) = %.4f  (true = %.3f,  e_b = %+.4f)\n', b_hat(end), b, e_b(end));
    fprintf('  max|e_phi| = %.3e,  final|e_phi| = %.3e\n\n', max(abs(e_phi)), abs(e_phi(end)));
end


function dx = rhs(t, x, J, k, b, g1, g2, u, d)
    phi   = x(1);  phid  = x(2);
    phi_h = x(3);  phid_h = x(4);
    k_hat = x(5);  b_hat  = x(6);

    ut = u(t);
    dt = d(t);

    phi_ddot = (-k*phid + b*ut + dt) / J;

    e = J*phi_ddot + k_hat*phid - b_hat*ut;

    % For the plots
    phi_hddot = (-k_hat*phid_h + b_hat*ut) / J;

    dk_hat = -g1 * e * phid;
    db_hat =  g2 * e * ut;

    dx = [phid;
          phi_ddot;
          phid_h;
          phi_hddot;
          dk_hat;
          db_hat];
end