% Lab02 - Part 1: Gradient method estimation of k, b
% Plant: J*phi_ddot = -k*phi_dot + b*u + d
% Linear parametric model: y = theta*^T * psi
%   y = J*phi_ddot, psi = [phi_dot; u], theta* = [-k; b]
% Mixed (series-parallel) identifier; parameter updates via gradient law.
% A separate parallel reconstruction of phi_hat is integrated for plotting.

clear; close all; clc;

%% Parameters
J  = 0.025;
k  = 0.2;        % true (in [0.1, 0.5])
b  = 1.0;        % true (in [0.5, 1.5])

gamma1 = 200;    % adaptation gain for k_hat
gamma2 = 200;    % adaptation gain for b_hat

T   = 60;
dt  = 1e-3;
tspan = 0:dt:T;

u = @(t) 0.25 * sin(0.5*pi*t);

%% Run both cases
run_case(J, k, b, gamma1, gamma2, u, @(t) 0,                 tspan, '(α) d(t) = 0',           '_a');
run_case(J, k, b, gamma1, gamma2, u, @(t) 0.02*sin(2*t),     tspan, '(β) d(t) = 0.02 sin(2t)', '_b');

%% ---- helpers ----
function run_case(J, k, b, g1, g2, u, d, tspan, tag, suffix)
    % State: x = [phi; phi_dot; phi_hat; phi_hat_dot; k_hat; b_hat]
    x0 = [0; 0; 0; 0; 0.3; 0.8];   % init k_hat=0.3, b_hat=0.8 (away from zero)
    f  = @(t, x) rhs(t, x, J, k, b, g1, g2, u, d);
    [tt, X] = ode45(f, tspan, x0);

    phi     = X(:,1);  phid    = X(:,2);
    phi_h   = X(:,3);  phid_h  = X(:,4);
    k_hat   = X(:,5);  b_hat   = X(:,6);

    e_phi = phi - phi_h;
    e_k   = k_hat - k;
    e_b   = b_hat - b;

    % --- phi vs phi_hat ---
    fig1 = figure('Name',['phi' suffix],'Position',[100 100 760 380]);
    plot(tt, phi,   'Color',[0.22 0.45 0.70],'LineWidth',1.5); hold on;
    plot(tt, phi_h, 'Color',[0.85 0.33 0.10],'LineWidth',1.2,'LineStyle','--');
    grid on; set(gca,'GridAlpha',0.3,'GridLineStyle','--','Box','on','FontSize',10);
    xlabel('Time [s]'); ylabel('$\phi,\ \hat{\phi}$ [rad]','Interpreter','latex');
    title(['Gradient estimator — ' tag]);
    legend({'$\phi(t)$','$\hat{\phi}(t)$'},'Interpreter','latex','Location','best');

    % --- errors ---
    fig2 = figure('Name',['err' suffix],'Position',[100 100 760 520]);
    subplot(3,1,1);
    plot(tt, e_phi, 'Color',[0.22 0.45 0.70],'LineWidth',1.2); grid on;
    ylabel('e_{\phi}(t)'); title(['\bfEstimation errors — ' tag],'Interpreter','tex');
    set(gca,'GridAlpha',0.3,'Box','on');
    subplot(3,1,2);
    plot(tt, e_k, 'Color',[0.85 0.33 0.10],'LineWidth',1.2); grid on;
    ylabel('e_k(t)'); set(gca,'GridAlpha',0.3,'Box','on');
    subplot(3,1,3);
    plot(tt, e_b, 'Color',[0.10 0.60 0.30],'LineWidth',1.2); grid on;
    ylabel('e_b(t)'); xlabel('Time [s]'); set(gca,'GridAlpha',0.3,'Box','on');

    fprintf('--- %s ---\n', tag);
    fprintf('  k_hat(T)=%.4f (true %.3f, err=%.3e)\n', k_hat(end), k, e_k(end));
    fprintf('  b_hat(T)=%.4f (true %.3f, err=%.3e)\n', b_hat(end), b, e_b(end));
    fprintf('  max|e_phi|=%.3e, final|e_phi|=%.3e\n\n', max(abs(e_phi)), abs(e_phi(end)));
end

function dx = rhs(t, x, J, k, b, g1, g2, u, d)
    phi     = x(1); phid    = x(2);
    phi_h   = x(3); phid_h  = x(4);
    k_hat   = x(5); b_hat   = x(6);

    ut = u(t); dt = d(t);

    % plant
    phi_ddot = (-k*phid + b*ut + dt) / J;

    % mixed-structure gradient error: e = J*phi_ddot + k_hat*phi_dot - b_hat*u
    e = J*phi_ddot + k_hat*phid - b_hat*ut;

    % parallel reconstruction of phi_hat (for plotting only)
    phi_hddot = (-k_hat*phid_h + b_hat*ut) / J;

    % gradient updates
    dk_hat = -g1 * e * phid;
    db_hat =  g2 * e * ut;

    dx = [phid;
          phi_ddot;
          phid_h;
          phi_hddot;
          dk_hat;
          db_hat];
end
