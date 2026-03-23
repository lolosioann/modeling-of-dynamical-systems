% Part 2: Least Squares Parameter Estimation
% Mass-spring-damper: m*xdd + c*xd + k*x = u(t)

% --- Ground truth & simulation ---
m_true = 0.8;  k_true = 10;  c_true = 0.3;
dt = 1e-4;     T  = 20;

u_fn = @(t) 20 * sin(15 * t);
f    = @(t, x) [x(2); -(k_true/m_true)*x(1) - (c_true/m_true)*x(2) + u_fn(t)/m_true];

tspan = 0:dt:T;
[t_sim, x_sim] = ode45(f, tspan, [0; 0]);
x_sim = x_sim';   % 2 x M  (row 1: position, row 2: velocity)

% =========================================================
% Helper: filter-based LS (Section 2.3.1)
% Filter both sides with 1/Lambda(s), Lambda(s) = s + lambda
% Uses bilinear (Tustin) discretisation.
% Since xd is measurable: (1/Lambda)*xdd = (s/Lambda)*xd
% Returns theta = [m; c; k]
% =========================================================
function theta = run_ls_filter(x_s, xd_s, u_s, Ts, lambda)
    % Bilinear coefficients for H1(z) = 1/(s+lambda)
    b1 = Ts * [1,  1];
    a1 = [2 + lambda*Ts,  lambda*Ts - 2];
    % Bilinear coefficients for Hs(z) = s/(s+lambda)
    bs = 2  * [1, -1];
    % (same denominator a1)

    % Filtered signals
    z     = filter(b1, a1, u_s);    % (1/Lambda) * u
    zeta1 = filter(bs, a1, xd_s);  % (s/Lambda) * xd  [replaces xdd/Lambda]
    zeta2 = filter(b1, a1, xd_s);  % (1/Lambda) * xd
    zeta3 = filter(b1, a1, x_s);   % (1/Lambda) * x

    % Skip filter settling transient (~5 time constants = 5/lambda seconds)
    skip = ceil(5 / (lambda * Ts));
    ii   = (skip + 1) : length(z);
    if length(ii) < 10
        theta = nan(3, 1);
        return;
    end

    Z     = [zeta1(ii)', zeta2(ii)', zeta3(ii)'];
    z_vec = z(ii)';
    theta = (Z' * Z) \ (Z' * z_vec);
end

% =========================================================
% (α) Ts = 0.05 s
% =========================================================
lambda = 10;   % filter pole: Lambda(s) = s + lambda
Ts  = 0.05;
idx = 1 : round(Ts/dt) : length(t_sim);

x_s  = x_sim(1, idx);
xd_s = x_sim(2, idx);
u_s  = u_fn(t_sim(idx)');

theta0 = run_ls_filter(x_s, xd_s, u_s, Ts, lambda);
m_hat = theta0(1);  c_hat = theta0(2);  k_hat = theta0(3);

fprintf('--- (α) Ts = %.3f s ---\n', Ts);
fprintf('         m         c         k\n');
fprintf('True:  %8.4f  %8.4f  %8.4f\n', m_true, c_true, k_true);
fprintf('LS:    %8.4f  %8.4f  %8.4f\n', m_hat, c_hat, k_hat);
fprintf('err%%:  %8.4f  %8.4f  %8.4f\n\n', ...
    abs(m_hat-m_true)/m_true*100, ...
    abs(c_hat-c_true)/c_true*100, ...
    abs(k_hat-k_true)/k_true*100);

% Re-simulate with estimated params
f_hat = @(t, x) [x(2); -(k_hat/m_hat)*x(1) - (c_hat/m_hat)*x(2) + u_fn(t)/m_hat];
[~, x_hat_sim] = ode45(f_hat, tspan, [0; 0]);
x_true = x_sim(1, :)';
x_hat  = x_hat_sim(:, 1);
e_x    = x_hat - x_true;

W = 720;  H = 420;

fig1 = figure(1);
set(fig1, 'Position', [100 600 W H]);

subplot(2,1,1);
hold on;
plot(t_sim, x_true, 'Color', [0.22 0.45 0.70], 'LineWidth', 1.5, 'DisplayName', 'x(t)');
plot(t_sim, x_hat,  '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'DisplayName', 'x\^(t)');
ylabel('x(t)  [m]', 'FontSize', 11);
title('\bfΑπόκριση: πραγματική vs εκτιμώμενο μοντέλο', 'FontSize', 12, 'Interpreter', 'tex');
legend('x(t)', 'x\^(t)', 'Location', 'northeast', 'FontSize', 10);
grid on;
set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--', 'FontSize', 10, 'Box', 'on', 'XLim', [0 T]);
hold off;

subplot(2,1,2);
plot(t_sim, e_x, 'Color', [0.47 0.67 0.19], 'LineWidth', 1.2);
xlabel('Time [s]', 'FontSize', 11);
ylabel('e_x(t)  [m]', 'FontSize', 11, 'Interpreter', 'tex');
title('\bfΣφάλμα απόκρισης  e_x(t) = x\^(t) - x(t)', 'FontSize', 12, 'Interpreter', 'tex');
grid on;
set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--', 'FontSize', 10, 'Box', 'on', 'XLim', [0 T]);
print('../outputs/alpha_response.svg', '-dsvg');

% =========================================================
% (β) Effect of sampling period
% =========================================================
Ts_vec = logspace(-3, 0, 80);
err_m  = nan(size(Ts_vec));
err_c  = nan(size(Ts_vec));
err_k  = nan(size(Ts_vec));

for j = 1:length(Ts_vec)
    Tsj  = Ts_vec(j);
    idxj = 1 : round(Tsj/dt) : length(t_sim);
    if length(idxj) < 5, continue; end

    x_j  = x_sim(1, idxj);
    xd_j = x_sim(2, idxj);
    u_j  = u_fn(t_sim(idxj)');

    try
        th = run_ls_filter(x_j, xd_j, u_j, Tsj, lambda);
        err_m(j) = abs(th(1) - m_true) / m_true * 100;
        err_c(j) = abs(th(2) - c_true) / c_true * 100;
        err_k(j) = abs(th(3) - k_true) / k_true * 100;
    catch
    end
end

fig2 = figure(2);
set(fig2, 'Position', [100 100 W H]);
hold on;
semilogx(Ts_vec, err_m, 'Color', [0.22 0.45 0.70], 'LineWidth', 1.8, 'DisplayName', '\epsilon_m');
semilogx(Ts_vec, err_c, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.8, 'DisplayName', '\epsilon_c');
semilogx(Ts_vec, err_k, 'Color', [0.47 0.67 0.19], 'LineWidth', 1.8, 'DisplayName', '\epsilon_k');
plot([0.05 0.05], [0 max([err_m, err_c, err_k])], '--k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
text(0.05, max([err_m, err_c, err_k])*0.92, '  T_s=0.05', 'FontSize', 9);
xlabel('T_s  [s]', 'FontSize', 11);
ylabel('\epsilon_p  [%]', 'FontSize', 11, 'Interpreter', 'tex');
title('\bfΣφάλμα εκτίμησης παραμέτρων vs T_s', 'FontSize', 12, 'Interpreter', 'tex');
legend('Location', 'northwest', 'FontSize', 10);
grid on;
set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--', 'FontSize', 10, 'Box', 'on');
hold off;
print('../outputs/beta_sampling.svg', '-dsvg');

% =========================================================
% (γ) White Gaussian noise — Monte Carlo
% =========================================================
sigma_x  = 0.05 * std(x_s);
sigma_xd = 0.05 * std(xd_s);
N_mc     = 500;

thetas_mc = zeros(3, N_mc);
for mc = 1:N_mc
    x_noisy  = x_s  + sigma_x  * randn(size(x_s));
    xd_noisy = xd_s + sigma_xd * randn(size(xd_s));
    thetas_mc(:, mc) = run_ls_filter(x_noisy, xd_noisy, u_s, Ts, lambda);
end

theta_mean = mean(thetas_mc, 2);
theta_std  = std(thetas_mc, 0, 2);
err_mc = abs(theta_mean - [m_true; c_true; k_true]) ./ [m_true; c_true; k_true] * 100;

fprintf('--- (γ) Monte Carlo (%d runs), sigma_x=%.2e  sigma_xd=%.2e ---\n', N_mc, sigma_x, sigma_xd);
fprintf('            m         c         k\n');
fprintf('True:    %8.4f  %8.4f  %8.4f\n', m_true,        c_true,        k_true);
fprintf('Clean:   %8.4f  %8.4f  %8.4f\n', theta0(1),     theta0(2),     theta0(3));
fprintf('MC mean: %8.4f  %8.4f  %8.4f\n', theta_mean(1), theta_mean(2), theta_mean(3));
fprintf('MC std:  %8.4f  %8.4f  %8.4f\n', theta_std(1),  theta_std(2),  theta_std(3));
fprintf('Err clean [%%]:  %.4f  %.4f  %.4f\n', ...
    abs(theta0(1)-m_true)/m_true*100, abs(theta0(2)-c_true)/c_true*100, abs(theta0(3)-k_true)/k_true*100);
fprintf('Err MC mean [%%]: %.4f  %.4f  %.4f\n', err_mc(1), err_mc(2), err_mc(3));

% Histogram of MC estimates
fig3 = figure(3);
set(fig3, 'Position', [100 100 W H]);
params     = {'m', 'c', 'k'};
true_vals  = [m_true, c_true, k_true];
clean_vals = [theta0(1), theta0(2), theta0(3)];
colors     = {[0.22 0.45 0.70], [0.85 0.33 0.10], [0.47 0.67 0.19]};

for p = 1:3
    subplot(1, 3, p);
    hist(thetas_mc(p,:), 30);
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', colors{p}, 'FaceAlpha', 0.75, 'EdgeColor', 'none');
    hold on;
    xline_val = true_vals(p);
    yl = ylim();
    plot([xline_val xline_val], yl, 'k-',  'LineWidth', 1.5, 'DisplayName', 'True');
    plot([clean_vals(p) clean_vals(p)], yl, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'DisplayName', 'Clean LS');
    xlabel(sprintf('\\hat{%s}', params{p}), 'FontSize', 11);
    ylabel('Count', 'FontSize', 10);
    title(sprintf('\\bf%s', params{p}), 'FontSize', 12, 'Interpreter', 'tex');
    if p == 1
        legend('True', 'Clean LS', 'Location', 'northwest', 'FontSize', 8);
    end
    grid on;
    set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--', 'FontSize', 9, 'Box', 'on');
    hold off;
end
print('../outputs/gamma_montecarlo.svg', '-dsvg');
