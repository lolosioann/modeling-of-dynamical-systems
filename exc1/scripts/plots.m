% Pole plot and impulse response for mass-spring-damper
% G(s) = 1 / (m*s^2 + c*s + k)

pkg load control;

m = 0.8;
k = 10.0;
c = 0.3;

num = 1;
den = [m, c, k];

sys = tf(num, den);
p = roots(den);
fprintf('Poles: %.4f %+.4fj,  %.4f %+.4fj\n', ...
	real(p(1)), imag(p(1)), real(p(2)), imag(p(2)));

W = 560;
H = 420;
xl = [-0.8 0.4];
yl = [-4.2  4.2];

% --- Pole plot ---
fig1 = figure(1);
set(fig1, 'Position', [100 100 W H]);

pzmap(sys);
ax1 = gca();
set(ax1, 'GridAlpha', 0.3, 'GridLineStyle', '--', ...
	'FontSize', 10, 'Box', 'on', 'LineWidth', 0.8, ...
	'XLim', xl, 'YLim', yl);
title('\bfPole–Zero Map', 'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'tex');
xlabel('Real', 'FontSize', 11);
ylabel('Imaginary', 'FontSize', 11);
print('../outputs/pole_plot.svg', '-dsvg');

% --- Impulse response ---
t = 0:1e-3:20;
[y, t_out] = impulse(sys, t);

fig2 = figure(2);
set(fig2, 'Position', [700 100 W H]);

ax2 = axes();
hold on;

fill([t_out; flipud(t_out)], [y; zeros(size(y))], ...
	[0.27 0.51 0.71], 'FaceAlpha', 0.13, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(t_out, y, 'Color', [0.22 0.45 0.70], 'LineWidth', 1.8);
plot(xlim, [0 0], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7, 'HandleVisibility', 'off');

xlabel('Time [s]', 'FontSize', 11);
ylabel('x(t)  [m / N\cdots]', 'FontSize', 11, 'Interpreter', 'tex');
title('\bfImpulse Response', 'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'tex');
grid on;
set(ax2, 'GridAlpha', 0.3, 'GridLineStyle', '--', ...
	'FontSize', 10, 'Box', 'on', 'LineWidth', 0.8, 'XLim', [0 20]);
hold off;
print('../outputs/impulse_response.svg', '-dsvg');
