% System Params
m=0.8;
k=10;
c=0.3;
dt=1.0e-4;
t=20;


% State-Space model + transfer function
A= [   0	1;
	-k/m -c/m];

B = [0 1/m];

C = [1 0];

% Simulation via ode45
u = @(t) 20 * sin(15 * t);
f = @(t, x) [x(2); -(k/m)*x(1) - (c/m)*x(2) + u(t)/m];

tspan = 0:dt:t;
x0 = [0; 0];
[t_out, x_out] = ode45(f, tspan, x0);
y_out = x_out(:, 1);

% Plot
fig = figure(1);
set(fig, 'Position', [100 100 720 380]);
hold on;
plot(t_out, y_out, 'Color', [0.22 0.45 0.70], 'LineWidth', 1.5);
xlabel('Time [s]', 'FontSize', 11);
ylabel('x(t)  [m]', 'FontSize', 11);
title('\bfSystem Response — u(t) = 20 sin(15t)', 'FontSize', 13, 'Interpreter', 'tex');
grid on;
set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--', 'FontSize', 10, 'Box', 'on');
xlim([0 t]);
hold off;
print('../outputs/sim_response.svg', '-dsvg');


