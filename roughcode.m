clear; clc; close all;

%% Constants Definition
d = 9.525; % [mm] truss tube OD
t = 1.587; % [mm] thickness
djj = 250; % [mm] joint-to-joint distance
I = 2.475e6; % [mm^4] moment of inertia
E = 69000; % [MPa] = N/mm^2
P = 222.4; % [N] applied load
L = 4000; % [mm] span length

%% Discretization
x = linspace(0,L,200);

%% Deflection function (piecewise for midspan point load)
v = zeros(size(x));
for k = 1:numel(x)
    if x(k) <= L/2
        v(k) = -(P*x(k)/(48*E*I))*(3*L^2 - 4*x(k)^2);
    else
        v(k) = -(P*(L-x(k))/(48*E*I))*(3*L^2 - 4*(L-x(k))^2);
    end
end

%% Plot
figure;
plot(x, v, 'b-', 'LineWidth',1.5);
xlabel('Beam length, x [mm]');
ylabel('Deflection v(x) [mm]');
title('Deflection of simply supported beam under midspan load');
grid on;

%% Show max deflection
[maxDefl, idx] = min(v); % note: deflection is downward, so min value
fprintf('Max deflection = %.4f mm at x = %.1f mm\n', maxDefl, x(idx));