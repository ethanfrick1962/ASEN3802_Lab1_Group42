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



%% Support reactions (symmetry)
RA = P/2; % [N]
RB = P/2; % [N]

%% Discretization
npts = 500;
x = linspace(0, L, npts);


%% Shear force diagram V(x)
V = zeros(size(x));
for k = 1:numel(x)
    if x(k) < L/2
        V(k) = RA;
    elseif x(k) > L/2
        V(k) = -RB;
    else
        V(k) = NaN; % break line at midspan jump
    end
end

%% Bending moment diagram M(x)
M = zeros(size(x));
for k = 1:numel(x)
    if x(k) <= L/2
        M(k) = RA * x(k);
    else
        M(k) = RA * x(k) - P * (x(k) - L/2);
    end
end

%% Axial force diagram N(x) (zero for vertical load)
N = zeros(size(x));

%% --- PLOTS ---
% 1. Reaction forces
figure;
bar([0 L],[RA RB],0.2);
title('Support Reactions');
xlabel('Location [mm]'); ylabel('Reaction force [N]');
set(gca,'XTick',[0 L],'XTickLabel',{'A1 (x=0)','A2 (x=L)'});
grid on;

% 2. Shear force
figure;
plot(x,V,'LineWidth',2);
title('Shear Force Diagram V(x)');
xlabel('x [mm]'); ylabel('Shear force [N]');
grid on;

% 3. Bending moment
figure;
plot(x,M, 'LineWidth', 2);
title('Bending Moment Diagram M(x)');
xlabel('x [mm');
ylabel('Moment [N*mm');
grid on;

% 5. Axial force
figure;
plot(x,N,'LineWidth',2);
title('Axial Force Diagram N(x)');
xlabel('x [mm]'); ylabel('Axial force [N]');
grid on;
