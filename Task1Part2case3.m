clear; clc; close all;

%% ---- Inputs ----
E = 69000;        % [N/mm^2] Young's modulus
I = 2.475e6;      % [mm^4] second moment of area
P = 222.4;        % [N] each point load
L = 4000;         % [mm] span
d = 500;          % [mm] distance from midspan to each load

% Load locations (symmetric about midspan)
a = L/2 - d;      % left load position
b = L/2 + d;      % right load position

%% ---- Reactions (static equilibrium) ----
RA = P;           % [N]
RB = P;           % [N]

%% ---- Discretization ----
npts = 1001;
x = linspace(0,L,npts);

%% ---- Shear force V(x) [N] ----
V = zeros(size(x));

% Left of first load
left1 = x < a;
V(left1) = RA;

% Between the two loads
mid = x >= a & x < b;
V(mid) = RA - P;

% Right of second load
right1 = x >= b;
V(right1) = RA - P - P;   % = -P

% Put discontinuities at load points
[~,ia] = min(abs(x - a));
[~,ib] = min(abs(x - b));
V([ia ib]) = NaN;

%% ---- Bending moment M(x) [N*mm] ----
M = zeros(size(x));

% Region 1: 0 <= x < a
reg1 = x <= a;
M(reg1) = RA .* x(reg1);

% Region 2: a <= x < b
reg2 = (x >= a & x <= b);
M(reg2) = RA .* x(reg2) - P .* (x(reg2) - a);

% Region 3: b <= x <= L
reg3 = x >= b;
M(reg3) = RA .* x(reg3) - P .* (x(reg3) - a) - P .* (x(reg3) - b);

% Max moment occurs at midspan (symmetry)
Mmax = M(x == L/2);

%% ---- Deflection v(x) [mm], downward negative ----
% For multiple point loads, closed-form gets messy. 
% We'll solve using direct integration of EI*v''=M numerically.

dx = x(2)-x(1);
theta = cumtrapz(x, M/(E*I));          % slope (integration of curvature)
v = cumtrapz(x, theta);                % deflection

% Boundary condition corrections (v(0)=0, v(L)=0)
v = v - (v(end)/L).*x;  

% Max deflection (symmetric → at midspan)
[~,iVmax] = min(v);
vmax = v(iVmax);
xvmax = x(iVmax);

%% ---- Axial force N(x) [N] ----
N = zeros(size(x));  % zero for vertical point loads

%% ---- Plots ----
% 1) Reactions
figure;
bar([0 L],[RA RB],0.2);
title('Support Reactions');
xlabel('Location [mm]'); ylabel('Reaction force [N]');
set(gca,'XTick',[0 L],'XTickLabel',{'A (x=0)','B (x=L)'}); grid on;

% 2) Shear force
figure;
plot(x,V,'LineWidth',2); hold on;
xline(a,'--'); xline(b,'--');
legend('V(x)','Loads','Location','best');
title('Shear Force Diagram V(x)');
xlabel('x [mm]'); ylabel('Shear force [N]'); grid on;

% 3) Bending moment
figure;
plot(x,M,'LineWidth',2); hold on;
plot(xvmax,M(iVmax),'ro','MarkerFaceColor','r');
title('Bending Moment Diagram M(x)');
xlabel('x [mm]'); ylabel('Moment [N·mm]'); grid on;

% 4) Deflection
figure;
plot(x,v,'LineWidth',2); hold on;
plot(xvmax,vmax,'ro','MarkerFaceColor','r');
title('Deflection v(x) (downward negative)');
xlabel('x [mm]'); ylabel('Deflection [mm]'); grid on;

% 5) Axial force
figure;
plot(x,N,'LineWidth',2);
title('Axial Force N(x)');
xlabel('x [mm]'); ylabel('Axial force [N]'); grid on;

%% ---- Console summary ----
fprintf('Reactions: RA = %.3f N, RB = %.3f N\n', RA, RB);
fprintf('M_max at midspan: %.1f N·mm\n', M(L/2 == x));
fprintf('v_max at x = %.1f mm: v_max = %.4f mm (down)\n', xvmax, vmax);