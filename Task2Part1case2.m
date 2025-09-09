clear;
clc;
close all;

%% ---- Inputs ----
E = 69000;        % [N/mm^2] Young's modulus
I = 2.475e6;      % [mm^4] second moment of area
P = 222.4;        % [N] point load
L = 4000;         % [mm] span
a = 3000;         % [mm] load location from left support (set a > L/2)

if ~(a > L/2 && a < L)
    error('Choose a load location a such that L/2 < a < L.');
end

%% ---- Reactions (static equilibrium) ----
RA = P*(L - a)/L;   % [N]
RB = P*a/L;         % [N]

%% ---- Discretization ----
npts = 1001;
x = linspace(0,L,npts);
left  = x <= a;
right = x >= a;

%% ---- Shear force V(x) [N] ----
V = nan(size(x));
V(left)  = RA;
V(right) = RA - P;

% Make the jump visible at x=a
[~,iMid] = min(abs(x - a));
V(iMid) = NaN;

%% ---- Bending moment M(x) [N*mm] ----
M = zeros(size(x));
M(left)  = RA .* x(left);
M(right) = RA .* x(right) - P .* (x(right) - a);

Mmax = RA*a;             % also equals RB*(L-a) = P*a*(L-a)/L
xMmax = a;

%% ---- Deflection v(x) [mm], downward negative
% Piecewise closed-form solution for simply supported beam, point load at a
% Derived from EI v'' = M with v(0)=v(L)=0 and continuity at x=a.
% Region 1: 0 <= x <= a
%   v1(x) = P*x*( -a*(2L^2 - 3La + a^2) + x^2*(L-a) ) / (6*E*I*L)
% Region 2: a <= x <= L
%   v2(x) = P*a*( L*(a^2 + 3x^2) - x*(2L^2 + a^2 + x^2) ) / (6*E*I*L)
v = zeros(size(x));
v(left)  = (P .* x(left) .* ( -a*(2*L^2 - 3*L*a + a^2) + x(left).^2.*(L - a) )) ...
           ./ (6*E*I*L);
v(right) = (P .* a .* ( L*(a^2 + 3*x(right).^2) - x(right).*(2*L^2 + a^2 + x(right).^2) )) ...
           ./ (6*E*I*L);

% Max deflection (numerical locate)
[~,iVmax] = min(v);  % most negative
vmax = v(iVmax);     % [mm]
xvmax = x(iVmax);    % [mm]

%% ---- Axial force N(x) [N] ----
N = zeros(size(x));  % zero for vertical point load on simply supported beam

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
xline(a,'--'); legend('V(x)','Load at x=a','Location','best');
title('Shear Force Diagram V(x)');
xlabel('x [mm]'); ylabel('Shear force [N]'); grid on;

% 3) Bending moment
figure;
plot(x,M,'LineWidth',2); hold on;
plot(xMmax,Mmax,'ro','MarkerFaceColor','r');
xline(a,'--');
text(xMmax,Mmax,sprintf('  M_{max}=%.1f N·mm',Mmax), ...
    'VerticalAlignment','bottom','HorizontalAlignment','left');
title('Bending Moment Diagram M(x)');
xlabel('x [mm]'); ylabel('Moment [N·mm]'); grid on;

% 4) Deflection
figure;
plot(x,v,'LineWidth',2); hold on;
plot(xvmax,vmax,'ro','MarkerFaceColor','r');
xline(a,'--');
text(xvmax,vmax,sprintf('  v_{max}=%.3f mm at x=%.1f mm',vmax,xvmax), ...
    'VerticalAlignment','bottom','HorizontalAlignment','left');
title('Deflection v(x) (downward negative)');
xlabel('x [mm]'); ylabel('Deflection [mm]'); grid on;

% 5) Axial force
figure;
plot(x,N,'LineWidth',2); xline(a,'--');
title('Axial Force N(x)'); xlabel('x [mm]'); ylabel('Axial force [N]'); grid on;

%% ---- Console summary ----
fprintf('Reactions: RA = %.3f N, RB = %.3f N\n', RA, RB);
fprintf('M_max at x = a = %.1f mm: M_max = %.1f N·mm\n', a, Mmax);
fprintf('v_max at x = %.1f mm: v_max = %.4f mm (down)\n', xvmax, vmax);