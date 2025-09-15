clear; clc; close all;

%% ---- Inputs ----
E = 69000;        % [N/mm^2]
I = 2.475e6;      % [mm^4]
L = 4000;         % [mm]
P = 222.4;        % [N]
a = 3000;         % [mm] load position for Case 2
d = 500;          % [mm] distance from midspan to each load (Case 3)

npts = 500;
x = linspace(0,L,npts);

%% ===== CASE 1: Midspan load =====
RA1 = P/2; RB1 = P/2;

V1 = zeros(size(x));
V1(x < L/2) = RA1;
V1(x > L/2) = -RB1;
V1(x==L/2) = NaN;

M1 = zeros(size(x));
M1(x <= L/2) = RA1 * x(x <= L/2);
M1(x > L/2)  = RA1 * x(x > L/2) - P*(x(x > L/2) - L/2);

v1 = zeros(size(x));
for k = 1:numel(x)
    if x(k) <= L/2
        v1(k) = -(P*x(k)/(48*E*I))*(3*L^2 - 4*x(k)^2);
    else
        v1(k) = -(P*(L-x(k))/(48*E*I))*(3*L^2 - 4*(L-x(k))^2);
    end
end

N1 = zeros(size(x));

%% ===== CASE 2: Off-center load =====
RA2 = P*(L-a)/L; RB2 = P*a/L;
left  = x <= a; right = x >= a;

V2 = nan(size(x));
V2(left)  = RA2;
V2(right) = RA2 - P;
[~,ia] = min(abs(x - a)); V2(ia) = NaN;

M2 = zeros(size(x));
M2(left)  = RA2 .* x(left);
M2(right) = RA2 .* x(right) - P .* (x(right) - a);

v2 = zeros(size(x));
v2(left)  = (P .* x(left) .* ( -a*(2*L^2 - 3*L*a + a^2) + x(left).^2.*(L - a) )) ...
           ./ (6*E*I*L);
v2(right) = (P .* a .* ( L*(a^2 + 3*x(right).^2) - x(right).*(2*L^2 + a^2 + x(right).^2) )) ...
           ./ (6*E*I*L);

N2 = zeros(size(x));

%% ===== CASE 3: Two symmetric loads =====
RA3 = P; RB3 = P;
aL = L/2 - d; bL = L/2 + d;

V3 = zeros(size(x));
V3(x < aL) = RA3;
V3(x >= aL & x < bL) = RA3 - P;
V3(x >= bL) = RA3 - 2*P;
[~,ia] = min(abs(x - aL)); [~,ib] = min(abs(x - bL));
V3([ia ib]) = NaN;

M3 = zeros(size(x));
M3(x <= aL) = RA3 .* x(x <= aL);
M3(x >= aL & x <= bL) = RA3 .* x(x >= aL & x <= bL) - P.*(x(x >= aL & x <= bL)-aL);
M3(x >= bL) = RA3 .* x(x >= bL) - P.*(x(x >= bL)-aL) - P.*(x(x >= bL)-bL);

theta = cumtrapz(x, M3/(E*I));
v3 = cumtrapz(x, theta);
v3 = v3 - (v3(end)/L).*x;   % enforce BCs

N3 = zeros(size(x));

%% ===== PLOTS =====

% 1) Support Reactions
figure;
bar([0 L],[RA1 RB1],0.3,'FaceAlpha',0.4); hold on;
bar([0 L],[RA2 RB2],0.3,'FaceAlpha',0.4);
bar([0 L],[RA3 RB3],0.3,'FaceAlpha',0.4);
legend('Case 1 Midspan','Case 2 Off-center','Case 3 Two Loads');
title('Support Reactions');
xlabel('Location [mm]'); ylabel('Reaction [N]');
set(gca,'XTick',[0 L],'XTickLabel',{'A (x=0)','B (x=L)'}); grid on;

% 2) Shear Force
figure;
plot(x,V1,'-','LineWidth',2); hold on;
plot(x,V2,'--','LineWidth',2);
plot(x,V3,':','LineWidth',2);
legend('Case 1','Case 2','Case 3');
title('Shear Force Diagram');
xlabel('x [mm]'); ylabel('V(x) [N]'); grid on;

% 3) Bending Moment
figure;
plot(x,M1,'-','LineWidth',2); hold on;
plot(x,M2,'--','LineWidth',2);
plot(x,M3,':','LineWidth',2);
legend('Case 1','Case 2','Case 3');
title('Bending Moment Diagram');
xlabel('x [mm]'); ylabel('M(x) [N·mm]'); grid on;

% 4) Axial Force
figure;
plot(x,N1,'-','LineWidth',2); hold on;
plot(x,N2,'--','LineWidth',2);
plot(x,N3,':','LineWidth',2);
legend('Case 1','Case 2','Case 3');
title('Axial Force Diagram');
xlabel('x [mm]'); ylabel('N(x) [N]'); grid on;

% 5) Deflection
figure;
plot(x,v1,'-','LineWidth',2); hold on;
plot(x,v2,'--','LineWidth',2);
plot(x,v3,':','LineWidth',2);
legend('Case 1','Case 2','Case 3');
title('Deflection Curve');
xlabel('x [mm]'); ylabel('v(x) [mm] (down)'); grid on;





% Linear Regression and Uncertainty Analysis

% Sweep load magnitudes
Pvals = linspace(50,1000,20); % [N]
defl_case1 = zeros(size(Pvals));
defl_case2 = zeros(size(Pvals));
defl_case3 = zeros(size(Pvals));

for i = 1:numel(Pvals)
    Ptest = Pvals(i);

    % ---- Case 1 Midspan ----
    vtemp = zeros(size(x));
    for k = 1:numel(x)
        if x(k) <= L/2
            vtemp(k) = -(Ptest*x(k)/(48*E*I))*(3*L^2 - 4*x(k)^2);
        else
            vtemp(k) = -(Ptest*(L-x(k))/(48*E*I))*(3*L^2 - 4*(L-x(k))^2);
        end
    end
    defl_case1(i) = min(vtemp);

    % ---- Case 2 Off-center ----
    RA = Ptest*(L-a)/L; RB = Ptest*a/L;
    vtemp = zeros(size(x));
    vtemp(x<=a) = (Ptest .* x(x<=a) .* ( -a*(2*L^2 - 3*L*a + a^2) + x(x<=a).^2.*(L - a) )) ...
                  ./ (6*E*I*L);
    vtemp(x>=a) = (Ptest .* a .* ( L*(a^2 + 3*x(x>=a).^2) - x(x>=a).*(2*L^2 + a^2 + x(x>=a).^2) )) ...
                  ./ (6*E*I*L);
    defl_case2(i) = min(vtemp);

    % ---- Case 3 Two symmetric ----
    RA = Ptest; RB = Ptest;
    Mtemp = zeros(size(x));
    Mtemp(x <= aL) = RA .* x(x <= aL);
    Mtemp(x >= aL & x <= bL) = RA .* x(x >= aL & x <= bL) - Ptest.*(x(x >= aL & x <= bL)-aL);
    Mtemp(x >= bL) = RA .* x(x >= bL) - Ptest.*(x(x >= bL)-aL) - Ptest.*(x(x >= bL)-bL);
    theta = cumtrapz(x, Mtemp/(E*I));
    vtemp = cumtrapz(x, theta);
    vtemp = vtemp - (vtemp(end)/L).*x;
    defl_case3(i) = min(vtemp);
end

% Perform linear regression (y = m*P + c)
[p1, S1] = polyfit(Pvals,defl_case1,1);
[p2, S2] = polyfit(Pvals,defl_case2,1);
[p3, S3] = polyfit(Pvals,defl_case3,1);

fit1 = polyval(p1,Pvals,S1);
fit2 = polyval(p2,Pvals,S2);
fit3 = polyval(p3,Pvals,S3);

% Confidence intervals (95%)
[fit1_ci,delta1] = polyval(p1,Pvals,S1);
[fit2_ci,delta2] = polyval(p2,Pvals,S2);
[fit3_ci,delta3] = polyval(p3,Pvals,S3);

%% Plot regression curves
figure;
errorbar(Pvals,defl_case1,delta1,'o'); hold on;
plot(Pvals,fit1,'-','LineWidth',1.5);
errorbar(Pvals,defl_case2,delta2,'s');
plot(Pvals,fit2,'-','LineWidth',1.5);
errorbar(Pvals,defl_case3,delta3,'^');
plot(Pvals,fit3,'-','LineWidth',1.5);

xlabel('Load P [N]'); ylabel('Max deflection [mm]');
legend('Case1 data','Case1 fit','Case2 data','Case2 fit','Case3 data','Case3 fit','Location','best');
title('Linear Regression of Deflection vs Load');
grid on;

%% Report regression slopes & R^2
R2_case1 = 1 - (S1.normr^2 / norm(defl_case1-mean(defl_case1))^2);
R2_case2 = 1 - (S2.normr^2 / norm(defl_case2-mean(defl_case2))^2);
R2_case3 = 1 - (S3.normr^2 / norm(defl_case3-mean(defl_case3))^2);

fprintf('\n--- Linear Regression Results ---\n');
fprintf('Case 1 slope = %.4e mm/N, R^2 = %.4f\n',p1(1),R2_case1);
fprintf('Case 2 slope = %.4e mm/N, R^2 = %.4f\n',p2(1),R2_case2);
fprintf('Case 3 slope = %.4e mm/N, R^2 = %.4f\n',p3(1),R2_case3);
