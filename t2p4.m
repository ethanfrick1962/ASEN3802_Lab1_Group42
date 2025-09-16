clc;clear;close all;
data = readmatrix('cases/Case_3.csv');
data = data(11:40,:);
v = -1.*data(:,6).*25.4; % midpoint displacement in mm
P = (data(:,1)./2).*4.44822162*1000; % total load over 2 converted to kg*mm/s^2
I = 2.475e6; % [mm^4] moment of inertia
inline = (data(:,5)+85.5).*4.44822162*1000; % inline loading converted to kg*mm/s^2
E = 69000.*1000; % [MPa] = N/mm^2
L = 4000; % [mm] span length

%a = roots([ones([62 1]) zeros([62 1]) ones(62,1).*(-3*(L^2)/4) (-v.*P./(6*E*I))]);
%ones()
%a = zeros(i,3);
%for i = 1:62

a = zeros(3,30);

for i = 1:30
    a(:,i) = roots([1 0 (-3*(L^2)/4) ((-v(i).*6.*E.*I)./P(i))]);
end

finala = mean(a(1,1:30))

a2 = (inline.*2.*I)./(39.6.*P.*250)
finala2 = mean(a2)


deflection = -P.*500*(3*4000^2-4*500^2)/(24*E*I);
deflection = deflection.*0.0393701;
internalStress = (P./1000).*500*250/(2*I);
internalForce = internalStress.*145.038*0.06138012;
measuredStress = inline./39.6;
measuredStress = measuredStress.*145.038*0.06138012;

figure()
subplot(1,2,1)
plot([0; data(:,1)],[0; data(:,6)])
hold on
plot([0; data(:,1)],[0; -deflection])
legend('Measured','Expected')
xlabel('Applied Load (lb)')
ylabel('Deflection at Midspan (in)')
title('Midspan Deflection vs Applied Load')

subplot(1,2,2)
plot([0; data(:,1)],[0; data(:,5)+85.5])
hold on
plot([0; data(:,1)],[0; internalForce])
legend('Measured','Expected')
xlabel('Applied Load (lb)')
ylabel('Internal Force (lb)')
title('Internal Stress vs Applied Load')

