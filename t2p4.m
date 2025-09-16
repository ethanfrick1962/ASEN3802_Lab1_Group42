clc;clear;close all;
data = readmatrix('cases/Case_3.csv');
data = data(11:40,:);
v = -1.*data(:,6).*25.4; % midpoint displacement in mm
P = (data(:,1)./2).*4.44822162*1000; % total load over 2 converted to kg*mm/s^2
I = 2.475e6; % [mm^4] moment of inertia
inline = (data(:,5)+85.5).*4.44822162*1000;
E = 69000; % [MPa] = N/mm^2
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
