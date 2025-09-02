%% Clear Workspace

clear; clc; close all;

%% Constants Definition

% d - Exterior diameter of truss tube [mm]
d = 9.525;

% t - Thickness of truss tube [mm]
t = 1.587;

% djj - Joint-to-joint distance [mm]
djj = 250;

% I - Moment of inertia of trusses from center [mm^4]
I = 2.475e6;

% E - Young's modulus for 6061-T6 aluminum [MPa]
E = 69000;

% P - Applied load on truss [N]
P = 222.4;

% L - Length of "beam" [mm]
L = 4000;

%% Defining Linspace

x = linspace(0, L, 100); % Define a linspace for the length of the beam

%% Plot: Part 1 Task 2 - Centered Load

v_centered = @(x) ...
    (1 / (E * I)) * ...
    (((P / 12) * (x .^ 3)) - (((P * (L ^ 2)) / (16)) .* x));

v_centered_solved = v_centered(x);

figure();
plot(x, v_centered_solved);

%% Plot: Part 1 Task 3 - Off-Centered Load



%% Plot: Part 1 Task 4 - 2 Symmetric Loads

% Note that 