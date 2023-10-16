% Run this script to run every function related to Part 2 of the Project

addpath('functions');
clearvars
clc
close all

% Choose initial conditions
x0 = [0 -5];

% Simulations
sliding_1(x0);
sliding_2(x0);
redesign_1(x0);
redesign_2(x0);

