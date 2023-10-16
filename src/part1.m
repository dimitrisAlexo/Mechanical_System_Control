% Run this script to run every function related to Part 1 of the Project

addpath('functions');
clearvars
clc
close all

% Choose initial conditions
x0 = [0 -5];

% Simulations
phase_portrait();
simulation(x0);

