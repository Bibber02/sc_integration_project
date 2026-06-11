clear;
clc;
close all;

h = 0.001;
x0 = [pi pi 0 0];
y0 = [x0(1) x0(2)];

linearize_nonlinear_plant;

% Temporarily overwrite until we figure out linearize plant again
load('linearized_plant_eq_pi_pi_ts_0.001.mat')

calc_lqr;
hwinit;
open_system("rotpen_LQG.slx");