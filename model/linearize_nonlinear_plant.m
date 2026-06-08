clear;
clc;
close all;

lin = linearize_rotpendulum();

A = lin.A;
B = lin.B;
C = lin.C;
D = lin.D;
Ad = lin.Ad;
Bd = lin.Bd;
Cd = lin.Cd;
Dd = lin.Dd;
sys_lin = lin.sys_lin;
sys_disc = lin.sys_disc;
Ts = lin.Ts;
x0 = lin.x0;
u0 = lin.u0;
p = lin.p;
