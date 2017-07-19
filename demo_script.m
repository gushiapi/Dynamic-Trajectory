%% This script shows an example of how to use optimalControlContinuous

clear
close all
clc

% load SC mat
load('SC.mat');

% normalize SC to ensure the convergence
N = length(A);
A = A/eigs(A,1);
A = A - eye(N);

% fix some parameters
rho = 100;  % balance between minimizing energy or minimizing distance from target state
T = 1;      % time to go from initial to target state
B = eye(N); % input matrix -- select specific columns to make those nodes control points
K = N-10;   % control set -- select specific 
KSet = sort(randperm(N,K)); % randon choose the control ndoes
BK = B(:,KSet);
% choose random initial and target states
x0 = randn(N,1);
xT = randn(N,1);

% run optimal control with all nodes
[x,u] = optimalControlContinuous(A,B,rho,x0,xT,T); % x is the state trajectory, u is the time-varying input

% run optimal control with control nodes
[xK,uK] = optimalControlContinuous(A,BK,rho,x0,xT,T); % x is the state trajectory, u is the time-varying input
