function [ x, u ] = optimalControlContinuous( A, B, rho, x0, xT, T )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solve the optimal control problem where 
% min_u \int_0^T ||x_T-x(t)||^2+\rho ||u(t)||^2 dt
% s.t. dx(t)/dt = A x(t) + B u(t), x(0) = x_0, x(T) = x_T
% 
% Parameters:
%       A: nxn, transition matrix for the linear dynamics, n is the number
%          of nodes
%       B: nxk, control matrix denoting the location of control nodes, n is
%          the number of nodes and k is the number of control regions
%       rho: weight for squared energy loss
%       x0: 1xn vector, the initial state when t = 0, i.e. x(0) = x0
%       xT: 1xn vector, the target state when t = T, i.e. x(T) = xT
%        T: control horizon, i.e. time allowed to perform the control
% Output:
%       x: nx(nStep+1), the optimal trajectory, n is the number of regions,
%          nStep is the number steps to divide the control horizon
%       u: nx(nStep+1), the optimal control input, n is the number of 
%          regions, nStep is the number steps to divide the control horizon
% Reference: 
%       [1] Gu, Shi, Richard F. Betzel, Marcelo G. Mattar, Matthew Cieslak, 
%       Philip R. Delio, Scott T. Grafton, Fabio Pasqualetti, and Danielle S. Bassett.
%       "Optimal trajectories of brain state transitions." NeuroImage 148 (2017): 305-317.
%       [2] Betzel, Richard F., Shi Gu, John D. Medaglia, Fabio Pasqualetti, and Danielle S. Bassett.
%       "Optimally controlling the human connectome: the role of network topology." Scientific reports 6 (2016).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set the paramters and calculate x,u with matrix operations
nStep = 1000;
p = size(A,1);

% redefine the variables and reduce the equation to stardard linear
% equations
I = eye(p);
A_tilde = [A, -1.0/(2*rho)*(B*B');-2*I, -A'];
b = [zeros(p);I]*2*xT;
b_tilde = A_tilde\b;
b1 = b_tilde(1:p); b2 = b_tilde(p+1:end);

% compute the constant vector c = [c1;c2]
c1 = x0+b1;
E = expm(-A_tilde*T);

% E = [E11, E12; E21, E22];
E11 = E(1:p,1:p); E12 = E(1:p,p+1:2*p);
E21 = E(p+1:2*p,1:p); E22 = E(p+1:2*p, p+1:2*p);

% calculate PT, c2
PT = E12\(c1- E11*xT- E11*b1-E12*b2);
c2 = E21*xT + E22*PT + E21*b1 + E22*b2;
P0 = c2 - b2;
c = [c1;c2];

% calculate x,P, u
xstar = zeros(2*p, nStep);

tt = expm(1/nStep*A_tilde);
ct = c;
for i = 1:nStep
    ct = tt*ct;
    xstar(:,i) = ct - b_tilde;
end
xstar0 = [x0;P0];
xstar = [xstar0, xstar];
x = xstar(1:p,:);
P = xstar(p+1:2*p,:);
u = -B'*P/(2*rho);
end