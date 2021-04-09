clear all
%% Load Real Model
load('barrassmodel.mat')


%% Load State Space Model
load('model_5_4_50Hz.mat')

f = 0.1;
T = 300;

Ts=stateSpaceModel.Ts;    

% State Space Matrices
A = stateSpaceModel.A;
B = stateSpaceModel.B;
C = stateSpaceModel.C;
D = 0;
G = ss(A,B,C,D,Ts);

%% Regulator design
Q = C'*C;
R = 80;
[K,S,E] = dlqr(A,B,Q,R);

% Computation fo the external input gain (disturbance/ref)

N = inv([A-eye(size(A)), B; C,0])*[zeros(size(A,1),1);1];
Nx = N(1:end-1,:);
Nu = N(end,:);
Nbar = Nu+K*Nx;

%% Current estimator design

% %Method A
% w2 = 100;
% QE = eye(size(A))*w2;
% RE = 1;
% G1 = eye(size(A)); %in the place of B

% Method B
QE=1;
RE=1;

[M,P,Z,EE] = dlqe(A,B,C,QE,RE);

PHIE = A-M*C*A; 
GAMMAE = B-M*C*B;

