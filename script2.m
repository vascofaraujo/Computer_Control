%% Load State Space Model
load('model_5_4_50Hz.mat')

Ts=stateSpaceModel.Ts;
        
% State Space Matrices
A = stateSpaceModel.A;
B = stateSpaceModel.B;
C = stateSpaceModel.C;
D = 0;
G = ss(A,B,C,D,Ts);

%% Regulator design
Q = C'*C;
R = 10;
[K,S,E] = dlqr(A,B,Q,R);

% Computation fo the external input gain (ref)
N = inv([A-eye(size(A)), B; C,0])*[zeros(size(A,1),1);1];
Nx = N(1:end-1,:);
Nu = N(end,:);
Nbar = Nu+K*Nx;

%% Rate of decay confirmation - LQR

% Get log of output
figure(10)
plot(log(abs(out.y.signals.values)));
title('Rate of decay R= 10- Log Scale')
% hold on
% Dx=200:Ts:1000;
% Dy=Mc*Dx + 4;
% plot(Dx,Dy)
% legend('Log of y','Aproximate slope')


% Points used to calculate slope M (  )
[X,Y] = ginput(2);
                   
%Calculate slope
Mc = (Y(2)-Y(1))/(X(2)-X(1));
%Mc=diff(Y)/diff(X);

%Get log of the closed loop system poles
lambda_max = log(abs(E(5,1)));% polo dominante 

%Confirmation 
confirmation_LQR=(Mc/lambda_max)*100;


%% Current estimator design

% % Method A
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

figure(6)
subplot(1,2,1)
zplane([],eig(PHIE))
subplot(1,2,2)
zplane([],EE)


%% Rate of decay confirmation - LQE

% Get log of estimation error
figure(11)
plot(log(abs(out.estimation_error.signals.values)));
title('Rate of decay RE = 0.01 - Log Scale')

% Points used to calculate slope M (  )
[X,Y] = ginput(2); 

%Calculate slope
Me = (Y(2)-Y(1))/(X(2)-X(1));

%Get log of the closed loop system poles
lambdaE_max = log(abs(EE(5,1)));% polo dominante pode nao ser o E(1,1)

%Confirmation
confirmation_LQE=(Me/lambdaE_max)*100;


%% Analysis of LQR

% Controlled System
lqr = ss((A-B*K),B,C,0,Ts);

%Transfer Function
lqr_tf = tf(lqr);

%Closed Loop Gain LQR
gain_lqr = sum(lqr_tf.Numerator{1, 1})/sum(lqr_tf.Denominator{1, 1});

% PZ Maps
figure(1)
pzplot(G,'r',lqr,'g')
legend('System','LQR')
title('Pole-Zero Map of the System and LQR Control R=40')

% Comparison between poles from dlqr and the CLoop poles of LQR
figure(2)
[NC_lqr,DC_lqr] = tfdata(lqr,'v');
subplot(1,2,1)
zplane(NC_lqr,DC_lqr)
ax = axis;
subplot(1,2,2)
zplane([],E)
axis(ax)

% Time response
figure(3)
subplot(1,2,1)
impulse(lqr)
subplot(1,2,2)
step(lqr)

% Frequency
figure(4)
bode(lqr)

% Symmetric Root Locus
[NumG,DenG] = tfdata(G,'v');
[NumG,DenG] = eqtflength(NumG,DenG);

SRL = tf(conv(NumG,fliplr(NumG)),conv(DenG,fliplr(DenG)));
p_srl = rlocus(SRL,1/R); %slr poles

figure(5)

subplot(1,2,2)
zplane([],p_srl)
title(['SRL poles for \rho = 1/R = ' num2str(1/R)])
ax = axis;

subplot(1,2,1)
rlocus(SRL)
hold on; zplane([],[]); hold off
axis(2*ax)

%% Analysis of the LQG

%State-Space Representation
lqg = ss([A -B*K; M*C*A PHIE-GAMMAE*K-M*C*B*K],[B; M*C*B+GAMMAE]*Nbar,[C zeros(size(C))],0,Ts);

%Transfer Function
lqg_tf=tf(lqg);

%Closed Loop Gain LQG
gain_lqg=sum(lqg_tf.Numerator{1, 1})/sum(lqg_tf.Denominator{1, 1});

%Gain and Phase Margins - Bode
figure(12)
margin(lqg_tf)
%Time Response






