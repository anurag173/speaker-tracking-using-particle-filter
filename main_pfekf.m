clear all;clc;
%% INITIALISATION AND PARAMETERS:
%x=[1.85:0.01:3.05; 1.85:0.01:3.05];                 % True trajectory of source / sequence of states
x=[1.85:0.01:3.05 ; [1.85*ones(1,60),2.5,3.05*ones(1,60)]]; % True trajectory of source / sequence of states
T = size(x,2);                                            % Number of time steps.
M=4;                                                % No. of microphone pairs
R = diag([1e-5 1e-5 ]);                             % Used to evaluate the likelihood
Q = diag([1e-7 1e-5]);                              % Used to evaluate the prior
P01 = 0.1;                                          % PF's Variance of x-ccordinate of each particle
P02 = 0.1;                                          % PF's Variance of y-ccordinate of each particle
N = 10;                                             % Number of particles.		    
Q_pfekf = 10*1e-5*eye(2);                           % EKF's process noise variance.
R_pfekf = 1e-6*eye(2);                              % EKF's measurement noise variance.
p(:,:,1)=[1,2.3;1,2.7];                             % positions of microphones
p(:,:,2)=[2.3,4;2.7,4];
p(:,:,3)=[4,2.7;4,2.3];
p(:,:,4)=[2.7,1;2.3,1];
track_x=zeros(1,T);
track_y=zeros(1,T);                                           
C = repmat(eye(2,2),[1,1,T]);                       % Jacobian matrix for EKF
%C=(M,2,T);
tau=ones(M,1,T);
xparticle_pfekf = ones(2,T,N);                      % Variable for storing particles at all times
Pparticle_pfekf = cell(N,1);                        % Particles for the covariance of x.
%% GENERATING PARTICLES AND THEIR (CO)VARIANCES (EACH PARTICLE IS 2x1)
for i=1:N,                                          
  xparticle_pfekf(1,1,i) =  sqrt(2.5) + randn(1,1); % for n particles at time t=1
  xparticle_pfekf(2,1,i) =  sqrt(2.5)+ randn(1,1);
  Pparticle_pfekf{i} = ones(2,2,T);
  for t=1:T,
    Pparticle_pfekf{i}(:,:,t)= diag([P01 P02]); 
  end;
end;
%% PLOT PARTICLES
a=xparticle_pfekf(:,1,:); b= a(:,:); plot(b(1,:),b(2,:),'o');
%% SOME DATA STRUCTURES FOR THE FILTER
xparticlePred_pfekf = ones(2,T,N);                  % One-step-ahead predicted values of the states.
PparticlePred_pfekf = Pparticle_pfekf;              % One-step-ahead predicted values of P.
yPred_pfekf = ones(2,T,N);                          % One-step-ahead predicted values of y.
w = ones(T,N);                                      % Importance weights.
muPred_pfekf = ones(2,T);                           % EKF O-s-a estimate of the mean of the states.
PPred_pfekf = ones(2,2);                            % EKF O-s-a estimate of the variance of the states.
mu_pfekf = ones(2,T,N);                             % EKF estimate of the mean of the states.
P_pfekf = ones(2,2,T);                              % EKF estimate of the variance of the states.
%% Main Particle filter loop
for t=2:T,    
  %fprintf('PF-EKF : t = %i / %i  \r \n',t,T);
  %x_est=mean(yPred_pfekf(:,t,:),3);
  %y(:,:,t)=C(:,:,t)*x_est;
  %for j=1:M
  %    C(j,:,t) = ((1/343)*((x_est - p(1,:,j)')/norm(x_est - p(1,:,j)') - (x_est - p(2,:,j)')/norm(x_est - p(2,:,j)')))';
  %    tau(j,1,t) = (1/343)*(norm(x_est - p(1,:,j)') - norm(x_est - p(2,:,j)'));
  %end
  for i=1:N,
    %% PREDICTION STEP:
    muPred_pfekf(:,t) = xparticle_pfekf(:,t-1,i);                                       % muPred_pfekf has no i dependence
    PPred_pfekf = Pparticle_pfekf{i}(:,:,t-1) + Q_pfekf; 
    %yPredTmp = feval('bshfun',muPred_pfekf(:,t),u(:,t),t);
    yPredTmp=C(:,:,t)*muPred_pfekf(:,t);
    %% APPLY THE EKF UPDATE EQUATIONS FOR EACH PARTICLE
    S = R_pfekf + C(:,:,t)*PPred_pfekf*C(:,:,t)';                                       % Innovations covariance.
    K = PPred_pfekf*C(:,:,t)'*inv(S);                                                   % Kalman gain.
    mu_pfekf(:,t,i) = muPred_pfekf(:,t) + K*(x(:,t)-yPredTmp);                          % Mean of proposal.
    P_pfekf(:,:,t) = PPred_pfekf - K*C(:,:,t)*PPred_pfekf;                              % Variance of proposal.
    xparticlePred_pfekf(:,t,i) = mu_pfekf(:,t,i) + sqrtm(P_pfekf(:,:,t))*randn(2,1);
    PparticlePred_pfekf{i}(:,:,t) = P_pfekf(:,:,t);
  end;
  %% EVALUATE IMPORTANCE WEIGHTS FOR PARTICLE FILTER
  for i=1:N,
    %yPred_pfekf(:,t,i) = feval('bshfun',xparticlePred_pfekf(:,t,i),u(:,t),t); 
    yPred_pfekf(:,t,i) = xparticlePred_pfekf(:,t,i);
    lik = exp(-0.5*(x(:,t)-yPred_pfekf(:,t,i))'*inv(R)*(x(:,t)-yPred_pfekf(:,t,i)) ) + 0.4;
    prior = exp(-0.5*(xparticlePred_pfekf(:,t,i)- xparticle_pfekf(:,t-1,i))'*inv(Q) * (xparticlePred_pfekf(:,t,i)-xparticle_pfekf(:,t-1,i) ))+ 1e-99;
    proposal = inv(sqrt(det(PparticlePred_pfekf{i}(:,:,t)))) * exp(-0.5*(xparticlePred_pfekf(:,t,i)-mu_pfekf(:,t,i))'*inv(PparticlePred_pfekf{i}(:,:,t)) * (xparticlePred_pfekf(:,t,i)-mu_pfekf(:,t,i)))+ 1e-99;
    w(t,i) = lik*prior/proposal;      
  end;  
  %% NORMALISATION OF WEIGHTS
  w(t,:) = w(t,:)./sum(w(t,:));                                                         % Normalise the weights.
  
  %% RESAMPLING USING IMPORTANCE FUNCTION
    outIndex = residualR(1:N,w(t,:)');                                                  % Residual resampling.
  
  xparticle_pfekf(:,t,:) = xparticlePred_pfekf(:,t,outIndex);                           % Keep particles with resampled indices.
  for i=1:N,
    Pparticle_pfekf{i} = PparticlePred_pfekf{outIndex(i)};  
  end;
end;   
%% EXTRACTING AND PLOTTING TRACKED VALUES
for t=1:T,
  track_x(t) = mean(yPred_pfekf(1,t,:));
  track_y(t) = mean(yPred_pfekf(2,t,:));  
end; 
figure(2)
plot(track_x(2:T),'r')
ylabel('x-coordinate')
hold on
plot(x(1,:),'g')
legend('Tracked x-ccordinate','Actual x-ccordinate')
figure(3)
plot(track_y(2:T),'r')
ylabel('y-coordinate')
hold on
plot(x(2,:),'g')
legend('Tracked y-ccordinate','Actual y-ccordinate')