close all
clearvars
sim  = 6000;
burnin = 3000;
N = 7;
M0 = mmread('../Data/M0.txt');
M1 = mmread('../Data/M1.txt');
M2 = mmread('../Data/M2.txt');
%
% creating the basic spatial matrices
%

obj = struct;
obj.M0 = M0;
obj.M1 = M1;
obj.M2 = M2;
obj.X = cell(N, 1);
obj.Y = cell(N, 1);

for i = 1:N
    obj.X{i} = load(['../data/ERA_',num2str(i),'.dat']);
    obj.Y{i} = load(['../data/Yts_',num2str(i),'.dat']);
end
options = optimset('Display','iter','MaxFunEvals',10000);
f = @(x) likelihood(x, obj);
theta0 = zeros(6,1);
%   setting up AMCMC for theta parameters
%%
L = cell(2,1);
L{2} = zeros(6,1);
L{1} = @(x) -likelihood(x, obj)  + sum(exp(x([1:2,4:5]))); %log correction
[AMCMC_obj]= AMCMC_MH_prec_init(L, 0.1, 0, 1);
AMCMC_obj{5}.burnin = 100;
AMCMC_obj{5}.batch = 5;
theta_sim = zeros(sim, 6);
for i=1:sim
    fprintf('i = %d\n', i);
   [AMCMC_obj]      = AMCMC_MH_prec_sample(AMCMC_obj);
   [AMCMC_obj]      = AMCMC_MH_RR(AMCMC_obj);
    theta_sim(i,:)  = AMCMC_obj{5}.X;
end
figure()
subplot(231)
plot(theta_sim(:,1))
title('tau')
subplot(232)
plot(theta_sim(:,2))
title('kappa')
subplot(233)
plot(theta_sim(:,3))
title('mu')
subplot(234)
plot(theta_sim(:,4))
title('tau_beta')
subplot(235)
plot(theta_sim(:,5))
title('kappa_beta')
subplot(236)
plot(theta_sim(:,6))
title('mu_beta')



predicit = cell(N,1);
 for j = 1:1
    obj_crossval = obj;
    obj_crossval.X = cell(N-1, 1);
    obj_crossval.Y = cell(N-1, 1);
    set_c = setdiff(1:N,j);
    for i = 1:(N-1)
        obj_crossval.X{i} = obj.X{set_c(i)};
        obj_crossval.Y{i} = obj.Y{set_c(i)};
    end
      predicit{j, 1} = zeros(length(obj_crossval.Y{j}),1);
      L{1} = @(x) -likelihood(x, obj_crossval) + sum(exp(x([1:2,4:5]))); %log correction
      [AMCMC_obj]= AMCMC_MH_prec_init(L, 0.1, 0, 1);
       AMCMC_obj{5}.burnin = 100;
       AMCMC_obj{5}.batch = 5;
       AMCMC_obj{5}.X = theta_sim(end,:)';
 for i=1:sim
    fprintf('j, i = %d ,%d\n', j, i);
    [AMCMC_obj]      = AMCMC_MH_prec_sample(AMCMC_obj);
    [AMCMC_obj]      = AMCMC_MH_RR(AMCMC_obj);
    theta_sim(i,:)  = AMCMC_obj{5}.X;
    
    tau        = exp(AMCMC_obj{5}.X(1));
    kappa      = exp(AMCMC_obj{5}.X(2));
    beta  = posterior_beta_sample(AMCMC_obj{5}.X, obj_crossval);
    Q_eps  = tau * (obj.M0 + kappa*obj.M1 + kappa^2*obj.M2 );	
    if i > burnin
        predicit{j, 1} =  predicit{j, 1} + exp(beta.*obj.X{j} + AMCMC_obj{5}.X(3) + chol(Q_eps)\randn(length(Q_eps),1));
    end
 end
end
mean((exp(obj.Y{1}) -predicit{1,1}/(sim- burnin)).^2)