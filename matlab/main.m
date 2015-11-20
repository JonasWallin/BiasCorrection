close all
clearvars
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
    obj.evalX{i} = load(['../data/ERA_',num2str(i),'.dat']);
    obj.X{i} = load(['../data/ERA_',num2str(i),'.dat']);
    obj.Y{i} = load(['../data/Yts_',num2str(i),'.dat']);
end
options = optimset('Display','iter','MaxFunEvals',10000);
f = @(x) likelihood(x, obj);
theta0 = zeros(6,1);
[theta, fval] = fminsearch(f, theta0, options);
[theta, fval] = fminsearch(f, theta, options);
[theta, fval] = fminunc(f, theta, options);
beta  = posterior_beta(theta, obj);
f_1 = @(x) likelihood_1(x, obj);

 theta0_1 = zeros(4,1);
 [theta_1, fval] = fminsearch(f_1, theta0_1, options);
 [theta_1, fval] = fminsearch(f_1, theta_1, options);
 [theta_1, fval] = fminunc(f_1, theta_1, options);
 beta_1  = posterior_beta_1(theta_1, obj);

 predicit = cell(N,4);
 options = optimset('MaxFunEvals',10000);
 for j = 1:N
    obj_crossval = obj;
    obj_crossval.X = cell(N-1, 1);
    obj_crossval.Y = cell(N-1, 1);
    set_c = setdiff(1:N,j);
    for i = 1:(N-1)
        obj_crossval.X{i} = obj.X{set_c(i)};
        obj_crossval.Y{i} = obj.Y{set_c(i)};
    end
    f = @(x) likelihood(x, obj_crossval);
    f_1 = @(x) likelihood_1(x, obj_crossval);

    [theta, fval] = fminunc(f, theta, options);
    [theta_1, fval] = fminunc(f_1, theta_1, options);
    
    beta  = posterior_beta(theta, obj_crossval);
    beta_1  = posterior_beta_1(theta_1, obj_crossval);
    
    predicit{j, 1} = beta;
    predicit{j, 2} = exp(obj.Y{j}) - exp(beta.*obj.evalX{j} + theta(3)) ;
    predicit{j, 3} = beta_1;
    predicit{j, 4} = exp(obj.Y{j}) - exp(beta_1.*obj.evalX{j} + theta_1(2)) ;
    
 end
 for j=1:N
     fprintf('mse_%d  = %.2f\n',j,mean(predicit{j, 2}.^2 ));
     
 end