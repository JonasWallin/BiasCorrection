close all
clearvars
N = 7;
M0 = mmread('../data/M0.txt');
M1 = mmread('../data/M1.txt');
M2 = mmread('../data/M2.txt');
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
    obj.evalX{i} = load(['../data/BCM',num2str(i),'.txt']);
    obj.X{i} = load(['../data/ERA',num2str(i),'.txt']);
    obj.Y{i} = load(['../data/Y',num2str(i),'.txt']);
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

%f_0 = @(x) likelihood_0(x, obj);
%theta0_0 = zeros(4,1);
%[theta_0, fval] = fminsearch(f_0, theta0_0, options);
%[theta_0, fval] = fminsearch(f_0, theta_0, options);
%[theta_0, fval] = fminunc(f_1, theta_0, options);
%beta_0  = posterior_beta_0(theta_0, obj);


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
    %f_0 = @(x) likelihood_0(x, obj_crossval);

    [theta, fval] = fminunc(f, theta, options);
    [theta_1, fval] = fminunc(f_1, theta_1, options);
    %[theta_0, fval] = fminunc(f_0, theta_0, options);

    beta  = posterior_beta(theta, obj_crossval);
    beta_1  = posterior_beta_1(theta_1, obj_crossval);
    beta_0  = posterior_beta_0(theta_0, obj_crossval);

    predicit{j, 1} = beta;
    predicit{j, 2} = exp(obj.Y{j}) - exp(beta.*obj.evalX{j} + theta(3)) ;
    predicit{j, 3} = beta_1;
    predicit{j, 4} = exp(obj.Y{j}) - exp(beta_1.*obj.evalX{j} + theta_1(2)) ;
    %predicit{j, 5} = beta_0;
    %predicit{j, 6} = exp(obj.Y{j}) - exp(beta_0.*obj.evalX{j} + theta_0(2)) ;
 end
 for j=1:N
     fprintf('mse_%d  = %.2f %.2f\n',j,mean(predicit{j, 2}.^2 ),mean(predicit{j, 4}.^2 ));
 end

  load ../data/ind.txt
mat = NaN*zeros(58,63);
figure(1)
subplot(231)
mat(ind) = exp(obj.Y{1});
imagesc(mat'); axis xy;colorbar
subplot(232)
mat(ind) = exp(beta.*obj.evalX{1} + theta(3));
imagesc(mat'); axis xy ;colorbar
subplot(233)
mat(ind) = predicit{1,2};
imagesc(mat'); axis xy ;colorbar

subplot(234)
mat(ind) = exp(obj.evalX{1});
imagesc(mat'); axis xy;colorbar
subplot(235)
mat(ind) = beta;
imagesc(mat'); axis xy ;colorbar
subplot(236)

