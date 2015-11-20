close all
clearvars

%%%%%%%%%%%%%%%%%%%%
N = 7;
smooth_error = 0;
smooth_beta = 1;

nu  = 1;
%%%%%%%%%%%%%%%%%%%%

obj = struct;
obj.X = cell(N, 1);
obj.Y = cell(N, 1);
obj.smooth_error = smooth_error;
obj.smooth_beta = smooth_beta;
obj.nu = nu;

load ../data/loc_x.txt;
load ../data/loc_y.txt;
loc = [loc_x loc_y];
obj.D = distancematrix(loc);
for i = 1:N
    obj.evalX{i} = load(['../data/BCM',num2str(i),'.txt']);
    obj.X{i} = load(['../data/ERA',num2str(i),'.txt']);
    obj.Y{i} = load(['../data/Y',num2str(i),'.txt']);
end
n_err = 1;
if smooth_error
  n_err = 2;
end
n_beta  =1;
if smooth_beta
  n_beta = 2;
end
options = optimset('Display','iter','MaxFunEvals',10000);
f = @(x) likelihood_cov(x, obj);
theta0 = zeros(n_err + n_beta +1,1);
[theta, fval] = fminsearch(f, theta0, options);
[theta, fval] = fminsearch(f, theta, options);
[theta, fval] = fminunc(f, theta, options);
beta  = posterior_beta_cov(theta, obj);

predict = cell(N,4);
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
  f = @(x) likelihood_cov(x, obj_crossval);

  [theta, fval] = fminunc(f, theta, options);
  beta  = posterior_beta_cov(theta, obj_crossval);

  predict{j, 1} = beta;
  predict{j, 2} = exp(obj.Y{j}) - exp(beta.*obj.evalX{j} + theta(3)) ;
end
for j=1:N
  fprintf('mse_%d  = %.2f \n',j,mean(predict{j, 2}.^2 ));
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

