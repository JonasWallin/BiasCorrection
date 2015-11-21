close all
clearvars

%%%%%%%%%%%%%%%%%%%%
N = 7;
smooth_error = 1;
smooth_beta = 0;
smooth_x = 1;
nu  = 1;
%%%%%%%%%%%%%%%%%%%%

obj = struct;
obj.X = cell(N, 1);
obj.Y = cell(N, 1);
obj.smooth_error = smooth_error;
obj.smooth_beta = smooth_beta;
obj.smooth_X = 0;%smooth_x;
obj.nu = nu;
%load ../data/reo.txt
%ireo = 1:777;
%ireo(reo) = 1:777
%obj.reo = ireo;
load ../data/loc_x.txt;
load ../data/loc_y.txt;
loc = [loc_x loc_y];
obj.D = distancematrix(loc);
for i = 1:N
  obj.Y{i} = load(['../data/Y',num2str(i),'.txt']);
  if smooth_x
    obj.evalX{i} = load(['../data/BCM',num2str(i),'.txt']);
    obj.X{i} = load(['../data/ERA',num2str(i),'.txt']);
  else
    obj.evalX{i} = load(['../data/BCMs',num2str(i),'.txt']);
    obj.X{i} = load(['../data/ERAs',num2str(i),'.txt']);
  end
end
n_err = 1;
if smooth_error
  n_err = 2;
end
n_beta  =1;
if smooth_beta
  n_beta = 2;
end

n_x = 0;
if smooth_x
 n_x = 4;
end
options = optimset('Display','iter','MaxFunEvals',10000);
f = @(x) likelihood_cov(x, obj);
theta0 = zeros(n_err + n_beta +1 + n_x,1);
[theta, fval] = fminsearch(f, theta0, options);
[theta, fval] = fminsearch(f, theta, options);
[theta, fval] = fminunc(f, theta, options);
beta  = posterior_beta_cov(theta, obj);

predict = cell(N,4);
options = optimset('MaxFunEvals',10000);
for j = 1%:N
  obj_crossval = obj;
  obj_crossval.X = cell(N-1, 1);
  obj_crossval.Y = cell(N-1, 1);
  set_c = setdiff(1:N,j);
  ym = [];
  for i = 1:(N-1)
    obj_crossval.X{i} = obj.X{set_c(i)};
    obj_crossval.Y{i} = obj.Y{set_c(i)};
    ym =[ym obj.Y{set_c(i)}];
  end
  f = @(x) likelihood_cov(x, obj_crossval);

  [theta, fval] = fminunc(f, theta, options);
  beta  = posterior_beta_cov(theta, obj_crossval);

  %if obj.smooth_X == 1
   % K_x = materncov(obj.D,1,exp(theta(end)),obj.nu);
   % K_x = bsxfun(@rdivide,K_x,sum(K_x,2));
  %else
    K_x = speye(777);
  %end
  predict{j, 1} = beta;
  predict{j, 2} = exp(obj.Y{j}) - exp(beta.*(K_x*obj.evalX{j}) + theta(end-obj.smooth_X)) ;
  predict{j,3} = exp(obj.Y{j}) - exp(mean(ym,2));
end
for j=1:N
  fprintf('mse_%d  = %.2f %.2f \n',j,mean(predict{j, 3}.^2 ), mean(predict{j, 2}.^2 ));
end

load ../data/ind.txt
mat = NaN*zeros(58,63);
figure(1)
subplot(231)
mat(ind) = exp(obj.Y{1});
imagesc(mat'); axis xy;colorbar
subplot(232)
mat(ind) = exp(predict{1,1}.*(K_x*obj.evalX{1}) + theta(end-obj.smooth_X));
imagesc(mat'); axis xy ;colorbar
subplot(233)
mat(ind) = predict{1,2};
imagesc(mat'); axis xy ;colorbar

subplot(234)
mat(ind) = exp(K_x*obj.evalX{1});
imagesc(mat'); axis xy;colorbar
subplot(235)
mat(ind) = beta;
imagesc(mat'); axis xy ;colorbar
subplot(236)

