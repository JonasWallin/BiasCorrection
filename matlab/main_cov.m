close all
clearvars

%%%%%%%%%%%%%%%%%%%%
N = 8;
smooth_error = 0;
smooth_beta = 1;
smooth_x = 1;
nu  = 1;
season = 1;
%%%%%%%%%%%%%%%%%%%%

obj = struct;
obj.X = cell(N, 1);
obj.Y = cell(N, 1);
obj.smooth_error = smooth_error;
obj.smooth_beta = smooth_beta;
obj.nu = nu;
load ../data/reo.txt
ireo = 1:777;
ireo(reo) = 1:777;
obj.reo = ireo;
load ../data/loc_x.txt;
load ../data/loc_y.txt;
load ../data/ind.txt
loc = [loc_x loc_y];
obj.D = distancematrix(loc);
for i = 1:N
obj.Y{i} = load(['../data/season',num2str(season),'/Y',num2str(i),'.txt']);
if smooth_x
obj.X{i} = load(['../data/season',num2str(season),'/ERA',num2str(i),'s.txt']);
obj.ERA{i} = load(['../data/season',num2str(season),'/ERA',num2str(i),'s.txt']);
obj.BCM{i} = load(['../data/season',num2str(season),'/BCM',num2str(i),'s.txt']);
else
obj.X{i} = load(['../data/season',num2str(season),'/ERA',num2str(i),'.txt']);
obj.ERA{i} = load(['../data/season',num2str(season),'/ERA',num2str(i),'.txt']);
obj.BCM{i} = load(['../data/season',num2str(season),'/BCM',num2str(i),'.txt']);
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

options = optimset('Display','iter','MaxFunEvals',10000);
f = @(x) likelihood_cov(x, obj);
theta0 = zeros(n_err + n_beta +1,1);
if smooth_error == 0 && smooth_beta == 1 && smooth_x == 1
  theta0 = [-3.1173; -0.0722; -4.5077; 2.3645; 0];
elseif smooth_error == 0 && smooth_beta == 1 && smooth_x == 0
  theta0 = [-2.5331; -0.3469; -4.6138; 2.0717; 0];
else
  theta0 = zeros(n_err + n_beta +2,1);
end
[theta, fval] = fminsearch(f, theta0, options);
[theta, fval] = fminsearch(f, theta, options);
[theta, fval] = fminunc(f, theta, options);
beta  = posterior_beta_cov(theta, obj);


options = optimset('MaxFunEvals',10000);
predict_ERA = cell(N,4);
mse_ERA = zeros(N,3);
coverage_ERA = zeros(N,3);
predict_BCM = cell(N,4);
mse_BCM = zeros(N,3);
coverage_BCM = zeros(N,3);

for j = 1:N
  obj_crossval = obj;
  obj_crossval.X = cell(N-1, 1);
  obj_crossval.Y = cell(N-1, 1);
  set_c = setdiff(1:N,j);
  ym = [];
  xx = [];
  for i = 1:(N-1)
    obj_crossval.X{i} = obj.X{set_c(i)};
    obj_crossval.Y{i} = obj.Y{set_c(i)};
    ym =[ym obj.Y{set_c(i)}];
    xx = [xx; ones(777,1) diag(obj.X{set_c(i)})];
  end

  %linear regression as comparisson:
  par_est = (xx'*xx)\(xx'*ym(:));
  alpha_est = par_est(1);
  beta_est = par_est(2:778);
  quantp_ERA = alpha_est + beta_est.*obj.ERA{j};
  quantp_BCM = alpha_est + beta_est.*obj.BCM{j};

  e_var = var(ym(:)-xx*par_est);
  par_var = e_var*diag(inv(xx'*xx));

  quant_var = par_var(1) + par_var(2:end).*obj.ERA{j}.^2 + e_var;
  coverage_ERA(j,2)= mean(abs((obj.Y{j} - quantp_ERA)./sqrt(quant_var)) < 1.96);
  quant_var = par_var(1) + par_var(2:end).*obj.BCM{j}.^2 + e_var;
  coverage_BCM(j,2)= mean(abs((obj.Y{j} - quantp_BCM)./sqrt(quant_var)) < 1.96);

  %%%%%%%%%%%%%%%%
  f = @(x) likelihood_cov(x, obj_crossval);

  [theta, fval] = fminunc(f, theta, options);
  [beta beta_var e_var] = posterior_beta_cov(theta, obj_crossval);

  %if obj.smooth_X == 1
   % K_x = materncov(obj.D,1,exp(theta(end)),obj.nu);
   % K_x = bsxfun(@rdivide,K_x,sum(K_x,2));
  %else
  %end
  quantps_ERA = beta.*obj.ERA{j} + theta(end-1);
  quantps_BCM = beta.*obj.BCM{j} + theta(end-1);
  quant_var = beta_var.*obj.ERA{j}.^2 + e_var;
  coverage_ERA(j,3)=mean(abs((obj.Y{j} - quantps_ERA)./sqrt(quant_var)) < 1.96);
  quant_var = beta_var.*obj.BCM{j}.^2 + e_var;
  coverage_BCM(j,3)=mean(abs((obj.Y{j} - quantps_BCM)./sqrt(quant_var)) < 1.96);

  predict_ERA{j, 1} = exp(mean(ym,2));
  predict_ERA{j, 2} = exp(quantp_ERA);
  predict_ERA{j, 3} = exp(quantps_ERA);
  predict_BCM{j, 1} = exp(mean(ym,2));
  predict_BCM{j, 2} = exp(quantp_BCM);
  predict_BCM{j, 3} = exp(quantps_BCM);

  errors_ERA{j, 1} = exp(obj.Y{j}) - exp(mean(ym,2));
  errors_ERA{j, 2} = exp(obj.Y{j}) - exp(quantp_ERA);
  errors_ERA{j, 3} = exp(obj.Y{j}) - exp(quantps_ERA);
  errors_BCM{j, 1} = exp(obj.Y{j}) - exp(mean(ym,2));
  errors_BCM{j, 2} = exp(obj.Y{j}) - exp(quantp_BCM);
  errors_BCM{j, 3} = exp(obj.Y{j}) - exp(quantps_BCM);

  mse_ERA(j,:)=[mean(errors_ERA{j, 1}.^2) mean(errors_ERA{j, 2}.^2) mean(errors_ERA{j, 3}.^2)];
  mse_BCM(j,:)=[mean(errors_BCM{j, 1}.^2) mean(errors_BCM{j, 2}.^2) mean(errors_BCM{j, 3}.^2)];

  fprintf('fold %d : ERA : %.2f %.2f (%.2f) %.2f (%.2f) ',j,mse_ERA(j,1),mse_ERA(j,2),coverage_ERA(j,2),mse_ERA(j,3),coverage_ERA(j,3))

  fprintf('BCM : %.2f %.2f (%.2f) %.2f (%.2f) \n',mse_BCM(j,1),mse_BCM(j,2),coverage_BCM(j,2),mse_BCM(j,3),coverage_BCM(j,3))


  mat = NaN*zeros(58,63);
  figure(j)
  subplot(231); mat(ind) = exp(obj.Y{j});
  imagesc(mat'); axis xy;colorbar
  subplot(232); mat(ind) = predict_BCM{j,2};
  imagesc(mat'); axis xy;colorbar
  subplot(233); mat(ind) = predict_BCM{j,3};
  imagesc(mat'); axis xy;colorbar
  subplot(234); mat(ind) = exp(obj.BCM{j});
  imagesc(mat'); axis xy;colorbar
  subplot(235); mat(ind) = errors_BCM{j,2};
  imagesc(mat'); axis xy;colorbar
  subplot(236); mat(ind) = errors_BCM{1,3};
  imagesc(mat'); axis xy;colorbar
end

for j=1:N
  fprintf('fold %d : ERA : %.2f %.2f (%.2f) %.2f (%.2f) ',j,mse_ERA(j,1),mse_ERA(j,2),coverage_ERA(j,2),mse_ERA(j,3),coverage_ERA(j,3))

  fprintf('BCM : %.2f %.2f (%.2f) %.2f (%.2f) \n',mse_BCM(j,1),mse_BCM(j,2),coverage_BCM(j,2),mse_BCM(j,3),coverage_BCM(j,3))

end


fprintf('mean   : ERA : %.2f %.2f (%.2f) %.2f (%.2f) ',mean(mse_ERA(:,1)),mean(mse_ERA(:,2)),mean(coverage_ERA(:,2)),mean(mse_ERA(:,3)),mean(coverage_ERA(:,3)))

fprintf('BCM : %.2f %.2f (%.2f) %.2f (%.2f) \n',mean(mse_BCM(:,1)),mean(mse_BCM(:,2)),mean(coverage_BCM(:,2)),mean(mse_BCM(:,3)),mean(coverage_BCM(:,3)))

for j = 1:N
  mat = NaN*zeros(58,63);
  figure(j)
  subplot(231); mat(ind) = exp(obj.Y{j});
  imagesc(mat'); axis xy;colorbar
  subplot(232); mat(ind) = predict_BCM{j,2};
  imagesc(mat'); axis xy;colorbar
  subplot(233); mat(ind) = predict_BCM{j,3};
  imagesc(mat'); axis xy;colorbar
  subplot(234); mat(ind) = exp(obj.BCM{j});
  imagesc(mat'); axis xy;colorbar
  subplot(235); mat(ind) = errors_BCM{j,2};
  imagesc(mat'); axis xy;colorbar
  subplot(236); mat(ind) = errors_BCM{1,3};
  imagesc(mat'); axis xy;colorbar
end
