close all
clearvars

%%%%%%%%%%%%%%%%%%%%
N = 8;
smooth_error = 1;
smooth_beta = 0;
smooth_x = 1;
season = 1;

%%%%%%%%%%%%%%%%%%%%
model_eb = 0;
model_b = 0;
model_e = 0;
if smooth_error == 1 && smooth_beta == 1
  model_eb = 1;
elseif smooth_beta == 1 && smooth_error == 0
  model_b = 1;
elseif smooth_beta == 0 && smooth_error == 1;
  model_e = 1;
end

%M0 = mmread('../data/M0.txt');
%M1 = mmread('../data/M1.txt');
%M2 = mmread('../data/M2.txt');

load ../data/reo.txt
ireo = 1:777;
ireo(reo) = 1:777;

load ../data/ind.txt
%use smooth_e and reo for best BCM results:
%M0 = M0(ireo,ireo);
%M1 = M1(ireo,ireo);
%M2 = M2(ireo,ireo);

m  =58;
n = 63;
mat = zeros(58,63);
mat(ind) = 1;
Y = inf*ones(58,63);
Y(ind) = 0;
newS = 58*63;
A = speye(newS);
A = A(ind,:);
%A = A(find(Y(mat(:)==1)<inf),:);

C = speye(m*n);
G1 = -spdiags(ones(n+1,2),[-1 1],n,n);
G2 = -spdiags(ones(m+1,2),[-1 1],m,m);
G = kron(G1,speye(m)) + kron(speye(n),G2);
G = G - spdiags(sum(G,2),0,m*n,m*n);
M0 = G'*G;
M1 = G + G';
M2 = C;
%A = speye(777);
obj = struct;
obj.A = A;
obj.M0 = M0;
obj.M1 = M1;
obj.M2 = M2;
obj.X = cell(N, 1);
obj.Y = cell(N, 1);

obj.smooth_error = smooth_error;
obj.smooth_beta = smooth_beta;
obj.reo = reo;

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

options = optimset('MaxFunEvals',10000);
if model_eb == 1
  f = @(x) likelihood(x, obj);
  theta = zeros(6,1);
elseif model_b == 1
  f = @(x) likelihood_1(x, obj);
  theta = zeros(4,1);
elseif model_e == 1
  f = @(x) likelihood_0(x, obj);
  theta = zeros(4,1);
end

[theta, fval] = fminsearch(f, theta, options);
[theta, fval] = fminsearch(f, theta, options);
[theta, fval] = fminunc(f, theta, options);

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

  if model_eb == 1
    f = @(x) likelihood(x, obj_crossval);
    [theta, fval] = fminunc(f, theta, options);
    mu = theta(3);
    [beta beta_var e_var] = posterior_beta(theta, obj_crossval);
  elseif model_b == 1
    f = @(x) likelihood_1(x, obj_crossval);
    [theta, fval] = fminunc(f, theta, options);
    mu = theta(2);
    [beta beta_var e_var] = posterior_beta_1(theta, obj_crossval);
  elseif model_e == 1
    f = @(x) likelihood_0(x, obj_crossval);
    [theta, fval] = fminunc(f, theta, options);
    mu = theta(3);
    [beta beta_var e_var] = posterior_beta_0(theta, obj_crossval);
  end

  quantps_ERA = beta.*obj.ERA{j} + mu;
  quantps_BCM = beta.*obj.BCM{j} + mu;
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
