function lik = likelihood_cov(p, obj)

n = length(obj.Y{1});
n_rep = length(obj.Y);

if obj.smooth_beta == 1
  n_beta = 2;
  Sigma_beta = materncov(obj.D, exp(p(1)),exp(p(2)),exp(p(end)));
else
  n_beta = 1;
  Sigma_beta = exp(p(1))*speye(n);
end
%Sigma_beta = Sigma_beta(obj.reo,obj.reo);
if obj.smooth_error == 1
  n_error = 2;
  Sigma_eps = materncov(obj.D,exp(p(n_beta+1)), exp(p(n_beta+2)),exp(p(end)));
else
  n_error = 1;
  Sigma_eps = exp(p(n_beta+1))*speye(n);
end
%Sigma_eps = Sigma_eps(obj.reo,obj.reo);
mu = p(n_beta+n_error+1);

%if obj.smooth_X == 1
%  n_x = 2;
%  K_x = materncov(obj.D,1,exp(p(n_beta+n_error+2)),obj.nu);
%  K_x = bsxfun(@rdivide,K_x,sum(K_x,2));
%else
%  K_x = speye(n);
%end
Xs = obj.X;


%prior on sigma2 for beta:
mu_sigma2 = var(obj.Y{1})/var(obj.X{1});
std_sigma = 0.01;
l = 0;% -(1/(2*std_sigma^2))*(exp(p(1))-mu_sigma2)^2;


for ii= 1:n_rep
  X = diag(Xs{ii});
  Sigma = X*Sigma_beta*X' + Sigma_eps;
  R = chol(Sigma);
  v = R\(obj.Y{ii}-mu);
  l = l - sum(log(diag(R))) - v'*v/2;
end
lik = -l;



