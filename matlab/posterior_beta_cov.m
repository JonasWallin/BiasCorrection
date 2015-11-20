function beta  = posterior_beta_cov(p, obj)

n = length(obj.Y{1});
n_rep = length(obj.Y);

if obj.smooth_beta == 1
  n_beta = 2;
  Sigma_beta = materncov(obj.D, exp(p(1)),exp(p(2)),obj.nu);
else
  n_beta = 1;
  Sigma_beta = exp(p(1))*speye(n);
end
if obj.smooth_error == 1
  n_error = 2;
  Sigma_eps = materncov(obj.D,exp(p(n_beta+2)), exp(p(n_beta+1)),obj.nu);
else
  n_error = 1;
  Sigma_eps = exp(p(n_beta+1))*speye(n);
end
mu = p(n_beta+n_error+1);

Q_E = inv(Sigma_eps);
Q_post = inv(Sigma_beta);
b = 0;

for ii=1:n_rep
  X = diag(obj.X{ii});
  Q_post = Q_post + X'*Q_E*X;
  b = b + X*Q_E*(obj.Y{ii}-mu);
end
beta  = Q_post\b;