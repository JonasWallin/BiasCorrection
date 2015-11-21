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
  Sigma_eps = materncov(obj.D,exp(p(n_beta+1)), exp(p(n_beta+2)),obj.nu);
else
  n_error = 1;
  Sigma_eps = exp(p(n_beta+1))*speye(n);
end
%Sigma_beta = Sigma_beta(obj.reo,obj.reo);
mu = p(n_beta+n_error+1);

%if obj.smooth_X == 1
%  n_x = 2;
%  K_x = materncov(obj.D,1,exp(p(n_beta+n_error+2)),obj.nu);
%  K_x = bsxfun(@rdivide,K_x,sum(K_x,2));
%else
%  K_x = speye(n);
%end
  if obj.smooth_X == 1
  Sigma_X = materncov(obj.D,exp(p(n_beta+n_error+2)),exp(p(n_beta+n_error+3)),obj.nu);
  sigma_xe = exp(p(n_beta+n_error+4));
  mux = p(n_beta+n_error+5);

  Q_XE = speye(n)/sigma_xe;
  Q_Xpost = inv(Sigma_X) + Q_XE;
  bx = 0;
  R = chol(Q_Xpost);
  Xs = cell(n_rep,1);
  for ii= 1:n_rep
    Xs{ii} = mux + Q_Xpost\(obj.X{ii}-mux)/sigma_xe;
  end
else
  Xs = obj.X;
end




Q_E = inv(Sigma_eps);
Q_post = inv(Sigma_beta);
b = 0;

for ii=1:n_rep
  X = diag(Xs{ii});
  Q_post = Q_post + X'*Q_E*X;
  b = b + X*Q_E*(obj.Y{ii}-mu);
end
beta  = Q_post\b;