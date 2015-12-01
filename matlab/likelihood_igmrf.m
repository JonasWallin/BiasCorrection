function lik = likelihood_1(theta, obj)


tau        = exp(theta(1));
mu         = theta(2);
tau_beta   = exp(theta(3));

n = length(obj.X{1});

Q_eps  = tau * speye(n);
Q_beta = tau_beta*obj.Qx;

[R_beta, p] = chol(Q_beta);
if p == 1; lik = inf;    return; end

Q_hat= Q_beta;
b = 0;
yQy = 0;

for i=1:length(obj.X)
  Ax    = sparse(1:n, 1:n, obj.X{i})*obj.A;
  AQ = Ax'*Q_eps;
  Q_hat = Q_hat + AQ * Ax;
  Y_mu  = (obj.Y{i} - mu);
  b     = b + AQ * Y_mu;
  yQy = yQy + (Y_mu' * (Q_eps * Y_mu))/2;
end

[R, p] = chol(Q_hat);
if p == 1; lik = inf;    return; end
v = R'\b;
lik = length(obj.X)*n*theta(1)/2 ...
    + length(obj.Qx-10)*theta(3)/2 - sum(log(diag(R))) ...
    + v'*v/2 - yQy;
%lik = lik - (theta(1)-10)^2/(2*0.0001);
lik = -lik;