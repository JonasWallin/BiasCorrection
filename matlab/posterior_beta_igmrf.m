function [beta beta_var e_var] = posterior_beta_igmrf(theta, obj)


tau        = exp(theta(1));
mu         = theta(2);
tau_beta   = exp(theta(3));

n = length(obj.X{1});

Q_eps  = tau * speye(n);
Q_beta = tau_beta*obj.Qx;

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

v = R'\b;
beta = obj.A*(R\v);

beta_var = obj.A*diag(Qinv(R));
e_var = ones(size(beta_var))/tau;

