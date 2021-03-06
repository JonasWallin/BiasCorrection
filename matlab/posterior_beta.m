function [beta beta_var e_var]  = posterior_beta(theta, obj)




tau        = exp(theta(1));
kappa      = exp(theta(2));
mu         = theta(3);
tau_beta   = exp(theta(4));
kappa_beta = exp(theta(5));
mu_b = theta(6);
Q_eps  = tau *       (obj.M0 + kappa*obj.M1      + kappa^2*obj.M2 );
Q_beta = tau_beta  * (obj.M0 + kappa_beta*obj.M1 + kappa_beta^2*obj.M2 );
%Q_beta = Q_beta(obj.reo,obj.reo);
%Q_beta = Q_beta(obj.reo,obj.reo);
%Q_eps = Q_eps(obj.reo,obj.reo);

mu_b = mu_b  *ones(length(Q_beta), 1);

b = Q_beta * mu_b;
Q_hat = Q_beta;
n = length(obj.X{1});
yQy = 0;
for i=1:length(obj.X)
    Ax    = sparse(1:n, 1:n, obj.X{i});
    AQ = Ax'*Q_eps;
    Q_hat = Q_hat + AQ * Ax;
    Y_mu  = (obj.Y{i} - mu);
    b     = b + AQ * Y_mu;
    yQy = yQy + (Y_mu' * (Q_eps * Y_mu))/2;
end

[R, p] = chol(Q_hat);

v = R'\b;
beta = R\v;

beta_var = diag(Qinv(R));
e_var = diag(Qinv(chol(Q_eps)));
