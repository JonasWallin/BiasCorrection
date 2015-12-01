function lik = likelihood_0(theta, obj)


tau        = exp(theta(1));
kappa      = exp(theta(2));
mu         = theta(3);
tau_beta   = exp(theta(4));

n = length(obj.X{1});

Q_eps  = tau * (obj.M0 + kappa*obj.M1 + kappa^2*obj.M2 );
Q_beta = tau_beta  * speye(n);

[R_eps, p] = chol(Q_eps);
if p == 1;    lik = inf;    return; end

b = 0;
Q_hat= Q_beta;

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
if p == 1;    lik = inf;    return; end
v = R'\b;
lik = length(obj.X) * sum(log(diag(R_eps))) ...
    + n*theta(4)/2 - sum(log(diag(R))) ...
    + v'*v/2 - yQy;
%lik = lik - kappa - tau - tau_beta - kappa_beta; %priors
lik = lik - (mu-3)^2/(2*0.05);
lik = -lik;
