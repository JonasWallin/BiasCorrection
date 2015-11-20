function lik = likelihood_1(theta, obj)


tau        = exp(theta(1));
mu         = theta(2);
tau_beta   = exp(theta(3));
kappa_beta = exp(theta(4));
mu_b = 0;%theta(6);
Q_eps  = tau * speye(length(obj.M0));	
Q_beta = tau_beta  * (obj.M0 + kappa_beta*obj.M1 + kappa_beta^2*obj.M2 );	

[R_eps, p] = chol(Q_eps);
if p == 1;    lik = inf;    return; end
[R_beta, p] = chol(Q_beta);
if p == 1;    lik = inf;    return; end

mu_b = mu_b  *ones(length(Q_beta), 1);

b = Q_beta * mu_b;
Q_hat= Q_beta;  
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
if p == 1;    lik = inf;    return; end
v = R'\b;
lik = length(obj.X) * sum(log(diag(R_eps))) ...
    + sum(log(diag(R_beta))) - sum(log(diag(R))) ...
    + v'*v/2 - yQy -mu_b'*Q_beta*mu_b/2;
lik = -lik;