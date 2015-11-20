function beta  = posterior_beta_sample(theta, obj)




tau        = exp(theta(1));
kappa      = exp(theta(2));
mu         = theta(3);
tau_beta   = exp(theta(4));
kappa_beta = exp(theta(5));
mu_b = 0;%theta(6);
Q_eps  = tau * (obj.M0 + kappa*obj.M1 + kappa^2*obj.M2 );	
Q_beta = tau_beta  * (obj.M0 + kappa_beta*obj.M1 + kappa_beta^2*obj.M2 );	

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
end 

[R, p] = chol(Q_hat);

v = R'\b;
beta = R\(v + randn(length(v),1));