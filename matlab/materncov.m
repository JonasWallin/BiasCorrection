function Sigma = materncov(D,sigma2,kappa,nu)

Sigma = sigma2*2^(1-nu)*((kappa*D).^nu).*besselk(nu,kappa*D)/gamma(nu);
Sigma(logical(eye(size(Sigma)))) = sigma2;
Sigma(isnan(Sigma)) = sigma2;