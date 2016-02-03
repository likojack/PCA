function theta=agauss_inv(EW,CW,p)
mu=norm(EW);
sigma=sqrt(CW);
X=mu./sigma;

f=@(theta,X) 1+((2.^(1/2).*X.*exp(-X.^2/2).*cos(theta).^2)./(4.*pi.^(1/2)) - (7186705221432913.*X.*exp(-X.^2/2).*cos(theta).^2) ...
    /36028797018963968 - (erf((2.^(1/2).*X.*cos(theta))./2).*exp((X.^2.*cos(theta).^2)./2 - X.^2./2).*cos(theta))./2 ...
    - (exp((X.^2.*cos(theta).^2)./2 - X.^2./2).*cos(theta))./2);

options = optimset('TolFun',1e-8,'TolX',1e-8);
theta=fminsearch(@(theta) abs(p-f(theta,X)),0,options);
