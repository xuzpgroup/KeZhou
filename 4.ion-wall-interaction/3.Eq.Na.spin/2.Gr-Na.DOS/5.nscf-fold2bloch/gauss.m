% -------------------------------------------------------------------------
function RES = gauss(X,x0,s) % bound. conditions
RES = 1/(s*sqrt(2*pi)) * exp(-((X-x0).^2)/(2*s^2));
end
% -------------------------------------------------------------------------

