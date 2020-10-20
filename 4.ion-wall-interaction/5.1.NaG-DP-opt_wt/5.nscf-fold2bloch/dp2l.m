% -------------------------------------------------------------------------
function RES = dp2l(X0,X1,X2) % distance from point {X0} to line {X1}-{X2}
% see http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
denom = X2 - X1;
denomabs = sqrt(dot(denom,denom));
if denomabs < eps
    display(X1); display(X2);
    error('X1 = X2');
end;
numer = cross( X0-X1 , X0-X2 );
numerabs = sqrt(dot(numer,numer));
RES = numerabs/denomabs;
end
% -------------------------------------------------------------------------
