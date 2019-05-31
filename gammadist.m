function y = gammadist(x,mu,s2)
%function y = gammadist(x,p,a)

a = s2/mu;
p = mu/a;

dx = x(2)-x(1);
F = gammainc( (x+dx/2)/a, p );
if size(x,1) > 1
    y = [ 0 ; diff(F)/dx ];
else
    y = [ 0 diff(F)/dx ];
end
