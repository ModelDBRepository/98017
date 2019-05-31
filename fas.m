function fi = fas( x, y, T )

% Givet en korrelationsfunktion med 0-punkten
% i mitten av x-sekvensen finner fas fasskillnaden
% mellan de korrelerade signalerna. Detta gors 
% genom att lokalisera tyngdpunkten i perioden
% narmast 0-punkten. For detta andamal maste denna
% periods granser finnas. Dessa definieras darfor 
% som minimum runt centralpiken.

ind = find( x>=0 & x<T );
x = x(ind);
y = y(ind);

e1 = sum( x.*y );
e2 = sum( x.*x.*y );
s2 = e2-e1*e1;
mu = e1
ind = 1;
for i = 2:length( x )
    x = [ x(end) x(1:end-1) ];
    e1 = sum( x.* y );
    e2 = sum( x.*x.*y );
    s22 = e2-e1*e1;
    if s22 < s2
        mu = e1 + x(i)
        s2 = s22;
    end
    
end
if mu > T/2
    mu = mu - T;
end

fi = mu;