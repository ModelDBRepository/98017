function [xc, x ] = crcorr( APs, t1, t2, netborder )

% Denna funktion gor ett populations-korskorrelogram 
dn = diff( netborder );
x = APs(:,1);
f0 = length( x ) / ((t2-t1)/1000*netborder(end));
Tend = 500;
APmx = binAPs( APs, 10/f0, t1, t2, netborder(end) );
xc = zeros( 2*Tend+1, length(netborder)*(length(netborder)-1)/2 );
for i = 1:length( netborder ) - 1 
  for j = 1:i
    ind = prod(1:(i-1))+j;
    if i == j
      a = sum( APmx(:,netborder(i)+1:netborder(i+1))' )';
      xc(:,ind) = xcorr( a/dn(i), Tend );
    else
      for k = netborder(i)+1:netborder(i+1)
	    for m = netborder(j)+1:netborder(j+1)
            a = xcorr( [APmx(:,k) APmx(:,m)], Tend );
            xc(:,ind) = xc(:,ind) + a(:,2); 
	    end
      end
      xc(:,ind) = xc(:,ind) / (dn(i)*dn(j));
    end
  end
end

x = -Tend:10/f0:Tend;