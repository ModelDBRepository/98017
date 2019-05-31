function [xc, -Tend:10/f0:Tend] = crcorr( APs, t1, t2, netborder )

% Denna funktion gor ett populations-korskorrelogram 

x = APs(:,1);
f0 = length( x ) / ((t2-t1)/1000*nCell(end));
Tend = 5*1000/f0;
APmx = binAPs( APs, 10/f0, t1, t2, nCell(end) );
xc = zeros( 2*Tend+1, 1 );
    if i == j
      a = sum( APmx(:,nCell(i)+1:nCell(i+1))' )';
      xc(:,ind) = xcorr( a/(nCell(i+1)-nCell(i)), Tend );
    else
      for k = 1:nCell(1)
	for m = nCell(1)+1:nCell(2)
	  xc(:,ind) = xc(:,ind) + xcorr( [APmx(:,k) APmx(:,m)], Tend );
	end
      end
      xc(:,ind) = xc(:,ind) / (nCell(1)*(nCell(2)-nCell(1)));
    end
  end
end

