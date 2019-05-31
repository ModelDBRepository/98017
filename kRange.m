function [kappa, tau] = kRange(APs, t1, t2, nCell)

% Denna funktion gor en spikkoherenskurva
% Ifall nCell ar ett tal gors en autokoherens
% Ifall nCell ar en vektor gors en korskoherens.
% Da korreleras celler 0 till nCell(1)-1 med celler
% nCell(1)-nCell(2)-1

APs = APs( find( (APs(:,1) >= t1) & (APs(:,1) < t2) ), : );
f0 = size( APs, 1 ) / ((t2-t1)/1000*nCell(end));
dtau = 1000*0.05/f0;
kappa = zeros( 41,1 );
tau = 0:dtau:2*1000/f0;
n = 2;
y = 0;
kappa(1) = 0;
while y < 1 & n<= length(tau)
	disp( ' ' )
	disp( strcat( 'kRange, #', int2str( n ) ) )
	APmx = binAPs( APs, tau(n), t1, t2, nCell(end) );
	y = coherence( APmx, nCell );
	kappa(n) = y;
	n = n + 1;
	disp( strcat( 'kappa = ', num2str(y) ) )
end
kappa( n:end ) = [];
tau( n:end ) = [];
