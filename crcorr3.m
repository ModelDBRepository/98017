function [xc, x ] = crcorr( APs, t1, t2, netborder )

% Denna funktion gor ett populations-korskorrelogram 
ind = find( (APs(:,1) >= t1) & (APs(:,1) < t2) );
APs = APs(ind,:);
dn = diff( netborder );
f0 = size( APs, 1 ) / ((t2-t1)/1000*netborder(end));
tau = 10/f0;
Tend = 100; % dvs 100 * tau ms bort
for i = 1:length(netborder)-1
    APs( find( (APs(:,2) >= netborder(i)) & (APs(:,2) < netborder(i+1)) ), 2 ) = i-1;
end
APmx = binAPs( APs, tau, t1, t2, length(netborder)-1 );
for i = 1:length(netborder)-1
    APmx(:,i) = APmx(:,i) / dn(i);
end
[xc, x] = xcorr( APmx, Tend );
x = x * tau;