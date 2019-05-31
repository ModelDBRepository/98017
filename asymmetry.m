function W = asymmetry(x,NMDA11,AMPA11,NMDA12,AMPA12)

% input
% x = n-by-1 vector with sorted percentage asymmetry values.
%     positive values indicate stronger connections
%     from module 1 to 2 than from 2 to 1.
% NMDA11 = internal NMDA strength for the symmetrical case
% AMPA11 = internal AMPA strength for the symmetrical case
% NMDA12 = inter-regional NMDA strength for the symmetrical case
% AMPA12 = inter-regional AMPA strength for the symmetrical case
% 
% output
% W = 8-by-n matrix where each column contains the 8 e-to-e
%     connection strength values for that asymmetry percentage:
%     [n11 a11 n22 a22 n12 a12 n21 a21]',
%     where n and a indicate NMDA and AMPA synapses, respectively,
%     and 12 indicates a connection from module 1 to module 2.

x = x(:)'; % Make sure x is horizontal
n = length(x);
if sum(x==0)>0
    x = [-fliplr(x(2:end)) 0 x(2:end)];
else
    x = [-fliplr(x(1:end)) x(1:end)];
end

% Inter-regional NMDA connection
yn=NMDA12*2*(0.5+0.5*x);
yn = round(1000*yn)/1000;
s_yn = yn + fliplr(yn);
ind7 = ceil(length(yn)/2):-1:1;
ind5 = floor(length(yn)/2+1):length(yn);
ind_more = find(s_yn(ind7) - max(s_yn) < 1e-10);
if length(ind_more) ~= length(ind7)
    yn(ind_more) = yn(ind_more)-0.001;
end

% Inter-regional AMPA connection
ya=AMPA12*2*(0.5+0.5*x);
ya = round(1000*ya)/1000;
s_ya = ya + fliplr(ya);
ind_more = find(s_ya(ind7) - max(s_ya) < 1e-10);
if length(ind_more) ~= length(ind7)
    ya(ind_more) = ya(ind_more)-0.001;
end

% Put it all together
W = zeros(8,n);
W(1,:) = NMDA11 + (yn(ind5)- NMDA12);
W(2,:) = AMPA11 + (ya(ind5)- AMPA12);
W(3,:) = NMDA11 - (yn(ind5)- NMDA12);
W(4,:) = AMPA11 - (ya(ind5)- AMPA12);
W(5,:) = yn(ind5);
W(6,:) = ya(ind5);
W(7,:) = yn(ind7);
W(8,:) = ya(ind7);
