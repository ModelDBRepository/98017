function [ bmp, s2, s2S ] = isbump( x, y, time, ind, window, Qt, s2S )

% isbump finds out whether there is a bump in a time interval.
% The criterion for having a bump is that the standard deviation
% surpasses an empirically found (but arbitrary) threshold of 1.5.
% Variance of activity is measured as the average variance in time
% windows of size window to eliminate effects of bump drift.
%
% x      : Vector with action potential time points
% y      : Vector with action potential angles, or whatever
% time   : 2-vector with beginning and end of measure
% ind    : vectors with cell indices of interest
% window : Time window in which variance is measured
% Qt     : Time window for spontaneous activity is 500 ms <= t < Qt

s2 = 0;

% Find only action potentials of interest
yind = [];
for i = 1:length( ind )
  yind = [ yind find( y == ind(i) ) ];
end
x = x(ind);
y = y(ind);

% Calculate standard deviation of spontaneous activity if none exists
if nargin < 7
    s2S = 0;
    timS = 500:window:Qt; 
    cnt = 0;
    for i = 1:length( timS ) -1
        indt = find( x>=timS(i) & x<timS(i+1) );
        if isempty( indt )
            cnt = cnt + 1;
        else
            xtmp = x( indt );
            ytmp = y( indt );
            rate = histc( ytmp, ind ) / (diff(timS(i:i+1))/1000); %Kolla om detta ar statistiskt riktigt
            tmp =  var( rate ) / mean( rate );
            s2S = s2S + tmp;
        end
    end
    cnt = length( timS )-1-cnt;
    if cnt > 0
        s2S = s2S / cnt;
    else
        s2S = 0
    end
end

time = time(1):window:time(2);
cnt = 0;
for i = 1:length( time ) -1
    indt = find( x>=time(i) & x<time(i+1) );
    if isempty( indt )
        cnt = cnt + 1;
    else
        ytmp = y(indt);
        xtmp = x(indt);
        rate = histc( ytmp, ind ) / (diff(time(i:i+1))/1000);
        tmp =  var( rate ) / mean( rate );
        s2 = s2 + tmp;
    end
end
cnt = length( time )-1-cnt;
if cnt > 0
    s2 = s2 / cnt;
else
    s2 = 0;
end
bmp = s2>1.5*s2S;  % Note, 1.5 is arbitrary border
