function [DTF,H]=DirTransFunc(AR,f,dt)

% AR{j}, j=1,...,p, where p is the order of the AR model, is the 
% stack of AR matrices in the AR model.
% Or AR is a d-by-p*d matrix, where d is the dimensionality of 
% the time series.
% f is the frequency of interest
% dt is the sampling time bin size in sec

z = exp(-2*pi*i*dt*f);

if iscell(AR)
    H = zeros(size(AR{1}));
    for m = 1:length(AR)
        H = H + AR{m}*z^(-m);
    end
else
    d = size(AR,1);
    H = zeros(d);
    for m = 1:size(AR,2)/d
        H = H + AR(:,d*(m-1)+1:d*m)*z^(-m);
    end
end

for m = 1:size(H,1)
    for n = 1:size(H,2)
        DTF(m,n) = abs(H(m,n))^2 / (H(m,:)*H(m,:)');
    end
end

