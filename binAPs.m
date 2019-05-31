function APmx = binAPs(APs,tau,t1,t2,nCell)

bins = ceil( ( t2 - t1 ) / tau );
APmx = zeros(bins,nCell);
x = APs(:,1);
y = APs(:,2);
siz = length(x);
for i = 1:siz
	if x(i) >= t1 & x(i) < t2 
		APmx(floor((x(i)-t1)/tau+1),y(i)+1) = APmx(floor((x(i)-t1)/tau+1),y(i)+1) + 1;
	end
end
APmx(bins+1:end,:) = [];

