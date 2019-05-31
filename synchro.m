function synchro(APs, netborder)

x = APs( :, 1 );
y = APs( :, 2 );
if nargin < 2
    netborder = [0 max( y )];
end
nTrain = length( netborder ) - 1;

% group spike trains
dt = 5;
mx = min( x );
Mx = max( x );
t = zeros( ceil((Mx - mx) / dt), nTrain );
for i = 1:length( x )
    ny = 0;
    j=1;
    while y(i) >= netborder(j)
        ny = ny + 1;
        j = j + 1;
    end
    t( ceil( (x(i)-mx)/dt+0.01 ), ny ) = t( ceil( (x(i)-mx)/dt+0.01 ), ny ) + 1;
end

N = 256;
for i = 1:1:nTrain        
    m = ceil( sqrt(nTrain));
    subplot( m, ceil( nTrain/m ), i )
    t(:,i) = t(:,i)-mean(t(:,i));
    y = fft(t(:,i),N);
    Py = y.*conj(y);
    x = 1000/dt*[0:N]/N;
    plot( x(1:N/2), Py(1:N/2) )
    set( gca, 'XLim', [ 3 100 ] )
end
