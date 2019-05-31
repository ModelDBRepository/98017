function kappaPlot( APs, netborder, t1, t2, autocross, deriv )

% autocross = 0 if autocoherence
% autocross = 1 if cross coherence
% autocross = 2 if cross correlation
% deriv = 0: coherence function plotted. 
% deriv = 1: its derivative plotted.

n = length( netborder ) - 1;
if autocross
    n1 = n;
    n2 = n;
else
    n1 = ceil( sqrt( n ) );
    n2 = ceil( n/n1 );
end
x = APs(:,1);
y = APs(:,2);
str = 'IE';

nCell = 0;
if autocross
    xx = [];
    yy = [];
    for i = 1:n
        if mod( i, 2 ) == 1  % I cells
    	    ind = find( ( y >= netborder(i) ) & ( y < netborder(i+1) ) );
	        xx = [ xx ; x(ind) ];
	        yy = [ yy ; y(ind)-netborder(i)+nCell(end) ];
            nCell = [ nCell nCell(end)+netborder(i+1)-netborder(i) ];
        else
            histx = [ netborder(i):netborder(i+1)-1 ];
            histy = histc(y, histx);
            f0 = mean(histy);
            hind = histx( find( histy>1.5*f0 ) ); % Find bump cells
            nCell = [ nCell nCell(end)+length( hind ) ];
            for k = 1:length(hind)
                ind = find( y == hind(k) );
                xx = [ xx ; x(ind) ];
                yy = [ yy ; (nCell(end-1)+k-1)*ones(length(ind),1) ];
            end
        end
    end
    if autocross == 2
        [y x] = crcorr3( [xx yy], t1, t2, nCell );
        for i = 1:length(netborder)-1
            for j = 1:i
                subplot( n1, n2, (i-1)*n1+j )
                plot( x, y(:,(i-1)*n1+j) )
                %fi = fas( x, y );
                title( sprintf( '%s-%d to %s-%d', str(mod(i-1,2)+1), ceil(i/2), str(mod(j-1,2)+1), ceil(j/2) ) )
                set( gca, 'YLim', [0 20] )
                set( gca, 'XLim', [min(x) max(x)] )
                grid on
                if i<length(netborder)-1
                    set( gca, 'XTickLabel', [] )
                else
                    xlabel( 'time(ms)' )
                end
            end
        end
    elseif autocross == 1
        for i = 2:length(netborder)-1
            for j = 1:i-1
                ind = find( ( ( yy >= nCell(j) ) & ( yy < nCell(j+1) ) ) );
                xxx = xx(ind);
                yyy = yy(ind) - nCell(j);
                ind = find( ( ( yy >= nCell(i) ) & ( yy < nCell(i+1) ) ) );
                xxx = [ xxx ; xx(ind) ];
                yyy = [ yyy ; yy(ind)-nCell(i)+nCell(j+1)-nCell(j) ];
            	[xxxx yyyy] = kRange( [xxx yyy], t1, t2, [ nCell(j+1)-nCell(j) nCell(i+1)-nCell(i)+nCell(j+1)-nCell(j) ] );
                if deriv
                    xxxx = xxxx(1:end-1) + (xxxx(2)-xxxx(1));
                    yyyy = diff( yyyy );
                end
              	subplot( n1-1, n2-1, (i-2)*(n1-1)+j )
   	            plot( x, y )
                title( sprintf( '%s-%d to %s-%d', str(mod(i-1,2)+1), ceil(i/2), str(mod(j-1,2)+1), ceil(j/2) ) )
                set( gca, 'XLim', [0 tau(end)] )
                if i == length(nCell)-1
                    xlabel( 'time (ms)' )
                end
            end
        end
    end
else
    for i = 1:n
        if mod( i, 2 ) == 1  % I cells
    	    ind = find( ( y >= netborder(i) ) & ( y < netborder(i+1) ) );
	        xx = x(ind);
	        yy = y(ind)-netborder(i);
            nCell = netborder(i+1)-netborder(i);
        else
            ind = [];
            xx = [];
            yy = [];
            histx = [ netborder(i):netborder(i+1)-1 ];
            histy = histc(y, histx);
            f0 = mean(histy);
            hind = histx( find( histy>1.5*f0 ) ); % Find bump cells
            nCell = length( hind )
            for j = 1:nCell
                ind = find( y == hind(j) );
                xx = [ xx ; x(ind) ];
                yy = [ yy ; (j-1)*ones(length(ind),1) ];
            end
        end
    
    	[kappa tau] = kRange( [xx yy], t1, t2, nCell );
	    subplot( n1, n2, i )
        if deriv
            dkappa = diff( kappa );
            dtau = tau(1:end-1) + (tau(2)-tau(1));
            plot( dtau, dkappa )
        else
        	plot( tau, kappa )
            set( gca, 'YLim', [0 1] )
        end
        title( sprintf( '%s-cells module %d', str(mod(i-1,2)+1), ceil(i/2) ) )
        set( gca, 'XLim', [0 tau(end)] )
        xlabel( 'time (ms)' )
    end
end