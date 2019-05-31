function kappa = coherence( APmx, nCell )

% Ifall nCell ar en vektor gors en korrelation mellan
% celler 0 till nCell(1)-1 och celler nCell(1) till nCell(2)-1

% Doing the autokorrelation
if length( nCell ) == 1
    siz = size( APmx, 2 );
    coMat = zeros( siz );
    for i = 1:siz
        for j = 1:i-1
            %mi = mean( APmx(:,i) );
            %mj = mean( APmx(:,j) );
            %si = sum( (APmx(:,i)-mi).*(APmx(:,i)-mi) );
            %sj = sum( (APmx(:,j)-mj).*(APmx(:,j)-mj) );
            %num = sum( (APmx(:,i)-mi).* (APmx(:,j)-mj) );
            %den = sqrt(si*sj);
            num = sum( APmx(:,i).*APmx(:,j) );
            den = sqrt( sum(APmx(:,i)) * sum(APmx(:,j)) );
            if( num == 0 | den == 0 )
                coMat(i,j) = 0;
            else
                coMat(i,j) = num / den;
            end
            coMat(j,i) = coMat(i,j);
        end
        coMat(i,i) = 1;
    end
    % Averaging
    kappa = 0;
    for i = 1:siz
        for j = 1:i
            kappa = kappa + coMat(i,j);
        end
    end
    kappa = kappa / (siz*(siz+1)/2);
else
    % Doing the cross correlation
    coMat = zeros( size( nCell(1), nCell(2)-nCell(1) ) );
    for i = 1:nCell(1)
        for j = 1:nCell(2)-nCell(1)
            si = sum( APmx(:,i) );
            sj = sum( APmx(:,nCell(1)+j) );
            mi = mean( APmx(:,i) );
            mj = mean( APmx(:,nCell(1)+j) );
	        num = sum( (APmx(:,i)-mi).* (APmx(:,nCell(1)+j)-mj) );
		    den = sqrt( mi * mj );
		    if( num == 0 | den == 0 )
		        coMat(i,j) = 0;
		    else
		        coMat(i,j) = num / den;
		    end
    	end
    end
    % Averaging
    kappa = 0;
    for i = 1:nCell(1)
        for j = 1:nCell(2)-nCell(1)
            kappa = kappa + coMat(i,j);
        end
    end
    kappa = kappa / (nCell(1)*(nCell(2)-nCell(1)));
end



