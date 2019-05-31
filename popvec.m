function [ popv, cut ] = popvec( x, y, ang, time )

% popvec calculates the population vector. To calculate the
% population vector in a ring, one first needs to find
% the angle section with the lowest variance. Then one calculates the
% mean position of the activity.
%
% x   : Action position times
% y   : Action position cell index
% ang : Matrix. left column: cell index, right colums: angle or
%       whichever value one wants to measure
% time: Vector of time intervals where the population vector is to
%       be calculated.
% popv: Population vector in radians
% cut: The angle of the cut of the ring that produces the smallest variance

popv = zeros( 1, length( time )-1 );
% for every time step, do
one = ones( size( y ) );
for i = 1:length( time )-1
  
    % only action potentials within the right time window
    ytmp = y( find( x>=time(i) & x < time(i+1) ) );
    
    % only action potentials from the right cells
    for j = 1:size( ang, 1 )
        cnt = length( find( ytmp == ang( j, 1 ) ) );
        dist(j) = cnt;
    end
    if sum( dist ) == 0 % no rate
        disp( 'Error: no rate' )
        popv = -1;
        cut = 0;
        return
    end
    dist = dist / sum( dist );
    dist = dist(:);
    
    % sort vector so that vectors can be rotated properly
    [ angl, ind ] = sort( ang(:,2) );
    dist = dist( ind );
    
    % find minimum variance
    s2 = sum( dist.*angl.*angl ) - sum( dist.*angl )*sum( dist.*angl );
    popv(i) = sum(angl.*dist);
    len = length( angl );
    cut = 1;
    for j = 1: len-1
        angl = [ angl(end) ; angl(1:end-1) ];
        s22 = sum( dist.*angl.*angl ) - sum( dist.*angl )*sum( dist.*angl );
    	% If variance smallest so far, calculate mean activity position.
        if s22 < s2
            s2 = s22;
            popv(i) = mod( (sum( angl.*dist )+j), len );
            cut = j+1;
        end
    end
    angl = [ angl(end) ; angl(1:end-1) ];
end
cut = angl( cut );