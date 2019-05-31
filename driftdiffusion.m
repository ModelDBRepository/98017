function [ dr, popv ] = driftdiffusion( x, y, ang, binw, time )

% drift calculates the drift of the population vector. Drift 
% is here defined as the standard deviation of diffusion process
% of the population vector of the bump. This function uses popvec.m,
% which needs to be in the same catalog.
%
% x    : vector of time points of action potentials
% y    : vector of cell indices of action potentials
% ang  : Matrix. left column: cell index, right colums: angle or
%        whichever value one wants to measure
% binw : time bin
% time : two-vector with start and stop times

tid = time(1):binw:time(2);
popv = popvec( x, y, ang, tid );
if popv == -1
    dr = -1;
    return
end
ddtpopv = ( diff( popv ) ./ (diff( tid(1:end-1) )/1000) );
wddt = ddtpopv.^2.*diff( tid(1:end-1) ) / (tid(end-1)-tid(1));
dr = sqrt( sum( wddt ));
%figure(6)
%subplot( 3,1,1 )
%plot( tid( 1:end-1), popv )
%subplot( 3,1,2 )
%plot( tid( 1:end-2 ), ddtpopv )
%subplot( 3,1,3 )
plot( tid( 1:end-2 ), wddt )
disp( dr )