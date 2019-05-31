% This script plots cell populations for showConnMulti.m

angle = 0:0.1:2*pi+0.1;
one = ones( size( angle ) );
for iiiii = 1:length( pop )
    hold on
    [ x y ] = pol2cart( angle, pop(iiiii).radius*one );
    x = (x+one*pop(iiiii).x);
    y = (y+one*pop(iiiii).y);
    plot( x, y, 'LineWidth', 2 )
    poptx(iiiii) = text( pop(iiiii).x, pop(iiiii).y, pop(iiiii).label, 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle' );
end
for iiiii = 1:length( pop )/2
    poptx(iiiii) = text( pop(2*iiiii).x, pop(2*iiiii-1).y/3, ['Area ' int2str(iiiii)], 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle' );
end
set( gca, 'YLim', [0 1] )
set( gca, 'XLim', [0 1] )
    