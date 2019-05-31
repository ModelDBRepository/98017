function checkrNMDA( directory, DoPlot, rows, cols )

% checkrNMDA displays how an iteration of rNMDA evolves
%
% directory: Path to the catalog which
%            contains the simulations to be plotted. Those
%            simulations should be located in directories found in
%            a subdirectory of that path named DATA.
% DoPlot: if DoPlot = 0 or is left out, the figure is plotted but not printed
%         if DoPlot = 1, the figure is both plotted and printed
% rows  : Number of rows on each page
% cols  : Number of columns on each page
%
% Change path names to suit your preferences. Places are marked % CHANGE %

if nargin < 1 | nargin > 4
    disp( 'usage1: checkrNMDA( directory ), ' )
    disp( 'usage2: checkrNMDA( directory, DoPlot ), ' )
    disp( 'usage3: checkrNMDA( directory, DoPlot ), ' )
    disp( 'usage4: checkrNMDA( directory, DoPlot, rows, cols ), ' )
    disp( 'Do "type checkrNMDA" for further information' )
    return
end
if nargin < 4
    rows = 4;
    cols = 3;
end
if nargin < 2
    DoPlot = 0;
end


thisdir = pwd;
d = getTreeSpec( directory )
d = d(2:end);

% plot the data in the catalogs in d
len = length( d );

% How many pages?
pages = floor( (len-1) / (rows*cols) ) + 1;
for i = 1: pages
    figure ( i )
    clf
    set( i,'Position', [ 550 * mod( ( i ), 2 ) 100 720 852 ] );
end

hand = zeros( 1, len ); % subplot handles
height = 0.93 / rows;
width = 0.90 / cols;
for ind = 1 : len
    page = floor( (ind-1) / ( rows*cols ) ) + 1;
    subp = ind - ( page-1 ) * rows * cols; 
    figure( page )
    if subp == 1
        clf
    end
    bottom = 1 - height * ( floor( (subp-1) / cols ) + 1 );
    left = 0.1 + mod( (subp-1), cols ) * width;
    hand( ind ) = subplot( 'Position', [ left bottom width-0.05 height-0.07 ]);
end

perc = 0:100;
for ind = 1 : len
    tit = sprintf( 'Iteration nr %d', ind );
    filename = strcat( d(ind).name, '/', 'rAMPA' );
    load( filename )
    filename = strcat( d(ind).name, '/', 'rNMDA' );
    load( filename )
    filename = strcat( d(ind).name, '/', 'relNMDAdoneSoFar' );
    load( filename )
  
    page = floor( (ind-1) / ( rows*cols ) ) + 1;
    subp = ind - ( page-1 ) * rows * cols;
    figure( page ) 
    subplot( hand( ind ) )
    
    % plot
    title( tit, 'Fontsize', 14 )
    hold on
    [ legh a b ] = plotyy( perc, rAMPA, perc, rNMDA );
    hold on
    c = plot( perc, relNMDAdoneSoFar(1:101)-0.5, 'r*' );
    h = [a;b;c];
    if subp == 1
        legend( h, 'gAMPA', 'gNMDA', 'Measured' );
    end
    set( gca, 'box', 'on' )
    
    set( gca, 'YLim', [0 1] )
    set( gca, 'FontSize', 14 )
    if mod((subp-1), cols )
        set( gca, 'YTickLabel', [] )
    end

end


% set titles and labels
for i = 1 : pages
    figure( i )
    thand( i ) = axes('Position',[0 0 1 1],'Visible','off');
    set(gcf,'CurrentAxes',thand( i ) )
    xl(i) = text(.5,.03,'rNMDA (%)','FontSize',14, 'HorizontalAlignment', 'center');
    yl(i) = text(.03, .5, 'Cell #', 'Rotation', 90, 'FontSize', 14 , 'VerticalAlignment', 'middle');
    tit(i) = text(.5, .99, sprintf( 'Relative NMDA and AMPA current contributions to total current' ), 'FontSize', 14, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

% adjust paper size for printing
if DoPlot
    for i = 1 : page
        figure( i )    
        set( gcf, 'PaperUnits', 'centimeters' );
        set( gcf, 'PaperType', 'A4');
        papersize = get( gcf, 'PaperSize' );
        left = 0.02;
        bottom = 0.02;
        width = papersize( 1 ) - 2 * left;
        height = papersize( 2 ) - 2 * left;
        myfiguresize = [ left, bottom, width, height ];
        set( gcf, 'PaperPosition', myfiguresize );
        print
    end
end
