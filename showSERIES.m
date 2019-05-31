function showSERIES( name, DoPlot, rows, cols )

% showSERIES displays a series of pages with simulations all
% belonging to a simulation series. By clicking on a graph, it can
% be enlarged. showSERIES should be in the same catalog as
% showSimulation.m, showBig.m, bumpanalysis.m
%
% name  : is either a path or a file name. 
%         1. When name is a path, it is a path to the catalog which
%            contains the simulations to be plotted. Those
%            simulations should be located in directories found in
%            a subdirectory of that path named DATA.
%         2. When name is a file name, it should be a file name
%            containing the name SERIES. That file contains the
%            simulations which are to be plotted
% DoPlot: if DoPlot = 0 or is left out, the figure is plotted but not printed
%         if DoPlot = 1, the figure is both plotted and printed
% rows  : Number of rows on each page
% cols  : Number of columns on each page
% bold  : 0 = only rastergram, 1 = bold + rastergram
%
% Change path names to suit your preferences. Places are marked % CHANGE %

if nargin < 1 | nargin > 4
    disp( 'usage1: showSERIES( name ), ' )
    disp( 'usage2: showSERIES( name, DoPlot ), ' )
    disp( 'usage3: showSERIES( name, DoPlot ), ' )
    disp( 'usage4: showSERIES( name, DoPlot, rows, cols ), ' )
    disp( 'usage5: showSERIES( name, DoPlot, rows, cols, bold ), ' )
    disp( 'Do "type showSERIES" for further information' )
    return
end
if nargin < 5
    bold = 0;
end
if nargin < 4
    rows = 4;
    cols = 4;
end
if nargin < 2
    DoPlot = 0;
end

thisdir = pwd;

% Is it a UNIX/Linux system or a Windows system
if( ispc )
    pathdelimiter = '\';
    rootsymbol = 'C:\';                 % CHANGE %
    if isempty( strmatch( rootsymbol, name ) )
        name = strcat( thisdir, pathdelimiter, name );
    end
    home = 'C:\Documents and Settings\freedin.admin\Desktop\Neuron\STANDARDFILER';  % CHANGE %
    fs = 10; %Fontsize
else
    pathdelimiter = '/';
    % Is name a relative path or file name?
    if ~( name(1) == '~' | name(1) == pathdelimiter ) 
        name = strcat( thisdir, pathdelimiter, name );
    end
    home = '/home/freedin/Neuron/STANDARDFILER/';  % CHANGE %
    rootsymbol = '/';                              % CHANGE %
    fs = 14; %Fontsize
end
    


if( name(1) == '~' )
    name = strcat( home, name(3:end) );
end

% Put all the directories containing simulations to be 
% displayed in one long list

% alt 1: name is a path
if isempty( findstr( name, 'SERIES' ) )
    name = strcat( name, pathdelimiter, 'DATA' );
    d = getTree( name );
    d(1) = [];
    if length( d ) == 0
        disp( 'No data in this catalog' )
        return
    end
else
% alt 2: name is a file named SERIES
    fid = fopen( name );
    if fid == -1
        disp( 'Error: Your file name is invalid' )
        return
    end
    d = [];
    [ txt, pos ] = readUntil( fid, rootsymbol );
    while pos >= 0
        if( txt(end) == '~' )
            txt = '~';
        else
            txt = '';
        end
        txt = strcat( txt, readUntil( fid, char(10) ) ); % char(10) = newline
        d = [ d ; struct( 'name', txt ) ];
        [ txt, pos ] = readUntil( fid, rootsymbol );
    end
    if length( d ) == 0
        disp( 'No data in this catalog' )
        return
    end
end

pos = findstr( name, home );
if pos
    titlename = strcat( '~', name( length( home ) + 1 : end ) );
else
    titlename = name;
end

% printTitleName is the file name but changed into TeX format,
% so that TeX characters like _ and \ are changed inte \_ and \\
if ispc
    printTitleName = [];
    ind = find( titlename == '\' );
    ind = [ 1 ind length( titlename ) ];
    for i = 1:length( ind )-1
        printTitleName = [ printTitleName titlename(ind(i):ind(i+1)) ];
    end
else
    printTitleName = titlename;
end
titlename = printTitleName;
ind = [ find( titlename == '_' | titlename == '^' ) length(titlename)+1 ];
printTitleName = titlename(1:ind(1)-1);
for i = 1:length( ind ) - 1
    printTitleName = [ printTitleName '\' titlename(ind(i):ind(i+1)-1) ];
end

% plot the data in the catalogs in d
len = length( d );

% How many pages?
pages = floor( (len*(bold+1)-1) / (rows*cols) ) + 1;
for i = 1: pages
    figure ( i )
    clf
    set( i,'Position', [ 300 * mod( ( i ), 2 ) 55 720 635 ] );
end

hand = zeros( 1, len ); % subplot handles
height = 0.93 / rows;
width = 0.90 / cols;
for ind = 1 : bold+1 : (bold+1)*len
    for i = 0:bold
        page = floor( (ind+i-1)*(bold+1) / ( rows*cols ) ) + 1;
        subp = ind+i - ( page-1 ) * rows * cols; 
        figure( page )
        if subp == 1
            clf
        end
        bottom = 1 - height * ( floor( (subp-1) / cols ) + 1 );
        left = 0.1 + mod( (subp-1), cols ) * width;
        hand( ind ) = subplot( 'Position', [ left bottom width-0.05 height-0.07 ]);
    end
end

for ind = 1 : len
    rastname = d(ind).name;
    pos = find( rastname == pathdelimiter );
    rastname = rastname( pos(end)+1:end );
    filename = strcat( d(ind).name, pathdelimiter, 'APs' );
    load( filename )
    x = APs(:,1);
    y = APs(:,2);
    filename = strcat( d(ind).name, pathdelimiter, 'Params' );
    load( filename )
    version = Params(3);
    if version == 6
        EXP_dt = 100;
        Params = [ Params(1:7) ; EXP_dt ; Params(8:end) ];
    end
    filename = strcat( d(ind).name, pathdelimiter, 'C' );
    load( filename )
    relNMDA = C( 32 );
    version = Params(1);
    filename = strcat( d(ind).name, pathdelimiter, 'Q' );
    load( filename )
	tStart = Params(5);
    filename = strcat( d(ind).name, pathdelimiter, 'X' );
    clear X
    if exist(filename) == 2
        load( filename )
        if length(X)>0
            if X(1) == 4 | X(1) == 5
                r_c = X(end)/(X(end)+X(end-1) );
                relNMDA = r_c;
            end      
        end
    end
  
    page = floor( (bold+1)*(ind-1) / ( rows*cols ) ) + 1;
    subp = ind*(bold+1) - ( page-1 ) * rows * cols;
    figure( page ) 
    subplot( hand( ind*(bold+1) ) )
    
    % specify call-back routine to be able to enlarge single simulations
    str = strcat( 'showBig(  ''',d(ind).name,''' )' );
    set( hand( ind ), 'ButtonDownFcn', str );
    
    
    % plot the rastergrams of the cell
    title( sprintf( '%s, %s %%', rastname, int2str(100*relNMDA) ), 'Fontsize', fs )
    hold on
    plot( x, y, '.', 'MarkerSize', 1 )
    set( gca, 'box', 'on' )
    
    % Plot text and draw lines in graph
    N = 0;
    nmod = Params(4);
    tstart = Params(5);
    tstop = Params(6);
    %tQ = Q(1);
    %xQ = find( x<tQ );
    for i = 1:nmod
        NI = Params(7+6*i);
        NE = Params(10+6*i);
        ipos(i) = N+NI;
        epos(i) = N+NE+NI;
        xpos(i) = N+NI;
        %fS(i) = length( find( y(xQ)>=NI ) )/(tQ/1000*NE);
        fS(i) = 0;
        fE(i) = length( find( y>=NI ) )/(tstop/1000*NE);
        fI(i) = length( find( y<NI ) )/(tstop/1000*NI);
        line( [0 tstop], [N+NI-0.5 N+NI-0.5], 'Color', 'k', 'LineWidth', 2 )
        if i>1
            line( [0 tstop], [N-0.5 N-0.5], 'Color', 'k', 'LineWidth', 2 )
        end
        N = N+NI+NE;
    end
            
    set( gca, 'YLim', [0 N] )
    set( gca, 'XLim', Params(5:6) )
    set( gca, 'FontSize', fs )
    for i = 1:length(ipos)
        tI(i) = text( 0.1*(tstop-tstart)+tstart, ipos(i), 'I', 'FontSize', fs+2, 'FontWeight', 'Bold', 'VerticalAlignment', 'top' ); 
        tE(i) = text( 0.1*(tstop-tstart)+tstart, epos(i), 'E', 'FontSize', fs+2, 'FontWeight', 'Bold', 'VerticalAlignment', 'top' );
        str = sprintf( 'fE: %.1f\nfS: %.1f\nfI: %.1f\n', fE(i), fS(i), fI(i) );
        tx(ind) = text( 0.1*(tstop-tstart)+tstart, xpos(i), str, 'FontSize', fs-2, 'VerticalAlignment', 'bottom' );
    end
    if mod((subp-1), cols )
        set( gca, 'YTickLabel', [] )
    end
    
    if bold
        page = floor( ((bold+1)*(ind-1)+1) / ( rows*cols ) ) + 1;
        subp = (ind*(bold+1)+1) - ( page-1 ) * rows * cols;
        figure( page ) 
        subplot( hand( ind*(bold+1)+1 ) )
        
        [ t, B ] = bold( APs, netborder );
        plot( t, B )
        leg = [];
        for i = 1:size( B, 2 )
            leg = [ leg ; strcat( 'Net', int2str( i ) ) ];
        end
        legend( leg )
        xlabel( 'Time (ms)', 'FontSize', fs )
        title( 'BOLD', 'FontSize', fs )

        % specify call-back routine to be able to enlarge single simulations
        str = strcat( 'showBold(  ''',d(ind).name,''' )' );
        set( hand( ind ), 'ButtonDownFcn', str );
    end

end

% set titles and labels
for i = 1 : pages
    figure( i )
    thand( i ) = axes('Position',[0 0 1 1],'Visible','off');
    set(gcf,'CurrentAxes',thand( i ) )
    xl(i) = text(.5,.03,'Time (ms)','FontSize',fs, 'HorizontalAlignment', 'center');
    yl(i) = text(.03, .5, 'Cell #', 'Rotation', 90, 'FontSize', fs, 'VerticalAlignment', 'middle');
    tit(i) = text(.5, .99, sprintf( 'Rastergrams in %s', printTitleName ), 'FontSize', fs, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

% adjust paper size for printing
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
    if DoPlot
            print
    end
end
