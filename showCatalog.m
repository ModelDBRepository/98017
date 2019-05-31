function showCatalog( name, savename, DoPlot, rows, cols, alt_title )

% showCatalog displays a series of pages with simulations all
% belonging to a simulation series. By clicking on a graph, it can
% be enlarged. showCatalog should be in the same catalog as
% showConnMulti.m, showBig.m, bumpanalysis.m
%
% name  : is either a path or a file name. 
%         1. When name is a path, it is a path to the catalog which
%            contains the simulations to be plotted. Those
%            simulations should be located in directories found in
%            a subdirectory of that path named DATA.
%         2. When name is a file name, it should be a file name
%            containing the name SERIES. That file contains the
%            simulations which are to be plotted
% savename: The name of the eps output file
% DoPlot: if DoPlot = 0 or is left out, the figure is plotted but not printed
%         if DoPlot = 1, the figure is both plotted and printed
% rows  : Number of rows on each page
% cols  : Number of columns on each page
% alt_title: An alternative title. These choices are possible:
%            1: the X->E connection strength
%
% Change path names to suit your preferences. Places are marked % CHANGE %
rratePPC = [];
rratePFC = [];
if nargin < 1 | nargin > 6
    disp( 'usage1: showCatalog( name ), ' )
    disp( 'usage2: showCatalog( name, savename ), ' )
    disp( 'usage3: showCatalog( name, savename, DoPlot ), ' )
    disp( 'usage4: showCatalog( name, savename, DoPlot, rows, cols ), ' )
    disp( 'usage5: showCatalog( name, savename, DoPlot, rows, cols, alt_title ), ' )
    disp( 'Do "help showCatalog" for further information' )
    return
end
if nargin < 6
    alt_title = -1;
end
if nargin < 5
    rows = 4;
    cols = 3;
elseif isempty(cols) || isempty(rows)
    rows = 4;
    cols = 3;
end    
if nargin < 3
    DoPlot = 0;
elseif isempty(DoPlot)
    DoPlot = 0;
end
if nargin < 2
    save = 0;
elseif length(savename) == 0
    save = 0;
else
    save = 1;
end
thisdir = pwd;
fs = 10; % Fontsize in plots

% Is name a relative path or file name?
if ispc
    home = 'C:\Documents and Settings\Fredrik Edin\Mina Dokument\Neuron\STANDARDFILER';  % CHANGE %    
    pathdelimiter = '\';
else
    home = '/E:';
    pathdelimiter = '/';
end
if ~( name(1:3) == 'C:\' )
    name = strcat( thisdir, pathdelimiter, name );
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
%     [ txt, pos ] = readUntil( fid, pathdelimiter );
    [ txt, pos ] = readUntil( fid, '/\' );
    while pos >= 0
        if isunix
            if( txt(end) == '~' )
                txt = '~';
            else
                txt = '';
            end
            txt = strcat( txt, readUntil( fid, char(10) ) ); % char(10) = newline
        else
            txt = strcat( 'C:', readUntil( fid, char(10) ) ); % char(10) = newline
        end
        d = [ d ; struct( 'name', txt ) ];
%         [ txt, pos ] = readUntil( fid, pathdelimiter );
        [ txt, pos ] = readUntil( fid, '/\' );
    end
    fclose(fid);
    if length( d ) == 0
        disp( 'No data in this catalog' )
        return
    end
end
if isunix
    str = '/afs/nada.kth.se/home/o/u1sxc4xo';    % CHANGE %
    pos = findstr( name, str );
    if pos
        titlename = strcat( '~', name( length( str ) + 1 : end ) );
    else
        titlename = name;
    end
end

% plot the data in the catalogs in d
len = length( d );

% How many pages?
pages = floor( (len-1) / (rows*cols) ) + 1;
for i = 1: pages
    figure ( i )
    clf
    if isunix
        set( i,'Position', [ 300 * mod( ( i ), 2 ) 100 720 659 ] );
    else
        set( i,'Position', [ 558 * mod( ( i ), 2 ) 35 720 685 ] );
    end
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

for ind = 1 : len
    rastname = d(ind).name;
    pos = find( rastname == '/' | rastname == '\' );
    rastname = rastname( pos(end)+1:end );
    filename = strcat( d(ind).name, pathdelimiter, 'APs.txt' );

    % This section is used to convert old files with no .txt extensions to
    % the new ones
    if exist(filename(1:end-4),'file') & exist(filename,'file')
        delete(filename(1:end-4))
    elseif exist(filename(1:end-4),'file')
        movefile(filename(1:end-4), filename);
    end
    
    % This is done because sometimes the simulations produces error outputs
    % in the data files which is removed here
    try
        load( filename )
    catch
        fid = fopen(filename);
        str = fread(fid);
        fclose(fid);
        errors = find(str > 57);
        linebreaks = find(str == 10);
        linebreaks_after_error = linebreaks(find(linebreaks>errors(end)));
        start = linebreaks_after_error(1) + 1;
        txtAPs = char(str(start:end)');
        APs = str2num(txtAPs);
        fid = fopen(filename,'w');
        fprintf(fid,'%s',txtAPs);
        fclose(fid)
    end
    x = APs(:,1);
    y = APs(:,2);
    rratePPC = [rratePFC length(find(x>4000 & y<160))];
    rratePFC = [rratePFC length(find(x>4000 & y>=160))];
    disp(['rate PPC: ' int2str(length(find(x>4000 & y<160)))])
    disp(['rate PFC: ' int2str(length(find(x>4000 & y>=160)))])
    filename = strcat( d(ind).name, pathdelimiter, 'Params.txt' );
    % This section is used to convert old files with no .txt extensions to
    % the new ones
    if exist(filename(1:end-4),'file') & exist(filename,'file')
        delete(filename(1:end-4))
    elseif exist(filename(1:end-4),'file')
        movefile(filename(1:end-4), filename);
    end
    load( filename )
    filename = strcat( d(ind).name, pathdelimiter, 'C.txt' );
    % This section is used to convert old files with no .txt extensions to
    % the new ones
    if exist(filename(1:end-4),'file') & exist(filename,'file')
        delete(filename(1:end-4))
    elseif exist(filename(1:end-4),'file')
        movefile(filename(1:end-4), filename);
    end
    load( filename )
    C = reshape(C,24,length(C)/24);
    if sum(C(12,find(C(1,:)==1 & C(2,:)>0))) > 0
        relNMDA = sum( C(12,find(C(1,:)==1 & C(2,:)>0)) .* C(8,find(C(1,:)==1 & C(2,:)>0))  )/sum(C(12,find(C(1,:)==1 & C(2,:)>0))); % of total synaptic weight
    else
        relNMDA = 0;
    end
    %relNMDA = C( 32 );
%     disp([2*C(12,9)/sum(C(12,[9 7])) 2*C(12,10)/sum(C(12,[8 10]))])  
%     disp([sum(C(12,[9 7])) sum(C(12,[8 10])) sum(C(12,[1 3])) sum(C(12,[2 4]))])  
%     disp([C(12,9) C(12,10)])  
    version = Params(1);
    if version > 1 
        filename = strcat( d(ind).name, pathdelimiter, 'Q.txt' );
        % This section is used to convert old files with no .txt extensions to
        % the new ones
        if exist(filename(1:end-4),'file') & exist(filename,'file')
            delete(filename(1:end-4))
        elseif exist(filename(1:end-4),'file')
            movefile(filename(1:end-4), filename);
        end
        load( filename )
	if version == 2
	    tStart = 0;
	    Icelltyp = 1;
	    tmp = [ Params(1:3) tStart Params(4:5) 100 100];
	    for i = 1:nmod
	        tmp = [ tmp Icelltyp Params(3+3*i:5+3*i) ];
	    end
	    Params = tmp;
            tmp = Params(1:8);
            for i = 1:Params(2)
                tmp = [ tmp Params(5+4*i) 1000 Params(7+4*i) Params(6+4*i) 1000 Params(8+4*i) ];
            end
            Params = tmp;
        elseif version > 2
	    tStart = Params(4);
            %filename = strcat( d(ind).name, '/', 'X' );
            %clear X
            %if exist(filename) == 2
            %    load( filename )
            %    if length(X)>0
            %        if X(1) == 4 | X(1) == 5
            %            r_c = X(end)/(X(end)+X(end-1) );
            %            relNMDA = r_c;
            %        end
            %    end
            %end      
        end
        if version == 3
            tmp = [ Params(1:6) ; 100 ; 100 ];
            for i = 1:Params(2)
                try
                    tmp = [ tmp ; Params(3+4*i) ; 1000 ; Params(5+4*i) ; Params(4+4*i) ; 1000 ; Params(6+4*i) ];
                catch
                    keyboard
                end
            end
            Params = tmp;
        elseif version == 4
            tmp = [ Params(1:6) ; 100 ; 100 ; Params(7:end) ];
            Params = tmp;
        end
    end
  
    page = floor( (ind-1) / ( rows*cols ) ) + 1;
    subp = ind - ( page-1 ) * rows * cols;
    figure( page ) 
    subplot( hand( ind ) )
    
    % specify call-back routine to be able to enlarge single simulations
    str = strcat( 'showBig(  ''',d(ind).name,''' )' );
    set( hand( ind ), 'ButtonDownFcn', str );

    % plot the rastergrams of the cell
    hold on
    plot( x, y, '.', 'MarkerSize', 1 )
    set( gca, 'box', 'on' )
    
    % Plot text and draw lines in graph
    
    % If net is the old 2-population version
    if version == 1
        NI = 2^( ceil( log2( max( y )/5 ) ) );
        NE = 4*NI;
        N = NI+NE;
        ipos = NI;
        epos = N;
        xpos = NI;
        tStop = 100 * ceil( max( x )/100 );
        tQ = 1000;
        xQ = find( x<tQ );
        line( [0 tStop], [NI-0.5 NI-0.5], 'Color', 'k', 'LineWidth', 2 )
	
    % else if net is the new multimodule version
    else
        N = 0;
        tStop = Params(5);
	if length(Q) > 0
            tQ = Q(1);
	else
	    tQ = Params(4);
	end
        %tQ = Q(1);
        %xQ = find( x<tQ );
        for i = 1:Params(2)
            NI = Params(5+6*i);
            NE = Params(8+6*i);
            ipos(i) = N+NI;
            epos(i) = N+NE+NI;
            xpos(i) = N+NI;
            line( [0 tStop], [N+NI-0.5 N+NI-0.5], 'Color', 'k', 'LineWidth', 2 )
            if i>1
                line( [0 tStop], [N-0.5 N-0.5], 'Color', 'k', 'LineWidth', 2 )
            end
            N = N+NI+NE;
        end
    end
        
    netborder = [ 0 ; cumsum(Params(11:3:end))];
    Efreq = [];
    for i = 2:2:length(netborder)
        Efreq = [Efreq length(find(y>=netborder(i) & y<netborder(i+1))) ...
		 / ((netborder(i+1)-netborder(i))*((tStop-tQ)/1000))];
    end
    Efreq = mean(Efreq);		 
    if alt_title == -1
        title( sprintf( '%s, %s%% f:%.2f Hz', rastname, int2str(100*relNMDA), Efreq ), 'Fontsize', fs )
    elseif alt_title == 1
        title( sprintf( '%s %s %s', 'g_{X\rightarrowE}:', num2str(C(12,2)), 'mS/cm^2' ), 'Fontsize', fs )
    end
    %title( sprintf( '%s, %s %%', rastname, int2str(100*relNMDA) ), 'Fontsize', fs )
    set( gca, 'YLim', [0 N] )
    set( gca, 'XLim', [tStart tStop] )
    set( gca, 'Fontsize', fs )
    if mod((subp-1), cols )
        set( gca, 'YTickLabel', [] )
    end

end

% set titles and labels
% for i = 1 : pages
%     figure( i )
    %thand( i ) = axes('Position',[0 0 1 1],'Visible','off');
    %set(gcf,'CurrentAxes',thand( i ) )
    %xl(i) = text(.5,.03,'Time (ms)','FontSize',14, 'HorizontalAlignment', 'center');
    %yl(i) = text(.03, .5, 'Cell #', 'Rotation', 90, 'Fontsize', fs , 'VerticalAlignment', 'middle');
    %tit(i) = text(.5, .99, sprintf( 'Rastergrams in %s', titlename ), 'Fontsize', fs, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
% end

% adjust paper size for printing
for i = 1 : pages
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
    if save
        print( i, '-depsc2', [savename '_p' int2str(i) '.eps'])
    end
end

disp(['rate PPC: ' int2str(mean(rratePPC)) ', rate PFC: ' int2str(mean(rratePFC))]);
