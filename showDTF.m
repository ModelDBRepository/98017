function showDTF( name, savename, DoPlot, rows, cols )

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
%
% Change path names to suit your preferences. Places are marked % CHANGE %

if nargin < 1 | nargin > 5
    disp( 'usage1: showCatalog( name ), ' )
    disp( 'usage2: showCatalog( name, savename ), ' )
    disp( 'usage3: showCatalog( name, savename, DoPlot ), ' )
    disp( 'usage4: showCatalog( name, savename, DoPlot, rows, cols ), ' )
    disp( 'Do "type showCatalog" for further information' )
    return
end
if nargin < 5
    rows = 4;
    cols = 3;
end
if nargin < 2
    save = 0;
elseif length(savename) == 0
    save = 0;
else
    save = 1;
end
if nargin < 3
    DoPlot = 0;
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
    for j = 0:3
        figure ( i+j*pages )
        clf
        if isunix
            set( i+j*pages,'Position', [ 300 * mod( ( i ), 2 ) 100 720 659 ] );
        else
            set( i+j*pages,'Position', [ 558 * mod( ( i ), 2 ) 35 720 685 ] );
        end
    end
end

hand = zeros( 1, 3*len ); % subplot handles
height = 0.93 / rows;
width = 0.90 / cols;
for ind = 1 : len
    for j = 0:3
        page = floor( (ind-1) / ( rows*cols ) ) + 1;
        subp = ind - ( page-1 ) * rows * cols; 
        figure( page + j*pages )
        if subp == 1
            clf
        end
        bottom = 1 - height * ( floor( (subp-1) / cols ) + 1 );
        left = 0.1 + mod( (subp-1), cols ) * width;
        hand( ind+j*len ) = subplot( 'Position', [ left bottom width-0.05 height-0.07 ]);
    end
end

for ind = 1 : len
    rastname = d(ind).name;
    pos = find( rastname == '/' | rastname == '\' );
    rastname = rastname( pos(end)+1:end );
    filename = strcat( d(ind).name, pathdelimiter, 'LFP_prox_dist.txt' );

    % This section is used to convert old files with no .txt extensions to
    % the new ones
    if exist(filename(1:end-4),'file') & exist(filename,'file')
        delete(filename(1:end-4))
    elseif exist(filename(1:end-4),'file')
        movefile(filename(1:end-4), filename);
    end
    LFP = load( filename );
    
    filename = strcat( d(ind).name, pathdelimiter, 'Params.txt' );
    % This section is used to convert old files with no .txt extensions to
    % the new ones
    if exist(filename(1:end-4),'file') & exist(filename,'file')
        delete(filename(1:end-4))
    elseif exist(filename(1:end-4),'file')
        movefile(filename(1:end-4), filename);
    end
    load( filename );
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
                tmp = [ tmp ; Params(3+4*i) ; 1000 ; Params(5+4*i) ; Params(4+4*i) ; 1000 ; Params(6+4*i) ];
            end
            Params = tmp;
        elseif version == 4
            tmp = [ Params(1:6) ; 100 ; 100 ; Params(7:end) ];
            Params = tmp;
        end
    end
  
    % extract the LFP
    if length(Q) > 0 %& Q(1,5) ~= 0
        % Normalize data to 0 mean and variance 1
%         tmp2 = LFP(Q(1,1)-tStart+1000:Q(1,1)-tStart+5000,2);
        tmp2 = LFP(Q(1,1)-tStart+1000:end,2);
        tmp2 = (tmp2 - mean(tmp2)) / std(tmp2);
        tmp3 = LFP(Q(1,1)-tStart+1000:end,3);
        tmp3 = (tmp3 - mean(tmp3)) / std(tmp3);
        % Store as LFPs
%         LFP_ppc = [LFP_ppc tmp2]; 
%         LFP_pfc = [LFP_pfc tmp3]; 
        LFP_ppc = tmp2; 
        LFP_pfc = tmp3; 
    end
    
    % Do power spectrum and directed transfer function
    dt = LFP(2,1)-LFP(1,1);
    fS = 1000/dt;
    fN = fS/2;
    
    [Pp f_sp] = psd(LFP_ppc,1000,fS);
    [Pf f_sp] = psd(LFP_pfc,1000,fS);
    
    f = 0:1:100;    
    [nD,D,order]=DTF([LFP_ppc LFP_pfc],fS,f);
    nD1to2 = reshape(nD(2,1,:),[],1);
    nD2to1 = reshape(nD(1,2,:),[],1);
    D1to2 = reshape(D(2,1,:),[],1);
    D2to1 = reshape(D(1,2,:),[],1);

    page = floor( (ind-1) / ( rows*cols ) ) + 1;
    subp = ind - ( page-1 ) * rows * cols;
    figure( page ) 
    subplot( hand( ind ) )
    hold on
    plot( f, nD1to2, 'r', f, nD2to1, 'k' )
    ylim([0 1])
    xlim([0 max(f)])
    set( gca, 'Fontsize', fs )
    set( gca, 'box', 'on' )
    if mod((subp-1), cols )
        set( gca, 'YTickLabel', [] )
    end
    title('nDTF')

    figure( page + pages ) 
    subplot( hand( ind+len ) )
    hold on
    plot( f, D1to2, 'r', f, D2to1, 'k' )
    ylim([0 1.1*max(max([D1to2 D2to1]))])
    xlim([0 max(f)])
    set( gca, 'Fontsize', fs )
    set( gca, 'box', 'on' )
%     if mod((subp-1), cols )
%         set( gca, 'YTickLabel', [] )
%     end
    title('DTF')
    
    figure( page + 2*pages ) 
    subplot( hand( ind+2*len ) )
    hold on
    plot( f_sp, Pp, 'r', f_sp, Pf, 'k' )
    ylim([0 1.1*(max(max([Pp Pf])))])
    xlim([0 max(f)])
    set( gca, 'Fontsize', fs )
    set( gca, 'box', 'on' )
%     if mod((subp-1), cols )
%         set( gca, 'YTickLabel', [] )
%     end
    title('Spectrum')
    
    
    % plot the rastergrams of the cell
    figure( page + 3*pages ) 
    subplot( hand( ind+3*len ) )
    hold on
    plot( x, y, '.', 'MarkerSize', 1 )
    set( gca, 'box', 'on' )
    
    % specify call-back routine to be able to enlarge single simulations
    str = strcat( 'showBig(  ''',d(ind).name,''' )' );
    set( gca, 'ButtonDownFcn', str );
    
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

    netborder = [ 0 ; cumsum(Params(11:3:end))];
%     title( sprintf( '%s, %s%% f:%.2f Hz', rastname, int2str(100*relNMDA), Efreq ), 'Fontsize', fs )
    title( 'Raster', 'Fontsize', fs )
    %title( sprintf( '%s, %s %%', rastname, int2str(100*relNMDA) ), 'Fontsize', fs )
    ylim( [0 N] )
    xlim( [tStart tStop] )
    set( gca, 'Fontsize', fs )
    if mod((subp-1), cols )
        set( gca, 'YTickLabel', [] )
    end
   
    if ind == 1
        leg(page) = legend(hand(ind),{'n1->2','n2->1'});
        leg(page+2*pages) = legend(hand(ind+len),{'1->2','2->1'});
        leg(page+2*pages) = legend(hand(ind+2*len),{'ppc','pfc'});
    end
    
    
%    title( sprintf( '%s, %s%% f:%.2f Hz', rastname, int2str(100*relNMDA), Efreq ), 'Fontsize', fs )

end

% set titles and labels
for i = 1 : pages
%     figure( i )
    %thand( i ) = axes('Position',[0 0 1 1],'Visible','off');
    %set(gcf,'CurrentAxes',thand( i ) )
    %xl(i) = text(.5,.03,'Time (ms)','FontSize',14, 'HorizontalAlignment', 'center');
    %yl(i) = text(.03, .5, 'Cell #', 'Rotation', 90, 'Fontsize', fs , 'VerticalAlignment', 'middle');
    %tit(i) = text(.5, .99, sprintf( 'Rastergrams in %s', titlename ), 'Fontsize', fs, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

% adjust paper size for printing
% for i = 1 : page
%     figure( i )    
%     set( gcf, 'PaperUnits', 'centimeters' );
%     set( gcf, 'PaperType', 'A4');
%     papersize = get( gcf, 'PaperSize' );
%     left = 0.02;
%     bottom = 0.02;
%     width = papersize( 1 ) - 2 * left;
%     height = papersize( 2 ) - 2 * left;
%     myfiguresize = [ left, bottom, width, height ];
%     set( gcf, 'PaperPosition', myfiguresize );
%     if DoPlot
%         print
%     end
%     if save
%         print( i, '-depsc2', [savename '_p' int2str(i) '.eps'])
%     end
% end
