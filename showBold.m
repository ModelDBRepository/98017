function showBold( name, choice, curr, ave, multipanel, fig0, filename )

% showBold displays a series of bold curves of simulations all
% belonging to a simulation series.
%
% name  : is either a path or a file name. 
%         1. When name is a path, it is a path to the catalog which
%            contains the simulations to be plotted. Those
%            simulations should be located in directories found in
%            a subdirectory of that path named DATA.
%         2. When name is a file name, it should be a file name
%            containing the name SERIES. That file contains the
%            simulations which are to be plotted
% choice: Type of bold signal. See bold.m for more
%            information. Default = 5
% curr  : If curr = 1, then synaptic currents are shown instead of
%            BOLD signal. Default = 0
% ave   : If ave = 1, then the average BOLD or current curve will
%            be plotted, otherwise single curves will be
%            plotted. Default = 0.
% multipanel: If multipanel = 1, then every signal will be plotted
%            in its own panel. Otherwise, all signals will be
%            plotted in the same panel. Default = 1.
% fig0  : The number of the figure where data should be plotted. Default = 1
%
% Change path names to suit your preferences. Places are marked % CHANGE %

if nargin < 1 | nargin > 7
    disp( 'usage1: showBold( name ), ' )
    disp( 'usage2: showBold( name, choice ), ' )
    disp( 'usage3: showBold( name, choice, curr ), ' )
    disp( 'usage4: showBold( name, choice, curr, ave ), ' )
    disp( 'usage5: showBold( name, choice, curr, ave, multipanel ), ' )
    disp( 'usage6: showBold( name, choice, curr, ave, multipanel, fig0), ' )
    disp( 'usage7: showBold( name, choice, curr, ave, multipanel, fig0, filename ), ' )
    disp( ' ' )
    disp( 'CHOICES:' )
    disp('1: Action Potentials from both inhibitory and excitatory populations')
    disp('2: Action Potentials from excitatory population only')
    disp('3: Sum of Excitatory Currents')
    disp('4: Sum of Excitatory Currents onto Excitatory Cells')
    disp('5: Sum of absolute values of all currents')
    disp('6: Sum of all currents')
    disp( ' ' )
    disp( 'Do "help showBold" for further information' )
    return
end
if  nargin < 2
    choice = 5;
end
if nargin < 3
    curr = 0;
end
if nargin < 4
    ave = 0;
end
if nargin < 5
    multipanel = 1;
end
if nargin < 6
    fig0 = 1;
end
if nargin < 7
    saving = 0;
else
    saving = 1;
end
fig0 = fig0-1;
smooth = 250;
rows = 4;
cols = 3;

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
if multipanel == 1
    pages = floor( (len-1) / (rows*cols) ) + 1;
    for i = 1: pages
        figure ( fig0+i )
        clf
        set( fig0+i,'Position', [ 300 * mod( ( i ), 2 ) 100 720 659 ] );
    end
else
    figure ( fig0+1 )
    clf
    set( fig0+1,'Position', [ 300 100 720 659 ] );
    set( gcf, 'Name', 'BOLD' )
    set( gcf, 'NumberTitle', 'off' )
end


% find number of subplots = no of modules in the nets
maxAntalModuler = 0;
for ind = 1 : len
    [ t_SERIE{ind} B_SERIE{ind} ] = bold( d(ind).name, choice, curr, smooth );
    maxAntalModuler = max( maxAntalModuler, length( B_SERIE{ind} ) );
end

for i = 1:maxAntalModuler
    leg_SERIE{i} = [];
end

% plot bold curve for every module in every simulation
% Go through all simulations
color = 'brgky';
type = {'-','-.','--'};
zr = '000000000000';
maxlen = length( int2str( len ) );
mB = 100000;
MB = 0;
uB = 0;
if ave
    aveB = zeros(length(B_SERIE{1}{i}),length(B_SERIE{1})+1);
end
for ind = 1 : len

    if curr
        t = t_SERIE{ind}/1000;
    else
        t = t_SERIE{ind};
    end
    B = B_SERIE{ind};
    % plot bold curves for every module in a simulation
    for i = 1:length( B )

        indstr = int2str( ind );
        indl = length( indstr );
        leg = leg_SERIE{i};
        leg = [ leg ; strcat( zr(1:maxlen-indl), indstr ) ];
        leg_SERIE{i} = leg;

        fignum = 1;
	if ave
	    aveB(:,i+1) = aveB(:,i+1) + B{i}'/len;
	    aveB(:,1) = t_SERIE{1}';
	else
	    if multipanel == 0
	        subplot(maxAntalModuler, 1, i);
                plot( t, B{i}, strcat( color(mod(ind-1,length(color))+1), type{mod(ind-1,length(type))+1} ) )
	    else
	        fignum = floor((ind-0.5)/(rows*cols)) + 1;
	        figure(fig0+fignum) 
	        subplot(rows, cols, ind-rows*cols*(fignum-1));
	        plot( t, B{i}, strcat( color(mod(i-1,length(color))+1), type{mod(i-1,length(type))+1} ) )
	        hold on
	        xlim([t(1) t(end)])
	    end
	end
	if curr
	    uB = max(uB, mean(B{i}));
	    mB = min(mB, min(B{i}));
	else
	    xlim([t(1) t(end)])
	    mB = min(mB, min(B{i}));
	    MB = max(MB, max(B{i}));
	end
	hold on
%        set( gca, 'XLim', [t(1) t(end)] )
        set( gca, 'FontSize', fs )
        if curr
	    if multipanel == 0
	        if i == 1
	            title( [ 'Bold curves in ' printTitleName ], 'FontSize', fs )
	        end
	    end
	else
	    if multipanel == 0
	        if i == 1
	            title( [ 'Bold curves in ' printTitleName ], 'FontSize', fs )
	        end
	    end
	    %ylabel( ['Module ' int2str( i ) ' BOLD signal (%)'], 'FontSize', fs )
	end
    end

    
%    set( gca, 'YLim', [0 N] )
%    set( gca, 'XLim', Params(5:6) )

end

if ave
    for i = 2:size(aveB,2)
        plot(t, aveB(:,i), strcat( color(mod(i,length(color))+1), type{mod(i,length(type))+1} ) )
    end
    xlim([t(1) t(end)])
    if curr == 0
        ylim([mB MB])
    else
        ylim([mB 1.5*uB])
    end
else
    for ind = 1:len
        fignum = floor((ind-0.5)/(rows*cols)) + 1;
        figure(fig0+fignum) 
        subplot(rows, cols, ind-rows*cols*(fignum-1));
        if curr == 0
            ylim([mB MB])
        else
            ylim([mB 1.5*uB])
	end
    end
end



% set titles and labels
if multipanel == 0
    figure( fig0+1 )
    subplot( maxAntalModuler, 1, maxAntalModuler )
    xlabel( 'time (s)', 'FontSize', fs )
    for i = 1 : maxAntalModuler
        subplot( maxAntalModuler, 1, i )
        legend( leg_SERIE{i} )
        set( gca, 'box', 'on' )
        set( gca, 'FontSize', fs )
    end
end


% adjust paper size for printing
for i = 1:fignum
    figure(fig0+fignum)
    set( gcf, 'PaperUnits', 'centimeters' );
    set( gcf, 'PaperType', 'A4');
    papersize = get( gcf, 'PaperSize' );
    left = 0.02;
    bottom = 0.02;
    width = papersize( 1 ) - 2 * left;
    height = papersize( 2 ) - 2 * left;
    myfiguresize = [ left, bottom, width, height ];
    set( gcf, 'PaperPosition', myfiguresize );
    if saving
        print( fig0+fignum, '-depsc2', [filename '_p' int2str(fignum) '.eps'] );
	if ave
	    save( [filename '.dat'], '-ascii', 'aveB' )
	end
    end
end

