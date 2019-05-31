function showSimulation( dirname, pr, wh, no, smooth )

% plots BUMP-data. For version 1-5 of working memory
% network simulator.
% dirname: The name of the simulation directory
% pr = 0 --> merely display figures
% pr = 1 --> print
% pr = 2 --> print and save to file
% pr = 3 --> save to file
% If wh = 0 or nothing  --> print only figures 5 and 6
%    wh = 1             --> print all figures
%    wh = 2             --> print only figure 7
% no = number of enlargements of activity. Default = 5
% smooth =  degree of smoothing of histogram. Default = 2
%
% Version 2.0
% Author: Fredrik Edin, 2004
% Address: freedin@nada.kth.se


if nargin < 2 | nargin > 5
    disp( 'usage1: showSimulation( dirname, pr ), ' )
    disp( 'pr = 0: Only show' )
    disp( 'pr = 1: Print figures 5 & 6' )
    disp( 'pr = 2: Print and save figures 5 & 6' )
    disp( 'pr = 3: Save figures 5 & 6' )
    disp( '  ' )
    disp( 'usage2: showSimulation( dirname, pr, wh ), ' )
    disp( 'wh = 0: As above' )
    disp( 'wh = 1: Also print and/or save figure 7' )
    disp( 'wh = 2: Print and/or save only figure 7' )
    disp( '  ' )
    disp( 'usage3: showSimulation( dirname, pr, wh, no ), ' )
    disp( 'no = number of enlargements of the rastergram.' )
    disp( '  ' )
    disp( 'usage4: showSimulation( dirname, pr, wh, no, smooth ), ' )
    disp( 'smooth = number of cells making up a point in the rastergram.' )
    return
end

if nargin < 3
    wh = 0;
end
if nargin < 4
    no = 5;
end
if nargin < 5
    smooth = 2;
end
screensize = get( 0, 'ScreenSize' );


% Determine whether to print or not. If yes, another window must be made
% so that both figures can be printed on the same page
P = pr>0 & wh<2; 
thisdir = pwd;
cd( dirname )

% load data files
load Params
version = Params(3);
if version == 6
    EXP_dt = 100;
    Params = [Params(1:7) EXP_dt Params(8:end) ];
end
load Connections
load C
C = [ reshape( C, 25, length(C)/25 ) ]';
load Q
nmod = Params(4);

% Name of simulation
fid = fopen( 'Parameters' );
fileName = pwd;
fileName( find( fileName == '.' ) ) = ',';

% printTitleName is the file name but changed into TeX format,
% so that TeX characters like _ and \ are changed inte \_ and \\
if ispc
    printFileName = [];
    ind = find( fileName == '\' );
    ind = [ 1 ind length( fileName ) ];
    for i = 1:length( ind )-1
        printFileName = [ printFileName fileName(ind(i):ind(i+1)) ];
    end
else
    printFileName = fileName;
end
fileName = printFileName;
ind = [ find( fileName == '_' | fileName == '^' ) length(fileName)+1 ];
printFileName = fileName(1:ind(1)-1);
for i = 1:length( ind ) - 1
    printFileName = [ printFileName '\' fileName(ind(i):ind(i+1)-1) ];
end


% Window placement
figs = get( 0, 'Children' );
for i = 1:length( figs ) 
    if figs(i) >= 5 & figs(i) <= 10
        close( figs(i) )
    end
end
if P  % Create figure for printing
    figure(8)
    clf
    set( 8, 'Position', [ 0 100 720 852 ] )
end
figure(5)
clf
set( 5,'Position', [ 0.01*screensize(3) 0.3*screensize(4) 0.5*screensize(3) 0.61*screensize(4) ] );
set( 5, 'Name', 'CONNECTIONS' )
set( 5, 'NumberTitle', 'off' )

% Network connection plot
nmod = Params(4);
minw = 0; % minw and maxw are the limits of the footprint diagram
maxw = 0;% They are scaled to the largest and smallest values of the footprint.
netborder = [ 0 ; Params( 13:3:end ) ];
netborder = cumsum( netborder );
xp = 0:1/(1+nmod):1;


% Create population structs (rings with I:s or E:s in connection diagram)
% struct variables:
% label: The cell type
% x, y: Position in the connection diagram
% radius: ring radius
% N: Number of cells in population
% ar: Arrows entering population
% in: The populations having incoming arrows to this population
rad = 0.1/nmod;
arrowind = 65; % Capital A
arrows = [];
diagram = [];
for i = 1:nmod
    if Params(5+6*i) == 1
        str = 'I-1';
    elseif Params(5+6*i) == 2
        str = 'I-IF';
    end
    pop(2*i-1) = struct( 'label', str, 'x', xp(i+1), 'y', 0.3, 'radius', rad, 'N', ...
                         Params(7+6*i), 'ar', [], 'in',[] );
    if Params(8+6*i) == 0 
        str = 'E-1';
    elseif Params(8+6*i) == 1 % Type of E-cell
        str = 'E-3';
    elseif Params(8+6*i) == 2 
        str = 'E-IF';
    end
    pop(2*i) = struct( 'label', str, 'x', xp(i+1), 'y', 0.7, 'radius', rad, 'N', ...
                       Params(10+6*i), 'ar', [], 'in', [] );
end


% Create arrows between populations. Every row in C contains a connection between
% cell populations = an arrow in the connection diagram. Every arrow is a
% struct with these variables
% pre: pre-population index. If prepopulation is 0 (external stimulation), pre is a negative number 
%        the positive value of which defines an angle
% post: the post population
% label: Every arrow has a label with either the connection weight (if
%           connection is flat), or a letter with an index to a plot showing the connection footprint.
% diag: The index of the diagram with the connection footprint, if there is one.
%
% If a connection is not flat, a diagram struct is created. That has a
% number of variables needed to plot the footprint.
for i = 1:size( C, 1 )
    postmod = C(i,1);
    premod = C(i,2);
    for j = 1:4
        pre = 2*(premod-1)+mod( j+1, 2 )+1; % Netborder indices of population
        post = 2*(postmod-1)+(j>2)+1;
        w = C(i,9+j);
        Jp = C(i,11+3*j);
        Jm = C(i,12+3*j);
        sgm = C(i,13+3*j);
        
        if w > 0 % there is a connection

            if ~pre % If pre is external, then it must define an angle
                pre = -(4/3-mod(post-1,2))*pi;
                %pre = -(3-mod(post-1,2))*pi/3;
            end
            % Determine labpos, position of arrow label. So that labels
            % will always be printed close to the arrowhead on the correct
            % side, 
            if pre<1
                labpos = mod( mod( post-1,2 ) - post>2, 2 ); % All external arrows to the left
            elseif pre == 3 | pre == 2
                labpos = 1;
            elseif pre == 1 | pre == 4   
                labpos = 0;
            end
            if mod( post, 2 ) == mod( pre, 2 )
                labpos = 1-labpos;
            end
            if pre == post
                labpos = -(5/3-mod(post-1,2))*pi;
            end
            
            % If connection footprint is not flat, display it. Otherwise
            % just show connection strength
            if Jp >= 0 & Jp ~= 1
                lab = char( arrowind );
                arrowind = arrowind + 1;
                ind = find( Connections(:,1) >= netborder( post ) & ...
                            Connections(:,1) < netborder( post+1 ) & ...
                            Connections(:,2) >= netborder(pre) & ...
                            Connections(:,2) < netborder(pre+1) );
                ind = ind(1:pop(pre).N);
                W = Connections(ind,4);
                ang = Connections(ind,3);
                % REMOVE paramind if it is not used
                diagram = [ diagram struct( 'label', lab, 'pre', pre, 'post', post, ...
                            'paramind', j, 'W', W, 'w', w, 'angle', ang, 'Jp', Jp, 'Jm', Jm', 'sgm', sgm ) ];
                minw = min( min( minw, W ) ); 
                maxw = max( max( maxw, W ) );
            else
                lab = sprintf( '%.2g', w );
            end
            % Label is always placed immediately after arrow head
            % Labpos = 0 --> left of arrow head
            %        = 1 --> right of arrow head
            % If pre and pos are equal, then labpos specifies the
            % incidence angle in radians-
            % A negative value of pre indicates the (positive) angle of incidence to 
            % post and is used when there exists no pre.             
            arrows = [ arrows struct( 'pre', pre, 'post', post, 'label', lab, 'labpos', labpos, 'diag', length( diagram ) ) ];
            if pre > 0
                pop(pre).ar = [ pop(pre).ar length(arrows) ];
            end
            pop(post).in = [pop(post).in arrows( end ).pre ];
        end
    end
end

% Plot network architecture figure
len = length( diagram ) + 1;
m = sum( C(:,3) > 0 ); % There exist delay distributions
tot = (3+ m + len); % Total number of subplots
hand = zeros( tot*(P+1) ); % Handles to subplots in the diagram

pos = get( gca, 'Position' ); 
hgt = pos(4);
btt = pos(2);
for h = 0:P
    figure(5+3*h)
    pos(4) = hgt*(1-0.7*h);
    pos(2) = btt+0.7*hgt*h;
    
    % create subplots in correct positions
    % hand(1): title
    % hand(2): connection diagram
    % hand(3-x): delay diagrams if there exists delay distributions
    % hand(x+1 alt 3): parameter text
    % the rest: footprint plots
    spac1 = 0.03; 
    spac2 = 0.06;
    spac3 = 0.05;
    spac4 = 0.6*h; % Height of evolution-of-firing-rate-figure 
    hand(1+tot*h) = axes( 'Position', [ 0, 0.95, 1, 0.03 ] );
    hand(2+tot*h) = subplot( 'Position', [ 0.25+spac1, 2*spac1+spac4*h, 0.5-2*spac1, 0.95-3*spac1-0.3*(m>0)-spac4 ] );
    p1 = 0.25+2*spac1; % Left edge of delay distribution diagram
    w = 0.25-3*spac1;
    for i = 1 : m
        hand(2+i+tot*h) = subplot( 'Position', [ p1, 0.95+(spac2-0.3)*(1-spac4), w, (0.3-2*spac1-spac2)*(1-spac4) ] );
        p1 = p1 + w + spac1;
    end
    n1 = floor( len / 2 );
    n = n1;
    for k = 0:1 % Footprint diagrams are created, first on the left side, then on the right
        for j = 1:n
           height = (0.9-spac3)/n;
           hand(m+2+j+k*n1+tot*h) = subplot( 'Position', [ 0.75*k+spac1, 0.95+(-j*height+spac2)*(1-spac4), 0.25-2*spac1, (height-2*spac2)*(1-spac4) ] );
        end
        n = len - n;
    end
   
    % Plot titel
    subplot( hand(1+tot*h) )
    t0 = text( 0.5, 0.5, printFileName, 'FontSize', 14, 'HorizontalAlignment', ...
	       'center' );
    set( hand(1+tot*h), 'Visible', 'off' )
    
    % Plot connection diagram
    cd( thisdir )
    subplot( hand(2+tot*h) )
    plotPop(pop)
    hold on
    plotArrows( arrows, pop )
    set( gca, 'XTick', [] )
    set( gca, 'YTick', [] )
    if h
        axis off
    end
    
    % Plot delay distributions
    k = 0;
    My = 0;
    for i = 1:size( C, 1 )
        dist = C( i, 3 ); % type of distribution
        mu = C( i, 4 ); % mean delay
        if dist > 0
            k = k + 1;
            axes( hand(2+k+tot*h) )
            xlabel( 'tid (s)' )
            str = sprintf( 'Delay modul %d - %d', C(i,2), C(i,1) );
            title( str )
        end
        if dist == 1 % Delta distribution
            line( [mu mu], [0 1] )
        elseif dist == 2 % Gamma distribution
            s2 = C(i, 5); % variance of delay distribution
            xm = mu+3*sqrt( s2 );
            dx = xm/1000;
            x = dx:dx:xm;
            y = gammadist(x,mu,s2);
            My = max( max(y), My );
            plot( x, y ) 
            set( gca, 'XLim', [0 2*mu ] );
            
            % Plot mean and standard deviation of distributions
            e1 = sum(x.*y)*dx;
            e2 = sum(x.*x.*y)*dx;
            s = sqrt(e2-e1*e1);
            line( e1*[1 1], My/4*[0.7 1.3], 'Color', 'r' )
            line( [e1-s e1+s], My/4*[1 1], 'Color', 'r' )
            line( (e1-s)*[1 1], My/4*[1.1 0.9], 'Color', 'r' )
            line( (e1+s)*[1 1], My/4*[1.1 0.9], 'Color', 'r' )
        end
        if k == 1
            ylabel( 'sannolikhetstathet' )
        else
            set( gca, 'YTickLabel', [] )
        end
    end
    if My > 0
        for i = 1:m
            set( hand(2+i+tot*h), 'YLim', [0 1.5*My] )
        end
    end
    
    % Plot other parameters
    tstart = Params(5);
    tstop = Params(6);
    dt = Params(7);
    subplot( hand( m+3+h*tot ) )
    str = sprintf( 'Parametrar:\n %s\ntstop: %s\ndt: %s\n', int2str( Params(5) ), ...
        int2str( Params(6) ), num2str( Params(7) ) );
    midx = get( gca, 'XLim' );
    midy = get( gca, 'YLim' );
    midx = (midx(1)+midx(2))/2;
    midy = (midy(1)+midy(2))/2;
    tx1 = text( midx, midy, str, 'HorizontalAlignment', 'center', ...
		'VerticalAlignment', 'middle' );
    set( gca, 'Visible', 'off' )
    set( gca, 'FontSize', 14 )
    
    % Plot connection footprint curves
    for j = 1:len-1
        subplot( hand (m+3+j+tot*h)  )
        plot( [0 2*pi], [1 1]*mean(diagram(j).W), 'r--' )
        hold on
        N = pop( diagram(j).pre ).N;
        plot( diagram(j).angle, diagram(j).W );
        set( gca, 'XLim', [0 2*pi] )
        set( gca, 'YLim', 1.1*[minw maxw] )
        set( gca, 'YTickLabel', [] )
        str = sprintf( 'w: %s\nJp: %s\nJm: %s\n%s: %s', num2str( diagram(j).w ), num2str( diagram(j).Jp ), ...
                       num2str( diagram(j).Jm ), '\sigma', num2str( diagram(j).sgm ) );
        tx(j) = text( 0.05*pi, minw+1.05*(maxw-minw), str, 'VerticalAlignment', 'top' );
        tx(len-1+j) = text( 1.8*pi, minw+1.05*(maxw-minw), diagram(j).label, 'VerticalAlignment', ...
            'top', 'HorizontalAlignment', 'right', 'FontSize', 14, 'LineWidth', 2 );
        if j == n1-1 | j == len-1 
            xlabel( 'Vinkel (rad)' )
        end
    end
end


% Plot network activity figure
cd( dirname )
figure(6)
clf
set( 6,'Position', [ 0.5*screensize(3) 0.47*screensize(4) 0.5*screensize(3) 0.40*screensize(4) ] );
set( 6, 'Name', 'RASTER' )
set( 6, 'NumberTitle', 'off' )
pos = get( gca, 'Position' );
hgt = pos(4);
btt = pos(2);
w = pos(3);
pos(3) = w*0.7;
for i = 0:P
    pos(4) = hgt*(1-(1-spac4*0.9)*i);
    figure( 6+2*i )
    hand(1+10*i) = subplot( 'Position', pos );
end
load APs
x = APs(:,1);
y = APs(:,2);
if exist( 'Q' ) == 1 & size( Q,2 ) == 5
    [ tQ p ] = min( Q(:,1) );
else
    tQ = tstop;
end

% Plot data
pos(1) = pos(1)+0.7*w;
pos(3) = 0.3*w;
for i = 0:P
    figure( 6+2*i )
    % Plot rastergram
    subplot( hand( 1+10*i ) )
    plot( x, y, 'o', 'MarkerSize', 2 )
    grid on
    ylabel( 'cell #', 'FontSize', 14 )
    title( 'Rastergram for E- och I-cellspopulationerna', 'FontSize', 14 )
    xlabel( 'tid (ms)', 'FontSize', 14 )
    set( gca, 'FontSize', 14 )
    set( gca, 'XLim', [tstart tstop] )
    box on
    
    % plot cue
    if exist( 'Q' ) == 1 & size( Q,2 ) == 5
        ln = [];
        ln = [ ln ; line( [ Q(:,1), Q(:,1)+Q(:,2) ]', [ Q(:,3), Q(:,3) ]') ];
        ln = [ ln ; line( [ Q(:,1), Q(:,1)+Q(:,2) ]', [ Q(:,3)+Q(:,4), Q(:,3)+Q(:,4) ]') ];
        ln = [ ln ; line( [ Q(:,1), Q(:,1) ]', [ Q(:,3), Q(:,3)+Q(:,4) ]') ];
        ln = [ ln ; line( [ Q(:,1)+Q(:,2), Q(:,1)+Q(:,2) ]', [ Q(:,3), Q(:,3)+Q(:,4) ]' ) ];
        for j = 1:size( Q, 1 )
            txt = sprintf( '%s uA/mc2', num2str( Q(j,5) ) );
            text( Q(j,1), Q(j,3), txt, 'HorizontalAlignment', 'right' )
        end
        set( ln, 'Color', 'r' )
        set( ln , 'LineStyle', '-.' )
    end
    
    % create histogram subplot
    pos(4) = hgt*(1-(1-spac4*0.9)*i);
    hand(2+10*i) = subplot( 'Position', pos );
    hold on
    grid on
    set( hand(2+10*i), 'XTick', [0:2.5:65] )
    set( hand(2+10*i), 'XTickLabel', ['  ';'  ';' 5';'  ';'10';'  ';'15';'  ';'20';'  ';'25';'  ';'30';'  ';'35';'  ';'40';'  ';'45';'  ';'50';'  ';'55';'  ';'60';'  ';'65'] )
    set( hand(2+10*i), 'YAxisLocation', 'right' )
    set( gca, 'FontSize', 14 )
    xlabel( 'f (Hz)' )
    title( 'Histogram' )

    % Continue plotting rastergram and histogram
    NN = 0;
    for k = 1:nmod
      % plot lines in rastergram as well as text
      subplot( hand( 1+10*i ) )
      NI = pop(2*k-1).N;
      NE = pop(2*k).N;
      N = NE+NI;
      NN = NN + N;
      xx = find( y>=NN-N & y<NN );
      yy = y(xx);
      xx = x(xx);
      fS = length( find( yy>=NN-NE & yy<NN & xx< tQ & xx > tstart) )/(NE*(tQ-tstart)/1000);
      if exist( 'Q' ) & size(Q,2) == 5
          fD = length( find( yy>=NN-NE & yy<NN & xx>tQ + Q(p,2) ) )/(NE*(tstop-tQ+Q(p,2))/1000);
      else
          fD = 0;
      end
      fI = length( find( yy<NN-NE ) )/(NI*(tstop-tstart)/1000);
      line( [0 tstop], [NN-NE-0.5 NN-NE-0.5], 'Color', 'k', 'LineWidth', 2 )
      if NN > N
          line( [0 tstop], [NN-N-0.5 NN-N-0.5], 'Color', 'k', 'LineWidth', 2 )
      end
      tI = text( 0.1*tstop, NN-NE-0.05*NI, 'I', 'FontSize', 16, 'FontWeight', 'Bold', 'VerticalAlignment', 'top' ); 
      tE = text( 0.1*tstop, NN-0.05*NI, 'E', 'FontSize', 16, ...
		 'FontWeight', 'Bold', 'VerticalAlignment', 'top' );
      
      % Histogram
      % Show connection curve E->E
      % Find the cell with maximum activity to center connection curve on
      % that one.
      subplot( hand(2+10*i) )
      if version == 0
          ind = 1;
          xC = Ee0N( :,1 );
          yC = Ee0N( :,2 );
      else
          tmp = find( pop( 2*k ).in == 2*k ); % Find E->E curve
          if ~isempty( tmp )  
              ind(k) = tmp(1);
              arr = pop( 2*k ).ar;
              xC = diagram( arrows( arr( ind(k) ) ).diag ).angle;
              yC = diagram( arrows( arr( ind(k) ) ).diag ).W;
          else
              xC = [0 2*pi];
              yC = [0 0];
          end
      end
      hind = NN-N-1:NN;
      h = histc( y, hind );
      h=h(2:end-1)/(tstop/1000);
      hind = hind(2:end-1);

      mm = find( h(NI+1:end) == max(h(NI+1:end) ) );
      mm = mm(1); % Index of cell with maximum activity
      M = max( 20, 10*ceil( max( h )/10 ) );
      xC = xC*(NE/N)/(2*pi/N); %Adjust width of curve
      xC = xC-NE/2+mm+NN-NE-1; % Midpoint position of connection curve at maximum activity (NN-NE-1)
      ind1 = find( xC<NN-NE );
      xC( ind1 ) = xC( ind1 ) + NE;
      ind2 = find( xC>=NN );
      xC( ind2 ) = xC( ind2 ) - NE; 
      if max( yC ) > 0 
          fac = 1/max(yC)*h(floor(mm+NI));
          yC = yC*fac;
      end
      [ xC ind ] = sort( xC );
      yC = yC(ind);
    
      %plot connection curve and show mean value of that curve 
      plot( h(1:NI), hind(1:NI) )
      plot( h(NI+1:end), hind(NI+1:end) )
      plot( yC, xC, '--' )
      meanW = mean( yC );
      plot( [meanW meanW], [NN-NE NN-1], 'r--' )

      line( [0 65], [NN-NE-0.5 NN-NE-0.5], 'Color', 'k', 'LineWidth', 2 )
      if NN > N
          line( [0 65], [NN-N-0.5 NN-N-0.5], 'Color', 'k', 'LineWidth', 2 )
      end
      str = sprintf( 'fS: %.1f\nfD: %.1f\nfI: %.1f\n', fS, fD, fI );
      tx = text( 0.6*M, NN-0.05*NI, str, 'FontSize', 14, 'VerticalAlignment', 'top' );
    end
    box on
    line( [0.01*M 0.01*M], [0 NN], 'Color', 'k', 'LineWidth', 2 )
    set( gca, 'YLim', [0 NN] )
    set( hand(2+10*i), 'XLim', [0 M] )
    set( hand( 1+10*i ), 'YLim', [0 NN] )
end

% DENNA SKA VERKLIGEN GORAS OM
% Utveckling av aktivitet i figur 7, visad i en serie histogram
figure(7)
clf
set( 7,'Position', [ 0.5*screensize(3) 0.05*screensize(4) 0.5*screensize(3) 0.43*screensize(4) ] );
set( 7, 'Name', 'ACTIVITY EVOLUTION' )
set( 7, 'NumberTitle', 'off' )
t0 = [ tstart : (tstop-tstart)/no : tstop ];
tlen = diff( t0 );
t0 = sort( [ t0 t0(2:end-1) ] );
smooth = 4; % To get a smoother histogram
color = 'ybrkg';
M = 0;
h = subplot( 'Position', [ 0.13 0.2 0.75 0.67 ] );
for i = 1:no
    %histogram
    t = find( x>=t0(2*i-1) & x<t0(2*i) );
    hind = 0:smooth:netborder(end);
    if length( y(t) ) > 0
        h = histc( y(t), hind );
        h=h(2:end-1)/(tlen(i)*smooth/1000);
        M = max( M, max( h ) );
        hind = hind(2:end-1);
        hold on
        hand7(i) = plot( hind(1:netborder(end)/smooth-1), h(1:netborder(end)/smooth-1), 'Color', color(mod(i,length(color))+1) );
        grid on
    end
end
for i = 2:length( netborder )-1
    line( [netborder(i)-0.5 netborder(i)-0.5], [0 65], 'Color', 'k', 'LineWidth', 2 )
end
box on

%set( gca, 'YTick', [10 20 30 40 50] )
%set( gca, 'YTickLabel', ['  ';'10';'20';'30';'40';'50'] )
set( gca, 'XLim', [0 netborder(end)] )
set( gca, 'YLim', [0 1.3*M] )
%set( gca, 'XTick', [0:20:NN] )
%set( gca, 'XTickLabel', [] )
%set( gca, 'YAxisLocation', 'top' )
set( gca, 'FontSize', 14 )
t = text( 1.1*netborder(end), 0.5*(1.3*M), printFileName, 'Rotation', 90, 'FontSize', 14, 'HorizontalAlignment', 'center' );
title( 'Utveckling av nataktiviteten', 'FontSize', 14 )
xlabel( 'Cell index' )
ylabel( 'Frekvens (Hz)' )
strmat = [];
len = length( [ int2str( t0(end) ) int2str( t0(end-1) ) ] ) + 3;
padding = '        ';
for i = 1:2:length( t0 )
    str = [int2str(t0(i)) '-' int2str(t0(i+1)) 'ms' ];
    str = [ padding( 1:len-length(str) ) str ];
    strmat = [ strmat ; str ];
end
l = legend( hand7, strmat );

if P
    name = strcat( fileName, 'Page1.eps' ); 
    figure(8)
    set( 8, 'PaperUnits', 'centimeters' );
    set( 8, 'PaperType', 'A4');
    papersize = get( 8, 'PaperSize' );
    marg = 0.03;
    left = marg;
    width = papersize( 1 ) - 2 * left;
    bottom = marg;
    height = papersize( 2 ) - 2 * marg;
    myfiguresize = [ left, bottom, width, height ];
    set( 8, 'PaperPosition', myfiguresize );
    if pr < 3
        print
    end
    if pr > 1
        print( '-depsc2', name )
    end
    close( 8 )
end

if pr & wh > 1
    name = strcat( fileName, 'Page1.eps' ); 
    figure(7)
    set( 7, 'PaperUnits', 'centimeters' );
    set( 7, 'PaperType', 'A4');
    set( 7,'PaperOrientation','landscape');
    papersize = get( 7, 'PaperSize' );
    marg = 0.05;
    left = marg;
    width = papersize( 1 ) - 2 * left;
    bottom = marg;
    height = papersize( 2 ) - 2 * marg;
    myfiguresize = [ left, bottom, width, height ];
    set( 7, 'PaperPosition', myfiguresize );
    if pr < 3
        print
    end
    if pr > 1
        print( '-depsc2', name )
    end
end

cd( thisdir )

% Make figures lie on top of eachother
% in the right order
figure(6)
figure(5)

men = uimenu( 'Label', 'F&unctions' );
uimenu( men, 'Label', '&Power spectrum', 'Callback', 'menuCall(''corr'')');
uimenu( men, 'Label', 'Network &autocoherence', 'Callback', 'menuCall(''acoher'')');
uimenu( men, 'Label', 'Network cross co&herence', 'Callback', 'menuCall(''crcoher'')');
uimenu( men, 'Label', 'Network cross &correlation', 'Callback', 'menuCall(''crcorr'')');
uimenu( men, 'Label', '&Advanced', 'Callback', 'menuCall(''advanced'')');
uimenu( men, 'Label', '&Bold', 'Callback', 'menuCall(''bold'')');
% Since this is a function, while the menus will be called from outside
% the function, we can only pass variables to the function called from the 
% menu by saving these variables in a file called menuvar.mat. Please feel
% free to add any number of desired variables.
save menuvar APs netborder tstart tstop Q dirname