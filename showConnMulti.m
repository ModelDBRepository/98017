function showConnMulti( dirname, pr, wh, no, smooth )

% plots BUMP-data. 
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
    disp( 'usage1: showConn( dirname, pr ), ' )
    disp( 'pr = 0: Only show' )
    disp( 'pr = 1: Print figures 5 & 6' )
    disp( 'pr = 2: Print and save figures 5 & 6' )
    disp( 'pr = 3: Save figures 5 & 6' )
    disp( '  ' )
    disp( 'usage2: showConn( dirname, pr, wh ), ' )
    disp( 'wh = 0: As above' )
    disp( 'wh = 1: Also print and/or save figure 7' )
    disp( 'wh = 2: Print and/or save only figure 7' )
    disp( '  ' )
    disp( 'usage3: showConn( dirname, pr, wh, no ), ' )
    disp( 'no = number of enlargements of the rastergram.' )
    disp( '  ' )
    disp( 'usage4: showConn( dirname, pr, wh, no, smooth ), ' )
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

P = pr>0 & wh<2; % Determine whether to print or not. If yes, another window must be made

thisdir = pwd;
cd( dirname )

% load data files
clear Ee0N
if exist( 'Ee0N', 'file' ) % 2 population network
    version = 1;
    nmod = 1;
    NI = 32;
    NE = 128;
    load Ee0N
    load Ie0N
    load Ei0G
    load Ii0G
    one1 = ones( size( Ee0N,1 ), 1 );
    one2 = ones( size( Ii0G,1 ), 1 );
    Connections = [ one1*(NI+NE/2) (NI:NI+NE-1)' Ee0N(:,1) Ee0N(:,2) ; ...
                    one1*(NI/2) (NI:NI+NE-1)' Ie0N(:,1) Ie0N(:,2) ; ...
                    one2*(NI+NE/2) (0:NI-1)' Ei0G(:,1) Ei0G(:,2) ; ...
                    one2*(NI/2) (0:NI-1)' Ii0G(:,1) Ii0G(:,2) ];

    % To translate the parameters of the old network into the same form as
    % those of the new network
    filename = 'Params.txt';
    if exist(filename(1:end-4),'file') & exist(filename,'file')
        delete(filename(1:end-4))
    elseif exist(filename(1:end-4),'file')
        movefile(filename(1:end-4), filename);
    end
    if exist( filename, 'file' )
        load( filename )
    end
    filename = 'Q.txt';
    if exist(filename(1:end-4),'file') & exist(filename,'file')
        delete(filename(1:end-4))
    elseif exist(filename(1:end-4),'file')
        movefile(filename(1:end-4), filename);
    end
    if exist( filename, 'file' )
        load( filename )
    elseif exist( 'Params', 'var' )
        Q = [ Params( 28:29 ) ; Params( 24:26 ) ];
    else
        Q = [];
    end
    if exist( 'Params', 'var' )
        C(1,:) = [ 1 ; 1 ; 0 ; 0 ; 0 ; Params( 24:-1:23 ) ; 1 ; 0 ; Params( 8:-1:7 ) ; ...
                   Params( 5:-1:4 ) ; Params( 13:15 ) ; Params( 19:21 ) ; Params( 16:18 ) ; ...
                   Params( 10:12 ) ]';
        C(2,:) = [ 1 ; 0 ; 0 ; 0 ; 0 ; -1 ; 2 ; 0 ; 0 ; 0 ; Params(9) ; 0 ; Params(6) ;...
                   ; -ones(12,1) ]';
        Params = [ 1 ; 1 ; 3 ; -1 ; 0.02 ; 100 ; 100 ; Params( 3:-1:1 ) ];
    end
    tStart = 0;
else % The multi-module network
    filenames = {'Connections.txt','Params.txt','C.txt','Q.txt'};
    for i = 1:length(filenames)
        if exist(filenames{i}(1:end-4),'file') & exist(filenames{i},'file')
            delete(filenames{i}(1:end-4))
        elseif exist(filenames{i}(1:end-4),'file')
            movefile(filenames{i}(1:end-4), filenames{i});
        end
        if exist( filenames{i}, 'file' )
            load( filenames{i} )
        end
    end
    version = Params(1);
    if version < 5
        C = [reshape( C, 24, length( C )/24 )]';
        C = [ C(:,1:8) zeros(size(C,1),1) C(:,9:end) ];
    else
        C = [reshape( C, 25, length( C )/25 )]';
    end
    nmod = Params(2);
    if version == 2
      tStart = 0;
      Icelltyp = 1;
      tmp = [ Params(1:3) tStart Params(4:5) 100 100 ];
      for i = 1:nmod
    	  tmp = [ tmp Icelltyp 1000 Params(4+3*i) Params(3+3*i) 1000 Params(5+3*i) ];
      end
      Params = tmp;
    elseif version == 3
      tStart = Params(4);
      tmp = [ Params(1:6) ; 100 ; 100 ];
      for i = 1:nmod
          tmp = [ tmp ; Params(3+4*i) ; 1000 ; Params(5+4*i) ];
          tmp = [ tmp ; Params(4+4*i) ; 1000 ; Params(6+4*i) ];
      end
      Params = tmp;
    elseif version == 4
      tStart = Params(4);
      Params = [ Params(1:6) ; 100 ; 100 ; Params(7:end) ];
    elseif version == 5
      tStart = Params(4);
    end
end

% Name of simulationen
fid = fopen( 'Parameters' );
fileName = pwd;
f = find( fileName == '/' | fileName == '\' );
fileName = fileName( f(end)+1:end );
fileName( find( fileName == '.' ) ) = ',';

% Font size
fs = 10;

% Window placement
figs = get( 0, 'Children' );
for i = 1:length( figs ) 
    if figs(i) >= 5 & figs(i) <= 10
        close( figs(i) )
    end
end
if P
    figure(8)
    clf
    set( 8, 'Position', [ 0 100 720 852 ] )
end
figure(5)
clf
set( 5,'Position', [ 5 365 705 355 ] );


if ~exist( 'Params', 'var' ) % In this case, don't plot parameters, just show the 
                         % parameter file on the matlab prompt
    txt = char( fread( fid ) )';
    text( 0, 1, txt, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top' )
    set( gca, 'Visible', 'off' )
    version = 0;
else

% Network connection plot
rad = 0.1/nmod;
arrowind = 65; % Capital A
arrows = [];
diagram = [];

nmod = Params(2);
minw = 0;
maxw = 0; % To determine axis limits
netborder = [ 0 ; Params( 11:3:end ) ];
netborder = cumsum( netborder );
xp = 0:1/(1+nmod):1;
% Create populations (rings with I:s or E:s)
for i = 1:nmod
    if Params(3+6*i) == 1 % Type of I-cell
        str = 'I-1';
    elseif Params(3+6*i) == 2
        str = 'I-IF';
    end
    pop(2*i-1) = struct( 'label', str, 'x', xp(i+1), 'y', 0.3, 'radius', rad, 'N', ...
                         Params(5+6*i), 'ar', [], 'in',[] );
    if Params(6+6*i) == 0 
        str = 'E-1';
    elseif Params(6+6*i) == 1 % Type of E-cell
        str = 'E-3';
    elseif Params(6+6*i) == 2 
        str = 'E-IF';
    end
    pop(2*i) = struct( 'label', str, 'x', xp(i+1), 'y', 0.7, 'radius', rad, 'N', ...
                       Params(8+6*i), 'ar', [], 'in', [] );
end

% Create arrows between populations
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
            % Determine labpos, position of arrow label
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
                % Since there are several connections between the same two
                % populations, the strength of this particular connection
                % needs to be found
                nConnBtwPop = length(ind)/(pop(pre).N);
                connS = sum(reshape(Connections(ind,4), [], nConnBtwPop));
                AMPAfactor = (1+1.5*(C(i,8)==0 & mod(pre,2)==0));
                [m,i] = min(abs(connS - w*AMPAfactor));
                ind = ind(pop(pre).N*(i(1)-1)+1:pop(pre).N*i(1));
                W = Connections(ind,4)/AMPAfactor;
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
            % A negative value of pre indicates the (positive) angle of incitence to 
            % post and is used when there exists no pre.             
            arrows = [ arrows struct( 'pre', pre, 'post', post, 'label', lab, 'labpos', labpos, 'diag', length( diagram ) ) ];
            if pre > 0
                pop(pre).ar = [ pop(pre).ar length(arrows) ];
            end
            pop(post).in = [pop(post).in arrows( end ).pre ];
        end
    end
end

% Plot figure 5
len = length( diagram ) + 1;
m = sum( C(:,3) > 0 ); % There exist delay distributions
tot = (3+ m + len);
hand = zeros( tot*(P+1) ); 

pos = get( gca, 'Position' ); 
hgt = pos(4);
btt = pos(2);
for h = 0:P
    figure(5+3*h)
    pos(4) = hgt*(1-0.7*h);
    pos(2) = btt+0.7*hgt*h;
    
    % create handles
    % hand(1): title
    % hand(2): connection diagram
    % hand(3-x): delay diagrams if there exists delay distributions
    % hand(x+1 alt 3): parameter text
    % the rest: footprint plots
    spac1 = 0.03;
    spac2 = 0.01;
    spac3 = 0.05;
    spac4 = 0.6*h;
    hand(1+tot*h) = axes( 'Position', [ 0, 0.95, 1, 0.03 ] );
    hand(2+tot*h) = subplot( 'Position', [ 0.25+spac1, 2*spac1+spac4*h, 0.5-2*spac1, 0.95-3*spac1-0.3*(m>0)-spac4 ] );
    p1 = 0.25+2*spac1;
    w = 0.4-spac1;
    if m
        hand(3+tot*h) = subplot( 'Position', [ p1, 0.95+(spac2-0.3)*(1-spac4), w, (0.3-2*spac1-spac2)*(1-spac4) ] );
        p1 = p1 + w + spac1;
    end
    n1 = floor( len / 2 );
    n = n1;
    for k = 0:1 % first left subplot, then right
        if n > 0
	    height = (1.5-spac3)/n;
            for j = 1:n
                height = (0.9-spac3)/n;
                hand(m+2+j+k*n1+tot*h) = subplot( 'Position', [ 0.75*k+spac1, 0.9+(-j*height+spac2)*(1-spac4), 0.25-2*spac1, (height-2*spac2)*(1-spac4) ] );
            end
	end
        n = len - n;
    end
   
    % Plot titel
    subplot( hand(1+tot*h) )
    t0 = text( 0.5, 0.5, fileName, 'FontSize', fs, 'HorizontalAlignment', ...
	       'center' );
    set( hand(1+tot*h), 'Visible', 'off' )
    % Plot connection diagram
    cd( thisdir )
    subplot( hand(2+tot*h) )
    plotPop
    hold on
    plotAr
    set( gca, 'XTick', [] )
    set( gca, 'YTick', [] )
    if h
        axis off
    end
    
    % Plot delay distributions
    if m
        axes( hand(3+tot*h) )
        color = 'brgykm';
        coli = 1;
        labstrind = 1;
        for i = 1:size( C, 1 )
            dist = C( i, 3 );
            mu = C( i, 4 );
            if dist == 1
                line( [mu mu], [0 1], 'Color', color(coli) )
                coli = mod( coli, length(color) ) + 1;
            elseif dist == 2
                s2 = C(i, 5);
                xm = mu+3*sqrt( s2 );
                dx = xm/1000;
                x = dx:dx:xm;
        		y = gammadist(x,mu,s2);
                plot( x, y, color(coli) )
                coli = mod( coli, length(color) ) + 1;
            
        		% Plot mean and standard deviation of distributions
        		My = ylim;
        		My = My(end);
        		e1 = sum(x.*y)*dx;
        		e2 = sum(x.*x.*y)*dx;
        		s = sqrt(e2-e1*e1);
        		line( e1*[1 1], My/4*[0.7 1.3], 'Color', 'r' )
        		line( [e1-s e1+s], My/4*[1 1], 'Color', 'r' )
        		line( (e1-s)*[1 1], My/4*[1.1 0.9], 'Color', 'r' )
        		line( (e1+s)*[1 1], My/4*[1.1 0.9], 'Color', 'r' )
		
    	    end
            hold on
            if dist > 0
                syn = '';
                if C(i,10) > 0
                    syn = 'Ii,';
                end
                if C(i,11) > 0
                    syn = [syn 'Ie,'];
                end
                if C(i,12) > 0
                    syn = [syn 'Ei,'];
                end
                if C(i,13) > 0
                    syn = [syn 'Ee,'];
                end
                syn = syn(1:end-1);
                labstr{labstrind} = sprintf( '%d%s%d,comp%d%s', C(i,2), '\rightarrow', C(i,1), C(i,7), syn );
                labstrind = labstrind + 1;
            end
        end
        l = legend( gca, labstr, -1 );
        xlabel( 'tid (ms)' )
        ylabel( 'sannolikhetstathet' )
        title( 'Delayfordelningar' )
    end
    
    % Plot other parameters
    subplot( hand( m+3+h*tot ) )
    str = sprintf( 'Parametrar:\ntStart: %s\ntStop: %s\ndt: %s\n', int2str( Params(4) ), ...
        int2str( Params(5) ), num2str( Params(6) ) );
    midx = get( gca, 'XLim' );
    midy = get( gca, 'YLim' );
    midx = (midx(1)+midx(2))/2;
    midy = (midy(1)+midy(2))/2;
    tx1 = text( midx, midy, str, 'HorizontalAlignment', 'center', ...
		'VerticalAlignment', 'middle' );
    set( gca, 'Visible', 'off' )
    set( gca, 'fontsize', fs )
    
    % Plot connection curves
    for j = 1:len-1
        subplot( hand (m+3+j+tot*h)  )
        plot( [0 2*pi], [1 1]*mean(diagram(j).W), 'r--' )
        hold on
        N = pop( diagram(j).pre ).N;
        plot( diagram(j).angle, diagram(j).W );
        set( gca, 'XLim', [0 2*pi], 'XTick', [] )
        set( gca, 'YLim', 1.1*[minw maxw] )
        set( gca, 'YTickLabel', [] )
        str = sprintf( 'w: %s\nJp: %s\nJm: %s\n%s: %s', num2str( diagram(j).w ), num2str( diagram(j).Jp ), ...
                       num2str( diagram(j).Jm ), '\sigma', num2str( diagram(j).sgm ) );
        tx(j) = text( 0.05*pi, minw+1.05*(maxw-minw), str, 'VerticalAlignment', 'top', 'FontSize', 7 );
        tx(len-1+j) = text( 1.8*pi, minw+1.05*(maxw-minw), diagram(j).label, 'VerticalAlignment', ...
            'top', 'HorizontalAlignment', 'right', 'FontSize', fs, 'LineWidth', 2 );
        if j == n1-1 | j == len-1 
            xlabel( 'Vinkel (rad)' )
        end
    end
end
end

% Plot figure 6
cd( dirname )
figure(6)
clf
set( 6,'Position', [ 719 70 560 650 ] );
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
load APs.txt
x = APs(:,1);
y = APs(:,2);
if exist( 'Q', 'var' ) & size( Q,2 ) == 5
    [ tQ p ] = min( Q(:,1) );
elseif exist( 'Q' )
    tQ = 0;
    Q = [];
else
    tQ = 3000;
    Q = [];
end
tstart = Params(4);
tStop = Params(5);

% Plot data
pos(1) = pos(1)+0.7*w;
pos(3) = 0.3*w;
for i = 0:P
    figure( 6+2*i )
    % Plot rastergram
    subplot( hand( 1+10*i ) )
    plot( x, y, 'o', 'MarkerSize', 2 )
    grid on
    ylabel( 'cell #', 'FontSize', fs )
    title( 'Rastergram for E- och I-cellspopulationerna', 'FontSize', fs )
    xlabel( 'tid (ms)', 'fontsize', fs )
    set( gca, 'fontsize', fs )
    set( gca, 'XLim', [tstart tStop] )
    box on
    
    % plot queue
    if exist( 'Q', 'var' ) & size( Q,2 ) == 5
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
    
    % Plot histogram
    pos(4) = hgt*(1-(1-spac4*0.9)*i);
    hand(2+10*i) = subplot( 'Position', pos );
    hold on
    grid on
    set( hand(2+10*i), 'XTick', [0:2.5:65] )
    set( hand(2+10*i), 'XTickLabel', ['  ';'  ';' 5';'  ';'10';'  ';'15';'  ';'20';'  ';'25';'  ';'30';'  ';'35';'  ';'40';'  ';'45';'  ';'50';'  ';'55';'  ';'60';'  ';'65'] )
    set( hand(2+10*i), 'YAxisLocation', 'right' )
    set( gca, 'fontsize', fs )
    xlabel( 'f (Hz)' )
    title( 'Histogram' )

    % Continue plotting rastergram and histogram
    NN = 0;
    for k = 1:nmod
      % rastergram
      subplot( hand( 1+10*i ) )
      if version == 0
          NI = 32;
          NE = 128;
      else
          NI = pop(2*k-1).N;
          NE = pop(2*k).N;
      end
      N = NE+NI;
      NN = NN + N;
      xx = find( y>=NN-N & y<NN );
      yy = y(xx);
      xx = x(xx);
      fS = length( find( yy>=NN-NE & yy<NN & xx< tQ & xx > 500) )/(tStop/1000*NE);
      if exist( 'Q' ) & size(Q,2) == 5
          fD = length( find( yy>=NN-NE & yy<NN & xx>tQ + Q(p,2) ) )/(tStop/1000*NE);
      else
          fD = 0;
      end
      fI = length( find( yy<NN-NE ) )/(tStop/1000*NI);
      line( [0 tStop], [NN-NE-0.5 NN-NE-0.5], 'Color', 'k', 'LineWidth', 2 )
      if NN > N
          line( [0 tStop], [NN-N-0.5 NN-N-0.5], 'Color', 'k', 'LineWidth', 2 )
      end
      tI = text( 0.1*tStop, NN-NE-0.05*NI, 'I', 'FontSize', fs+2, 'FontWeight', 'Bold', 'VerticalAlignment', 'top' ); 
      tE = text( 0.1*tStop, NN-0.05*NI, 'E', 'FontSize', fs+2, ...
		 'FontWeight', 'Bold', 'VerticalAlignment', 'top' );
      tA = text( 0.9*tStop, NN-0.05*NI, ['Area ' int2str(k)], 'FontSize', fs+2, ...
		 'FontWeight', 'Bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'Right' );
      
      % Histogram
      % Show connection curve E->E
      % Find the cell with maximum activity
      subplot( hand(2+10*i) )
      if version == 0
          ind = 1;
          xC = Ee0N( :,1 );
          yC = Ee0N( :,2 );
      else
          tmp = max(find( pop( 2*k ).in == 2*k ));
          if isempty( tmp ) 
              xC = [ 0 2*pi ];
              yC = [ 0 0 ];
          else
              ind(k) = tmp;
              arr = pop( 2*k ).ar;
              xC = diagram( arrows( arr( ind(k) ) ).diag ).angle;
              yC = diagram( arrows( arr( ind(k) ) ).diag ).W;
          end
      end
      hind = NN-N-1:NN;
      h = histc( y, hind );
      h=h(2:end-1)/(tStop/1000);
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
      tx = text( 0.6*M, NN-0.05*NI, str, 'fontsize', fs, 'VerticalAlignment', 'top' );
    end
    box on
    line( [0.01*M 0.01*M], [0 NN], 'Color', 'k', 'LineWidth', 2 )
    set( gca, 'YLim', [0 NN] )
    set( hand(2+10*i), 'XLim', [0 M] )
    set( hand( 1+10*i ), 'YLim', [0 NN] )
end

% Uppforstoring i fig7
%figure(7)
%clf
%set( 7,'Position', [ 7 91 700 467 ] );
%t0 = [ tstart : (tStop-tstart)/no : tStop ];
%tlen = diff( t0 );
%t0 = sort( [ t0 t0(2:end-1) ] );
    
% subplots
% $$$ pos = get( gca, 'Position' );
% $$$ pos(4) = pos(4)*0.95;
% $$$ width = pos(3);
% $$$ spacing = 0.05;
% $$$ panelwidth = ( width - (no-1)*spacing ) / no;
% $$$ 
% $$$ cpos = pos;
% $$$ cpos(3) = 0.7*panelwidth;
% $$$ hpos = pos;
% $$$ hpos(3) = 0.3*panelwidth;
% $$$ hpos(1) = cpos(1)+cpos(3);
% $$$  
% $$$ for i = 1:no
% $$$   % rastergram
% $$$   hand2(2*i-1) = subplot( 'Position', cpos );
% $$$   t = find( x>=t0(2*i-1) & x<t0(2*i) );
% $$$   plot( x(t), y(t), 'o', 'MarkerSize', 2 )
% $$$   grid on
% $$$   line( [t0(2*i) t0(2*i)], [0 NN-1], 'Color', 'k', 'LineWidth', 2 )
% $$$   set( gca, 'fontsize', fs )
% $$$   set( gca, 'XTick', [t0(2*i-1) t0(2*i) ] )
% $$$   set( gca, 'XLim', [t0(2*i-1) t0(2*i) ] )
% $$$   set( gca, 'YLim', [ 0 NN ] )
% $$$   if i>1
% $$$     set( gca, 'YTickLabel', [] )
% $$$   else 
% $$$     ylabel( 'cell #', 'fontsize', fs )
% $$$   end
% $$$   NN = 0;
% $$$   for k = 1:nmod
% $$$     % rastergram
% $$$     subplot( hand2(2*i-1) )
% $$$     if version == 0
% $$$         NI = 32;
% $$$         NE = 128;
% $$$     else
% $$$         NI = pop(2*k-1).N;
% $$$         NE = pop(2*k).N;
% $$$     end
% $$$     N = NE+NI;
% $$$     NN = NN + N;
% $$$     line( [t0(2*i-1) t0(2*i)], [NN-NE-0.5 NN-NE-0.5], 'Color', 'k', ...
% $$$ 	  'LineWidth', 2 )
% $$$     if NN > N
% $$$         line( [0 tStop], [NN-N-0.5 NN-N-0.5], 'Color', 'k', 'LineWidth', 2 )
% $$$     end
% $$$     fE = length( find( y(t)>=NN-NE & y(t)<NN ) )/(tlen(i)/1000*NE);
% $$$     fI = length( find( y(t)<NN-NE & y(t)>= NN-N ) )/(tlen(i)/1000*NI);
% $$$     tI(i) = text( 0.1*tlen(i)+t0(2*i-1), NN-NE-0.05*NI, 'I', 'FontSize', 16, 'FontWeight', 'Bold', 'VerticalAlignment', 'top' ); 
% $$$     tE(i) = text( 0.1*tlen(i)+t0(2*i-1), NN-0.05*NI, 'E', 'FontSize', ...
% $$$ 		  16, 'FontWeight', 'Bold', 'VerticalAlignment', 'top' ); 
% $$$     str = sprintf( 'f\nfE: %.1f\nfI: %.1f\n', fE, fI );
% $$$     %tx = text( 0.6*M, NN-0.05*NI, str, 'fontsize', fs, 'VerticalAlignment', 'top' );
% $$$ 
% $$$     %histogram
% $$$     if length( y(t) ) > 0
% $$$         hand2(2*i) = subplot( 'Position', hpos );
% $$$         hind = NN-N-1:smooth:NN;
% $$$         h = histc( y(t), hind );
% $$$         h=h(2:end-1)/(tlen(i)*smooth/1000);
% $$$         M = 10*ceil( max(h)/10 );
% $$$         hind = hind(2:end-1);
% $$$         hold on
% $$$         plot( h(1:NI/smooth), hind(1:NI/smooth) )
% $$$         plot( h(NI/smooth+1:end), hind(NI/smooth+1:end) )
% $$$         grid on
% $$$         line( [0 65], [NN-NE-0.5 NN-NE-0.5], 'Color', 'k', 'LineWidth', 2 )
% $$$         if NN > N
% $$$             line( [0 65], [NN-N-0.5 NN-N-0.5], 'Color', 'k', 'LineWidth', 2 )
% $$$         end
% $$$     end
% $$$ 
% $$$   end
% $$$   box on
% $$$ 
% $$$   % Visa kopplingskurvan
% $$$   % finn maxaktivitets-cellen
% $$$   set( hand2(2*i), 'XTick', [10 20] )
% $$$   set( hand2(2*i), 'XTickLabel', ['  ';'20'] )
% $$$   set( hand2(2*i), 'YLim', [0 NN] )
% $$$   set( hand2(2*i), 'XLim', [0 M] )
% $$$   set( hand2(2*i), 'YTick', [0:20:NN] )
% $$$   set( hand2(2*i), 'YTickLabel', [] )
% $$$   set( hand2(2*i), 'XAxisLocation', 'top' )
% $$$   set( gca, 'fontsize', fs )
% $$$   
% $$$   cpos(1)=cpos(1)+panelwidth+spacing;
% $$$   hpos(1)=hpos(1)+panelwidth+spacing;
% $$$   box on
% $$$ end
% $$$ cpos(1) = cpos(1)-spacing/2;
% $$$ cpos(3) = 0.02;
% $$$ hand2(2*(no+1)) = axes( 'Position', cpos );
% $$$ t = text( 0.5, 0.5, fileName, 'Rotation', 90, 'fontsize', fs, 'HorizontalAlignment', 'center' );
% $$$ axis off
% $$$ 
% $$$ a = axes;
% $$$ pos = get( a, 'Position' );
% $$$ delete( a )
% $$$ a = axes( 'Position', pos );
% $$$ title( 'Forstoring av enskilda tidpunkter', 'fontsize', fs )
% $$$ axis off
% $$$ pos(2) = 0.2*pos(2);
% $$$ pos(4) = pos(2)+0.1;
% $$$ a = axes( 'Position', pos );
% $$$ t = text( 0.5, 0, 'tid (ms) samt frekvens (Hz)', 'fontsize', fs );
% $$$ set( t, 'HorizontalAlignment', 'center' )
% $$$ set( t, 'VerticalAlignment', 'bottom' )
% $$$ axis off;
% $$$ 
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

figure(5)
men = uimenu( 'Label', 'F&unctions' );
uimenu( men, 'Label', '&Power spectrum', 'Callback', 'menuCall(''corr'')');
uimenu( men, 'Label', 'Network &autocoherence', 'Callback', 'menuCall(''acoher'')');
uimenu( men, 'Label', 'Network cross co&herence', 'Callback', 'menuCall(''crcoher'')');
uimenu( men, 'Label', 'Network cross &correlation', 'Callback', 'menuCall(''crcorr'')');
uimenu( men, 'Label', '&Advanced', 'Callback', 'menuCall(''advanced'')');

figure(6)
% Since this is a function, while the menus will be called from outside
% the function, we can only pass variables to the function called from the 
% menu by saving these variables in a file called menuvar.mat. Please feel
% free to add any number of desired variables.
save menuvar APs netborder tstart tStop Q
