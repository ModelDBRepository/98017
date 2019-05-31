function netChar3D( filename, DoPrint, lab, meas, rows, cols )

% netChar3D takes a file name beginning with netCharSERIE. This
% file contains names of catalogues with simulation series. These
% series can have arbitrary names as long as they contain the word
% SERIES. Every catalogue is represented by a four-dimensional
% index based on the following criteria:
% drift, measured as the standard deviation of a diffusion process
%        of the population vector
% distractibility, measured as total current needed to extinguish
% the bump or move it more than halfway towards the distracting
% stimulus
% variance of the spontaneous activity
% variance of the delay activity

% remember to change path names to your destination directory.
% those places are marked by % CHANGE %

% filename: name of the netCharSERIE file
% DoPrint = 0, only display
% DoPrint = 1, print
% lab = 0 or nothing: Incrementing label
% lab = 1: Specified label
% meas: A vector with the indices to be calculated. If the vector is
%           empty, all measures will be calculated.
% rows: Number of rasterplot rows per page
% cols: Number of rasterplot columns per page

% Author: Fredrik Edin, 2004
% email: freedin@nada.kth.se

if nargin < 1 | nargin > 6
    disp( 'usage 1: netChar3D( filename, DoPrint )' )
    disp( 'filename: Name of the netCharSERIE file' )
    disp( 'DoPrint = 0: Display only' )
    disp( 'DoPrint = 1: Print and display' )
    disp( ' ' )
    disp( 'usage 2: netChar3D( filename, DoPrint, lab )' )
    disp( 'lab = 0: Incrementing integers as labels' )
    disp( 'lab = 1: Prespecified labels from filename' )
    disp( ' ' )
    disp( 'usage 3: netChar3D( filename, DoPrint, lab )' )
    disp( 'meas: A vector with indices to be calculated' )
    disp( ' ' )
    disp( 'usage 4: netChar3D( filename, DoPrint, lab, rows, cols )' )
    disp( 'rows: Number of rasterplot rows per page' )
    disp( 'cols: Number of rasterplot columns per page' )
    disp( ' ' )
    return
end
if nargin < 4
    meas = [ 1:6 ];
elseif isempty( meas )
    meas = [ 1:6 ];
end
if nargin < 6
    rows = 4; % Rows of plots per page
    cols = 3;
end

homedir = pwd;
filesdir = '/afs/nada.kth.se/home/o/u1sxc4xo/Neuron/Program/STANDARDFILER'; % CHANGE %
if nargin < 3
    lab = 0
end
if lab
    [ filer, beg, fin, labv ] = textread( filename, '%s %d %d %s', 'commentstyle', 'shell' );
else
    [ filer, beg, fin ] = textread( filename, '%s %d %d', 'commentstyle', 'shell' );
end
len = length( filer );

% How many pages? Position of windows, subplots, etc.
pages = 2 * ( floor( (len-1) / (rows*cols) ) + 1 );
for i = 1 : pages/2
    figure ( i )
    clf
    set( i,'Position', [ 550 * mod( ( i ), 2 ) 100 720 852 ] );
end
figure( pages+10 )
clf
handl = zeros( 1, len*2 ); % subplot handles
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
    handl( ind ) = subplot( 'Position', [ left bottom width-0.05 height-0.07 ] );
    figure( page+pages/2 )
    if subp == 1
        clf
    end
    handl( ind+len ) = subplot( 'Position', [ left bottom width-0.05 height-0.07 ] );
end


% Below are the different measures for determining the performance of each
% network. A file with filename containing the string netCharSERIE contains
% data for the different points, see that file for more information.
% Several simulations can be used to characterize a single simulation, and
% which simulations is indicated in the entCharSERIE file.
INDEX = [];
errv = [];
labt = [];
images = [];
for kat = 1:length( filer )
    savefile = filer{kat};
    last = find( savefile == '/' );
    savedir = savefile(1:last(end));
    cd( savedir )
    savefile = strcat( 'netCharINDEX-to-', savefile(last(end)+1:end),':',int2str(beg(kat)),'-',int2str(fin(kat)),'.mat');
    if exist( savefile ) > 1
        % Network already characterized, no need to recalculate
        load( savefile );
        INDEX = [ INDEX NetCharINDEX ]
        errv = [ errv NetCharINDEX.err ];
        binwfinal = NetCharINDEX.binwf;
        binw = NetCharINDEX.binw;
        images = [ images NetCharINDEX.im ];
        drfinal = NetCharINDEX.drf;
        dr = NetCharINDEX.dr;
    else

        % Parse text from simulation series file.
        serie = [];
        txt = [ textread( filer{kat}, '%s', 'commentstyle', 'shell', 'delimiter', '\n' ) ]';
        ind = 1;
        while ind <= length( txt )
            posb = find( txt{ind} == '/' );
            if isempty( posb )
                txt(ind) = [];
            else
                ind = ind+1;
            end
        end
        if isempty( txt )
            disp( 'Error: Your file does not contain any series' )
            return
        end 
        if( beg( kat ) == 0 | fin(kat) == 0 )
            beg(kat) = 1;
            fin(kat) = length( txt );
        end    
        for ind = beg(kat):fin( kat )
            tx = txt{ ind };
            posb = find( tx == '/' );
            pose = length( tx );
            serie = [ serie ; struct( 'name', tx(posb(1):pose ) ) ];
            disp( [ tx(posb(1):pose ) filer{kat} ] )
        end
        if length( serie ) == 0
            disp( 'No data in this catalog' )
            return
        end
    
        delta = [];
        ampwidth = [];   
        distS = 0;
        s2S = zeros( length( serie ), 1 );
        for fil = 1:length( serie )
            simdir = serie(fil).name;
            cd( simdir )
            
            if fil == 1
                % Index no 1: Drift as a diffusion process
                % The population vector is measured for time interval
                % t(i)-t(i+1). Drift is measured as mean square drift. To
                % find the correct time bin, see where drift levels off and
                % pick that value.
                load APs
                x = APs(:,1);
                y = APs(:,2);
                load Q
                load Params
		if Params(1) == 2
		  tStart = 0;
		  Icelltyp = 1;
		  tmp = [ Params(1:3) tStart Params(4:5) ];
		  for i = 1:nmod
		    tmp = [ tmp Icelltyp Params(3+3*i:5+3*i) ];
		  end
		  Params = tmp;
		elseif version == 3
		  tStart = Params(4);
		end
		cd( filesdir )
                begin = Q(1)+Q(2);
                finish = Params(5);
                NN = 0;
                ind = [];
                ang = [];
                for j = 10:4:length( Params ) 
                    N = sum(Params(j-1:j));
                    ind = [ ind'  [ NN+Params(j-1) : NN+N-1 ] ]';
                    ang = [ ang ; [ 0:2*pi/Params(j):2*pi*(1-1/Params(j))]' ];
                    NN = NN + N;
                end
                ang = [ ind ang ];
                dr = [];
                b = 1:11;
                binw =2.^b;
                figure(pages+2)
                clf
                for bb = 1:length(b)
    		        subplot( length(b), 1, bb )
                    dr = [ dr driftdiffusion( x, y, ang, binw(bb), [ begin finish ] ) ];
                end
                ind = find( dr == -1 );
                dr(ind) = [];
                binw(ind) = [];
                ddr = diff( dr );
                rddr = abs( ddr./dr(1:end-1) );
                ind = find( rddr == min( rddr ) );
                binwfinal = binw( ind );
                if isempty( dr )
                    drfinal = -1;
                else
                    drfinal = dr(ind);
                end
                
                % Index no 5: Drift measured by linear regression. Since we're measuring 
                % from a ring, we need to first find a proper cut in the ring
                yang = zeros( size( y ) );
                for j = 1:size( ang, 1 )
                    angind = find( y == ang( j, 1 ) ); 
                    if ~isempty( angind )
                        yang( angind ) = ang( j, 2 );
                    end
                end
                [ tmp, cut ] = popvec(  x, y, ang,  [ begin finish ] );
                tmp = find( yang > cut );
                yang( tmp ) = yang( tmp ) - 2 * pi; 
                A = [ ones( size( x ) ), x ];
                linrg = A\yang;
                
                
                % Index no 2 and no 3: Variance of spontaneous activity and
                % bump activity
                window = 500;
                [  tmp, s2, vars ] = isbump( x, y, [ begin finish ], ind, window, Q(1) );
                s2S(fil) = vars;

                
                
            else
                %Distractability
                %1. All simulations in a series
                %    For every simulation, do
                %2. Measure population vector last 500ms before distracter
                %3. Add distracter 180 degrees away from bump
                %4. During analysis, check popvec 250 ms before cue. Chect
                %    popvec one more time 1 second later. Store difference delta.
                %    Find simulation where smallest current makes delta >
                %    90 degrees. This is index of distractability.
                load APs
                load Params
		if Params(1) == 2
		  tStart = 0;
		  Icelltyp = 1;
		  tmp = [ Params(1:3) tStart Params(4:5) ];
		  for i = 1:nmod
		    tmp = [ tmp Icelltyp Params(3+3*i:5+3*i) ];
		  end
		  Params = tmp;
		elseif version == 3
		  tStart = Params(4);
		end
                load X
                cd( filesdir )
                finish = Params(5);
                for j = 10:4:length( Params ) 
                    N = sum(Params(j-1:j));
                    ind = [ ind ; [ NN +Params(j-1) : NN + N - 1 ]' ];
                    ang = [ ang ; [ 0:2*pi/Params(j):2*pi(1-1/Params(j))]' ];
                    NN = NN + N;
                end
                x = APs(:,1);
                y = APs(:,2);
                isAlive = 1;
                window = 500;
                tid = [ X(2)-window X(2) ];
                [  isAlive, tmp2, vars ] = isbump( x, y, tid, ind, window, Q(1) )
                s2S(fil) = vars;
                if isAlive
                    window = 500;
                    tid = [ X(2)-window X(2) X(2)+X(3)+1010-window X(2)+X(3)+1010 ]';
                    [ distIsAlive, tmp2, vars ] = isbump( x, y, tid(3:4), ind, window, Q(1), vars )
                    if distIsAlive
                        popv = zeros( 2, 1 );
                        popv(1) = popvec( x, y, ang, tid(1:2) );
                        popv(2) = popvec( x, y, ang, tid(3:4) );
                        delta = [ delta min( diff( popv ), 2*pi-diff(popv) ) ];
                        ampwidth = [ ampwidth X(5)*X(4) ];
                    else
                        delta = [ delta pi ];
                        ampwidth = [ ampwidth X(5)*X(4) ];
                    end
               end
            end
        end
        
        % Find weakest effective distracter. Err is one if either all
        % distracters are unable to distract or all distracters are able.
        % Then no threshold is found.
        dind = find( delta >= pi/2 );
	    dind2 = find( delta < pi/2 );
        err = 0;
        if isempty( dind ) | isempty( dind2 ) 
            if isempty( delta )
                dind = -1;
            elseif isempty( dind )
                dind = length( delta(end) );
            end
            err = 1;
            images = [ images 1 ];
       else
            images = [ images dind(1) ];
        end
        errv = [ errv err ];
        if dind(1)>0
            ind2 = ampwidth(dind(1));
        else
            ind2 = -1;
        end
        
        % Putting together index.
        ind1 = drfinal; % Drift as diffusion process
        ind3 = mean( s2S ); % Spontaneous activity
        ind4 = s2; % Bump activity
        ind5 = linrg(2); % Drift using linear regression
        position = [ ind1 ind2 ind3 ind4 ind5 ];
        num = kat;
        % Save index for this nat
        cd( savedir )
        if lab
            INDEX = [ INDEX struct( 'num', num, 'pos', position, 'im', images(kat), 'serie', serie, 'err', err, 'binwf', binwfinal, 'binw', binw, 'drf', drfinal, 'dr', dr, 'lab', labv( kat ) ) ];
        else
            INDEX = [ INDEX struct( 'num', num, 'pos', position, 'im', images(kat), 'serie', serie, 'err', err, 'binwf', binwfinal, 'binw', binw, 'drf', drfinal, 'dr', dr, 'lab', int2str( num ) ) ];
        end
        NetCharINDEX = INDEX(kat);
        save( savefile, 'NetCharINDEX' );
        im = images(kat);
    end
    
    % Plotting
    % Choose a simulation with a threshold distracter, if there is one. 
    % Add index with a number in a 3D-plot. Save both pictures in this
    % folder.
    cd( savedir )
    cd( INDEX(kat).serie(INDEX(kat).im).name )
    load Q
    load Params
    load APs
    x = APs(:,1);
    y = APs(:,2);
    page = floor( (kat-1) / ( rows*cols ) ) + 1;
    subp = kat - ( page-1 ) * rows * cols;
    figure( page ) 
    subplot( handl( kat ) )
    % specify call-back routine so that one can click on a picture to
    % get a big display. 
    str = strcat( 'showBig(  ''serie(ind).name,'' )' );
    set( handl( kat ), 'ButtonDownFcn', str );

    % plot the rastergrams of the cell
    str = sprintf( 'Net # %s', INDEX(kat).lab );
    title( str, 'Fontsize', 14 )
    hold on
    plot( x, y, '.', 'MarkerSize', 1 )
    set( gca, 'box', 'on' )
    if Params(1) == 1
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
        fS = length( find( y(xQ)>=NI ) )/(tQ/1000*NE); 
        fE = length( find( y>=NI ) )/(tStop/1000*NE);
        fI = length( find( y<NI ) )/(tStop/1000*NI);
    else
        if Params(1) == 2
	  tStart = 0;
	  Icelltyp = 1;
	  tmp = [ Params(1:3) tStart Params(4:5) ];
	  for i = 1:nmod
	    tmp = [ tmp Icelltyp Params(3+3*i:5+3*i) ];
	  end
	  Params = tmp;
	elseif version == 3
	  tStart = Params(4);
	end
        N = 0;
        tStop = Params(5);
        tQ = Q(1);
        xQ = find( x<tQ );
        for i = 1:Params(2)
            NI = Params(5+4*i);
            NE = Params(6+4*i);
            ipos(i) = N+NI;
            epos(i) = N+NE+NI;
            xpos(i) = N+NI;
            fS(i) = length( find( y(xQ)>=NI ) )/(tQ/1000*NE); 
            fE(i) = length( find( y>=NI ) )/(tStop/1000*NE);
            fI(i) = length( find( y<NI ) )/(tStop/1000*NI);
            line( [0 tStop], [N+NI-0.5 N+NI-0.5], 'Color', 'k', 'LineWidth', 2 )
            if i>1
                line( [0 tStop], [N-0.5 N-0.5], 'Color', 'k', 'LineWidth', 2 )
            end
            N = N+NI+NE;
        end
    end
        
    set( gca, 'YLim', [0 N] )
    set( gca, 'FontSize', 14 )
    for i = 1:length(ipos)
        tI(i) = text( 0.1*tStop, ipos(i), 'I', 'FontSize', 16, 'FontWeight', 'Bold', 'VerticalAlignment', 'top' ); 
        tE(i) = text( 0.1*tStop, epos(i), 'E', 'FontSize', 16, 'FontWeight', 'Bold', 'VerticalAlignment', 'top' );
        str = sprintf( 'fE: %.1f\nfS: %.1f\nfI: %.1f\nbin: %d\n', fE(i), fS(i), fI(i), binwfinal );
        tx(ind) = text( 0.1*tStop, xpos(i), str, 'FontSize', 12, 'VerticalAlignment', 'bottom' );
    end
    if mod((subp-1), cols )
        set( gca, 'YTickLabel', [] )
    end
    
    % plot drift as function of step size
    if ~isempty( dr )
        page = pages/2+floor( (kat-1) / (rows*cols) ) + 1
        figure( page )
        subplot( handl(kat+len) )
        semilogx( binw, dr )
        str = sprintf( 'Net # %s', INDEX(kat).lab );
        title( str, 'Fontsize', 14 )
        set( gca, 'YLim', [ 0 max( dr ) ] )
    end
    
    % plot the 3D-plot. Choose any index you wish
    figure( pages+10 )
    px = INDEX(kat).pos(1);
    py = INDEX(kat).pos(4);
    pz = INDEX(kat).pos(3);
    if errv(kat)
        plot3( px, py, pz,'ro', 'Markersize', 5 );
    else
        plot3( px, py, pz, 'x', 'MarkerSize', 5); %fixa label
    end
    hold on
    grid on
    title( 'Scatterplot av 1-modulsnatverks performance' )
    xlabel( 'drift (celler/sek)' )
    ylabel( 'varians av delayaktivitet (Hz2)' )
    zlabel( 'varians av spontanaktivitet (Hz2)' )
    labt = [ labt text( px, py, pz, INDEX(kat).lab) ]; 
    clear savedir savefile
end


% Print or save pictures
if DoPrint
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
        print
    end
    figure( pages+10 )
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


   cd( homedir )