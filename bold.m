function [ t_sig, BOLD ] = bold( dirname, choice )
% This function computes the BOLD signal as a convolusion of the underlying
% neural signal and a hemodynamic response function. The hemodynamic
% response function is calculated using Karl Fristons spm_hrf. Five
% possible underlying neural activities are possible to use. These are shown
% below.
%
% This function is tightly coupled to the WM-Net and it is hard to use it
% using other signals.
%
% Fredrik Edin, 2004
% freedin@nada.kth.se


thisdir = pwd;
cd( dirname )
fs = 14; % Font size
load Q
load Params 
cd( thisdir )

% establish netborder
netborder = [0];
for i = 13:3:length(Params)
    netborder = [netborder Params(i) ];
end
netborder = cumsum( netborder );
t1 = Params( 5 );
t2 = Params( 6 );


% load underlying signal into variable signal: The BOLD signal is caused by
% one of four possible underlying signals:
% 1: Action Potentials from both inhibitory and excitatory populations
% 2: Action Potentials from excitatory population only
% 3: Sum of Excitatory Currents
% 4: Sum of Excitatory Currents onto Excitatory Cells
% 5: Sum of absolute values of all currents
% 6: Sum of all currents
% N: Number of signals
switch choice
    case {1,2}
        cd( dirname )
        load APs
        cd( thisdir )
        N = (length( netborder )-1)/2;
        dt_sig = 50; % Find dt as smallest time difference in vector
        t_sig = [ dt_sig*ceil(Params( 5 )/dt_sig)+dt_sig/2:dt_sig:Params( 6 ) ]';
        for i = 1:N
            if choice == 1
                ind1 = 2*i-1;
                ind2 = 2*i+1;
                titlestring = sprintf( 'BOLD signal from both \ninhibitory and excitatory populations' );
            elseif choice == 2
                ind1 = 2*i;
                ind2 = 2*i+1;
                titlestring = sprintf( 'BOLD signal calculated \nfrom excitatory population only' );
            end
            x = APs(:,1);
            y = APs(:,2);
            ind = find( y >= netborder(ind1) & y<netborder(ind2) );
            s = x(ind);
            signal{i} = binAPs( [ s zeros( size( s ) ) ], dt_sig, t1, t2, 1 );
            listofsignals = [];
        end
   case {3,4,5,6}
        cd( dirname )
        d = dir( 'CURRENT*-params' );
        cd( thisdir )
        N = 1;
        % 3: Sum of Excitatory Currents
        % 4: Sum of Excitatory Currents onto Excitatory Cells
        % 5: Sum of absolute values of all currents
        % 6: Sum of all currents
        if choice == 3
            allowedSynapses = [ 1 2 3 5 6 7 8 ];
            nmod = 2; % 2 = connections onto both excitatory and inhibitory cells
            titstr = 'excitatory currents';
        elseif choice == 4
            allowedSynapses = [ 1 2 3 5 6 7 8 ];
            nmod = 1; % 1 = connections onto excitatory cells only 
            titstr = 'abs( E \rightarrow E currents )';
        elseif choice == 5
            allowedSynapses = [ 0 1 2 3 4 5 6 7 8 ];
            nmod = 2; % 2 = connections onto both excitatory and inhibitory cells
            titstr = 'abs( all currents )';
        elseif choice == 6
            allowedSynapses = [ 0 1 2 3 4 5 6 7 8 ];
            nmod = 2; % 2 = connections onto both excitatory and inhibitory cells
            titstr = 'all currents';
        end
        listofsignals = [];
        for j = 2:2:length( netborder )
            signal{j/2} = [];
            for i = 1:length( d )
                cd( dirname )
                s = load( d(i).name );
                if s(1) < j & s(1) >= j-nmod & sum( find( allowedSynapses == s(2) ) ) > 0
                    % The file CURRENTx-Params indicates that this current
                    % should be included. Then load CURRENTx and add to
                    % signal.
                    str = d(i).name;
                    sig = load( str( 1:8 ) ); 
                    t_sig = sig(:,1);
                    dt_sig = round( mean( diff( t_sig ) ) );
                    if choice == 3 | choice == 4
                        sig = -sig(:,2); % The minus sign is to make excitatory currents positive
                        sig( find( sig < 0 ) ) = 0;
                    elseif choice == 5
                        sig = abs( sig(:,2) );
                    elseif choice == 6
                        sig = -sig(:,2); % The minus sign is to make excitatory currents positive
                    end
                    if isempty( signal{j/2} )
                        signal{j/2} = sig;
                    else
                        signal{j/2} = signal{j/2} + sig;
                    end
                    listofsignals = [ listofsignals ; s ];
                end
                cd( thisdir )
            end
            titlestring = sprintf( 'BOLD signal from \n%s', titstr );
        end
end

% Make sure that signal vector and hrf vector have equal time resolution
% t1 = Params( 5 );
% t2 = Params( 6 );
% if dt_sig > dt_hrf
%     t_hrf = t1:dt_hrf:t2;
%     t_sig = t1:dt_sig:t2;
%     for i = 1:length( signal )
%         signal{i} = interp1( t_sig, signal{i}, t_hrf )
%     end
%     t_sig = t_hrf;
% elseif dt_sig < dt_hrf
%     t_sig = hrf(1):dt_sig:hrf(end);
%     hrf = interp1( t_hrf, hrf, t_sig );
% end

% Compute hrf function
cd( thisdir )
[ hrf, p, t_hrf ] = spm_hrf( dt_sig/1000 ); % 1000 factors needed to switch between
t_hrf = t_hrf*1000;                                   % ms and s.


% Calculate BOLD signal
minQ = min( Q(:,1) );
Qind = find( t_sig < minQ );
zero = zeros( t_hrf(end)/dt_sig, 1 );
t_sig = [ t_sig ; [t_sig(end)+dt_sig:dt_sig:t_sig(end)+length(zero)*dt_sig]' ];
for i = 1:N
    sig = signal{i};
    mean_sig = mean( sig( Qind ) )
    sig = ( sig - mean_sig ) / mean_sig; % First subtract mean from signal

    % Pad signal vector with trailing zeros with having the length of the hrf
    % curve. Must be done after mean signal has been subtracted.
    sig = [ sig ; zero ];
    
    B = conv( sig, hrf ); % Then convolve signal
    BOLD{i} = B(1:end-length( hrf )+1,:);
end

% Convert time from ms to s
t_sig = t_sig / 1000;

% plot BOLD signal
if nargout == 0
    figure( 13 )
    set( gcf, 'Name', 'BOLD' )
    set( gcf, 'NumberTitle', 'off' )
    a = subplot(1,1,1);
    pos = get( a, 'Position' );
    delete(a)
    tit_h = 0.05;
    h = (pos(4)-(N)*tit_h)/N;
    for i = 1:N
        offset_y = pos(2)+pos(4)-(h+tit_h)*i;
        thispos = [ pos(1), offset_y, pos(3), h ];
        subplot( 'Position', thispos )
        plot( t_sig, BOLD{i} )
        legend( ['Module ' int2str(i)] )
        set( gca, 'XLim', [t_sig(1) t_sig(end)] )
        set( gca, 'FontSize', 14 )
        if i == 1
            title( titlestring, 'FontSize', fs )
        end
        if i == N
            xlabel( 'time (s)', 'FontSize', fs )
        end
        if i == ceil( N/2 )
            ylabel( 'BOLD signal (%)', 'FontSize', fs )
        end
    end
end
