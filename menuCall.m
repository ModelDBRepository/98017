function menuCall( option )

load menuvar
8;
switch option
    case( 'corr' )
        figure( 9 )
        synchro( APs, netborder )
    case( 'acoher' )
        figure( 10 )
        if size( Q, 2 ) == 5 & size( Q, 1 ) > 1 & Q(1,2)+Q(1,1) < Q(2,1)
    	    kappaPlot( APs, netborder, Q(1,1)+Q(1,2), Q(2,1), 0, 0 )
        elseif size( Q, 2 ) == 5
            kappaPlot( APs, netborder, Q(1)+Q(2), tStop, 0, 0 )
        else
            kappaPlot( APs, netborder, tstart, tStop, 0, 0 )
        end
    case( 'crcoher' )
        figure( 12 )
        if size( Q, 2 ) == 5 & size( Q, 1 ) > 1 & Q(1,1)+Q(1,2)< Q(1,1)
            kappaPlot( APs, netborder, Q(1,1)+Q(1,2), Q(2,1), 1, 1 )
        elseif size( Q, 2 ) == 5
            kappaPlot( APs, netborder, Q(1)+Q(2), tStop, 1, 1 )
        else
            kappaPlot( APs, netborder, tstart, tStop, 1, 1 )
        end
    case( 'crcorr' )
        figure( 12 )
        if size( Q, 2 ) == 5 & size( Q, 1 ) > 1 & Q(1,1)+Q(1,2) < Q(2,1)
            kappaPlot( APs, netborder, Q(1,1)+Q(1,2), Q(2,1), 2 )
        elseif size( Q, 2 ) == 5
            kappaPlot( APs, netborder, Q(1)+Q(2), tStop, 2 )
        else
            kappaPlot( APs, netborder, tstart, tStop, 2 )
        end
    case( 'advanced' )
        dlgstr = {'Power spectrum', 'Network autocoherence', 'Network cross coherence', 'Network cross correlation' };
        str = {'corr', 'acoher', 'crcoher', 'crcorr'};
        [s,v] = listdlg( 'PromptString', 'Synchronization Analysis', 'SelectionMode', 'single', 'ListString', dlgstr );
    	answer = inputdlg( {'start time', 'stop time'}, 'Please fill in start and stop times', 1 );
        t1 = eval(answer{1});
        t2 = eval(answer{2});
	    if ~isempty( answer )
            Q = [ t1, 0, 0, 0, 0 ; t2, 0, 0, 0, 0 ];
            for i = 1:length(s)
                option = str{s(i)};
                save menuvar APs netborder tstart tStop Q
                menuCall( option )
            end
        end
    case( 'bold' )
        dlgstr = {'Action potentials', 'E Cell Action Potentials','Excitatory currents', 'abs( sum of currents )', 'sum of currents'};
        [s, v] = listdlg( 'PromptString', 'The BOLD activity determined by ...', 'SelectionMode', 'single', 'ListString', dlgstr );
        bold(dirname, s);
end
