function bumpanalysis( rows, cols )

% bumpanalysis reads simulation series from the file AllSeries 
% in the path (which you need to set below). It creates a menu 
% of simulation series, the rastergrams of which are displayed
% by clickning with the mouse on the menu. The program uses
% showCatalog, which needs to be in the same catalog.
% 
% rows: Number of colums per page 
% cols: Number of rows per page

if nargin < 2
    rows = 4;
    cols = 3;
end

buf = textread( '/afs/nada.kth.se/home/o/u1sxc4xo/Neuron/Program/STANDARDFILER/AllSeries', '%s', 'whitespace', '\n' );
buf{ length( buf ) + 1 } = 'Avsluta';
k = 1;
pos = get( 0, 'ScreenSize' );

while( ~strcmp( buf(k), 'Avsluta' ) )
    k = menu( 'Vilken katalog vill du lasa?', buf );
    close all

    if k < length( buf )
        showCatalog( buf{ k }, 0, rows, cols );
    end
end
