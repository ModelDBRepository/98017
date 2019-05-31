function showBig( catalogName )

% showBig is intended to do an enlargement of a simulation if the
% user clicks with the mouse on the smaller picture of the
% simulation. Used together with bumpanalysis.m and showCatalog.m.


thisdir = pwd;
% Display on screen if left button is pressed
disp( get( gcf, 'SelectionType' ) )
if strcmp( get( gcf, 'SelectionType' ), 'normal' )
    showConnMulti( catalogName, 0 )
    disp( 'DISPLAYING DATA' )

% Print if right button is pressed
elseif strcmp( get( gcf, 'SelectionType' ), 'alt' )
    showConnMulti( catalogName, 1 )
    disp( 'PRINTING DATA' )

% Save to file if middle button or both left and right mouse
% buttons are pressed
elseif strcmp( get( gcf, 'SelectionType' ), 'extend' )
    showConnMulti( catalogName, 3 )
    disp( 'SAVING DATA' )
end
cd( thisdir );
