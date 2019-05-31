function [ txt, pos ] = readUntil( fi, str )

% readUntil reads from ascii file fi until it finds the pattern str.
% It positions the pointer at the beginning of the pattern and
% returns the string txt containing all characters between the original
% position in the file and the new position. It also returns the position 
% of the pointer pos. If no pattern is found, readUntil positions the 
% pointer where it was at the call of the function and returns txt = '', 
% pos = -1.

txt = '';
pos = ftell( fi );
buf = char( fread( fi )' );
p = [];
for i = 1:length(str)
    p = [p findstr( buf, str(i) )];
end
if isempty( p )
    fseek( fi, pos, 'bof' );
    pos = -1;
else
    pos = pos + p(1) - 1;
    fseek( fi, pos, 'bof' );
    txt = buf( 1:p(1)-1 );
end
txt;