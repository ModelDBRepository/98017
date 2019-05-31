function dirs = findSims(key, pos, start)

% This function goes through all directories to find simulations with
% parameters at index pos matching those in key. Uses function
% getTreeSpec.m, which must be in the same directory.
%
% key is a vector containing a pattern of numerical values
% pos is a vector of indices
% start is the directory where you want to start your search
%
% Author: Fredrik Edin, 2004
% Address: freedin@nada.kth.se

thisdir = pwd;
if nargin < 3
    start = '~/Neuron/Program/MySim';
end
if isempty( find( start == '~' ) )
    start = strcat( thisdir, '/', start );
end
homepw = '~/Neuron/Program/STANDARDFILER/'; % The directory
                                            % where this file exists

% Find all directories in which to search
cd( homepw )
tree = getTreeSpec( start );

% Find directories where pattern matches
dirs = [];
for i = 1:length( tree )
    cd( tree(i).name )
    if exist( 'Params' ) > 1
        load( 'Params' )
        stammer = 1;
        for j = 1:length( pos )
            if key( j ) ~= Params( pos( j ) )
                stammer = 0;
                break
            end
        end
        if stammer
            dirs = [ dirs ; tree(i) ];
        end
        clear Params
    end
end

% return to where you started
cd( thisdir )
