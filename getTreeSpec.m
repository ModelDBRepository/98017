function dirs = getTreeSpec(root)

% Finds the directory structure beginning with root and places
% all the directories in an array. Directories beginning with '.' 
% and directories with the name 'i686' are ignored
%
% Author: Fredrik Edin, 2004
% Address: freedin@nada.kth.se
%

thdir = pwd;
if ispc
    pathdelimiter = '\';
else
    pathdelimiter = '/';
end
if nargin > 0
    cd( root );
end
dirs = dir;
dirs(2) = [];
i = 2;
while i <= length( dirs )
    if ~(dirs(i).isdir) | dirs(i).name(1)=='.' | strcmp( dirs(i).name, 'i686' ) 
        dirs(i)=[];
    else
        i = i + 1;
    end
end
pw = pwd;
for i = 1:length( dirs )
    dirs(i).name = strcat( pw, pathdelimiter, dirs(i).name );
end

i = 2;
while i <= length( dirs )
    cd( dirs(i).name )
    
    % New dir
    tmp = dir;
    tmp(1:2) = [];
    j = 1;
    while j <= length( tmp )
        if ~tmp(j).isdir | tmp(j).name(1)=='.' | strcmp( tmp(j).name, 'i686' ) 
            tmp(j)=[];
        else
            tmp(j).name = strcat( dirs(i).name, pathdelimiter, tmp(j).name );
            j = j + 1;
        end
    end
    dirs = [ dirs ; tmp ];
    i = i + 1;
end

cd( thdir )