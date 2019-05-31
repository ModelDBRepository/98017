function compareCatalogs(homecat)

% Put the contents of remotecat in a file called Catalogs.txt with 
% each entry in the catalog on a separate row.
% Then the program will print all the entries that are the same in the two
% catalogs so that one avoids overwriting.

remote = textread('Catalogs.txt','%s');
home = dir(homecat);
home(1:2) = [];
home = struct2cell(home);
home = home(1,:);
intsect = intersect(home,remote);
disp(char(intsect))

