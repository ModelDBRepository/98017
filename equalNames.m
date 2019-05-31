function equalNames(filename)

% Copy the names of all files to a file. This is done in putty/dayhoff
% by ls, then marking the text in the terminal window, then opening emacs
% in the terminal window and right clicking to paste the text into the 
% document. Next run this document to see whether any filenames are the
% same.

fid = fopen(filename);
str = char(fread(fid))';
i = 1;
while length(str) > 0
    [token, str] = strtok(str);
    filelist{i} = token;
    i = i + 1;
end
if isempty(filelist{end})
    filelist(end) = [];
end

filelist = sort(filelist);
for i = 1:length(filelist)
    for j = i+1:length(filelist)
        if strcmp(filelist{i}, filelist{j})
            disp(['strings same: a)' filelist{i} ', b)' filelist{j}])
        end
    end
end

disp('Strings compared')
for i = 1:length(filelist)
    disp(filelist{i})
end
