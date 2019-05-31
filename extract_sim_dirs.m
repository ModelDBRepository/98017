function cellarray = extract_sim_dirs(SERIESfile)
% EXTRACT_SIM_DIRS takes SERIES-file named SERIESfile and finds the directories
%   for the simulations in that file, extracts them and returns them in the
%   cell array named cellarray

fid = fopen(SERIESfile);
str = char([fread(fid)]');
fclose(fid);
startind = findstr(str,'C:');
endind = find((str(1:end-7)=='/' | str(1:end-7)=='\') & str(8:end)=='.')+14;
cellarray = cell(length(startind),1);
for i = 1:length(cellarray)
    cellarray{i} = str(startind(i):endind(i));
end

