%
% nets_cellfind - find a string in an array of cells
% Steve Smith 2015
%
% cellindex = nets_cellfind(cellarray,string); 
%

function [cellindex] = nets_cellfind(cellarray,string);

cellindex=[]; % default return if nothing found

for i=1:length(cellarray)
  if strcmp(cellarray{i},string)==1
    cellindex=[cellindex i];
  end
end

if isempty(cellindex)
  for i=1:length(cellarray)
    if strncmp(cellarray{i},string,length(string))==1
      cellindex=[cellindex i];
    end
  end
end

