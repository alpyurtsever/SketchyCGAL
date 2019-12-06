function [A, B] = qapread(fname)

FID = fopen(fname, 'r');
if (FID == -1); error('File cannot be opened.'); end

TLINE = strtrim(fgetl(FID));
DIMENSION = sscanf(TLINE, '%d');

A = zeros([DIMENSION,DIMENSION]);
for INDEX = 1:DIMENSION
    ROW = [];
    while numel(ROW) < DIMENSION
        TLINE = strtrim(fgetl(FID));
        ROW = [ROW; sscanf(TLINE, '%f')]; %#ok
    end
    A(INDEX,:) = ROW';
end
 
B = zeros([DIMENSION,DIMENSION]);
for INDEX = 1:DIMENSION
    ROW = [];
    while numel(ROW) < DIMENSION
        TLINE = strtrim(fgetl(FID));
        ROW = [ROW; sscanf(TLINE, '%f')]; %#ok
    end
    B(INDEX,:) = ROW';
end

fclose(FID);

end
%% Last edit: Alp Yurtsever - December 05, 2019