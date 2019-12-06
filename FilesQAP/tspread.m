function [A, B, info] = tspread(fname)

info = struct();

FID = fopen(fname, 'r');
if (FID == -1); error('File cannot be opened.'); end

while true
    
    TLINE = strtrim(fgetl(FID));
    if strcmp('EOF', TLINE) || ~ischar(TLINE), break; end
    
    [TOKEN, REMAIN] = strtok(TLINE, ':');
    TOKEN = strtrim(TOKEN);
    REMAIN = strtrim(REMAIN(2:end));
    
    switch TOKEN
        case 'NAME'
            NAME = REMAIN;
            info.NAME = NAME;
            
        case 'COMMENT'
            COMMENT = REMAIN;
            info.COMMENT = COMMENT;
            
        case 'TYPE'
            TYPE = REMAIN;
            info.TYPE = TYPE;
            
        case 'DIMENSION'
            DIMENSION = sscanf(REMAIN, '%d');
            info.DIMENSION = DIMENSION;
            
        case 'EDGE_WEIGHT_TYPE'
            EDGE_WEIGHT_TYPE = REMAIN;
            info.EDGE_WEIGHT_TYPE = EDGE_WEIGHT_TYPE;
            
        case 'EDGE_WEIGHT_FORMAT'
            EDGE_WEIGHT_FORMAT = REMAIN;
            info.EDGE_WEIGHT_FORMAT = EDGE_WEIGHT_FORMAT;
            
        case 'NODE_COORD_SECTION'
            for INDEX = 1:DIMENSION
                DATA = sscanf(fgetl(FID), '%f');
                NODE_COORD_SECTION(INDEX, :) = DATA(2:end)'; %#ok
            end
            info.NODE_COORD_SECTION = NODE_COORD_SECTION;
            
        case 'EDGE_WEIGHT_SECTION'
            switch EDGE_WEIGHT_FORMAT
                case 'FULL_MATRIX'
                    A = zeros([DIMENSION,DIMENSION]);
                    for INDEX = 1:DIMENSION
                        DATA = sscanf(fgetl(FID), '%f');
                        A(INDEX, :) = DATA';
                    end
                case 'UPPER_COL'
                    A = zeros(DIMENSION,DIMENSION);
                    for SUBJ = 2:DIMENSION
                        DATA = sscanf(fgetl(FID), '%f');
                        A(1:(SUBJ-1),SUBJ) = DATA;
                        A(SUBJ,1:(SUBJ-1)) = DATA';
                    end
                case 'UPPER_ROW'
                    A = zeros(DIMENSION,DIMENSION);
                    for SUBI = 1:(DIMENSION-1)
                        DATA = sscanf(fgetl(FID), '%f');
                        A(SUBI,(SUBI+1):end) = DATA';
                        A((SUBI+1):end,SUBI) = DATA;
                    end
                case 'LOWER_COL'
                    A = zeros(DIMENSION,DIMENSION);
                    for SUBJ = 1:(DIMENSION-1)
                        DATA = sscanf(fgetl(FID), '%f');
                        A((SUBJ+1):end,SUBJ) = DATA;
                        A(SUBJ,(SUBJ+1):end) = DATA';
                    end
                case 'LOWER_ROW'
                    A = zeros(DIMENSION,DIMENSION);
                    for SUBI = 2:DIMENSION
                        DATA = sscanf(fgetl(FID), '%f');
                        A(SUBI,1:(SUBI-1)) = DATA';
                        A(1:(SUBI-1),SUBI) = DATA;
                    end
                case 'LOWER_DIAG_COL'
                    DATA = [];
                    A = nan(DIMENSION,DIMENSION);
                    for SUBJ = 1:DIMENSION
                        for SUBI = SUBJ:DIMENSION
                            if isempty(DATA)
                                DATA = sscanf(fgetl(FID), '%f');
                            end
                            A(SUBI,SUBJ) = DATA(1);
                            A(SUBJ,SUBI) = DATA(1);
                            DATA(1) = [];
                        end
                    end
                case 'LOWER_DIAG_ROW'
                    DATA = [];
                    A = nan(DIMENSION,DIMENSION);
                    for SUBI = 1:DIMENSION
                        for SUBJ = 1:SUBI
                            if isempty(DATA)
                                DATA = sscanf(fgetl(FID), '%f');
                            end
                            A(SUBI,SUBJ) = DATA(1);
                            A(SUBJ,SUBI) = DATA(1);
                            DATA(1) = [];
                        end
                    end
                case 'UPPER_DIAG_COL'
                    DATA = [];
                    A = nan(DIMENSION,DIMENSION);
                    for SUBJ = 1:DIMENSION
                        for SUBI = 1:SUBJ
                            if isempty(DATA)
                                DATA = sscanf(fgetl(FID), '%f');
                            end
                            A(SUBI,SUBJ) = DATA(1);
                            A(SUBJ,SUBI) = DATA(1);
                            DATA(1) = [];
                        end
                    end
                case 'UPPER_DIAG_ROW'
                    DATA = [];
                    A = nan(DIMENSION,DIMENSION);
                    for SUBI = 1:DIMENSION
                        for SUBJ = SUBI:DIMENSION
                            if isempty(DATA)
                                DATA = sscanf(fgetl(FID), '%f');
                            end
                            A(SUBI,SUBJ) = DATA(1);
                            A(SUBJ,SUBI) = DATA(1);
                            DATA(1) = [];
                        end
                    end
                otherwise
                    error('Unknown EDGE_WEIGHT_FORMAT.');
            end
                    
        case 'DISPLAY_DATA_TYPE'
            DISPLAY_DATA_TYPE = REMAIN;
            info.DISPLAY_DATA_TYPE = DISPLAY_DATA_TYPE;
            
        case 'DISPLAY_DATA_SECTION'
            for INDEX = 1:DIMENSION
                DATA = sscanf(fgetl(FID), '%f');
                DISPLAY_DATA_SECTION(INDEX, :) = DATA(2:end)'; %#ok
            end
            info.DISPLAY_DATA_SECTION = DISPLAY_DATA_SECTION;
            
        otherwise
            warning(['Unknown option ', TOKEN, '.']);
    end
end
fclose(FID);

% We handle different edge weight types below
if exist('NODE_COORD_SECTION','var')
    switch EDGE_WEIGHT_TYPE
        case 'EXPLICIT'
            % DONE ALREADY IN "EDGE_WEIGHT_SECTION" PART ABOVE
        case {'EUC_2D','EUC_3D'}
            A = round(squareform(pdist(NODE_COORD_SECTION,'euclidean')));
        case {'MAN_2D','MAN_3D'}
            A = round(squareform(pdist(NODE_COORD_SECTION,'cityblock')));
        case {'MAX_2D','MAX_3D'}
            A = round(squareform(pdist(NODE_COORD_SECTION,'chebychev')));
        case 'CEIL_2D'
            A = ceil(squareform(pdist(NODE_COORD_SECTION,'euclidean')));
        case 'GEO'
            A = squareform(pdist(NODE_COORD_SECTION,@geodistance));
        case 'ATT'
            A = squareform(pdist(NODE_COORD_SECTION,@attdistance));
        case {'XRAY1','XRAY2'}
            error('Edge weight types XRAY1 and XRAY2 are not implemented.');
            % to be implemented later
        otherwise
            error('Unknown edge weight type.');
    end
end

% Create B (scaled canonical tour matrix)
n = size(A,1);
B = spdiags(ones(n,1),1,n,n);
B(1,n) = 1;
B = 0.5*(B+B');

end

function ARCLEN = geodistance(x,X)

PI = 3.141592;

degs = round(x(1));
mins = x(1) - degs;
xLat = PI * (degs + 5 * mins / 3 ) / 180;

degs = round(x(2));
mins = x(2) - degs;
xLon = PI * (degs + 5 * mins / 3 ) / 180;

degs = round(X(:,1));
mins = X(:,1) - degs;
XLat = PI * (degs + 5 * mins / 3 ) / 180;

degs = round(X(:,2));
mins = X(:,2) - degs;
XLon = PI * (degs + 5 * mins / 3 ) / 180;

RRR = 6378.388;
q1 = cos( xLon - XLon );
q2 = cos( xLat - XLat );
q3 = cos( xLat + XLat );
ARCLEN = ceil( RRR * acos( 0.5*((1+q1).*q2 - (1-q1).*q3) ));

end

function dij = attdistance(x,X)
% x is 1 by 2 coordiantes of one point
% X is N by 2 coordiantes of N points

xdiff = x(1) - X(:,1);
ydiff = x(2) - X(:,2);
rij = sqrt( (xdiff.^2 + ydiff.^2)/10 );
dij = round(rij);
dij(dij<rij) = dij(dij<rij) + 1;

end
%% Last edit: Alp Yurtsever - December 05, 2019