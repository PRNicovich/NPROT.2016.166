function ClusterList = dbscanViaELKI(PathTo1txtFile, Epsilon, MinPts, varargin)

% Run DBSCAN algorithm through the ELKI command line interface
% Download ELKI from :
% https://elki-project.github.io/releases/
% Input path to CSV file from ThunderSTORM, 
% Function returns clustering of each point in an M x 1 vector

% Input path information for elki.jar, input data, and output folder
% pathToELKIjar = 'E:\MATLAB\ELKI\elki.jar';
pathToELKIjar = 'C:\Users\Rusty Nicovich\Documents\MATLAB\ELKI\elki-bundle-0.7.5.jar';

if exist(pathToELKIjar, 'file')
	% Continue
else
    error('Path to ELKI file invalid.\n');
end

if ischar(PathTo1txtFile)
    DataFile = PathTo1txtFile;
    OutPathBase = fileparts(PathTo1txtFile);
elseif isnumeric(PathTo1txtFile)
    % Write to temp file
    OutPathBase = pwd;
    [~, tmpName] = fileparts(tempname);
    DataFile = strcat(OutPathBase, '\', [tmpName, '.dat']);
    dlmwrite(DataFile, PathTo1txtFile, 'delimiter', '\t', 'newline', 'pc')
    
%     figure(3)
%     plot3(PathTo1txtFile(:,5), PathTo1txtFile(:,6), PathTo1txtFile(:,7), '.');
end

[~, tmpName] = fileparts(tempname);
OutPath = strcat(OutPathBase, [tmpName, '.clst']);

% Assemble API call string

callString = {'java -jar $JARPATH$ KDDCLIApplication -dbc.in "$DATAPATH$"', ...
    '-string.comment ^\s*(\D).*$', ... % Any row beginning with a letter plus or minus some white space, not a number, is ignored
    '-parser.labelIndices $LABELINDICES$', ... % Specific for 1.txt file format.
    '-algorithm clustering.DBSCAN', ...
    '-dbscan.epsilon $EPSILON$',...
    '-dbscan.minpts $MINPTS$', ...
    '-resulthandler ClusteringVectorDumper -clustering.output "$CLUSTERVECTORFILE$"'};

% 
if (nargin == 3) || ((nargin == 4) && strcmp(varargin{1}, 'Zeiss'))
    % default, go with a 1.txt file format
    if size(PathTo1txtFile, 2) == 3
    % 
        labelIndices = '4';

    else
        labelIndices = '0,1,2,3,6,7,8,9,10,11,12,13'; % Specific for 1.txt file format.  Position data cols 5-6 (4-6 zero-indexed)
    end
elseif (nargin == 4) && ~strcmp(varargin{1}, 'Zeiss')
    if strcmp(varargin{1}, 'ThunderSTORM')
        if size(PathTo1txtFile, 2) == 3
            labelIndices = '4';
        else
            labelIndices = '0,3,4,5,6,7,8'; % Specific for ThunderSTORM CSV file format.  Position data cols 2-3 (1-2 zero-indexed)
        end
    elseif strcmp(varargin{1}, 'XY');
        labelIndices = '2,3,4'; % Data input is 2 column matrix.  labelIndices is placeholder.
    end
end

% Use strrep to replace variable place-holders with variable values
executeString = strjoin(callString, ' ');
executeString = strrep(executeString, '$JARPATH$', sprintf('"%s"', pathToELKIjar));
executeString = strrep(executeString, '$LABELINDICES$', labelIndices);
executeString = strrep(executeString, '$DATAPATH$', DataFile);
executeString = strrep(executeString, '$CLUSTERVECTORFILE$', OutPath);
executeString = strrep(executeString, '$EPSILON$', num2str(Epsilon));
executeString = strrep(executeString, '$MINPTS$', num2str(MinPts));

% disp(executeString)

% Call formatted string
% Remove '-echo' to keep things from displaying in comnand window

if ispc
    [~,~] = dos(executeString, '-echo');
elseif ismac
    [~,~] = system(executeString,'-echo');
end


% Read in temp file
fid = fopen(OutPath, 'r');
clusterResults = fread(fid, 'uint8=>char');
fclose(fid);
delete(OutPath);

if isnumeric(PathTo1txtFile)
    delete(DataFile);
end

clusterResults = strsplit(clusterResults', ' ');
clusterResults(~cellfun(@all, isstrprop(clusterResults, 'digit'))) = [];
ClusterList = str2double(clusterResults);
    
    
