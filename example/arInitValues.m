%%%% Read a txt file and execute it ... in the setup

% Open the text file in read mode
fileID = fopen('neu_initValues.txt', 'r');

% Check whether the file was opened successfully
if fileID == -1
    error('The file could not be opened.');
end

% Read the file content line by line
fileContent = {};
while ~feof(fileID)
    line = fgetl(fileID);
    fileContent = [fileContent, line];
end

% Close the file
fclose(fileID);

% Execute the read code
for i = 1:length(fileContent)
    try
        eval(fileContent{i});
    catch
        fprintf('Error executing line %d\n', i);
    end
end
