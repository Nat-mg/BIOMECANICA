function Data= dataImport (filedir)
% Import data from the txt file generated in Nexus Vicon
% Usefull for a sampling frequency of 250Hz for kinematic data
% Modifided from the Auto-generated Scrip by MATLAB on 2022/08 
%           by Claudia Romero 
%==========================================================================
% INPUT:
%--------------------------------------------------------------------------
% file = location and name of the csv file
%==========================================================================
% OUTPUT:
%--------------------------------------------------------------------------   
%  Data = Table with the captures data
%==========================================================================

%% Set up the Import Options and import the data
try 
    imported_data = string(readcell(filedir));
catch
    disp('Wrong file extension')
    Data = [];
    return
end

% Identify position of the data
x_pos = sum(contains(imported_data,'X'),2);
[~,x_pos] = max(x_pos);

row_labels = x_pos - 1;
labels = imported_data(row_labels,:);
column_label = ~ismissing(labels);
column_label = find(column_label);
labels = extractAfter(labels,':');
labels = replace(labels,' ','_');
% Add Time, Frames and Sub Frames
labels(1:column_label(1)-1) = imported_data(row_labels+1,1:column_label(1)-1);
for i = column_label(1):3:length(labels)
    labels(i+1) = strcat(labels(i),'y');
    labels(i+2) = strcat(labels(i),'z');
    labels(i) = strcat(labels(i),'x');
end

row_data = [x_pos + 2, length(imported_data)];
data = double(imported_data(row_data(1):row_data(2),:));

Data = array2table(data);
Data.Properties.VariableNames = labels;
end %function