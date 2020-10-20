function jsonString = jsonFixReplace(jsonString, str2Remove, str2Replace)
rem_ind_start = num2cell(strfind(jsonString, str2Remove)); % Find start string to remove
index2rm = cell2mat(cellfun(@(x) x:x+length(str2Remove)-1,rem_ind_start,'UniformOutput',false)); % Create index of char to remove
jsonString(index2rm) = []; % Remove temporary vector
if ~isempty(rem_ind_start)
    jsonString = [jsonString(1:rem_ind_start{1}-1), str2Replace, jsonString(rem_ind_start{1}:end)];
end
end