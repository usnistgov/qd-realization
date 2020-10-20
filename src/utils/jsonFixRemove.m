function jsonString = jsonFixRemove(jsonString, str2Remove)
rem_ind_start = num2cell(strfind(jsonString, str2Remove)); % Find start string to remove
index2rm = cell2mat(cellfun(@(x) x:x+length(str2Remove)-1,rem_ind_start,'UniformOutput',false)); % Create index of char to remove
jsonString(index2rm) = []; % Remove temporary vector

end