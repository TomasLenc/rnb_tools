function new_labels = rm_spaces(old_labels)
% input is a cell of strings 
% output is a cell where ANY whitespace has been removed from each string 

new_labels = old_labels; 

for i=1:length(old_labels)
     new_labels{i} = old_labels{i}(find(~isspace(old_labels{i}))); 
end