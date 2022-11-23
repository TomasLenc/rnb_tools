function x = ensure_row(x)

if isempty(x)
    return
elseif isrow(x)
    return
elseif iscolumn(x)
    x = x';  
    return
else
    error('cannot convert to column...too many dimensions...'); 
end
