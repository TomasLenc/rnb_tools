function x = ensure_col(x)

if isempty(x)
    return
elseif iscolumn(x)
    return
elseif isrow(x)
    x = x';  
    return
else
    error('cannot convert to column...too many dimensions...'); 
end
