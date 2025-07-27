function mkdir2(fpath)
% makes a new directory if it doesn't exist yet  

if ~isdir(fpath)
    mkdir(fpath)
end
