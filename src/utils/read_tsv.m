function tbl=read_tsv(fname)

tbl = readtable(fname, 'Delimiter', '\t', 'FileType', 'text');                          
