function [outvar] = load_StrArray(filename,filepath,inclusion_fields,fields_idx,cells_idx)
% loads cell array varialbes from .mat file into structure array (See save_StrArray.m to save
% structures in this format) while allowing the specification of which
% fields/varialbes to includeor exclude as well as which cells from each cell array 
% (cell identifier variables are always loaded and do not need to be specified).
% To load all variables and all cells use load_StrArray(...,...,'IncludeField','all',[]) 
% and to load only the cell identifier varilabes use load_StrArray(...,...,'ExcludeField','all',[]) 
% Order of fields will be lost, except for the cell identifier variables (eg. animal name, cell name ect.).

Key_variables  = {'animal_name', 'patching_date', 'experimentator', 'slice_nr', 'cellname', ...
    'eye_inj_ord', 'brain_contra_ipsi', 'hemisphere', 'MD'}; % these are always loaded, and place at the top of the strucutre

if isempty(filepath)
    filepath = cd;
end
cd(filepath)
temp = whos('-file',[filename]);
varibale_list = {temp.name}';

if strcmp(inclusion_fields,'IncludeField') && any(contains(fields_idx,'all'))
    fields = cat(2,Key_variables,varibale_list(~ismember(varibale_list,Key_variables))');
    
elseif strcmp(inclusion_fields,'ExcludeField') && any(contains(fields_idx,'all'))
    fields = Key_variables;
    
elseif strcmp(inclusion_fields,'IncludeField')
    fields = cat(2,Key_variables,fields_idx);
    
elseif strcmp(inclusion_fields,'ExcludeField')
    fields = cat(2,Key_variables,varibale_list(~ismember(varibale_list,Key_variables))');
    fields(ismember(fields,idfields_idxx)) = [];
end

if nargin == 5 && ~isempty(cells_idx)
    file_info = matfile(filename);
    for i = 1:length(fields)
        eval([fields{i} ...
            ' = file_info.' fields{i} '(1,[' num2str(cells_idx) ']);'])
    end
else
    load([filepath '\' filename],'-mat',fields{:})
end

for i = 1:length(fields)
    for ii = 1:length(eval(fields{i}))
        if iscell(eval(fields{i}))
            eval(['outvar(ii).' fields{i} '=' fields{i} '{ii};']);
        elseif isnumeric(eval(fields{i}))
            eval(['outvar(ii).' fields{i} '=' fields{i} '(ii);']);
        end
    end
end

end
