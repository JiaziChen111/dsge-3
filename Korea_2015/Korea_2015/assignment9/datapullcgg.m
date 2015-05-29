load data;

endvars = {'da' 'pie' 'r' 'rstar' 'tau' 'x' 'dy'}; 
 
%endvars = sort(endvars);

for i=1:length(endvars)
    eval(strcat(endvars{i}, '=','sim(',num2str(i),',:)'';'));
end
