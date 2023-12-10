function [v,c] = voronoielimnation(idx,v,c)
v(idx,:)=[];
eliminated_cells = 0;
for j=1:length(c)
    if ~isempty(find(c{j-eliminated_cells}==idx,1))
        c(j-eliminated_cells)=[];
        eliminated_cells=eliminated_cells+1;
    else
        c{j-eliminated_cells}(c{j-eliminated_cells}>idx)=c{j-eliminated_cells}(c{j-eliminated_cells}>idx)-1;
    end
end
end