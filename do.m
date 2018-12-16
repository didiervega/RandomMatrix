% EULER module load mpt
classdef do
    %UNTITLED Summary of this class goes here
    
    properties

    end
    
    methods (Static = true)    
	
        function nameS = getWithoutExt(f)
            
             %nC = find(isstrprop(f, 'punct')) -1; 
             nC = findstr(f, '.') -1;
             last = size(nC,2);
             nameS = sprintf('%s', f(1:nC(last)));
             if last > 1
                nameS(nC(1:(last-1))+1) = '';
             end
        end	
            
    end
    
    
end

