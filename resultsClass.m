classdef resultsClass < handle & dynamicprops
    
    properties
        sys
        mod
    end
    
    methods
        function obj = resultsClass()
            obj.sys = array2table(zeros(0,2), 'VariableNames',{'t', 'sys'});
            obj.mod = array2table(zeros(0,2), 'VariableNames',{'t', 'mod'});
        end
        
        function add(obj, t, field, val)
            try
                obj.(field) = [obj.(field); {t,val}];
%                 obj.(field){end+1} = val;
            catch
                obj.addprop(field);
                obj.(field) = array2table(zeros(0,2), 'VariableNames',{'t', (field)});
                
                    fprintf(['Logger: new field added to the logger object: ', field,'\n']);
                
%                 obj.fields{end+1} = field;
%                 obj.(field) = {};
%                 obj.(field){end+1} = val;
                obj.(field) = [obj.(field); {t,val}];
            end
            
        end
        
    end
end

