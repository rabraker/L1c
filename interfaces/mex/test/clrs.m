classdef clrs
    
    methods (Static)
        function status = in_terminal()
           status = ~isempty(javachk('jvm'));
        end
        
        function clr = red()
            if clrs.in_terminal()
                clr='[0;31m';
            else
                clr = '';
            end
        end
        function clr = green()
            if clrs.in_terminal()
               clr = '[0;32m';
            else
                clr = '';
            end
        end
        function clr = std()
            if clrs.in_terminal()
                clr = '[m'; 
            else
                clr='';
            end
          end
          
          function str = fail_str(str)
            str = [clrs.red, str, clrs.std];
          end
          function str = pass_str(str)
            str = [clrs.green, str, clrs.std];
          end
          
    end
    
    
end

