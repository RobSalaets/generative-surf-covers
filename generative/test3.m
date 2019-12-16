classdef test3 < handle
   properties
      ev
   end
   methods
      function obj = test3(val)
         if nargin > 0
            obj.ev = val;
         end
      end
   end
end