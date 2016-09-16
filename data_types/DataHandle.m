classdef DataHandle < handle
    properties (SetAccess = private, GetAccess = public)
        data_;
    end % properties
    methods 
        function obj = DataHandle (data)
            obj.data_ = [];
            if (nargin > 0) ; obj.data_ = data ; end % if 
        end % constructor
        
%         function setDataField (obj, fieldname, newValue)
%             warning ('This function is not tested!');
%             
%             obj.data_.(fieldname) = newValue;
%         end % function 
    end % methods
end % class