classdef StopWatch < handle
   properties (SetAccess = private)
      avgTime_;
      nMeasurements_;
      timer_;
      infoStr_;
   end
   methods (Access = public)
      function obj = StopWatch (infoStr, avgTime, nMeasurements)
          %% STOPWATCH Constructor to create new stop-watch
          obj.infoStr_       = ''; 
          obj.avgTime_       = 0;
          obj.nMeasurements_ = 0;
          
          if (nargin == 1)
            obj.infoStr_ = infoStr;
          end % if
          if (nargin > 1)
            obj.avgTime_ = avgTime;
            obj.nMeasurements_ = nMeasurements;
          end % if
      end 
      
      function start (obj)
         obj.timer_ = tic; 
      end
      
      function avgTime = getAverageTime (obj)      
         if (nargout == 0)
             fprintf ('%s: %.3fs (avg.)\n', obj.infoStr_, obj.avgTime_);
         else
             avgTime = obj.avgTime_;
         end % if
      end
      
      function duration = stop (obj)
         duration           = toc (obj.timer_);
         tmp_               = obj.avgTime_ * obj.nMeasurements_;
         tmp_               = tmp_ + duration;
         obj.nMeasurements_ = obj.nMeasurements_ + 1;
         obj.avgTime_       = tmp_ / obj.nMeasurements_;
      end     
      
      function showAvgTime (obj)
        fprintf ('%.3fs --- %s\n', obj.avgTime_, obj.infoStr_);
      end % function
      
%       function obj = plus (lhs, rhs)
%           infoStr       = sprintf ('"%s" + "%s"', lhs.infoStr_, rhs.infoStr_);
%           avgTime       = lhs.avgTime_ + rhs.avgTime_;
%           nMeasurements = 1;
%           obj = StopWatch (infoStr, avgTime, nMeasurements);
%       end % if
   end
end