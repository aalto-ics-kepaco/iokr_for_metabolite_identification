classdef StopWatch2 < handle
    properties (SetAccess = private)
        avgTime_;
        nMeasurements_;
        timer_;
        infoStr_;
        clockType_;  
        currentTime_; 
        isPaused_;
    end % properties
    
    properties (Access = private)
        f_getDuration_; 
        f_startClock_;
    end % properties
   
   
    methods (Access = public)
    function obj = StopWatch2 (infoStr, clockType)
    %% STOPWATCH Constructor to create new stop-watch
        obj.reset();
        if (nargin < 1)
            infoStr = ''; 
        end % if
        if (nargin < 2)
            clockType = 'WALL';
        end % if

        obj.infoStr_   = infoStr; 
        obj.clockType_ = clockType;

        switch (upper (obj.clockType_))
            case 'WALL'
                obj.f_startClock_  = 'obj.timer_ = tic';
                obj.f_getDuration_ = 'toc (obj.timer_)';
            case 'CPU'
                obj.f_startClock_  = 'obj.timer_ = cputime';
                obj.f_getDuration_ = 'cputime - obj.timer_';
            otherwise
                error ('StopWatch2:StopWatch2:InvalidArgument', 'Invalid clock-type: %s', ...
                    obj.clockType_);
        end % switch
    end % function
      
    function start (obj)
    %% START starts the stopwatch
        feval (obj.f_startClock_);
    end
    
    function duration = stop (obj)
    %% STOP stops the stopwatch and averages the measured time
        duration           = feval (obj.f_getDuration_);
        tmp_               = obj.avgTime_ * obj.nMeasurements_;
        tmp_               = tmp_ + duration;
        obj.nMeasurements_ = obj.nMeasurements_ + 1;
        obj.avgTime_       = tmp_ / obj.nMeasurements_;
    end % function
    
    function pause (obj)
        
    end % function
    
    function reset (obj)
    %% RESET resets the stopwatch to zero.
        obj.avgTime_       = 0;
        obj.nMeasurements_ = 0;
        obj.currentTime_   = 0;
        obj.isPaused_      = false;
    end % function 
      
      
      

      
      function avgTime = getAverageTime (obj)      
         if (nargout == 0)
             fprintf ('%s: %.3fs (avg.)\n', obj.infoStr_, obj.avgTime_);
         else
             avgTime = obj.avgTime_;
         end % if
      end % function
      
  
      
      function showAvgTime (obj)
        fprintf ('%.3fs --- %s\n', obj.avgTime_, obj.infoStr_);
      end % function

    end % methods
end