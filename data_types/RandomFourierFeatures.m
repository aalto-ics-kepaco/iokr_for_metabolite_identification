classdef RandomFourierFeatures < handle
    properties (Access = private)
        randomMatrix_; 
        D_;
        d_;
        gamma_;
    end % private properties
    
    methods (Access = private)
        function Psi = getRandomFourierFeatures_ (obj, Y, gamma)
            W   = sqrt (2 * gamma) * obj.randomMatrix_;
            Psi = sqrt (1 / ceil (obj.D_ / 2)) * [cos(W*Y); sin(W*Y)];
        end % function
    end % private methods
    
    methods (Access = public)
        function obj = RandomFourierFeatures (d, D, gamma)
        %% RANDOMFOURIERFEATURES Constructor of the class
        %   obj = RANDOMFOURIERFEATURES (d, D) creates an object which can
        %   be used to approximate the feature-vectors of a gaussian kernel
        %   using random fourier features with dimension D. The random
        %   matrix is initiated and has size D x d, with d beeing the
        %   dimension of the original feature vectors. 
        %   NOTE: The random fourier feature dimension will be ceil(D/2)*2.
        %
        %   obj = RANDOMFOURIERFEATURES (d, D, gamma) see above.
        %   Additionally the gamma parameter for the gaussian to
        %   approximate is set. 
        %   NOTE: If 'isempty (gamma)' is true, than this value is ignored,
        %   e.g. if only the random seed should be set. 
        %
        %   obj = RANDOMFOURIERFEATURES (..., randomSeed) same as above,
        %   but the random matrix is initialized using a pre-defined random
        %   seed for 'randn'. 
            if (nargin < 2)
                error ('RandomFourierFeatures:RandomFourierFeatures:InvalidArgument', ...
                    'Not enough input arguments.');
            end % if
            
            if (nargin < 3) || (isempty (gamma))
                obj.gamma_ = [];
            else
                obj.gamma_ = gamma;    
            end % if

            if (~ all ([d, D] > 0))
                error ('RandomFourierFeatures:RandomFourierFeatures:InvalidArgument', ...
                    'Input and approximated feature dimension must be greater zero.');
            else
                obj.d_ = d;
                obj.D_ = D;
            end % if
            
            obj.resetRandomMatrix();
        end % function
        
        function resetRandomMatrix (obj)
        %% RESETRANDOMMATRIX re-initializes the random matrix
        %   RESETRANDOMMATRIX (obj) a new random matrix is created and
        %   replaces the one in the object obj.
        %
        %   RESETRANDOMMATRIX (..., randomSeed) same as above, but a
        %   pre-defined random seed is used.
            obj.randomMatrix_ = randn (ceil (obj.D_ / 2), obj.d_);
        end % function 
        
        function Psi = getRandomFourierFeatures (obj, Y, gamma)
        %% GETRANDOMFOURIERFEATURES random fourier features for the input
        %   Psi = GETRANDOMFOURIERFEATURES (Y) returns the random fourier
        %   features for the input Y (d x m, with m being the number of 
        %   examples). Psi is of dimension D x m. 
        %
        %   Psi = GETRANDOMFOURIERFEATURES (Y, gamma) see above. Addionally
        %   the gamma parameter of the object can be overwritten.
            if (nargin < 2) % NOTE: 'obj' is a paramter
                error ('RandomFourierFeatures:getRandomFourierFeatures:InvalidArgument', ...
                    'Not enought input arguments.');
            end % if
            
            if (size (Y, 1) ~= obj.d_)
                error ('RandomFourierFeatures:getRandomFourierFeatures:InvalidArgument', ...
                    'Input feature dimension does not match: d = %d but obj.d_ = %d.', ...
                    size (Y, 1), obj.d_);
            end % if
            if (nargin < 3) 
                if (isempty (obj.gamma_))
                    error ('RandomFourierFeatures:getRandomFourierFeatures:InvalidArgument', ...
                        'Not enough input arguments. No gamma-value is predefined so it needs to be provided.');
                end % if
                
                Psi = obj.getRandomFourierFeatures_ (Y, obj.gamma_);
            else
                Psi = obj.getRandomFourierFeatures_ (Y, gamma);
            end % if
        end % function
    end % public methods
end % class