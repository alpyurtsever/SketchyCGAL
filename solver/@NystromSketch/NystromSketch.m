classdef NystromSketch < matlab.mixin.SetGet
    %NYSTROMSKETCH implements a class definition for the sketching method
    %described [TYUC2017Nys].
    %
    %[TYUC2017Nys] J.A. Tropp, A. Yurtsever, M. Udell and V. Cevher. Fixed-
    %Rank Approximation of a Positive-Semidefinite Matrix from Streaming
    %Data. In Proc. 31st Conference on Neural Information Processing Systems
    %(NIPS), Long Beach, CA, USA, December 2017.
    %
    %Coded by: Alp Yurtsever
    %Ecole Polytechnique Federale de Lausanne, Switzerland.
    %Laboratory for Information and Inference Systems, LIONS.
    %contact: alp.yurtsever@epfl.ch
    %Created: April 12, 2017
    %Last modified: November 26, 2019
    %
    %Nys???SKETCHv1.0
    %Copyright (C) 2017 Laboratory for Information and Inference Systems
    %(LIONS), Ecole Polytechnique Federale de Lausanne, Switzerland.
    %This code is a part of Nys???SKETCH toolbox.
    %Please read COPYRIGHT before using this file.
    %
    % THIS IS A MINIMAL COPY OF NYSTROM SKETCH FROM OUR "SKETCH" TOOLBOX.
    % FOR A MORE DETAILED IMPLEMENTATION WITH ADDITIONAL OPTIONS, VISIT
    % https://github.com/alpyurtsever/SKETCH
    
    
    %% properties
    properties (Access = private)
        Omega     % (n x k) dimensional test matrix for the range of A (std Gaussian + orthonormalization)
        S         % (n x k) dimensional range sketch
    end
    
    %% methods
    methods
        
        %% Constructor
        function obj = NystromSketch(n, R, field)
            if R > n, error('Sketch-size cannot be larger than the problem size.'); end
            if strcmp(field,'real')
                obj.Omega = randn(n,R);
            elseif strcmp(field,'complex')
                obj.Omega = randn(n,R) + 1i*randn(n,R);
            else
                error('Input ''field'' should be ''real'' or ''complex''.')
            end
            obj.S = zeros(n,R);
        end
        
        %% Reconstruct
        function [ U, Delta ] = Reconstruct( obj )
            S = obj.S;
            n = size(obj.S,1);
            % nu = eps*norm(Y);
            sigma = sqrt(n)*eps*max(vecnorm(S));
            S = S + sigma*obj.Omega;
            B = obj.Omega' * S;
            B = 0.5*(B+B');
            if ~any(B(:)), Delta = 0; U = zeros(n,1); 
            else
                C = chol(B);
                [U, Sigma, ~] = svd( S / C, 'econ' );
                Delta = max( Sigma.^2 - sigma*eye(size(Sigma)), 0 );
            end
            if nargout == 1; U = U*Delta*U'; end
        end
        
        %% Rank One Update
        function obj = RankOneUpdate( obj, v, eta )
            obj.S = (1-eta)*obj.S + eta*(v*(v'*obj.Omega));
        end
        
        %% Property set methods
        function obj = set.S(obj,value)
            if isequal(size(value), size(obj.S)) || isempty(obj.S)
                obj.S = value;
            else
                error('Size of input does not match with sketch size.')
            end
        end
        
    end
end

