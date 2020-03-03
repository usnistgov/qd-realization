function [x] = normalRandomGenerator(mu,sigma, varargin)
%NORMALRANDOMGENERATOR generates a random number using a normal distribution
%
%    [x] = NORMALRANDOMGENERATOR(mu,sigma)
% 		mu    - mean distribution
% 		sigma - dev std distribution
%
%    [x] = NORMALRANDOMGENERATOR(mu,sigma, N)
%       generates N numbers

%% Varargin processing 
if isempty(varargin)
    realizations = 1;
else
    realizations =varargin{1};
end
    
%% Output
x = (randn(realizations,1))*sigma + mu;

end