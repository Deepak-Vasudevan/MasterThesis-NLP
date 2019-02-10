function [a,b] = HeadLoss(rC,qMax,qMin)
%% Diurnal profile generation
% 
% Syntax: 
%   [a,b] = HeadLoss(rCoeff,qMax,qMin)
% 
% Description:
%   This function returns a random diurnal profile determined by the
%   category that is mapped to the corresponding node
% 
% Input Argument:
%   rCoeff       - Resistance Coefficient of the pipe
%   qMax         - Max flow defined by the diameter of the pipe
%   qMin         - Min flow of the pipe
%
% Output Argument:
%   [a,b]        - Minima of the integral square error over [qMin,qMax]

%% Input Options
if nargin < 2 || isempty(qMin)
	qMin = 0;
end

%% Function Code

% Method : 1 Direct Implementation

% Check if second differential is greater than zero (determine if the
% first differential at local extrema point is indeed a minima)

D = (qMax^8 + qMin^8 - 16*(qMax^5*qMin^3 + qMax^3*qMin^5))/60;

fn = @(q)2.*q.^4;
I = integral(fn,qMin,qMax);

% Find the minima if exists

if I > 0 && D > 0
        A = [4*qMax 5; 3*qMax 4];
        b = [4.122*rC*qMax^0.852; 3.1153*rC*qMax^0.852];
        sol = A\b;
        
        a = sol(1);
        b = sol(2);
end

% % Method : 2  Using OPTI Toolbox
% 
% fun = @(w)(integral(@(q)((w(1)*q.^2 + w(2).*q - rC.*q.^1.852).^2),qMin,qMax));
% 
% % Initial value
% w0 = [1;1];
% 
% % Optimization function
% Opt = opti('fun',fun,'x0',w0);
% 
% % Solve for minima
% [w,~,~,~] = solve(Opt);
% 
% a = w(1);
% b = w(2);
