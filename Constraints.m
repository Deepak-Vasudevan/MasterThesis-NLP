function nlcons = Constraints(x,np,nn,nv,k,A12,A10,h0,e,a0,b0,M,d,solver,method,t)
%% Nonlinear optimization - Constraint Definition
% 
% Syntax: 
%   nlcons = Constraints(x,np,nn,nv,k,A12,A10,h0,e,a0,b0,M,d)
% 
% Description:
%   This function returns a vector of constarint values the objective
%   function needs to satisfy for the optimization process. The constraints
%   are calculated based on the input values of the network.
% 
% Input Argument:
%   x       - Input vector that needs to satisfy the constraints 
%   np      - Number of pipes in the network
%   nn      - Number of nodes/junctions in the network
%   nv      - Number of pressure values to be present in the network
%   k       - Time samples
%   A12     - Pipe-junction incidence matrix
%   A10     - Pipe-source incidence matrix
%   h0      - Source head value
%   e       - Junction elevation vector
%   a0      - Roughness coefficient factor from quardatic model for headloss 
%   b0      - Roughness coefficient factor from quardatic model for headloss
%   M       - Positive matrix for valve placement and control
%   d       - Demand vector based on the time instance
%   solver  - Optimization algorithm
%   method  - Regularization method
%   t       - Regularization factor
%
% Output Argument:
%   nlcons  - Constraint vector output

%% Input Options

if nargin < 16 || isempty(t)
	t = 0;
end

if nargin < 15 || isempty(method)
	method = 'Relaxation';
end

if nargin < 14 || isempty(solver)
	solver = 'BONMIN';
end


%% Function Code

nl =(1:k)';

pidx = [(nl-1)*nn+1 (nl-1)*nn+nn];
qidx = [k*nn+((nl-1)*2*np+1) k*nn+((nl-1)*2*np+2*np)];
vidx = k*(nn +2*np)+1:k*(nn +2*np)+2*np;

C1 = zeros(2*np,k);
C2 = zeros(2*np,k);
C3 = zeros(nn,k);

for nl = 1:k

    C1((nl-1)*2*np+1:(nl-1)*2*np+2*np) = diag(x(qidx(nl,1):qidx(nl,2)))*( - A12*x(pidx(nl,1):pidx(nl,2)) - A12*e - A10*h0...
        +(a0.*x(qidx(nl,1):qidx(nl,2)).^2+b0.*x(qidx(nl,1):qidx(nl,2))));
    
    C2((nl-1)*2*np+1:(nl-1)*2*np+2*np) = - A12*x(pidx(nl,1):pidx(nl,2)) - A12*e - A10*h0...
        -(a0.*x(qidx(nl,1):qidx(nl,2)).^2+b0.*x(qidx(nl,1):qidx(nl,2))) - M*x(vidx);
    
    C3((nl-1)*nn+1:(nl-1)*nn+nn) = A12'*x(qidx(nl,1):qidx(nl,2)) - d(:,nl);
    
    %C3((nl-1)*nn+1:(nl-1)*nn+nn) = A12'*x(qidx(nl,1):qidx(nl,2)) - d(:,nl) - 0.5*c.*(A12'*l).*(x(pidx(nl,1):pidx(nl,2)).^alpha);
end

    C5 = sum(x(vidx));
    
    switch solver
        case 'BONMIN'
            C4 = x(vidx(1:np)) + x(vidx(np+1:end));
            nlcons = [C1;C2;C3;C4;C5];
        case 'IPOPT'
            switch method
                case 'Penalty'
                    nlcons = [C1;C2;C3;C5];
                case 'Relaxation'
                    C4 = sum(x(vidx(1:np)).*x(vidx(np+1:end))) - t;
                    nlcons = [C1;C2;C3;C4;C5];
            end
    end