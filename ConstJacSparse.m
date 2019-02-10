function sjac = ConstJacSparse(x,np,nn,k,A12,A10,h0,e,a0,b0,M,solver,method)
%% Nonlinear optimization - Sparse Constraint Jacobian
% 
% Syntax: 
%   sjac = ConstJacSparse(x,np,nn,k,A12,A10,h0,e,a0,b0,M)
% 
% Description:
%   This function returns the jacobian of the constraint vector in a sparse form 
% 
% Input Argument:
%   x       - Input vector that needs to satisfy the constraints 
%   np      - Number of pipes in the network
%   nn      - Number of nodes/junctions in the network
%   k       - Time samples
%   A12     - Pipe-junction incidence matrix
%   A10     - Pipe-source incidence matrix
%   h0      - Source head value
%   e       - Junction elevation vector
%   a0      - Roughness coefficient factor from quardatic model for headloss 
%   b0      - Roughness coefficient factor from quardatic model for headloss
%   M       - Positive matrix for valve placement and control
%   solver  - Optimization algorithm
%   method  - Regularization method
%
% Output Argument:
%   sjac    - Jacobian matrix of the constraint vector

%% Input Options
if nargin < 13 || isempty(method)
	method = 'Relaxation';
end

if nargin < 12 || isempty(solver)
	solver = 'BONMIN';
end


%% Function Code

nl =(1:k)';

pidx = [(nl-1)*nn+1 (nl-1)*nn+nn];
qidx = [k*nn+((nl-1)*2*np+1) k*nn+((nl-1)*2*np+2*np)];
vidx = k*(nn +2*np)+1:k*(nn +2*np)+2*np;

o = numel(find(A12~=0));
colC = nn+4*np;

    switch solver
        case 'BONMIN'
            ntriplets = 3*k*(o + 2*np) + 4*np;
            rowC = k*(4*np+nn)+np+1;
        case 'IPOPT'
            switch method
                case 'Penalty'
                    ntriplets = 3*k*(o + 2*np) + 2*np;
                    rowC = k*(4*np+nn)+1;
                case 'Relaxation'
                    ntriplets = 3*k*(o + 2*np) + 4*np;
                    rowC = k*(4*np+nn)+2;
            end
    end

C = zeros(ntriplets,3);

idx=0;

for nl = 1:k
    
    % Constraint: 1 
    
    posI = (nl-1)*2*np;
    
    Cp = diag(x(qidx(nl,1):qidx(nl,2)))*(-A12);  % if any of the expected locations have zeros, sparse will be an issue.
    CpS = sparse(Cp);
    [CI,CJ,CX] = find(CpS);
        C(idx+1:idx+o,1) = posI + CI;
        C(idx+1:idx+o,2) = CJ;
        C(idx+1:idx+o,3) = CX;
    idx = idx+o;
    
    posJ = nn;
    Cq = - A12*x(pidx(nl,1):pidx(nl,2)) - A12*e - A10*h0...
        +(3*a0.*x(qidx(nl,1):qidx(nl,2)).^2+2*b0.*x(qidx(nl,1):qidx(nl,2)));
%     Cq = zeros(2*np,1);
        C(idx+1:idx+2*np,1) = (posI+1:posI+2*np)';
        C(idx+1:idx+2*np,2) = (posJ+1:posJ+2*np)';
        C(idx+1:idx+2*np,3) = Cq;
    idx = idx+2*np;
    
    % Constraint: 2  
    
    posI = k*2*np +(nl-1)*2*np;
    
    [CI,CJ,CX] = find(A12);
        C(idx+1:idx+o,1) =  posI + CI;
        C(idx+1:idx+o,2) =  CJ;
        C(idx+1:idx+o,3) = -CX;
    idx = idx+o;
    
    posJ = nn;
    Cq = -(2*a0.*x(qidx(nl,1):qidx(nl,2))+b0);
        C(idx+1:idx+2*np,1) = (posI+1:posI+2*np)';
        C(idx+1:idx+2*np,2) = (posJ+1:posJ+2*np)';
        C(idx+1:idx+2*np,3) = Cq;
    idx = idx+2*np;
    
    posJ = nn+2*np;
        C(idx+1:idx+2*np,1) = (posI+1:posI+2*np)';
        C(idx+1:idx+2*np,2) = (posJ+1:posJ+2*np)';
        C(idx+1:idx+2*np,3) = - M(1:2*np+1:end)';
    idx = idx+2*np;
    
    % Constraint: 3
    
    posI = k*4*np +(nl-1)*nn;

%     Cp = -0.5*c.*(A12'*l).*(alpha*x(pidx(nl,1):pidx(nl,2)).^(alpha-1));
%         C(idx+1:idx+2*np,1) = (posI+1:posI+nn)';
%         C(idx+1:idx+2*np,2) = (1:nn)';
%         C(idx+1:idx+2*np,3) = Cp(1:nn+1:end)';
%     idx = idx+nn;

    posJ = nn;
    [CI,CJ,CX] = find(A12');
        C(idx+1:idx+o,1) = posI + CI;
        C(idx+1:idx+o,2) = posJ + CJ;
        C(idx+1:idx+o,3) = CX;
    idx = idx+o;
    
end

    % Constraint: 4
    
    posI = k*(4*np+nn);
    posJ = nn+2*np;
    
    switch solver
        case 'BONMIN'
            C(idx+1:idx+2*np,1) = [posI+1:posI+np posI+1:posI+np]';
            C(idx+1:idx+2*np,2) = [posJ+1:posJ+np posJ+np+1:posJ+2*np]';
            C(idx+1:idx+2*np,3) = ones(2*np,1);
        idx = idx+2*np;
            posI = k*(4*np+nn)+np;
            posJ = nn+2*np;
        case 'IPOPT'
            if strcmp(method,'Relaxation')==1
                C(idx+1:idx+2*np,1) = (posI+1)*ones(2*np,1);
                C(idx+1:idx+2*np,2) = (posJ+1:posJ+2*np)';
                C(idx+1:idx+2*np,3) = [x(vidx(np+1:end));x(vidx(1:np))] ;
            idx = idx+2*np;
                posI = k*(4*np+nn)+1;
                posJ = nn+2*np;
            end
    end
    
    % Constraint: 5
    
    idx = idx+1;
        C(idx:end,1) = (posI+1)*ones(2*np,1);
        C(idx:end,2) = (posJ+1:posJ+2*np)';
        C(idx:end,3) = ones(2*np,1);
        
	% C = sortrows(C,[2 1]);
        
I = C(:,1);
J = C(:,2);
X = C(:,3);

sjac = sparse(I,J,X,rowC,colC);
