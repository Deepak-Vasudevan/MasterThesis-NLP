function sjacstr = JacobianStr(np,nn,k,A12,solver,method)
%% Nonlinear optimization - Sparse Constraint Jacobian structure
% 
% Syntax: 
%   sjacstr = JacobianStr(np,nn,k,A12)
% 
% Description:
%   This function returns the jacobian of the constraint vector in a sparse form 
% 
% Input Argument: 
%   np      - Number of pipes in the network
%   nn      - Number of nodes/junctions in the network
%   k       - Time samples
%   A12     - Pipe-junction incidence matrix
%   solver  - Optimization algorithm
%   method  - Regularization method
%
% Output Argument:
%   sjacstr - Structural matrix of the Constraint Jacobian matrix

%% Input Options
if nargin < 6 || isempty(method)
	method = 'Relaxation';
end

if nargin < 5 || isempty(solver)
	solver = 'BONMIN';
end

%% Function Code
N = k*(nn+2*np)+2*np;

C1 = zeros(2*np*k,N);
C2 = zeros(2*np*k,N);
C3 = zeros(nn*k,N);

I = A12~=0;
      
C1(:,1:nn) = repmat(I,k,1);                                         
C1(:,nn+1:nn+2*np) = repmat(diag(ones(2*np,1)),k,1);                 
    
C2(:,1:nn) = repmat(I,k,1);                                         
C2(:,nn+1:nn+4*np) = repmat(diag(ones(2*np,1)),k,2);                 
     
% C3(:,1:nn) = repmat(diag(ones(nn,1)),k,1);                        
C3(:,nn+1:nn+2*np) = repmat(I',k,1);                                

C5 = [zeros(1,nn+2*np) ones(1,2*np)];                                

    switch solver
        case 'BONMIN'
            C4 = [zeros(np,nn+2*np) diag(ones(np,1)) diag(ones(np,1))];          
            sjacstr = [C1;C2;C3;C4;C5];
        case 'IPOPT'
            switch method
                case 'Penalty'
                    sjacstr = [C1;C2;C3;C5];
                case 'Relaxation'
                    C4 = [zeros(1,nn+2*np) ones(1,2*np)];                       
                    sjacstr = [C1;C2;C3;C4;C5];
            end
    end
    
% sjacstr = sparse(C2);
