function sjacstr = LeakJacobianStr(np,nn,nv,l,k,A12,A13)
%% Nonlinear optimization - Sparse Constraint Jacobian structure
% 
% Syntax: 
%   sjacstr = LeakJacobianStr(np,nn,l,k,A12)
% 
% Description:
%   This function returns the jacobian of the constraint vector in a sparse form 
% 
% Input Argument: 
%   np      - Number of pipes in the network
%   nn      - Number of nodes/junctions in the network
%   l       - Leakage vector
%   k       - Time samples
%   A12     - Pipe-junction incidence matrix
%
% Output Argument:
%   sjacstr - Structural matrix of the Constraint Jacobian matrix

%% Function Code
N = nn+2*np+nv;

C1 = zeros(2*np*k,N);
C2 = zeros(2*np*k,N);
C3 = zeros(nn*k,N);

I = A12~=0;

A12f = -1*A12;
A12f(A12f<0) = 0;           % Start node incidence matrix 

t = logical(A12f'*l);
      
C1(:,1:nn) = repmat(I,k,1);                                         
C1(:,nn+1:nn+2*np) = repmat(diag(ones(2*np,1)),k,1);                
C1(:,nn+2*np+1:nn+2*np+nv) = repmat(A13,k,1);
    
C2(:,1:nn) = repmat(I,k,1);                                         
C2(:,nn+1:nn+2*np) = repmat(diag(ones(2*np,1)),k,1);                
C1(:,nn+2*np+1:nn+2*np+nv) = repmat(A13,k,1);
     
C3(:,1:nn) = repmat(diag(t),k,1);                                   
C3(:,nn+1:nn+2*np) = repmat(I',k,1);                                

sjacstr = [C1;C2;C3];