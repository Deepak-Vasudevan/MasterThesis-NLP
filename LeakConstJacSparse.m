function sjac = LeakConstJacSparse(x,np,nn,nv,k,A12,A10,A13,h0,e,rC,l,pLen,alpha,c)
%% Nonlinear optimization - Sparse Constraint Jacobian
% 
% Syntax: 
%   sjac = LeakConstJacSparse(x,np,nn,k,A12,A10,h0,e,rC,l,pLen,alpha,c)
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
%   rC      - Roughness coefficient factor
%   l       - Leakage location vector
%   pLen    - Leakage pipe lengths   
%   alpha   - Leakage exponent
%   c       - Leakage coefficient
%
% Output Argument:
%   sjac    - Jacobian matrix of the constraint vector

%% Input Options
if nargin < 13 || isempty(c)
	c = 0.7;
end

if nargin < 12 || isempty(alpha)
	alpha = 1;
end


%% Function Code

nl =(1:k)';

if alpha~=1  && alpha~=2
    alpha = 1;
end

lc = numel(find(l == 1));

pidx = [(nl-1)*nn+1 (nl-1)*nn+nn];
qidx = [k*nn+((nl-1)*2*np+1) k*nn+((nl-1)*2*np+2*np)];

A12f = -1*A12;
A12f(A12f<0) = 0;           % Start node incidence matrix 

o = numel(find(A12~=0));
colC = nn+2*np+nv;
rowC = k*(4*np+nn);

ntriplets = k*(3*o + 4*np + lc);
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
        +(2.852.*rC.*x(qidx(nl,1):qidx(nl,2)).^1.852);
%     Cq = zeros(2*np,1);
        C(idx+1:idx+2*np,1) = (posI+1:posI+2*np)';
        C(idx+1:idx+2*np,2) = (posJ+1:posJ+2*np)';
        C(idx+1:idx+2*np,3) = Cq;
    idx = idx+2*np;
    
    posJ = nn+2*np;
    Cv = diag(x(qidx(nl,1):qidx(nl,2)))*(-A13);
    CvS = sparse(Cv);
    [CI,CJ,CX] = find(CvS);
        C(idx+1:idx+nv,1) = posI + CI;
        C(idx+1:idx+nv,2) = posJ + CJ;
        C(idx+1:idx+nv,3) = CX;
    idx = idx+nv;
    
    % Constraint: 2  
    
    posI = k*2*np +(nl-1)*2*np;
    
    [CI,CJ,CX] = find(A12);
        C(idx+1:idx+o,1) =  posI + CI;
        C(idx+1:idx+o,2) =  CJ;
        C(idx+1:idx+o,3) = -CX;
    idx = idx+o;
    
    posJ = nn;
    Cq = -1.852.*rC.*x(qidx(nl,1):qidx(nl,2)).^0.852;
        C(idx+1:idx+2*np,1) = (posI+1:posI+2*np)';
        C(idx+1:idx+2*np,2) = (posJ+1:posJ+2*np)';
        C(idx+1:idx+2*np,3) = Cq;
    idx = idx+2*np;
    
    posJ = nn + 2*np;
    [CI,CJ,CX] = find(A13);
        C(idx+1:idx+nv,1) =  posI + CI;
        C(idx+1:idx+nv,2) =  posJ + CJ;
        C(idx+1:idx+nv,3) = -CX;
    idx = idx+nv;
    
    % Constraint: 3
    
    posI = k*4*np +(nl-1)*nn;
    posJ = nn;
    
    if alpha == 1
        Cp = -c.*diag(A12f'*(l.*pLen));
        Cp = sparse(Cp);
        [CI,CJ,CX] = find(Cp);
            C(idx+1:idx+lc,1) = posI + CI;
            C(idx+1:idx+lc,2) = CJ;
            C(idx+1:idx+lc,3) = CX;
        idx = idx+lc;

        Cq = A12' + A12f'*diag((c*0.5*1.852)*(l.*pLen).*(rC.*x(qidx(nl,1):qidx(nl,2)).^0.852));
        %Cq = sparse(Cq);
        [CI,CJ,CX] = find(Cq);
            C(idx+1:idx+o,1) = posI + CI;
            C(idx+1:idx+o,2) = posJ + CJ;
            C(idx+1:idx+o,3) = CX;
        idx = idx+o;
        
%     elseif alpha == 2
%         Cp = - c.*(A12f'*(l.*pLen)).*(2*x(pidx(nl,1):pidx(nl,2)) + A12f'*(rC.*x(qidx(nl,1):qidx(nl,2)).^1.852));
%         Cp = sparse(Cp);
%         [CI,CJ,CX] = find(Cp);
%             C(idx+1:idx+lc,1) = posI + CI;
%             C(idx+1:idx+lc,2) = CJ;
%             C(idx+1:idx+lc,3) = CX;
%         idx = idx+lc;
% 
%         Cq = A12' - c.*(A12f'*(l.*pLen)).*((1.852*rC*x(pidx(nl,1):pidx(nl,2))*x(qidx(nl,1):qidx(nl,2)).^0.852)...
%             +(0.926*rC^2*x(qidx(nl,1):qidx(nl,2)).^2.704));
%         %Cq = sparse(Cq);
%         [CI,CJ,CX] = find(Cq);
%             C(idx+1:idx+o,1) = posI + CI;
%             C(idx+1:idx+o,2) = posJ + CJ;
%             C(idx+1:idx+o,3) = CX;
%         idx = idx+o;
        
    end
    
end
        
I = C(:,1);
J = C(:,2);
X = C(:,3);

sjac = sparse(I,J,X,rowC,colC);

