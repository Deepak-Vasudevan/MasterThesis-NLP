function nlcons = LeakConstraints(x,np,nn,nv,k,A12,A10,A13,h0,e,rC,d,l,pLen,alpha,c)
%% Nonlinear optimization - Leakage Constraint Definition
% 
% Syntax: 
%   nlcons = LeakConstraints(x,np,nn,k,A12,A10,h0,e,rC,M,d,v,l,pLen,alpha,c)
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
%   k       - Time samples
%   A12     - Pipe-junction incidence matrix
%   A10     - Pipe-source incidence matrix
%   h0      - Source head value
%   e       - Junction elevation vector
%   rC      - Roughness coefficient factor
%   M       - Positive matrix for valve placement and control
%   d       - Demand vector based on the time instance
%   v       - Pressure valve location vector
%   l       - Leakage location vector
%   pLen    - Leakage pipe lengths     
%   alpha   - Leakage exponent
%   c       - Leakage coefficient
%
% Output Argument:
%   nlcons  - Constraint vector output

%% Input Options

if nargin < 16 || isempty(c)
	c = 0.7;
end

if nargin < 15 || isempty(alpha)
	alpha = 1;
end

%% Function Code

if alpha~=1  && alpha~=2
    disp(num2str(alpha));
    alpha = 1;
    disp ('The value of alpha has to be 1 or 2 for correct implementation. Proceeding with default value (alpha = 1).')
end

nl =(1:k)';

pidx = [(nl-1)*nn+1 (nl-1)*nn+nn];
qidx = [k*nn+((nl-1)*2*np+1) k*nn+((nl-1)*2*np+2*np)];
vidx = k*(nn +2*np)+1:k*(nn +2*np)+nv;

A12f = -1*A12;
A12f(A12f<0) = 0;           % Start node incidence matrix 

C1 = zeros(2*np*k,1);
C2 = zeros(2*np*k,1);
C3 = zeros(nn*k,1);

for nl = 1:k
    
    hf = rC.*x(qidx(nl,1):qidx(nl,2)).^1.852;
    
    C1((nl-1)*2*np+1:(nl-1)*2*np+2*np) = diag(x(qidx(nl,1):qidx(nl,2)))*(- A12*x(pidx(nl,1):pidx(nl,2)) - A12*e - A10*h0 + hf - A13*x(vidx));
    
    C2((nl-1)*2*np+1:(nl-1)*2*np+2*np) = - A12*x(pidx(nl,1):pidx(nl,2)) - A12*e - A10*h0 - hf - A13*x(vidx);
    
    if alpha == 1
        C3((nl-1)*nn+1:(nl-1)*nn+nn) = A12'*x(qidx(nl,1):qidx(nl,2)) - d(:,nl) ...
            - c.*(A12f'*(l.*pLen)).*(x(pidx(nl,1):pidx(nl,2)) - 0.5*A12f'*(hf.*l));
        
%     elseif alpha == 2
%         C3((nl-1)*nn+1:(nl-1)*nn+nn) = A12'*x(qidx(nl,1):qidx(nl,2)) - d(:,nl) - c.*(A12f'*(l.*pLen)).*(x(pidx(nl,1):pidx(nl,2))...
%             + 0.5*A12f'*hf).^alpha;
    end

end
    
    nlcons = [C1;C2;C3];