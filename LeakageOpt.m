function [x,fval,exitflag,info] = LeakageOpt(Pipes,Junctions,np,nn,nv,k,A12,A10,A13,h0,e,d,Wts,l,alpha,c)
%% Optimized network pressure setting
% 
% Syntax: 
%   [x,fval,exitflag,info] = LeakageOpt(Pipes,Junctions,np,nn,k,A12,A10,h0,e,M,d,Wts,l)
% 
% Description:
%   This function calculates and returns the optimal settings for the
%   pressure valves along with the minimized vector and value
% 
% Input Argument:
%   Pipes       - Pipe cluster details (table)
%   Junctions   - Junction cluster details (table)
%   np          - Number of pipes in the network
%   nn          - Number of nodes/junctions in the network
%   k           - Time samples
%   A12         - Pipe-junction incidence matrix
%   A10         - Pipe-source incidence matrix
%   h0          - Source head value
%   e           - Junction elevation vector
%   M           - Positive matrix for valve placement and control
%   d           - Demand vector based on the time instance
%   Wts         - Optimization weights (normalized/non-normalized)
%   x           - Optimized output vector obtained from valve location 
%   l           - Leakage location vector
%   alpha       - Leakage exponent
%   c           - Leakage coefficient

%
% Output Argument:
%   x           - Minimized output vector
%   fval        - Minimized function value
%   exitflag    - Optimization status
%   info        - Optimization details

%% Input Options
if nargin < 16 || isempty(c)
	c = 0.7;
end

if nargin < 15 || isempty(alpha)
	alpha = 1;
end

if nargin < 14 || isempty(l)
	l = zeros(2*np,1);
end

%% Function Code
    
% Initial Input
   
    p0 = repmat(30+randi(30,nn,1),k,1);
    q0 = repmat(0.00001*randi(50,np,1),2*k,1);
    eta0 = zeros(nv,1);
    x =[p0;q0;eta0];
    % f = 1;
    
    pLen = repmat(Pipes.L,2,1).*l.*0.5;
    rC = repmat(Pipes.rC,2,1); 
    
    Wts = [Wts;zeros(nv,1)];
    
    %x0 =[ones(nn,1);ones(2*np,1)];
    
% Boundaries
    
    qlb = 0.0001*ones(2*np*k,1);
    qub = repmat(Pipes.qMax,2*k,1);

    plb = repmat(Junctions.MinPre,k,1); 
    pubidx = logical((abs(A12))'*l);
    
    pub = (1-pubidx).*(80*ones(nn*k,1)) + repmat(pubidx.*(Junctions.MinPre+Junctions.Z+1),k,1);
    % pub = (1-pubidx).*(80*ones(nn*k,1)) + repmat(pubidx.*(Junctions.MinPre+Junctions.Z),k,1) + f*pubidx.*(abs(A12)'*hl);

    lb = [plb;qlb;zeros(nv,1)];
    ub = [pub;qub;ones(nv,1)];

    
    disp ('NLP: IPOPT')

    % Objective function
        fun = @(x) sum(x.*Wts);

    % Gradient
        grad = @(x) [repmat(Junctions.nodeWt,k,1);zeros(2*np*k,1);zeros(nv,1)]';

    % Constraints
    nlcon = @(x) LeakConstraints(x,np,nn,nv,k,A12,A10,A13,h0,e,rC,d,l,pLen,alpha,c);
    nlrhs = [zeros(2*np*k,1);zeros(2*np*k,1);zeros(nn*k,1)];
    nle = [ones(2*np*k,1);-1*ones(2*np*k,1);zeros(nn*k,1)];

    % Jacobian
    jac = @(x) LeakConstJacSparse(x,np,nn,nv,k,A12,A10,A13,h0,e,rC,l,pLen,alpha,c);
    jacstr = @() LeakJacobianStr(np,nn,nv,l,k,A12,A13);

    % Create OPTI Object
    opts = optiset('solver','IPOPT','maxiter',1500);
    Opt = opti('fun',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'jac',jac,'jacstr',jacstr,...
               'options',opts,'ndec',length(x));

    % NLP Optimization
    [x,fval,exitflag,info] = solve(Opt,x);
            
    switch exitflag
        case 1
            disp('Optimal Solution');
        case 0
            disp('Iteration / Function Evaluation / Time Limit Reached');
        case -1
            disp('Infeasible Problem');
        case -2
            disp('Unbounded / Solver Error');
        case -3
            disp('Solver Specific Error');
        case -5
            disp('User exited');
    end
    