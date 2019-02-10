function [x,fval,exitflag,info] = ValveOpt(Pipes,Junctions,np,nn,nv,k,A12,A10,h0,e,M,d,Wts,solver,method)
%% Smart valve location
% 
% Syntax: 
%   [x,fval,exitflag,info] = ValveOpt(Pipes,Junctions,nv,solver,method)
% 
% Description:
%   This function calculates and returns the optimal location for the
%   pressure valves along with the minimized vector and value
% 
% Input Argument:
%   Pipes       - Pipe cluster details (table)
%   Junctions   - Junction cluster details (table)
%   np          - Number of pipes in the network
%   nn          - Number of nodes/junctions in the network
%   nv          - Number of pressure values to be present in the network
%   k           - Time samples
%   A12         - Pipe-junction incidence matrix
%   A10         - Pipe-source incidence matrix
%   h0          - Source head value
%   e           - Junction elevation vector
%   M           - Positive matrix for valve placement and control
%   d           - Demand vector based on the time instance
%   Wts         - Optimization weights (normalized/non-normalized) 
%   solver      - Nonlinear optimization method
%   method      - Reformulation method (penality or relaxation)

%
% Output Argument:
%   x           - Minimized output vector
%   fval        - Minimized function value
%   exitflag    - Optimization status
%   info        - Optimization details

%% Input Options
if nargin < 15 || isempty(method)
	method = 'Relaxation';
end

if nargin < 14 || isempty(solver)
	solver = 'BONMIN';
end

%% Function Code

    a0 = repmat(Pipes.a0,2,1);
    b0 = repmat(Pipes.b0,2,1);
    vidx = k*(nn +2*np)+1:k*(nn +2*np)+2*np;

    rho = 0;                                % Penalty factor
    t = 1;                                  % Relaxation factor

% Initial Guess
%     p0 = repmat(Junctions.MinPre+Junctions.Z,f,1); % 
    p0 = repmat(30+randi(30,nn,1),k,1);
    q0 = repmat(0.00001*randi(50,np,1),2*k,1);
    v0 = zeros(2*np,1);
    x =[p0;q0;v0];                          %x0 =[ones(nn,1);ones(2*np,1);zeros(2*np,1)];
 
    lb = [5+repmat(Junctions.Z,k,1);0.0001*ones(2*np*k,1);zeros(2*np,1)];
    ub = [80*ones(nn*k,1);repmat(Pipes.qMax,2*k,1);ones(2*np,1)];

%%
tic;
    switch solver
        case 'BONMIN'
            disp ('NLP: BONMIN')
            
            % Objective function
                fun = @(x) sum(x.*Wts);
                
            % Gradient
                grad = @(x) [repmat(Junctions.nodeWt,k,1);zeros(2*np*k,1);zeros(2*np,1)]';
                
            % Constraints
            nlcon = @(x) Constraints(x,np,nn,nv,k,A12,A10,h0,e,a0,b0,M,d,solver);
            nlrhs = [zeros(2*np*k,1);zeros(2*np*k,1);zeros(nn*k,1);ones(np,1);nv];
            nle = [ones(2*np*k,1);-1*ones(2*np*k,1);zeros(nn*k,1);-1*ones(np,1);0];

            % Jacobian
            jac = @(x) ConstJacSparse(x,np,nn,k,A12,A10,h0,e,a0,b0,M,solver);
            jacstr = @() JacobianStr(np,nn,k,A12,solver);

            % Integer Constraints
            xtype = [repmat('C',1,k*(nn+2*np)) repmat('B',1,2*np)];

            % Create OPTI Object
            % opts = optiset('solver',solver,'maxiter',1500,'display','final');
            opts = optiset('solver',solver,'maxiter',1500);
            Opt = opti('fun',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'jac',jac,'jacstr',jacstr,...
                       'xtype',xtype,'options',opts,'ndec',length(x));
                   
            % MINLP Optimization
            [x,fval,exitflag,info] = solve(Opt,x);
            
        case 'IPOPT'
            switch method
                case 'Penalty'
                    disp ('NLP: IPOPT - Penalty Method')
                    
                    % Objective function
                        % fun = @(x) sum(x.*Wts) + rho*sum(x(vidx(1:np)).*(1-x(vidx(np+1:end)))); % Forcing complementarity into function
                        fun = @(x) sum(x.*Wts) + rho*sum(x(vidx(1:np)).*x(vidx(np+1:end))); % Complementary considered in output vector
                        
                    % Gradient
                        % grad = @(x) [repmat(Junctions.nodeWt,k,1);zeros(2*np*k,1);rho*[(1-x(vidx(np+1:end)));-x(vidx(1:np))]]'; % Type: 1
                        grad = @(x) [repmat(Junctions.nodeWt,k,1);zeros(2*np*k,1);rho*[x(vidx(np+1:end));x(vidx(1:np))]]'; % Type: 2
                        
                    % Constraints
                    nlcon = @(x) Constraints(x,np,nn,nv,k,A12,A10,h0,e,a0,b0,M,d,solver,method,t);
                    nlrhs = [zeros(2*np*k,1);zeros(2*np*k,1);zeros(nn*k,1);nv];
                    nle = [ones(2*np*k,1);-1*ones(2*np*k,1);zeros(nn*k,1);0];
                        
                    % Jacobian
                    jac = @(x) ConstJacSparse(x,np,nn,k,A12,A10,h0,e,a0,b0,M,solver,method);
                    jacstr = @() JacobianStr(np,nn,k,A12,solver,method);
                    
                    % Initial Optimization
                    % opts = optiset('solver',solver,'maxiter',1500,'display','final');
                    opts = optiset('solver',solver,'maxiter',1500);
                    Opt = opti('fun',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'jac',jac,'jacstr',jacstr,...
                               'options',opts,'ndec',length(x));
                    [x,fval,exitflag,info] = solve(Opt,x);
                    
                    % Modified NLP Optimization
                    pBeta = 1.1;
                    pAlpha = 0.01;
                    pEta = 10e-6;
                    it = 0;

                    rho = pAlpha*fval;
                    
                    while max(min(x(vidx),ones(2*np,1)-x(vidx))) > pEta

                        fun = @(x) sum(x.*Wts) + rho*sum(x(vidx(1:np)).*x(vidx(np+1:end)));
                        grad = @(x) [repmat(Junctions.nodeWt,k,1);zeros(2*np*k,1);rho*[x(vidx(np+1:end));x(vidx(1:np))]]';
                        OptItr = opti('fun',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'jac',jac,'jacstr',jacstr,...
                               'options',opts,'ndec',length(x)); 

                        [x,fval,exitflag,info] = solve(OptItr,x);
                        rho = pBeta*rho;

                        it = it + 1;
                        if it > 50
                            break;
                        end
                    end
                           
                case 'Relaxation'
                    disp ('NLP: IPOPT - Relaxation Method')
                    
                    % Objective function
                        fun = @(x) sum(x.*Wts);
                        
                    % Gradient
                        grad = @(x) [repmat(Junctions.nodeWt,k,1);zeros(2*np*k,1);zeros(2*np,1)]';
                        
                    % Constraints
                    nlcon = @(x) Constraints(x,np,nn,nv,k,A12,A10,h0,e,a0,b0,M,d,solver,method,t);
                    nlrhs = [zeros(2*np*k,1);zeros(2*np*k,1);zeros(nn*k,1);0;nv];
                    nle = [ones(2*np*k,1);-1*ones(2*np*k,1);zeros(nn*k,1);-1;0];
                    
                    % Jacobian
                    jac = @(x) ConstJacSparse(x,np,nn,k,A12,A10,h0,e,a0,b0,M,solver,method);
                    jacstr = @() JacobianStr(np,nn,k,A12,solver,method);

                    % Initial Optimization
                    % opts = optiset('solver',solver,'maxiter',1500,'display','final');
                    opts = optiset('solver',solver,'maxiter',1500);
                    Opt = opti('fun',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'jac',jac,'jacstr',jacstr,...
                               'options',opts,'ndec',length(x)); 
                    [x,fval,exitflag,info] = solve(Opt,x);
                    
                    % Modified NLP Optimization       
                    tmin = 10e-15;
                    rEta = 10e-6;
                    rBeta = 10e-4;

                    while (max(min(x(vidx),ones(2*np,1)-x(vidx))) > rEta && t > tmin)
                        nlcon = @(x) Constraints(x,np,nn,nv,k,A12,A10,h0,e,a0,b0,M,d,solver,method,t);
                        OptItr = opti('fun',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'jac',jac,'jacstr',jacstr,...
                               'options',opts,'ndec',length(x)); 

                        [x,fval,exitflag,info] = solve(OptItr,x);
                        t = rBeta*t;
                    end  
            end
         % Cleaning up
         [~,vPos] = sort(x(vidx),'descend');
     
         x(vidx(vPos(1:nv))) = 1;
         x(vidx(vPos(nv+1:end))) = 0;
    end
toc;    
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
    