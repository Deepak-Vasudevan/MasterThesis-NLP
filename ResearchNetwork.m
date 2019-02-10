%% Smart Water Networks - Optimization using MINLP

clear; 
clc;
close all;

orig_state = warning;
warning('off','all');

%% Input data and files used for processing

disp('Read input files');
tic

% Case Study : System comprising 22 nodes, 3 reservoirs, and 37 links

opts = detectImportOptions('ResearchData\CaseStudy.xlsx','Sheet','Junction_Data');
Junctions = readtable('ResearchData\CaseStudy.xlsx',opts);

opts = detectImportOptions('ResearchData\CaseStudy.xlsx','Sheet','Pipe_Data');
Pipes = readtable('ResearchData\CaseStudy.xlsx',opts);

opts = detectImportOptions('ResearchData\CaseStudy.xlsx','Sheet','Source_Data');
Source = readtable('ResearchData\CaseStudy.xlsx',opts);

opts = detectImportOptions('ResearchData\CaseStudy.xlsx','Sheet','Demand_Pattern');
Demand = readtable('ResearchData\CaseStudy.xlsx',opts);

toc

%% Frictional Headloss - Quadratic model

disp('Quadratic Model - Frictional Headloss');
tic
qMin = 0;

for i = 1:height(Pipes)
    Pipes.qMax(i) = (pi*(Pipes.Diameter(i))^2)/4;    % ******** Q = pi*D^2/4 *****
    Pipes.rC(i) = 10.670*Pipes.L(i)/((Pipes.Diameter(i))^4.871 * Pipes.RCoeff(i)^1.852);
    
    [Pipes.a0(i),Pipes.b0(i)] = HeadLoss(Pipes.rC(i),Pipes.qMax(i),qMin);
    
    % Initialization of the boundary condition for leaking pipes
    % (Min pressure of node endings define boundary condition)
    
    fnIdx = find(strcmp(Junctions.MUID,Pipes.FROMNODE(i))==1);
    tnIdx = find(strcmp(Junctions.MUID,Pipes.TONODE(i))==1);
    
    if isempty(fnIdx) == 0 && isempty(tnIdx) == 0
        Pipes.FnMinPre(i) = Junctions.MinPre(fnIdx);
        Pipes.TnMinPre(i) = Junctions.MinPre(tnIdx);
    elseif isempty(fnIdx) == 0 && isempty(tnIdx) == 1
        Pipes.FnMinPre(i) = Junctions.MinPre(fnIdx);
    elseif isempty(fnIdx) == 1 && isempty(tnIdx) == 0
        Pipes.TnMinPre(i) = Junctions.MinPre(tnIdx);
    end
       
    if strcmp(Pipes.FROMNODE(i),'wNode_23') == 1
        Pipes.FnMinPre(i) = 24.66;
    elseif strcmp(Pipes.FROMNODE(i),'wNode_24') == 1
        Pipes.FnMinPre(i) = 24.60;
    elseif strcmp(Pipes.FROMNODE(i),'wNode_25') == 1
        Pipes.FnMinPre(i) = 24.55;
    end
    
    if strcmp(Pipes.TONODE(i),'wNode_23') == 1
        Pipes.TnMinPre(i) = 24.66;
    elseif strcmp(Pipes.TONODE(i),'wNode_24') == 1
        Pipes.TnMinPre(i) = 24.60;
    elseif strcmp(Pipes.TONODE(i),'wNode_25') == 1
        Pipes.TnMinPre(i) = 24.55;
    end

end

clear qMin;
toc

%% Clustering Junction Nodes

% Clustering Method (applicable only for large scale networks, else cluster count can be 1)
    % 1. Random clustering based on the network centroid
    % 2. Clustering based on the location

disp('Clustering Junction Nodes');
tic

% Cluster count
n = 1;
    
mode = 2;

[Junctions,Pipes] = Clustering(n,Junctions,Pipes,mode);

toc

%% Otimization Step 1 : Location of the control valves

disp('Setup - Control Valve location');

% Solvers: IPOPT, BONMIN
% Methods: Penalty, Relaxation

tic;
solver = 'BONMIN';
method = 'Relaxation';

% Initialization: Optimize over individual clusters (Choose cluster or loop over entire network) 

%     ClSel = 1;                                              
% 
%     ClFList = (PipeCluster.FNcl==ClSel);
%     ClTlist = (PipeCluster.TNcl==ClSel);
%     ClList = [PipeCluster(ClFList,:);PipeCluster(ClTlist,:)];
%     [~,Plist,~] = unique(ClList.MUID,'stable');
%     Pipes = ClList(Plist,:);
% 
%     ClList = JunctionCluster.cl==ClSel;
%     Junctions = JunctionCluster(ClList,:);
% 
%     BNodecount = sum(Pipes.BV==1);
% 
%     BoundNodes = zeros(BNodecount,1);
%     id=1;
% 
%     for i = 1:height(Pipes)
%         if Pipes.FNcl(i)~=ClSel
%             BoundNodes(id) = find(strcmp(JunctionCluster.MUID,Pipes.FROMNODE(i))==1);
%             id=id+1;
%         elseif Pipes.TNcl(i)~=ClSel
%             BoundNodes(id) = find(strcmp(JunctionCluster.MUID,Pipes.TONODE(i))==1);
%             id=id+1;
%         end
%     end
% 
%     Junctions = [Junctions;JunctionCluster(BoundNodes,:)];

% ************ Process variables and vectors ******************************

    n0 = 1;                                                 % Number of sources
    nn = height(Junctions);                                 % Number of junctions
    np = height(Pipes);                                     % Number of pipes
    nv = 7;                                                 % Number of control valves

    k = 24;                                                 % Time duration
    s = 24;                                                 % Time sampling
    f = k/s;                                                % Time multipliers
    h0 = [54.66;54.60;54.50];                               % Reservoir head
    % d = table2array(repmat(Junctions(:,'Demand'),1,24));    % Junction demand
    p_ub = 80;                                              % Pressure upper bound
    p_lb = 30;                                              % Pressure lower bound
    d = Junctions.Demand*Demand.RelativeDemand';
    e = Junctions.Elev;                                     % Junction elevation

% ************ Incidence Matrices ***************************************** 

    A12 = zeros(2*np,nn);
    A10 = zeros(2*np,n0);

    M = zeros(2*np,2*np);                                   % Mij = H^m - Hj^l 

    for i = 1:np
        fnIdx = find(strcmp(Junctions.MUID,Pipes.FROMNODE(i))==1);
        tnIdx = find(strcmp(Junctions.MUID,Pipes.TONODE(i))==1);

        M(i,i) = p_ub - (p_lb + Pipes.TnMinPre(i));
        M(i+np,i+np) = p_ub - (p_lb + Pipes.FnMinPre(i));

        if isempty(fnIdx) == 0
            A12(i,fnIdx)=-1;
        end

        if isempty(tnIdx) == 0
            A12(i,tnIdx)=1;
        end    

        if strcmp(Pipes.FROMNODE(i),'wNode_23') == 1
            A10(i,1)=-1;
        elseif strcmp(Pipes.FROMNODE(i),'wNode_24') == 1 
            A10(i,2)=-1;
        elseif strcmp(Pipes.FROMNODE(i),'wNode_25') == 1
            A10(i,3)=-1;
        end

        if strcmp(Pipes.TONODE(i),'wNode_23') == 1
            A10(i,1)=1;
        elseif strcmp(Pipes.TONODE(i),'wNode_24') == 1
            A10(i,2)=1;
        elseif strcmp(Pipes.TONODE(i),'wNode_25') == 1
            A10(i,3)=1;
        end    
    end
    A12(np+1:end,:) = -1*A12(1:np,:);
    A10(np+1:end,:) = -1*A10(1:np,:);

    A12 = sparse(A12);
    A10 = sparse(A10);
    M = sparse(M);
    
   % M = 0.1*M;

% ************ Node weights ***********************************************

    fNode = grpstats(Pipes,'FROMNODE','sum','DataVars',{'L'});
    fNode.Properties.VariableNames = {'Node','GroupCount','L'};
    fNode.Properties.RowNames = {};
    tNode = grpstats(Pipes,'TONODE','sum','DataVars',{'L'});
    tNode.Properties.VariableNames = {'Node','GroupCount','L'};
    fNode.Properties.RowNames = {};

    Nodes = [fNode;tNode];
    NodeWeights = grpstats(Nodes,'Node','sum','DataVars',{'L','GroupCount'});
    NodeWeights.GroupCount=[];
    NodeWeights.Properties.VariableNames = {'MUID','nodeWt','linkIdx'};
    NodeWeights.nodeWt = NodeWeights.nodeWt/2;

    Junctions = innerjoin(Junctions,NodeWeights,'Keys','MUID');

    ClusterWts = grpstats(Junctions,'cl','sum','DataVars',{'nodeWt'});
    TotalWt = sum(ClusterWts.sum_nodeWt);

    vecLength = f*(2*np+nn)+2*np;

    x = zeros(vecLength,nv);
    fval = zeros(f,nv);
    exitflag = zeros(f,nv);
    infoArray = arrayfun(@(w) struct,f:nv,'UniformOutput',0);

    Wts = [repmat(Junctions.nodeWt,f,1);zeros(f*2*np,1);zeros(2*np,1)];
   
    Wts = Wts/TotalWt;                            % Normalization Factor
 toc
%% Optimal Valve Count
 
    for vl = 1:nv
        disp(['********************************* Iteration: ' num2str(vl) ' *****************************************']);
        for dm = 1:24
        [x(:,vl),fval(dm,vl),exitflag(dm,vl),infoArray{dm,vl}] = ValveOpt(Pipes,Junctions,np,nn,vl-1,f,A12,A10,h0,e,M,d(:,dm),Wts,solver,method);
        
        end
        disp('***************************************************************************************');
    end

%% Otimization Step 2 : Optimizing pressure setting of the control valves under leakage conditions

disp('Optimization - Control Valve leakage settings');

 nv = 6;
 % Note: Input the optimal valve count from Step 1 based on relative
 % improvement to the objective function observed from the function plot
 % - Automate later! 
 
% ************ Leakage Parameters *****************************************

    alpha = 1;                              % Leakage exponent
    c = 0.7;                                % Leakage coefficient
    
    l = zeros(np,1);                        % Leakage vector
    
    lcount = randi(3);                      % Leakage count
    lpos = randi(np,lcount,1);              % Leakage pointer    
    l(lpos) = 1;
    
    for i = 1:np
        if ~contains(Pipes.FROMNODE(i),'wNode_') == 0 && l(i) == 1
            l(i) = 0;
        end
        
        if ~contains(Pipes.TONODE(i),'wNode_') == 0 && l(i) == 1
            l(i) = 0;
        end
    end
    
    l = repmat(l,2,1);                      % Leakage location vector
    
	A13 = zeros(2*np,nv);
    valvLoc = find(x(end -2*np+1:end,nv+1)>0);
    for i=1:nv
        A13(valvLoc(i),i)=1;
    end
    A13 = sparse(A13);   

    disp('********************************* Valve Settings **************************************');
        [lx,lfval,lexitflag,linfoArray] = LeakageOpt(Pipes,Junctions,np,nn,nv,f,A12,A10,A13,h0,e,d,Wts(1:f*(nn+2*np)),l,alpha,c);
    disp('***************************************************************************************');
    
warning(orig_state);