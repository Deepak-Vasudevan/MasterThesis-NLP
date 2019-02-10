function [Junctions,Pipes] = Clustering(n,Junctions,Pipes,mode)
%% Clustering
% 
% Syntax: 
%   [Junctions,Pipes] = Clustering(n,Junctions,Pipes,mode)
% 
% Description:
%   This function reduces the network size in consideration by clustering
%   the junction nodes based on their location or distance from the network 
%   centroid usually defind by the input method
% 
% Input Argument:
%   n                   - Cluster count
%   Junctions           - Junction details (table)
%   Pipes               - Pipe details (table)
%   mode                - Clustering method (location based clustering or centroid based random clustering)
%
% Output Argument:
%   [Junctions,Pipes]   - Junctions and Pipes updated with cluster information

%% Input Options
if nargin < 4 || isempty(mode)
	mode = 2;
end

%% Function Code
%% Clustering Junction Nodes - Random uniform distribution based on distance from the centroid

switch mode
    case 1
disp('Clustering Junction Nodes based on centroid distance');

% Obtain the centroid of the distribution network
JnData = table2array(Junctions(:,{'OBJECTID','X','Y'}));
[~,C] = kmeans(JnData(:,2:3),n);

% Find the distance of each junction from the centroid
for i =1:height(Junctions)
    Junctions.Dist(i) = sqrt((Junctions.X(i) - C(1))^2 + (Junctions.Y(i) - C(2))^2);
end

% Sort Junction nodes by distance from the centroid
Junctions = sortrows(Junctions,'Dist');

% Assign Junction nodes to random clusters
cl = randi(n,height(Junctions),1);

for i = 1:height(Junctions)
    Junctions.cl(i) = cl(i);
end

t = grpstats(Junctions,'cl','mean','DataVars',{'Dist'});
lgnd=[];

figure(1);
hold on;
for i=1:n
    idx = cl==i;
    plot(Junctions.X(idx),Junctions.Y(idx),'color',rand(1,3),'marker','.','LineStyle','none');
    lgnd = [lgnd ; 'Cluster ' num2str(i,'%02.f') ': ' num2str(t.GroupCount(i),'%03.f')];
    text(Junctions.X(idx),Junctions.Y(idx),num2str(Junctions.cl(idx)));
end
title('Clustering');
legend(lgnd);
grid on;
hold off;


%% Clustering Junction Nodes - Location based grouping
    case 2
disp('Clustering Junction Nodes based on location');

% Cluster based on k-means cluster groups
JnData = table2array(Junctions(:,{'OBJECTID','X','Y'}));
[cl,~] = kmeans(JnData(:,2:3),n);

for i = 1:height(Junctions)
    Junctions.cl(i) = cl(i);
end

t = grpstats(Junctions,'cl','mean','DataVars',{'Elev'});
lgnd=[];

figure(1);
hold on;
for i=1:n
    idx = cl==i;
    plot(Junctions.X(idx),Junctions.Y(idx),'color',rand(1,3),'marker','.','LineStyle','none');
    hold on;
    lgnd = [lgnd ; 'Cluster ' num2str(i,'%02.f') ': ' num2str(t.GroupCount(i),'%03.f')];
end
title('Clustering');
legend(lgnd);
grid on;
hold off;

% Boundary nodes/valve location for each cluster

for i = 1:height(Pipes)
    fnIdx = find(strcmp(Junctions.MUID,Pipes.FROMNODE(i))==1);
    tnIdx = find(strcmp(Junctions.MUID,Pipes.TONODE(i))==1);
    
    if isempty(fnIdx) == 0 && isempty(tnIdx) == 0
        Pipes.FNcl(i) = Junctions.cl(fnIdx);
        Pipes.TNcl(i) = Junctions.cl(tnIdx);       
    elseif isempty(fnIdx) == 0 && isempty(tnIdx) == 1
        Pipes.FNcl(i) = Junctions.cl(fnIdx);
        Pipes.TNcl(i) = Junctions.cl(fnIdx);  
    elseif isempty(fnIdx) == 1 && isempty(tnIdx) == 0
        Pipes.FNcl(i) = Junctions.cl(tnIdx);
        Pipes.TNcl(i) = Junctions.cl(tnIdx);  
    end
    
    if Pipes.FNcl(i) ~= Pipes.TNcl(i)
        Pipes.BV(i) = 1;
    else
        Pipes.BV(i) = 0;
    end
end
end