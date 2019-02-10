%% Plots and Graphs

%% Demand Curve

figure(1);
hold on;
plot(Demand.RelativeDemand);
xlabel('Time Period');
ylabel('Relative Demand');
grid on;

%% AZP for different valves - extended period

figure(2); 
hold on ; 
lgnd= [] ;
for i = 1:nv
    plot(fval(:,i),'color',rand(1,3));
    lgnd = [lgnd ; 'Valve Count : ' num2str(i,'%2d')];
end
title('Extended duration optimization');
xlabel('Time duration');
ylabel('Average Network Pressure (m)');
legend(lgnd);
grid on;
hold off;

%% Optimal Valve count

figure(3); 
hold on ; grid on;
plot(fval(1,:),'-o');
title('Optimal valve count');
xlabel('Number of Valves');
ylabel('Average Network Pressure (m)');
grid on;
hold off;

%% Basic Network

opts = detectImportOptions('ResearchData\mw_RESJunction_Initial.xlsx','Sheet','mw_RESJunction');
InputData = readtable('ResearchData\mw_RESJunction_Initial.xlsx',opts);

Pressure = InputData(:,{'Pressure','Pressure_Min','Pressure_Max','Pressure_Avg','TStep'});
%Pressure.Pr_Min = bsxfun(@times,Pressure.Pressure_Min,repmat(Junctions.nodeWt,25,1));
DemandVal = InputData(:,{'Demand','Demand_Min','Demand_Max','Demand_Avg','TStep'});
RowList = DemandVal.Demand == 0;
DemandVal(RowList,:)= [];

S1PressureStat = grpstats(Pressure,'TStep','mean');
S1PressureStat.Properties.VariableNames = {'TStep','GroupPrCount','Pressure','Pressure_Min','Pressure_Max','Pressure_Avg'};
S1DemandStat = grpstats(DemandVal,'TStep','mean');
S1DemandStat.Properties.VariableNames = {'TStep','GroupDemCount','Demand','Demand_Min','Demand_Max','Demand_Avg'};

func =  @mean;
DerivedPressure = varfun(func,S1PressureStat(:,4:6));
DerivedPressure.Properties.VariableNames={'Minimum','Maximum','Average'} ;
DerivedDemand = varfun(func,S1DemandStat(:,4:6));
DerivedDemand.Properties.VariableNames={'Minimum','Maximum','Average'} ;

VarNames = table({'Demand','Pressure'}','VariableNames',{'VarName'});
DNetStat1 = horzcat(VarNames,vertcat(DerivedDemand,DerivedPressure));

%% Leakage Network

opts = detectImportOptions('ResearchData\mw_RESJunction_Leakage.xlsx','Sheet','mw_RESJunction');
InputData = readtable('ResearchData\mw_RESJunction_Leakage.xlsx',opts);

LeakID = strcmp(InputData.MUID,'L_01')==1;
LeakData = InputData(LeakID,:);

Pressure = InputData(:,{'Pressure','Pressure_Min','Pressure_Max','Pressure_Avg','TStep'});
DemandVal = InputData(:,{'Demand','Demand_Min','Demand_Max','Demand_Avg','TStep'});
RowList = DemandVal.Demand == 0;
DemandVal(RowList,:)= [];

S2PressureStat = grpstats(Pressure,'TStep','mean');
S2PressureStat.Properties.VariableNames = {'TStep','GroupPrCount','Pressure','Pressure_Min','Pressure_Max','Pressure_Avg'};
S2DemandStat = grpstats(DemandVal,'TStep','mean');
S2DemandStat.Properties.VariableNames = {'TStep','GroupDemCount','Demand','Demand_Min','Demand_Max','Demand_Avg'};

func =  @mean;
DerivedPressure = varfun(func,S2PressureStat(:,4:6));
DerivedPressure.Properties.VariableNames={'Minimum','Maximum','Average'} ;
DerivedDemand = varfun(func,S2DemandStat(:,4:6));
DerivedDemand.Properties.VariableNames={'Minimum','Maximum','Average'} ;

VarNames = table({'Demand','Pressure'}','VariableNames',{'VarName'});
DNetStat2 = horzcat(VarNames,vertcat(DerivedDemand,DerivedPressure));

%% Optimized Network

opts = detectImportOptions('ResearchData\mw_RESJunction_Optimal.xlsx','Sheet','mw_RESJunction');
InputData = readtable('ResearchData\mw_RESJunction_Optimal.xlsx',opts);

Pressure = InputData(:,{'Pressure','Pressure_Min','Pressure_Max','Pressure_Avg','TStep'});
DemandVal = InputData(:,{'Demand','Demand_Min','Demand_Max','Demand_Avg','TStep'});
RowList = DemandVal.Demand == 0;
DemandVal(RowList,:)= [];

S3PressureStat = grpstats(Pressure,'TStep','mean');
S3PressureStat.Properties.VariableNames = {'TStep','GroupPrCount','Pressure','Pressure_Min','Pressure_Max','Pressure_Avg'};
S3DemandStat = grpstats(DemandVal,'TStep','mean');
S3DemandStat.Properties.VariableNames = {'TStep','GroupDemCount','Demand','Demand_Min','Demand_Max','Demand_Avg'};

func =  @mean;
DerivedPressure = varfun(func,S3PressureStat(:,4:6));
DerivedPressure.Properties.VariableNames={'Minimum','Maximum','Average'} ;
DerivedDemand = varfun(func,S3DemandStat(:,4:6));
DerivedDemand.Properties.VariableNames={'Minimum','Maximum','Average'} ;

VarNames = table({'Demand','Pressure'}','VariableNames',{'VarName'});
DNetStat3 = horzcat(VarNames,vertcat(DerivedDemand,DerivedPressure));

%% Optimized Leakage Network

opts = detectImportOptions('ResearchData\mw_RESJunction_LeakOptimal.xlsx','Sheet','mw_RESJunction');
InputData = readtable('ResearchData\mw_RESJunction_LeakOptimal.xlsx',opts);

LeakID = strcmp(InputData.MUID,'L_01')==1;
OptLeakData = InputData(LeakID,:);

Pressure = InputData(:,{'Pressure','Pressure_Min','Pressure_Max','Pressure_Avg','TStep'});
DemandVal = InputData(:,{'Demand','Demand_Min','Demand_Max','Demand_Avg','TStep'});
RowList = DemandVal.Demand == 0;
DemandVal(RowList,:)= [];

S4PressureStat = grpstats(Pressure,'TStep','mean');
S4PressureStat.Properties.VariableNames = {'TStep','GroupPrCount','Pressure','Pressure_Min','Pressure_Max','Pressure_Avg'};
S4DemandStat = grpstats(DemandVal,'TStep','mean');
S4DemandStat.Properties.VariableNames = {'TStep','GroupDemCount','Demand','Demand_Min','Demand_Max','Demand_Avg'};

func =  @mean;
DerivedPressure = varfun(func,S4PressureStat(:,4:6));
DerivedPressure.Properties.VariableNames={'Minimum','Maximum','Average'} ;
DerivedDemand = varfun(func,S4DemandStat(:,4:6));
DerivedDemand.Properties.VariableNames={'Minimum','Maximum','Average'} ;

VarNames = table({'Demand','Pressure'}','VariableNames',{'VarName'});
DNetStat4 = horzcat(VarNames,vertcat(DerivedDemand,DerivedPressure));

%% Model Plots

figure(4);
plot(0:24,S1DemandStat.Demand,'r-^');
grid on;hold on;
plot(0:24,S2DemandStat.Demand,'g--o');
plot(0:24,S3DemandStat.Demand,'b:.');
plot(0:24,S4DemandStat.Demand,'c-.*');
hold off;
xlabel('Time (hr)');
xticks(0:2:24);
xtickformat(['00:' num2str('%02d')]);
legend('Before optimization - benchmark network','Before optimization - network with leakage',...
    'After optimization - benchmark network','After optimization - network with leakage');
ylabel('Average Demand (l/s)');

figure(5);
plot(0:24,S1PressureStat.Pressure,'r-^');
grid on; hold on;
plot(0:24,S2PressureStat.Pressure,'g--o');
plot(0:24,S3PressureStat.Pressure,'b:.');
plot(0:24,S4PressureStat.Pressure,'c-.*');
hold off;
xlabel('Time (hr)');
xticks(0:2:24);
xtickformat(['00:' num2str('%02d')]);
legend('Before optimization - benchmark network','Before optimization - network with leakage',...
    'After optimization - benchmark network','After optimization - network with leakage');
ylabel('Average Network Pressure (m)');

figure(6);
plot(0:24,LeakData.Demand,'r-');
grid on;hold on;
plot(0:24,OptLeakData.Demand,'g-');
hold off;
xlabel('Time (hr)');
xticks(0:2:24);
xtickformat(['00:' num2str('%02d')]);
legend('Leakage flow before optimization','Leakage flow after optimization');
ylabel('Flow (l/s)');
