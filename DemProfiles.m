%% Smart Water Network - Scenarios
% Creates pseudo random demand profiles based on available consumer data

GenDem = load('ResearchData\DemPattern.txt');
GenDem = GenDem';

opts = detectImportOptions('Forbrukning.xls','Sheet','Forbrukning');
opts.SelectedVariableNames = {'OBJECTID','Tjanstenr','Arsforbruk','CONSCODE'};
Forbruk = readtable('Forbrukning.xls',opts);
Forbruk.ProfType = zeros(size(Forbruk,1),1); 

ProfCat = struct;
ProfCat(1).P1 = {'220','320','321','325'};
ProfCat(1).P21 = {'820','824','825','826','827','828'};
ProfCat(1).P22 = {'430','800','823'};
ProfCat(1).P3 = {'221','433','498'};
ProfCat(1).P4 = {'X'};

DiProfile = struct;

%% Demand Curve Generation - General Pattern

disp('Curve Generation');
tic

for i = 1:size(Forbruk,1)
    if sum(contains(ProfCat(1).P1,Forbruk.CONSCODE(i)))==1
        DiProfile(i).ObjectId = Forbruk.OBJECTID(i); 
        DiProfile(i).Category = Forbruk.CONSCODE(i);
        DiProfile(i).BaseProfile = DiurnalProfile(1);
        DiProfile(i).Classify = 'Regular';
    elseif sum(contains(ProfCat(1).P21,Forbruk.CONSCODE(i)))==1
        DiProfile(i).ObjectId = Forbruk.OBJECTID(i); 
        DiProfile(i).Category = Forbruk.CONSCODE(i);
        DiProfile(i).BaseProfile = DiurnalProfile(2,'Business');
        DiProfile(i).Classify = 'Business';
    elseif sum(contains(ProfCat(1).P22,Forbruk.CONSCODE(i)))==1
        DiProfile(i).ObjectId = Forbruk.OBJECTID(i); 
        DiProfile(i).Category = Forbruk.CONSCODE(i);
        DiProfile(i).BaseProfile = DiurnalProfile(2,'Regular');
        DiProfile(i).Classify = 'Regular';
    elseif sum(contains(ProfCat(1).P3,Forbruk.CONSCODE(i)))==1
        DiProfile(i).ObjectId = Forbruk.OBJECTID(i); 
        DiProfile(i).Category = Forbruk.CONSCODE(i);
        DiProfile(i).BaseProfile = DiurnalProfile(3,'Flat');
        DiProfile(i).Classify = 'Industrial/Flat';
    elseif strcmp(ProfCat(1).P4,Forbruk.CONSCODE{i})==1
        DiProfile(i).ObjectId = Forbruk.OBJECTID(i); 
        DiProfile(i).Category = Forbruk.CONSCODE(i);
        DiProfile(i).BaseProfile = DiurnalProfile(4,'Special');
        DiProfile(i).Classify = 'Special';
    end
end

toc

%% File Generation 

disp('File Generation');
tic;

if exist('ResearchData\DemandProfile.txt', 'file') == 2
    delete('DemandProfile.txt');
end

if exist('ResearchData\DemandMultiplier.txt', 'file') == 2
    delete('ResearchData\DemandMultiplier.txt');
end

if exist('ResearchData\CyclicPattern.txt', 'file') == 2
    delete('ResearchData\CyclicPattern.txt');
end

if exist('ResearchData\ProfilePattern.txt', 'file') == 2
    delete('ResearchData\ProfilePattern.txt');
end

eof = 0;
noWeeks = 1;

for i = 1:size(DiProfile,2)
    
    if i == size(DiProfile,2)
        eof=1;
    end
    
    % Basic Profiles (Demand with no leakage)
    
    Category = cell2mat(DiProfile(i).Category);
    FileGen('Basic',DiProfile(i).ObjectId,Category,DiProfile(i).BaseProfile,DiProfile(i).Classify,noWeeks,eof);
    
    SmoothProfile = smooth(DiProfile(i).BaseProfile,0.125);
    sftfact = randi(3)*ceil(rand()-1/2 *randi(5));
    DiProfile(i).WeekendProfile = circshift(SmoothProfile,sftfact);
    
    FileGen('Smooth',DiProfile(i).ObjectId,Category,SmoothProfile,DiProfile(i).Classify,noWeeks,eof);
end


% for i = 1:size(DiProfile,2)
%     
%     if i == size(DiProfile,2)
%         eof=1;
%     end
%     
%     % Basic Profiles (Demand with no leakage)
%     
%     Category = cell2mat(DiProfile(i).Category);
%     FileGen('WeekdayBasic',DiProfile(i).ObjectId,Category,DiProfile(i).BaseProfile,DiProfile(i).Classify,noWeeks,eof);
%     
%     WeekendProfile = smooth(DiProfile(i).BaseProfile,0.125);
%     DiProfile(i).WeekendProfile = circshift(WeekendProfile,2);
%     
%     FileGen('WeekendBasic',DiProfile(i).ObjectId,Category,WeekendProfile,DiProfile(i).Classify,noWeeks,eof);
% 
%     % Profiles with Usage Demand as 90%
%     
%     WeekdayUse = 0.9*DiProfile(i).BaseProfile;
%     DiProfile(i).WeekdayUse = WeekdayUse;
%     FileGen('WeekdayUse',DiProfile(i).ObjectId,Category,WeekdayUse,DiProfile(i).Classify,noWeeks,eof);
% 
%     WeekendUse = 0.9*DiProfile(i).WeekendProfile;
%     DiProfile(i).WeekendUse = WeekendUse;
%     FileGen('WeekendUse',DiProfile(i).ObjectId,Category,WeekendUse,DiProfile(i).Classify,noWeeks,eof);
%     
%     % Profiles with Leakage Demand as 10%
%     
%     WeekdayLeak = 0.1*DiProfile(i).BaseProfile;
%     DiProfile(i).WeekdayLeak = WeekdayLeak;
%     FileGen('WeekdayLeak',DiProfile(i).ObjectId,Category,WeekdayLeak,DiProfile(i).Classify,noWeeks,eof);
%     
%     WeekendLeak = 0.1*DiProfile(i).WeekendProfile;
%     DiProfile(i).WeekendLeak = WeekendLeak;
%     FileGen('WeekendLeak',DiProfile(i).ObjectId,Category,WeekendLeak,DiProfile(i).Classify,noWeeks,eof);
% end

toc;

% %% Leakage Nodes
% 
%  NodeCnt = 4;
%  Duration = [1 7 14 30 91 183 365];
%  
%  EmitterCoeff = table({'Round hole'},{'Plastic'},0.52,'VariableNames',{'LeakType','PipeMtrl','Coeff'});
%  EmitterCoeff(2,:) = table({'Longitudinal crack'},{'Plastic/Cement'},1.04);
%  EmitterCoeff(3,:) = table({'Circumferential crack'},{'Plastic'},0.45);
%  EmitterCoeff(4,:) = table({'Corrosion cluster'},{'Metal'},0.93);
 
 
 