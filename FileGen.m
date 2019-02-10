function FileGen(PClass,ObjID,DCat,DPat,DClass,noWeek,eof)
%% Generates text files in a particular format
%
% Syntax: 
%   FileGen(Pattern_class,ObjectID,Demand_Category,Demand_Pattern,Demand_Class,No_of_weeks,End_of_file)
% 
% Description:
%   This function generates text files with the required pattern data
% 
% Input Argument:
%   PClass  - Pattern classification type bassed on the day and catergory
%   ObjID   - Object associated to the pattern
%   DCat    - Profile code associated to the object
%   DPat    - Generated diurnal multipliers for the object
%   DClass  - Demand Classification based on generalized patterns
%   noWeek  - The number of days/weeks needed for simulation
%   eof     - Last iteration identifier for post processing

%% Input Options
    if nargin < 6 || isempty(noWeek)
        noWeek = 1;
    end
    
    if nargin < 7 || isempty(eof) 
        eof = 0;
    end

%% Function Code

switch PClass
    case 'WeekdayBasic'
        DPID = 'WDB_';
        dayId = 1; 
    case 'WeekendBasic'
        DPID = 'WEB_';
        dayId = 2;
    case 'WeekdayUse'
        DPID = 'WDU_';
        dayId = 1;
    case 'WeekendUse'
        DPID = 'WEU_';
        dayId = 2;
    case 'WeekdayLeak' 
        DPID = 'WDL_';
        dayId = 1;
    case 'WeekendLeak'
        DPID = 'WEL_';
        dayId = 2;
    case 'Basic'
        DPID = 'G_';
        dayId = 1; 
    case 'Smooth'
        DPID = 'S_';
        dayId = 1;         
end

    % Create Diurnal Profiles for the object
    fileName = 'DemandProfile.txt';
    fileID = fopen(fileName,'a');
    fprintf(fileID,'%s%d\t<NULL>\t%s\t%s\n',DPID,ObjID,DCat,DClass);
    fclose(fileID);
    
    % Map multipliers after smoothing to relevant profiles
    fileName = 'DemandMultiplier.txt';
    fileID = fopen(fileName,'a');
    seq=0;
    for i = 1:24
        fprintf(fileID,'%s%d\t%d\t30-12-1899 %02d:00:00\t30-12-1899 %02d:00:00\t%5.3f\n',DPID,ObjID,seq,i-1,i,round(DPat(i),3));
        seq = seq+1;
    end
    fclose(fileID);
    
    if eof == 1
        fid  = fopen(fileName,'r');
        f=fread(fid,'*char')';
        fclose(fid);
        f = strrep(f,'24:00:00','00:00:00');
        fid  = fopen(fileName,'w');
        fprintf(fid,'%s',f);
        fclose(fid); 
    end
    
    DaySeq = 1:7*noWeek;
    profName = extractAfter(DPID,2);
            
    if dayId == 1
    % Create Cyclic Patterns for the diurnal profiles (only created once where the demand is set for weekdays)        
        fileName = 'CyclicPattern.txt';
        fileID = fopen(fileName,'a');
        fprintf(fileID,'%s%d\t<NULL>\t<NULL>\t<NULL>\tDiurnal\t<NULL>\t<NULL>\t<NULL>\t<NULL>\t<NULL>\t0\n',profName,ObjID);
        fclose(fileID);
        
        profSeq = find(mod(DaySeq,7)~=6 & mod(DaySeq,7)~=0);
        
    elseif dayId == 2
        profSeq = find(mod(DaySeq,7)==6 | mod(DaySeq,7)==0);
    end
    
    % Map Diurnal profiles to the patterns
    fileName = 'ProfilePattern.txt';
    fileID = fopen(fileName,'a');
    for i = 1:length(profSeq)
        fprintf(fileID,'INSERT INTO mw_DPProfileD ([ProfileID],[Sqn],[PatternID]) VALUES (''%s%d'',%d,''%s%d'');\n'...
            ,profName,ObjID,profSeq(i),DPID,ObjID);
    end
    if dayId == 2
        fprintf(fileID,'\n');
    end
    fclose(fileID);