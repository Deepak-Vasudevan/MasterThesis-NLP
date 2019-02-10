function pat = DiurnalProfile(categ,opt)
%% Diurnal profile generation
% 
% Syntax: 
%   pat = DiurnalProfile(category, option)
% 
% Description:
%   This function returns a random diurnal profile determined by the
%   category that is mapped to the corresponding node
% 
% Input Argument:
%   categ       - Classification category of the demand curve
%   opt         - Additional modifications to profile
%
% Output Argument:
%   pat         - the random generated pattern accoring to the category

%% Input Options
if nargin < 2 || isempty(opt)
	opt = 'Regular';
end

%% Function Code

% Variable declarations

tc = 3600*24;                    % 24 hour samples in seconds for better randomization
t = 0:tc-1;

switch categ
    
    % Residential profiles - Two distinctive peaks generally at 07:00 and
    % 19:00; sharp patterns with variations created by phase shifts
    case 1
        phi = 3600*round(randn);
        pN = 2;
        Ac = rand*ceil(t);
        
        C1 = -cos(pN*(2*pi*t/tc+phi));
        C1 = C1 + abs(max(C1));

        C2 = -Ac.*cos(pN*(2*pi*t/tc+phi));
        C2 = zscore(C2);
        C2 = C2 + 0.5*randn(size(t));

        C = C1 + C2;

        DemPat = abs(downsample(C,3600));
        DemPat = smooth(DemPat,0.125);
        if min(DemPat) <0
            DemPat = DemPat + abs(min(DemPat));
        end
        pat = DemPat*24/sum(DemPat); 
        
    % Business profiles - Broad spread of demand multipliers during day and
    % little or no (fixed work hours - 08:00 to 18:00) demand during nights. 
    % No phase shifts applied. 
    case 2
        phi = 0;
        pN = 2;    
        Ac = 3 + randi(4);
        win = [zeros(6,1);ones(12,1);zeros(6,1)]';
        
        C1 = -Ac*cos((2*pi*t/tc+phi)) - Ac/2*cos((pN*2*pi*t/tc+phi));
        C1 = C1 + abs(max(C1));

        C2 = -(1-rand*cos(pN*(2*pi*t/tc+phi))).*cos(2*pi*t/24);
        C2 = C2 + randn(size(t));

        C = C1 + C2;

        DemPat = abs(downsample(C,3600));
        DemPat = smooth(DemPat,0.125);
        if min(DemPat) <0
            DemPat = DemPat + abs(min(DemPat));
        end
        
        if strcmp(opt,'Regular')==1
            pat = DemPat*24/sum(DemPat); 
        elseif strcmp(opt,'Business')==1
            DemPat = DemPat.*win;
            pat = DemPat*24/sum(DemPat); 
        end
        
    % Industrial profiles - Flat or consistant demand with small intermittant
    % peaks. No phase shifts or demand spreading.    
    case 3
        phi = 0;
        pN = 3;
        Ac = randi(5);
        
        C = -Ac*cos(pN*(2*pi*t/tc+phi));
        C = C + 0.5*randn(size(t));
        C = C + abs(max(C));

        DemPat = abs(downsample(C,3600));
        idx = DemPat<= 0.7*max(DemPat);
        DemPat(idx)=0.7*max(DemPat);
        DemPat = smooth(DemPat,0.125);
        if min(DemPat) <0
            DemPat = DemPat + abs(min(DemPat));
        end
        
        if strcmp(opt,'Regular')==1
            pat = DemPat*24/sum(DemPat); 
        elseif strcmp(opt,'Flat')==1
            pat = DemPat/sum(DemPat); 
        end
        
    % Restauraunt profiles - More distinct eak during the second half of 
    % the day. No major phase shifts or demand spread.    
    case 4
        phi = 0;
        pN = 2;
        Ac = 3 + randi(4);
        win = [zeros(9,1);ones(13,1);zeros(5,1)]';
        
        C1 = -Ac*sin((2*pN*pi*t/tc+phi)) - Ac/2*cos((pN*4*pi*t/tc+phi));
        C1 = C1 + abs(max(C1));

        C2 = -(1-rand*cos(pN*(2*pi*t/tc+phi))).*cos(2*pi*t/24);
        C2 = C2 + randn(size(t));

        C = C1 + C2;

        DemPat = downsample(C,3600);
        DemPat = smooth(DemPat,0.125);
        if min(DemPat) <0
            DemPat = DemPat + abs(min(DemPat));
        end
        DemPat = DemPat.*win;
        
        if strcmp(opt,'Regular')==1
            pat = DemPat*24/sum(DemPat); 
        elseif strcmp(opt,'Special')==1
            DemPat = DemPat.*win;
            pat = DemPat*24/sum(DemPat); 
        end
end
