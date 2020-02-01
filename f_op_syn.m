function [osn,ost,os,ostn] = Nop_syn(ts1, ts2, dim, step)

% % % % % % % % % test for ordinal pattern correlation (FOR TESTING PURPOSES!!!) %%%%%%%%%%%%%
% % % % % % % clearvars;clc;close all;
% % % % % % % load /Users/johannm/Desktop/Academics/Data_Academics/Univ9large.mat;
% % % % % % % 
% % % % % % % x = Univ9large(2,:,1);
% % % % % % % y = Univ9large(3,:,1);
% % % % % % % ts1 = x;    ts2 = y;
% % % % % % % % ts1 = x(randperm(numel(x)));    ts2 = y(randperm(numel(y)));
% % % % % % % dim = 6;
% % % % % % % step = 'cons'; %'cons'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% function here %%%%%%%%%%% 
D = dim;                               %Dimension of the Ordinal Patterns (OrdPat)

ts_1 = ts1;                             %Time series X(t)
ts_2 = ts2;                             %Time series Y(t)

window = step;
M = numel(ts_1);                             %Number of samples of a time series

%%%%%%%% testing inputs
% if numel(num2str(M)) <= numel(num2str(factorial(D))) %comparing order of magnitude M>>D!
%     error('Remember: M >> D!, try another D')
% end
% if numel(ts_1) ~= numel(ts_2)
%     error('Times series must have the same length!')
% end
%% Chosing the symbol sequence length. (hence of the ordinal patterns too) due to sliding or consecutive windowing
if strcmp(window, 'slid')
    L = (M - D) + 1;    %Symbol sequence Lenght |S(t)|. Number of possible D-based partitions with 1 sample step (sliding windowing)
elseif strcmp(window, 'cons')
    L = floor(M/D);     %Symbol sequence Lenght |S(t)|. Number of possible D-based consecutive partitions (consecutive windowing)
    % vectores inicio y final para cortar la ts
    iniX = 1:D:M;
    final = iniX + D-1;
end
%% Building up the ordinal patterns
XDpartitions = zeros(L, D);                  %All original D-based partitions from the X(t) (just for comparing issue)
S1 = zeros(L, D);                             %S1(t) symbol sequence. Ordinal Patterns whose appear in my time series

YDpartitions = zeros(L, D);                  %All original D-based partitions from the Y(t) (just for comparing issue)
S2 = zeros(L, D);                             %S2(t) symbol sequence. Ordinal Patterns whose appear in my time series

% To calculate OS over one vector per channel (that is, all time series as
% one vector, and then dot product between them), take norm as:
% norm(0:D-1);
if D <= 9
    permutacions = perms(0 : D - 1);             %All D! possible Ordinal Patterns based on dimension D
    norm_OP = norm(permutacions(1,:));               %for normalisation of ordinal patterns vectors
else 
    norm_OP = norm(0:D-1);
end

is = 0;  %auxiliar index for final sample (sliding window)
for op = 1 : L                                  %along all posible partitions (D-dependency)
    
    if strcmp(window, 'slid')                       %   construyendo los patrones de orden con ventana MOVIBLE
        finalid = D + is;                           %cutoff sample for each partition (for the 1-1 sample sliding window)
        
        XDpartitions(op, : ) = ts_1(op : finalid);  %cutting up my X(t) (sample by sample) and storing it
        YDpartitions(op, : ) = ts_2(op : finalid);  %cutting up my Y(t) (sample by sample) and storing it
        
        is = is + 1;                                %increasing my auxiliar index (for the 1-1 sample sliding window)
    elseif strcmp(window, 'cons')                   %   construyendo los patrones de orden con ventana CONSECUTIVA
        XDpartitions(op, : ) = ts_1(iniX(op) : final(op));   %cutting up my X(t) (sample by sample) and storing it
        YDpartitions(op, : ) = ts_2(iniX(op) : final(op));   %cutting up my Y(t) (sample by sample) and storing it
    end
    
    %%%%% This part is for to achieve the BP Method........................
    %%%%% The Sympbol Sequence S(t) is based on the BP codification rule
    %%%%% Ordinal Patterns (OP) belonging to S(t) is a D-dimensional vector
    %%%%% (0,1,...,D-1) obtained by comparing the vectors(rows) from
    %%%%% original D-based partitions "XDpartitions" from the X(t).
    
    [~, id1] = sort(XDpartitions(op, : ), 'ascend'); %original indexes of XDpartition values after being sorted
    ordinalP1 = zeros(1, D);                         %one specific OP from previous XDpartition based on D
    ordp_bit1 = 0;                                   %the first bit in OP is always zero
    
    [~, id2] = sort(YDpartitions(op, : ), 'ascend'); %original indexes of YDpartition values after being sorted
    ordinalP2 = zeros(1, D);                         %one specific OP from previous YDpartition based on D
    ordp_bit2 = 0;                                   %the first bit in OP is always zero
    
    for idx = 1 : D                                 %along all elements in "id" vector
        ordinalP1(id1(idx)) = ordp_bit1;               %creating the OP according to values in each partition of X(t)
        ordp_bit1 = ordp_bit1 + 1;                    %increasing OP bit until (D - 1)
        
        ordinalP2(id2(idx)) = ordp_bit2;               %creating the OP according to values in each partition of Y(t)
        ordp_bit2 = ordp_bit2 + 1;                    %increasing OP bit until (D - 1)
    end
    S1(op, :) = ordinalP1;                    %storing the OP that appears in X(t) sample by sample
    S2(op, :) = ordinalP2;                    %storing the OP that appears in Y(t) sample by sample
end
S1_norm = S1/norm_OP;
S2_norm = S2/norm_OP;

ost = dot(S1_norm,S2_norm,2); % To get a time series value of OS (floor(M/D) values)
os = mean(ost); % To get a single value of OS

min = dot((0:D-1),(D-1:-1:0))/norm(0:D-1)^2; % Minimum possible value

osn = 2*(((os-min)/(1-min))-0.5);% Normalize to get a value between -1 and 1
ostn = 2*(((ost-min)/(1-min))-0.5); 

% %% Box randomizing 
% 
% S1_norm_rand = S1_norm(randperm(length(S1_norm)), : );
% S2_norm_rand = S2_norm(randperm(length(S2_norm)), : );
% 
% % One Value per inner product (floor(M/D) values)
% 
% ost_rand = dot(S1_norm_rand, S2_norm_rand, 2);
% 
% % Single Value of Randomized OS:
% 
% os_rand = mean(ost_rand); %VALOR DE LA ORDINAL SYNCHRONISATION RANDOM CAJAS
% 
% ost_nrand = 2*(((ost_rand-min)/(1-min))-0.5);% Rand between -1 and 1
% os_nrand = 2*(((os_rand-min)/(1-min))-0.5);% Rand between -1 and 1