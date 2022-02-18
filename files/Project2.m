
%% Project: Detect Arbitrage from SPX Option Prices
% Author: Travis Choy
% Date Created: 10/6/2021
% Date Modified: 10/7/2021
% Class: FNCE 4820

%% Housekeeping

close all
clear
clc

%% Example Using Linprog Function

% c = [-1; -1/3]; % f in linprog
% b = [2; 1; 2; 1; -1; 2];
% A = [1 1; 1 1/4; 1 -1; -1/4 -1; -1 -1; -1 1];
% 
% [x fval] = linprog(c, A, b); % linprog always minimizes (to maximize use -c) using dual simplex method
% clc
% 
% Aeq = [1 1/4]; % Aeq and beq account for equality constraints...
% beq = [1/2]; % in this case, we're saying x1 + x2/4 = 1/2
% 
% [x fval] = linprog(c, A, b, Aeq, beq);
% clc
% 
% lb = [-1, -0.5]; % sets upper and lower bounds of xi
% ub = [1.5, 1.25]; % here we're saying -1 <= x1 <= 1.5 and -0.5 <= x2 <= 1.25
% 
% [x fval] = linprog(c, A, b, Aeq, beq, lb, ub);
% 
% disp(x)
% disp(fval)

%% Load in Data

[today, rawData, fileName, path] = getRawData();
[cleanData] = parseCBOEOptionData(rawData);
[mat, strikes, callBid, callAsk, putBid, putAsk] = singleMaturity(cleanData);

%% Type-A Arbitrage

A = zeros(length(strikes) + 2,length(strikes) * 4);
for i = 1:length(strikes)+2 % gives us i rows
    for j = 1:length(strikes)*2 % gives us j columns/2
        j = j*2;
        if j >= length(strikes)*2+1 % RHS, upper triangle
            if j/2-length(strikes) >= i
                if i == 1
                    A(i,j-1) = 0 - strikes(j/2-length(strikes));
                    A(i,j) = strikes(j/2-length(strikes));
                else
                    A(i,j-1) = strikes(i-1) - strikes(j/2-length(strikes));
                    A(i,j) = -strikes(i-1) + strikes(j/2-length(strikes));
                end
            end
        elseif j <= length(strikes)*2 % LHS, lower triangle
            if i > 1 & i < length(strikes)+2
                if i > j/2
                    A(i,j-1) = strikes(j/2) - strikes(i-1);
                    A(i,j) = strikes(i-1) - strikes(j/2);
                end
            elseif i == length(strikes)+2
                A(i,j-1) = -1;
                A(i,j) = 1;
            end
        end
    end
end
A = -A;
    
b = zeros(length(strikes)+2,1);

c = zeros(length(strikes) * 4,1);
for i = 1:length(strikes)
    c(i*2-1) = -callBid(i);
    c(i*2) = callAsk(i);
    c(i*2-1+length(strikes)*2) = -putBid(i);
    c(i*2+length(strikes)*2) = putAsk(i);
end

test = false;

if test == true
    temp = c(5);
    c(5) = c(6);
    c(6) = temp;
    clear temp
end

Aeq = zeros(1,length(strikes)*4); % [0x1 + 0x2 + ... ] = beq
beq = [0]; % Aeq = 0

ub = Inf(1,length(strikes)*4);
lb = zeros(1,length(strikes)*4);

[x, fval] = linprog(c, A, b, Aeq, beq, lb, ub);
%disp(x)
disp(fval)

if fval == 0
    fprintf('The optimal solution is zero. There is no type-A arbitrage.\n')
else
    fprintf('The optimal solution is greater than zero; there is type-A arbitrage.\n')
end

%% Type-B Arbitrage

cDual = -transpose(b); % They're all zeros but whatever
%b = zeros(length(strikes)*4,1);
ADual = -transpose(A);
bDual = transpose(c);

AeqDual = zeros(1,length(strikes)+2); % [0x1 + 0x2 + ... ] = beq
beqDual = [0]; % Aeq = 0

ubDual = Inf(1,length(strikes)+2);
lbDual = zeros(1,length(strikes)+2);

[y, fvalDual] = linprog(cDual, ADual, bDual, AeqDual, beqDual, lbDual, ubDual);
%disp(y)
disp(fvalDual)

fprintf('If there is an optimal solution, there is no type-B arbitrage.\n')

%% Closing Postion

r = 0.0025;
t = caldays(between(today, mat, 'days')) / 365;
c(end+1) = -exp(-r * t); % Short the zero-coupon bond
c(end+1) = exp(-r * t); % Long the zero-coupon bond

A(:,end+1) = 1;
A(:,end+1) = -1;
A(end,(end-1):end) = 0;

% We need a position to close. We were never given one so let's just
% pretend to close a random one:

position = 23; % position 23 would be our 12th ((23+1)/2) call bid (because it's even and in the first half of the data)

for i = 1:2
    ub(end+1) = inf;
    lb(end+1) = 0;
    Aeq(end+1) = 0;
end

ub(position) = 1; % We're going to bound the 5th x variable in the linprog function such that it is <= 1 and >= 1 that way it has to equal 1. Thus, we have a position in it.
lb(position) = 1;

[x, fval] = linprog(c, A, b, Aeq, beq, lb, ub);
%disp(x)
%disp(fval)

if x(end-1) > 0
    fprintf('You need to enter into %.2f zero-coupon bonds in the short position.\n', x(end-1))
end
if x(end) > 0
    fprintf('You need to enter into %.2f zero-coupon bonds in the long position.\n', x(end))
elseif x(end-1) <= 0
    fprintf("You don't need to enter into any zero-coupon bonds, the optimal solution is minimized without them.\n")
end

fprintf("Our solution is the minimized bid/ask spread. Final cost is %.2f.\n", fval)

%% Functions

function [today, rawData,fileName ,path] = getRawData()
% Author: Daniel Brown

%getRawData This function retrieves data from a
%   csv file downloaded from the CBOE. It does not do anything to the data
%   but creates a table
[fileName,path] = uigetfile('*.csv');
%path = '/Users/travischoy/Desktop/';
%fileName  = 'SPX Options ver 1.csv'; % SPX Options 2Feb21 335pm EST.csv
% Smaller sample size: SPX Options ver 1.csv

pathFile = strcat(path,fileName);
rawData = readtable(pathFile); % Previously had readtable(pathFile,'NumHeaderLines',3) but function wasn't working so I changed it. I don't know what 'NumHeaderLines',3 was supposed to do but...

% cannot think of how to retrieve the first three rows using readtable so
% reading the table again and throwing it away.  The file isn't so big that
% this is not an issue. Also, the file does not parse well so I have to
% jump through hoops to extract the actual date

fid = fopen(pathFile);

% get the third line which has the date
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);

commas = strfind(tline,',');
offset = commas(1,2);
today = datetime(extractBetween(tline,7,offset-1));

fclose(fid);
%dummyData = readtable(pathFile,'NumHeaderLines',1);
%rawDate = dummyData.Var1(2);
%commas = strfind(rawDate,',');
%offset = commas{1,1}(2);
%today = datetime(extractBetween(rawDate,7,offset-1));

end

function [cleanData] = parseCBOEOptionData(rawData)
% Author: Daniel Brown

%parseCBOEOptionData
%   This function takes raw data from the CBOE and creates a table
%   containing the maturities, strikes, call bid, call ask, put bid and put
%   ask prices.

ExpirationDate = datetime(rawData.ExpirationDate,'InputFormat','eee MMM dd yyyy');
IsMonthly = (~contains(rawData.Calls,"W"));
Strike = rawData.Strike;
CallBid = rawData.Bid;
CallAsk = rawData.Ask;
PutBid = rawData.Bid_1;
PutAsk = rawData.Ask_1;

cleanData = table(ExpirationDate,IsMonthly,Strike,CallBid,CallAsk,PutBid,PutAsk);
end

function [mat, strikes,callBid,callAsk,putBid,putAsk] = singleMaturity(cleanData)
% Author: Daniel Brown

%singleMaturity retrieves a set of options for a single maturity
%   The CBOE data set has weekly and monthly call and put option prices for
%   various maturities. This function retrieves a single set of prices

uniqueMats = unique(cleanData.ExpirationDate);
idxMat = listdlg('ListString',uniqueMats,'SelectionMode','Single','PromptString','Choose an Expiration Date');
mat = uniqueMats(idxMat);
matRows = (cleanData.ExpirationDate==mat);
matData = cleanData(matRows,:);

if (all(matData.IsMonthly))
    data = matData;
    disp ("Monthly Data only");
elseif (all(~matData.IsMonthly))
    data = matData;
    disp ("Weekly Data only");
else
    idx = listdlg('ListString',{'Weekly', 'Monthly'},'SelectionMode','Single','PromptString','Check Weekly or Monthly Options for Arbitrage?');
    if idx==1
        rows = (~matData.IsMonthly);
        disp("Using Weekly Data")
    else
        rows = (matData.IsMonthly);
        disp ("Using Monthly Data")
    end    
    data = matData(rows,:);
end

strikes = data.Strike;
callBid = data.CallBid;
callAsk = data.CallAsk;
putBid = data.PutBid;
putAsk = data.PutAsk;

end




