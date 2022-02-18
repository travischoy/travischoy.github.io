
%% Final Project
% Travis Choy
% FNCE 4820

%% Housekeeping

clear
close all
clc

%% Part 1

DowJones = xlsread('DowJones.xlsx','Historical Prices');
DowJones = DowJones(:,2:end);
DowSize = size(DowJones);

DowReturns = zeros(DowSize(1)-1,DowSize(2));

for j = 1:DowSize(2) % for each company
    for i = 1:(DowSize(1)-1) % for each date
        DowReturns(i,j) = log( DowJones(i+1,j) / DowJones(i,j) ); 
    end
end

% Part 1.1 and 1.2

V = cov(DowReturns(1:ceil(size(DowReturns,1)*2/3),2:end)); % Sample covariance

benchmark = var(DowReturns(1:ceil(size(DowReturns,1)*2/3),1));
sample_ret= DowReturns(1:ceil(size(DowReturns,1)*2/3),:);
sample_cov = cov(sample_ret);
beta = sample_cov(:,1) / benchmark;
V_SingleFactor = eye(30);

for j = 1:DowSize(2) 
    for i = 1:DowSize(2) 
        V_SingleFactor(i,j) = beta(i) * beta(j) * benchmark;
    end
end

V_SingleFactor = V_SingleFactor(2:end,2:end); % Single factor covariance

% Part 1.3

% Long Only
   
F = zeros(29,1);
A = -eye(29,29); % all weights > 0
b = zeros(29,1);
Aeq = ones(1,29);
beq = 1;
x_sample = quadprog(V, F, A, b, Aeq, beq); clc
x_SingleFactor = quadprog(V_SingleFactor, F, A, b, Aeq, beq); clc
    
PortRet_LongOnly_13_Sample = zeros(floor((size(DowReturns,1)/3)),1);
PortRet_LongOnly_13_SingleFactor = zeros(floor((size(DowReturns,1)/3)),1);

for i = 1:floor(size(DowReturns,1)/3)
    
    PortRet_LongOnly_13_Sample(i,1) = DowReturns(i+40,2:30) * x_sample;
    PortRet_LongOnly_13_SingleFactor(i,1) = DowReturns(i+40,2:30) * x_SingleFactor;
    
end

% Long Short
   
F = zeros(29,1);
A = zeros(1,29);
b = 0;
Aeq = ones(1,29);
beq = 1;
x_sample = quadprog(V, F, A, b, Aeq, beq); clc
x_SingleFactor = quadprog(V_SingleFactor, F, A, b, Aeq, beq); clc
    
PortRet_LongShort_13_Sample = zeros(floor((size(DowReturns,1)/3)),1);
PortRet_LongShort_13_SingleFactor = zeros(floor((size(DowReturns,1)/3)),1);

for i = 1:floor(size(DowReturns,1)/3)
    
    PortRet_LongShort_13_Sample(i,1) = DowReturns(i+40,2:end) * x_sample;
    PortRet_LongShort_13_SingleFactor(i,1) = DowReturns(i+40,2:end) * x_SingleFactor;
    
end

% Part 1.4

% Long Only

PortValueLongOnly_14_Sample = zeros(length(PortRet_LongOnly_13_Sample)+1,1);
PortValueLongOnly_14_Sample(1) = 1000;
for i = 1:length(PortRet_LongOnly_13_Sample)
    PortValueLongOnly_14_Sample(i+1) = PortValueLongOnly_14_Sample(i) * (1 + PortRet_LongOnly_13_Sample(i));
end

PortValueLongOnly_14_SingleFactor = zeros(length(PortRet_LongOnly_13_SingleFactor)+1,1);
PortValueLongOnly_14_SingleFactor(1) = 1000;
for i = 1:length(PortRet_LongOnly_13_SingleFactor)
    PortValueLongOnly_14_SingleFactor(i+1) = PortValueLongOnly_14_SingleFactor(i) * (1 + PortRet_LongOnly_13_SingleFactor(i));
end

% Long Short

PortValueLongShort_14_Sample = zeros(length(PortRet_LongShort_13_Sample)+1,1);
PortValueLongShort_14_Sample(1) = 1000;
for i = 1:length(PortRet_LongShort_13_Sample)
    PortValueLongShort_14_Sample(i+1) = PortValueLongShort_14_Sample(i) * (1 + PortRet_LongShort_13_Sample(i));
end

PortValueLongShort_14_SingleFactor = zeros(length(PortRet_LongShort_13_SingleFactor)+1,1);
PortValueLongShort_14_SingleFactor(1) = 1000;
for i = 1:length(PortRet_LongShort_13_SingleFactor)
    PortValueLongShort_14_SingleFactor(i+1) = PortValueLongShort_14_SingleFactor(i) * (1 + PortRet_LongShort_13_SingleFactor(i));
end

% Graph

PortReturnsDJI = DowReturns(ceil((size(DowReturns,1)*2/3))+1:end,1);

PortValueDJI = zeros(length(PortReturnsDJI)+1,1);
PortValueDJI(1) = 1000;
for i = 1:length(PortReturnsDJI)
    PortValueDJI(i+1) = PortValueDJI(i) * (1 + PortReturnsDJI(i));
end

days = (1:floor(size(DowReturns,1)/3)+1);
days = days';

figure(1)
plot(days, PortValueLongOnly_14_Sample)
hold on
plot(days, PortValueLongOnly_14_SingleFactor)
plot(days, PortValueLongShort_14_Sample)
plot(days, PortValueLongShort_14_SingleFactor)
plot(days, PortValueDJI)
title('Part 1.4 Plot')
xlabel('Days')
ylabel('Portfolio Value [$]')
legend('Long Only - Sample Covariance','Long Only - Single Factor Covariance','Long Short - Sample Covariance','Long Short - Single Factor Covariance','DJI','Location','northwest')
hold off

% Part 1.5

% Returns
PortRet_LongOnly_15_Sample = zeros(floor((size(DowReturns,1)/3)),1);
PortRet_LongOnly_15_SingleFactor = zeros(floor((size(DowReturns,1)/3)),1);
PortRet_LongShort_15_Sample = zeros(floor((size(DowReturns,1)/3)),1);
PortRet_LongShort_15_SingleFactor = zeros(floor((size(DowReturns,1)/3)),1);
for i = 1:floor((size(DowReturns,1)/3))
   
    returns = DowReturns(i:i+floor((size(DowReturns,1)*2/3)),2:end);
    H = cov(returns);
    benchmark = var(DowReturns(i:i+floor(size(DowReturns,1)*2/3),1));
    sample_ret= DowReturns(i:i+floor(size(DowReturns,1)*2/3),:);
    sample_cov = cov(sample_ret);
    beta = sample_cov(:,1) / benchmark;
    V_SingleFactor = eye(30);
    for j = 1:DowSize(2) 
        for k = 1:DowSize(2) 
            V_SingleFactor(k,j) = beta(k) * beta(j) * benchmark;
        end
    end
    V_SingleFactor = V_SingleFactor(2:end,2:end);
    F = zeros(29,1);
    A_longonly = -eye(29,29); % all weights > 0
    A_longshort = zeros(1,29);
    b_longonly = zeros(29,1);
    b_longshort = 0;
    Aeq = ones(1,29);
    beq = 1;
    x_longonly_sample = quadprog(H, F, A_longonly, b_longonly, Aeq, beq); clc
    x_longonly_singlefactor = quadprog(V_SingleFactor, F, A_longonly, b_longonly, Aeq, beq); clc
    x_longshort_sample = quadprog(H, F, A_longshort, b_longshort, Aeq, beq); clc
    x_longshort_singlefactor = quadprog(V_SingleFactor, F, A_longshort, b_longshort, Aeq, beq); clc
    
    PortRet_LongOnly_15_Sample(i,1) = DowReturns(i+40,2:30) * x_longonly_sample;
    PortRet_LongOnly_15_SingleFactor(i,1) = DowReturns(i+40,2:30) * x_longonly_singlefactor;
    PortRet_LongShort_15_Sample(i,1) = DowReturns(i+40,2:30) * x_longshort_sample;
    PortRet_LongShort_15_SingleFactor(i,1) = DowReturns(i+40,2:30) * x_longshort_singlefactor;
    
end

% Long Only

PortValueLongOnly_15_Sample = zeros(length(PortRet_LongOnly_15_Sample)+1,1);
PortValueLongOnly_15_Sample(1) = 1000;
for i = 1:length(PortRet_LongOnly_15_Sample)
    PortValueLongOnly_15_Sample(i+1) = PortValueLongOnly_15_Sample(i) * (1 + PortRet_LongOnly_15_Sample(i));
end

PortValueLongOnly_15_SingleFactor = zeros(length(PortRet_LongOnly_15_SingleFactor)+1,1);
PortValueLongOnly_15_SingleFactor(1) = 1000;
for i = 1:length(PortRet_LongOnly_15_SingleFactor)
    PortValueLongOnly_15_SingleFactor(i+1) = PortValueLongOnly_15_SingleFactor(i) * (1 + PortRet_LongOnly_15_SingleFactor(i));
end

% Long Short

PortValueLongShort_15_Sample = zeros(length(PortRet_LongShort_15_Sample)+1,1);
PortValueLongShort_15_Sample(1) = 1000;
for i = 1:length(PortRet_LongShort_15_Sample)
    PortValueLongShort_15_Sample(i+1) = PortValueLongShort_15_Sample(i) * (1 + PortRet_LongShort_15_Sample(i));
end

PortValueLongShort_15_SingleFactor = zeros(length(PortRet_LongShort_15_SingleFactor)+1,1);
PortValueLongShort_15_SingleFactor(1) = 1000;
for i = 1:length(PortRet_LongShort_15_SingleFactor)
    PortValueLongShort_15_SingleFactor(i+1) = PortValueLongShort_15_SingleFactor(i) * (1 + PortRet_LongShort_15_SingleFactor(i));
end

% Graph

figure(2)
plot(days, PortValueLongOnly_15_Sample)
hold on
plot(days, PortValueLongOnly_15_SingleFactor)
plot(days, PortValueLongShort_15_Sample)
plot(days, PortValueLongShort_15_SingleFactor)
plot(days, PortValueDJI)
title('Part 1.5 Plot')
xlabel('Days')
ylabel('Portfolio Value [$]')
legend('Long Only - Sample Covariance','Long Only - Single Factor Covariance','Long Short - Sample Covariance','Long Short - Single Factor Covariance','DJI','Location','northwest')
hold off

%% Part 2

% Part 2.1 and 2.2

V2 = cov(DowReturns(1:ceil((size(DowReturns,1)*2/3)),2:end));

% Part 2.3

TargetReturns = xlsread('DowJones.xlsx','Target Returns');
TargetReturns = TargetReturns';
PortTargetRet = 0.12;

% Part 2.4

% Long Only

leverage = 0;

V4 = zeros(size(V)*2);
V4(1:size(V,1),1:size(V,2)) = V2;

A = [-TargetReturns, zeros(1,length(TargetReturns)); -eye(length(TargetReturns)), -eye(length(TargetReturns)); zeros(length(TargetReturns),length(TargetReturns)), -eye(length(TargetReturns)); zeros(1,length(TargetReturns)), ones(1,length(TargetReturns))];
b = [-PortTargetRet; zeros(size(A,1)-2,1); leverage];
c = zeros(length(TargetReturns)*2,1);
Aeq = ones(1, length(TargetReturns)*2);
beq = 1;

Weights_LongOnly = quadprog(V4, c, A, b, Aeq, beq); clc

% Long Short

leverage = inf;

A = [-TargetReturns, zeros(1,length(TargetReturns)); -eye(length(TargetReturns)), -eye(length(TargetReturns)); zeros(length(TargetReturns),length(TargetReturns)), -eye(length(TargetReturns)); zeros(1,length(TargetReturns)), ones(1,length(TargetReturns))];
b = [-PortTargetRet; zeros(size(A,1)-2,1); leverage];
c = zeros(length(TargetReturns)*2,1);
Aeq = ones(1, length(TargetReturns)*2);
beq = 1;

Weights_LongShort = quadprog(V4, c, A, b, Aeq, beq); clc

% 30/130

leverage = 0.3;

A = [-TargetReturns, zeros(1,length(TargetReturns)); -eye(length(TargetReturns)), -eye(length(TargetReturns)); zeros(length(TargetReturns),length(TargetReturns)), -eye(length(TargetReturns)); zeros(1,length(TargetReturns)), ones(1,length(TargetReturns))];
b = [-PortTargetRet; zeros(size(A,1)-2,1); leverage];
c = zeros(length(TargetReturns)*2,1);
Aeq = ones(1, length(TargetReturns)*2);
beq = 1;

Weights_30_130 = quadprog(V4, c, A, b, Aeq, beq); clc

% Part 2.5

% Long Only

PortReturnsLongOnly_25 = zeros(floor((size(DowReturns,1)/3)),1);
for i = 1:floor((size(DowReturns,1)/3))
    PortReturnsLongOnly_25(i) = DowReturns(ceil((size(DowReturns,1)/3))+i,2:size(DowReturns,2)) * Weights_LongOnly(1:length(Weights_LongOnly)/2);
end

PortValueLongOnly_25 = zeros(length(PortReturnsLongOnly_25)+1,1);
PortValueLongOnly_25(1) = 1000;
for i = 1:length(PortReturnsLongOnly_25)
    PortValueLongOnly_25(i+1) = PortValueLongOnly_25(i) * (1 + PortReturnsLongOnly_25(i));
end

% Long Short

PortReturnsLongShort_25 = zeros(floor((size(DowReturns,1)/3)),1);
for i = 1:floor((size(DowReturns,1)/3))
    PortReturnsLongShort_25(i) = DowReturns(ceil((size(DowReturns,1)/3))+i,2:size(DowReturns,2)) * Weights_LongShort(1:length(Weights_LongShort)/2);
end

PortValueLongShort_25 = zeros(length(PortReturnsLongShort_25)+1,1);
PortValueLongShort_25(1) = 1000;
for i = 1:length(PortReturnsLongShort_25)
    PortValueLongShort_25(i+1) = PortValueLongShort_25(i) * (1 + PortReturnsLongShort_25(i));
end

% 30/130

PortReturns30130_25 = zeros(floor((size(DowReturns,1)/3)),1);
for i = 1:floor((size(DowReturns,1)/3))
    PortReturns30130_25(i) = DowReturns(ceil((size(DowReturns,1)/3))+i,2:size(DowReturns,2)) * Weights_30_130(1:length(Weights_30_130)/2);
end

PortValue30130_25 = zeros(length(PortReturns30130_25)+1,1);
PortValue30130_25(1) = 1000;
for i = 1:length(PortReturns30130_25)
    PortValue30130_25(i+1) = PortValue30130_25(i) * (1 + PortReturns30130_25(i));
end

% Graph

figure(3)
plot(days, PortValueLongOnly_25)
hold on
plot(days, PortValueLongShort_25)
plot(days, PortValue30130_25)
plot(days, PortValueDJI)
title('Part 2.5 Plot')
xlabel('Days')
ylabel('Portfolio Value [$]')
legend('Long Only','Long Short','30/130','DJI','Location','northwest')
hold off

% Part 2.6

PortReturnsLongOnly_26 = zeros(floor((size(DowReturns,1)/3)),1);
PortReturnsLongShort_26 = PortReturnsLongOnly_26;
PortReturns30130_26 = PortReturnsLongOnly_26;
PortValueLongOnly_26 = zeros(floor((size(DowReturns,1)/3))+1,1);
PortValueLongOnly_26(1) = 1000;
PortValueLongShort_26 = PortValueLongOnly_26;
PortValue30130_26 = PortValueLongOnly_26;

A = [-TargetReturns, zeros(1,length(TargetReturns)); -eye(length(TargetReturns)), -eye(length(TargetReturns)); zeros(length(TargetReturns),length(TargetReturns)), -eye(length(TargetReturns)); zeros(1,length(TargetReturns)), ones(1,length(TargetReturns))];
c = zeros(length(TargetReturns)*2,1);
Aeq = ones(1, length(TargetReturns)*2);
beq = 1;
V4 = zeros(size(V)*2);

for i = 1:floor(size(DowReturns,1)/3)

    V2 = cov(DowReturns(i:i+floor((size(DowReturns,1)*2/3)),2:end));
    
    % Long Only
    V4(1:size(V,1),1:size(V,2)) = V2;
    leverage = 0;
    b = [-PortTargetRet; zeros(size(A,1)-2,1); leverage];
    Weights_LongOnly_26 = quadprog(V4, c, A, b, Aeq, beq); clc
    PortReturnsLongOnly_26(i) = DowReturns(ceil((size(DowReturns,1)/3))+i,2:size(DowReturns,2)) * Weights_LongOnly_26(1:length(Weights_LongOnly_26)/2);
    PortValueLongOnly_26(i+1) = PortValueLongOnly_26(i) * (1 + PortReturnsLongOnly_26(i));
    
    % Long Short
    V4(1:size(V,1),1:size(V,2)) = V2;
    leverage = inf;
    b = [-PortTargetRet; zeros(size(A,1)-2,1); leverage];
    Weights_LongShort_26 = quadprog(V4, c, A, b, Aeq, beq); clc
    PortReturnsLongShort_26(i) = DowReturns(ceil((size(DowReturns,1)/3))+i,2:size(DowReturns,2)) * Weights_LongShort_26(1:length(Weights_LongShort_26)/2);
    PortValueLongShort_26(i+1) = PortValueLongShort_26(i) * (1 + PortReturnsLongShort_26(i));
    
    % 30/130
    V4(1:size(V,1),1:size(V,2)) = V2;
    leverage = 0.3;
    b = [-PortTargetRet; zeros(size(A,1)-2,1); leverage];
    Weights_30130_26 = quadprog(V4, c, A, b, Aeq, beq); clc
    PortReturns30130_26(i) = DowReturns(ceil((size(DowReturns,1)/3))+i,2:size(DowReturns,2)) * Weights_30130_26(1:length(Weights_30130_26)/2);
    PortValue30130_26(i+1) = PortValue30130_26(i) * (1 + PortReturns30130_26(i));

end

figure(4)
plot(days, PortValueLongOnly_26)
hold on
plot(days, PortValueLongShort_26)
plot(days, PortValue30130_26)
plot(days, PortValueDJI)
title('Part 2.6 Plot')
xlabel('Days')
ylabel('Portfolio Value [$]')
legend('Long Only','Long Short','30/130','DJI','Location','northwest')
hold off
