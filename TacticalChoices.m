clear;
clc;
load PreliminaryAnalyses.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project part III - Tactical Choices w/MONTHLY DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rI = rm(:, 2:16);                          
[r, c] = size(rI); 
rB = rm(:, 1);  

%% 1.	Evaluation of asset moments

%% SAMPLE MOMENTS (PreliminaryAnalyses)
% MM = mean(rm(:, 2:16));
% MV = cov(rm(:, 2:16));

w = 60;

MMS = zeros(r-w, c);
MVS = zeros(r-w, c, c); 
% sample estimator
                                                    
for j = w:(r-1)                                                             % r-1 as we use data up to r-1 to allocate 1-step-ahead
    % sample estimator
    MMS(j+1-w, :) = mean(rI(j-w+1:j, :));                                       
    MVS(j+1-w, :, :) = cov(rI(j-w+1:j, :));                                      
end
        
%% EWMA METHOD 

% smoothing factor
lambda = 0.95;                                                                       
% pre-allocating for MEAN and COVARIANCE
EWMean = zeros(size(rm, 1)-w, size(rm, 2)-1); 
% drop w-1 to be coherent with the ROLLING METHOD approach
% we still have a FORWARD-LOOKING perspective
EWVar = zeros(size(rm, 1)-w, size(rm, 2)-1, size(rm, 2)-1); 
% this is a three-dimensional array as we have COVARIANCES
    
% initialization with sample moments of the first rolling sample
EWMean(1, :) = mean(rm(1:w, 2:size(rm, 2)));                                    
EWVar(1, :, :) = cov(rm(1:w, 2:size(rm, 2)));                                     
for j = (w+1):(size(rm, 1)-1)
    EWMean(j-w+1, :) = lambda*EWMean(j-w, :)+(1-lambda)*rm(j, 2:size(rm, 2));
    EWVar(j-w+1, :, :) = lambda*squeeze(EWVar(j-w, :, :))+(1-lambda)*(rm(j, 2:size(rm, 2))'*rm(j, 2:size(rm, 2)));
end

%% 2.	Evaluation of portfolio weights and realized returns

%% EWMA METHOD 

%--------------------------------------------------------------------------
% Portfolio 4 - NO SHORT SELLING and UPPER BOUND
%--------------------------------------------------------------------------

p4 = Portfolio;                       
p4 = p4.setAssetList(lab(2:16));
p4 = p4.setDefaultConstraints;
p4 = p4.setInitPort(0);                                                     % set initial portfolio to 0
I = zeros(size(rI, 2), 2);                                                  % variable to update initPort
                                                                            % the first column represents the GMV, the second column represents the MaxSharpe
p4 = p4.setBounds(0, 0.13);                                                 % upper bound at 13% for each index

% pre allocations
PortW_S4 = zeros(2, r-w, c);                                                % pre-allocation for WEIGHTS of GMV and MS over time 
PortRet_S4 = zeros(r-w, 2);                                                 % pre-allocation for RETURNS of GMV and MS over time

% loop for GMV and MS portfolios
for j = (w+1):r
    % set moments
    p4 = p4.setAssetMoments(EWMean(j-w,:), squeeze(EWVar(j-w, :, :)));
    % GMV
    p4 = p4.setInitPort(I(:, 1));
    pwgtGMV = p4.estimateFrontierLimits('Min');
    PortW_S4(1,j-w,:) = pwgtGMV';
    PortRet_S4(j-w,1) = rI(j, :)*(pwgtGMV);
    I(:, 1) = pwgtGMV;                             
    % MaxSharpe
    p4 = p4.setInitPort(I(:, 2));
    pwgtMS = p4.estimateFrontierLimits('Max');
    PortW_S4(2, j-w, :) = pwgtMS';
    PortRet_S4(j-w, 2) = rI(j, :)*(pwgtMS);
    I(:, 2) = pwgtMS;
end

%--------------------------------------------------------------------------
% Portfolio 5 - SHORT SELLING and GROUP CONSTRAINTS
%--------------------------------------------------------------------------

p5 = Portfolio(p4);                                                         % set initial portfolio object

p5 = p5.setBounds(-0.2);                                                    % LOWER BOUND => minimum investment for each index
GA = [1 1 0 0 1 1 0 0 1 0 0 0 0 1 1];                                       % CORE
GB = [0 0 0 0 0 0 1 1 0 0 1 1 0 0 0];                                       % GIIPS 
GC = [0 0 1 1 0 0 0 0 0 1 0 0 1 0 0];                                       % NORDIC
p5 = setGroups(p5, GA, 0.4, 1);                                             % BOUNDS for CORE countries
p5 = addGroups(p5, GB, 0, 0.4);                                             % BOUNDS for GIIPS countries 
p5 = addGroups(p5, GC, 0, 0.4);                                             % BOUNDS for NORDIC countries

% pre allocations
p5Wgt = zeros(2, r-w, c);                                                   % pre-allocation for WEIGHTS
p5Ret = zeros(r-w, 2);                                                      % pre-allocation for RETURNS
I = zeros(size(rI, 2), 2);                                                  % variable to update initPort

for j = (w+1):r
    % set asset moments
    p5 = p5.setAssetMoments(EWMean(j-w, :), squeeze(EWVar(j-w, :, :))); 
    % GMV Portfolio
    p5 = p5.setInitPort(I(:, 1));
    pwgtGMV = p5.estimateFrontierLimits('Min');
    PortW_S5(1, j-w, :) = pwgtGMV';
    PortRet_S5(j-w, 1) = rI(j, :)*(pwgtGMV);
    I(:, 1) = pwgtGMV;
    % Max Sharpe Portfolio
    p5 = p5.setInitPort(I(:, 2));
    pwgtMS = estimateMaxSharpeRatio(p5);
    PortW_S5(2, j-w, :) = pwgtMS';
    PortRet_S5(j-w, 2) = rI(j, :)*(pwgtMS);
    I(:, 2) = pwgtMS; 
end

%--------------------------------------------------------------------------
% Portfolio 6 - TURNOVER CONSTRAINT
%--------------------------------------------------------------------------

p6 = Portfolio(p5);
p6 = p6.setInitPort(0);

% pre-allocation for RETURNS and WEIGHTS
PortRet_S6 = zeros(r-w, 2);
PortW_S6 = zeros(2, r-w, c);
j = w+1;
p6 = p6.setAssetMoments(EWMean(j-w, :), squeeze(EWVar(j-w, :, :)));
% GMV
pwgtGMV = p6.estimateFrontierLimits('Min');
PortW_S6(1, j-w, :) = pwgtGMV';
PortRet_S6(j-w, 1) = rI(j, :)*(pwgtGMV);
GMV0W = pwgtGMV;
% MaxSharpe
pwgtMS = estimateMaxSharpeRatio(p6);
PortW_S6(2, j-w, :) = pwgtMS';
PortRet_S6(j-w, 2) = rI(j, :)*(pwgtMS);
MS0W = pwgtMS;

% portfolio Turnover 
% add turnover from the second iteration - otherwise the turnover will not work
p6 = Portfolio(p6, 'Turnover', 0.05);                                       % we cannot change more than 5% of the wealth allocated to each asset

for j = (w+2):r
    % set moments
    p6 = p6.setAssetMoments(EWMean(j-w, :), squeeze(EWVar(j-w, :, :)));
    % GMV
    p6 = Portfolio(p6, 'InitPort', GMV0W);                                  % set initial portfolio
    pwgtGMV = p6.estimateFrontierLimits('Min');
    PortW_S6(1, j-w, :) = pwgtGMV';
    PortRet_S6(j-w, 1) = rI(j, :)*(pwgtGMV);
    GMV0W = pwgtGMV;
    % MaxSharpe
    p6 = Portfolio(p6, 'InitPort', MS0W);                                   % set initial portfolio
    try
        pwgtMS = estimateMaxSharpeRatio(p6);
        PortW_S6(2, j-w, :) = pwgtMS';
        PortRet_S6(j-w, 2) = rI(j, :)*(pwgtMS);
        MS0W = pwgtMS;
    catch
        try
        p6 = Portfolio(p6, 'Turnover', 0.15);
        pwgtMS = estimateMaxSharpeRatio(p6);
        PortW_S6(2, j-w, :) = pwgtMS';
        PortRet_S6(j-w, 2) = rI(j, :)*(pwgtMS);
        MS0W = pwgtMS;
        p6 = Portfolio(p6, 'Turnover', 0.05);
        catch
        p6 = setTurnover(p6, []);
        pwgtMS = p6.estimateFrontierLimits('Max');
        PortW_S6(2, j-w, :) = pwgtMS';
        PortRet_S6(j-w, 2) = rI(j, :)*(pwgtMS);
        MS0W = pwgtMS;
        p6 = Portfolio(p6, 'Turnover', 0.05);
        end
    end
end

%% plot weights' evolution over time - EWMA METHOD
figure
set(gcf, 'Position', [100, 100, 1300, 600])
subplot(3, 2, 1)
a = area(datenum(D(w+1:r, :)), squeeze(PortW_S4(1, :, :)));                                  
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';  
ylim([0 1])
xlim([datenum(D(w+1, :)) datenum(D(r, :))]);
title('No Short Selling and Upper Bound - GMV')
tick3 = datenum(2016:2022, 1, 1);
set(gca, 'xtick', tick3);
datetick('x', 'yyyy', 'keepticks', 'keeplimits')
subplot(3, 2, 2)
a = area(datenum(D(w+1:r, :)),squeeze(PortW_S4(2, :, :)));                  % 2 for MS
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';  
ylim([0 1])
xlim([datenum(D(w+1, :)) datenum(D(r, :))]);
title('No Short Selling and Upper Bound - MS')
tick3 = datenum(2016:2022, 1, 1);
set(gca, 'xtick', tick3);
datetick('x', 'yyyy', 'keepticks', 'keeplimits')
subplot(3, 2, 3)
a = area(datenum(D(w+1:r, :)), squeeze(PortW_S5(1, :, :)));
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';  
datetick('x', 'mmm-yy')
xlim([datenum(D(w+1, :)) datenum(D(r, :))]);
title('Short selling and group constraints - GMV')
tick3 = datenum(2016:2022, 1, 1);
set(gca, 'xtick', tick3);
datetick('x', 'yyyy', 'keepticks', 'keeplimits')
subplot(3, 2, 4)
a = area(datenum(D(w+1:r, :)), squeeze(PortW_S5(2, :, :)));
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';  
xlim([datenum(D(w+1, :)) datenum(D(r, :))]);
title('Short selling and group constraints - MS')
tick3 = datenum(2016:2022, 1, 1);
set(gca, 'xtick', tick3);
datetick('x', 'yyyy', 'keepticks', 'keeplimits')
subplot(3, 2, 5)
a = area(datenum(D(w+1:r, :)), squeeze(PortW_S6(1, :, :)));
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';  
xlim([datenum(D(w+1, :)) datenum(D(r, :))]);
title('Turnover constraints - GMV')
tick3 = datenum(2016:2022, 1, 1);
set(gca, 'xtick', tick3);
datetick('x', 'yyyy', 'keepticks', 'keeplimits')
subplot(3, 2, 6)
a = area(datenum(D(w+1:r, :)), squeeze(PortW_S6(2, :, :)));
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';  
xlim([datenum(D(w+1, :)) datenum(D(r, :))]);
title('Turnover constraints - MS')
tick3 = datenum(2016:2022, 1, 1);
set(gca, 'xtick', tick3);
datetick('x', 'yyyy', 'keepticks', 'keeplimits')
legend(lab(2:16))
lgnd = legend('show');
lgnd.Position(1) = 0.03;
lgnd.Position(2) = 0.385;
movegui('center')
sgtitle('EWMA METHOD') 

%% SAMPLE METHOD

%--------------------------------------------------------------------------
% Portfolio 7 - NO SHORT SELLING and UPPER BOUND
%--------------------------------------------------------------------------

p7 = Portfolio;                       
p7 = p7.setAssetList(lab(2:16));                                            % assets' name
p7 = p7.setDefaultConstraints;                                              % NO SHORT SELLING CONSTRAINT
p7 = p7.setInitPort(0);                                                     % set initial portfolio to 0
I = zeros(size(rI, 2), 2);                                                  % variable to update initPort
                                                                            % the first column represents the GMV, the second column represents the Max Sharpe
p7 = p7.setBounds(0, 0.13);                                                 % upper bound at 13% for each index (1/15*2)

% pre allocations
PortW_S7 = zeros(2, r-w, c);                                                % pre-allocation for WEIGHTS of GMV and MS over time 
PortRet_S7 = zeros(r-w, 2);                                                 % pre-allocation for RETURNS of GMV and MS over time

%% loop for GMV and MS portfolios 

for j = (w+1):r
    % set moments
    p7 = p7.setAssetMoments(MMS(j-w, :), squeeze(MVS(j-w, :, :)));               
    % GMV
    p7 = p7.setInitPort(I(:, 1));
    pwgtGMV = p7.estimateFrontierLimits('Min');                               
    PortW_S7(1, j-w, :) = pwgtGMV';
    PortRet_S7(j-w, 1) = rI(j, :)*(pwgtGMV);
    I(:, 1) = pwgtGMV;                                                      % I is a variable created to show the WEIGHTS of GMV and MS portfolios, updating at the end of each strategy 
    % MaxSharpe
    p7 = p7.setInitPort(I(:, 2));
    pwgtMS = p7.estimateFrontierLimits('Max');                                
    PortW_S7(2, j-w, :) = pwgtMS';
    PortRet_S7(j-w, 2) = rI(j, :)*(pwgtMS);
    I(:,2) = pwgtMS;
end

%--------------------------------------------------------------------------
% Portfolio 8 - SHORT SELLING and GROUP CONSTRAINTS
%--------------------------------------------------------------------------

p8 = Portfolio(p7);                                                         % set initial portfolio object

p8 = p8.setBounds(-0.2);                                                    % LOWER BOUND => minimum investment for each index


GA = [1 1 0 0 1 1 0 0 1 0 0 0 0 1 1];                                       % CORE 
GB = [0 0 0 0 0 0 1 1 0 0 1 1 0 0 0];                                       % GIIPS 
GC = [0 0 1 1 0 0 0 0 0 1 0 0 1 0 0];                                       % NORDIC
p8 = setGroups(p8, GA, 0.4, 1);                                             % BOUNDS for CORE countries
p8 = addGroups(p8, GB, 0, 0.4);                                             % BOUNDS for GIIPS countries 
p8 = addGroups(p8, GC, 0, 0.4);                                             % BOUNDS for NORDIC countries

% pre allocations
p8Wgt=zeros(2,r-w,c);                                                       % pre-allocation for WEIGHTS
p8Ret=zeros(r-w,2);                                                         % pre-allocation for RETURNS
I=zeros(size(rI,2),2);                                                      % variable to update initPort

for j = (w+1):r
    % set asset moments
    p8 = p8.setAssetMoments(MMS(j-w, :), squeeze(MVS(j-w, :, :))); 
    % GMV Portfolio
    p8 = p8.setInitPort(I(:, 1));
    pwgtGMV = p8.estimateFrontierLimits('Min');
    PortW_S8(1, j-w, :) = pwgtGMV';
    PortRet_S8(j-w, 1) = rI(j, :)*(pwgtGMV);
    I(:, 1) = pwgtGMV;
    % Max Sharpe Portfolio
    p8 = p8.setInitPort(I(:, 2));
    pwgtMS = estimateMaxSharpeRatio(p8);
    PortW_S8(2, j-w, :) = pwgtMS';
    PortRet_S8(j-w, 2) = rI(j, :)*(pwgtMS);
    I(:, 2) = pwgtMS; 
end

%--------------------------------------------------------------------------
% Portfolio 9 - TURNOVER CONSTRAINT
%--------------------------------------------------------------------------

p9 = Portfolio(p8);
p9 = p9.setInitPort(0);

% pre-allocation for RETURNS and WEIGHTS
PortRet_S9 = zeros(r-w, 2);
PortW_S9 = zeros(2, r-w, c);
j=w+1;
p9 = p9.setAssetMoments(MMS(j-w, :), squeeze(MVS(j-w, :, :)));
% GMV
pwgtGMV = p9.estimateFrontierLimits('Min');
PortW_S9(1, j-w, :) = pwgtGMV';
PortRet_S9(j-w, 1) = rI(j, :)*(pwgtGMV);
GMV0W = pwgtGMV;
% MaxSharpe
pwgtMS = estimateMaxSharpeRatio(p9);
PortW_S9(2, j-w, :) = pwgtMS';
PortRet_S9(j-w, 2) = rI(j, :)*(pwgtMS);
MS0W = pwgtMS;

% add turnover from the second iteration - otherwise the turnover will not work
p9 = Portfolio(p9, 'Turnover', 0.05);                                       % we cannot change more than 5% of the wealth allocated to each asset

for j = (w+2):r
    % set moments
    p9 = p9.setAssetMoments(MMS(j-w, :), squeeze(MVS(j-w, :, :)));
    % GMV
    p9 = Portfolio(p9, 'InitPort', GMV0W);                                  % set initial portfolio
    pwgtGMV = p9.estimateFrontierLimits('Min');
    PortW_S9(1, j-w, :) = pwgtGMV';
    PortRet_S9(j-w, 1) = rI(j, :)*(pwgtGMV);
    GMV0W = pwgtGMV;
    % MaxSharpe
    p9 = Portfolio(p9, 'InitPort', MS0W);                                   % set initial portfolio
    try                                                                     
        pwgtMS = estimateMaxSharpeRatio(p9);
        PortW_S9(2, j-w, :) = pwgtMS';
        PortRet_S9(j-w, 2) = rI(j, :)*(pwgtMS);
        MS0W = pwgtMS;
    catch
        try
        p9 = Portfolio(p9, 'Turnover', 0.15);
        pwgtMS = estimateMaxSharpeRatio(p9);
        PortW_S9(2, j-w, :) = pwgtMS';
        PortRet_S9(j-w, 2) = rI(j, :)*(pwgtMS);
        MS0W = pwgtMS;
        p9 = Portfolio(p9, 'Turnover', 0.05);
        catch
        p9 = setTurnover(p9, []);
        pwgtMS = p9.estimateFrontierLimits('Max');
        PortW_S9(2, j-w, :) = pwgtMS';
        PortRet_S9(j-w, 2) = rI(j, :)*(pwgtMS);
        MS0W = pwgtMS;
        p9 = Portfolio(p9, 'Turnover', 0.05);
        end
    end
end

%% plot weights' evolution over time - SAMPLE MOMENTS
figure
set(gcf, 'Position', [100, 100, 1300, 600])
subplot(3, 2, 1)
a = area(datenum(D(w+1:r, :)), squeeze(PortW_S7(1, :, :)));                                  
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';  
ylim([0 1])
xlim([datenum(D(w+1, :)) datenum(D(r, :))]);
title('No Short Selling and Upper Bound - GMV')
tick3 = datenum(2016:2022, 1, 1);
set(gca, 'xtick', tick3);
datetick('x', 'yyyy', 'keepticks', 'keeplimits')
subplot(3, 2, 2)
a = area(datenum(D(w+1:r, :)), squeeze(PortW_S7(2, :, :)));                 % 2 for MS
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';  
ylim([0 1])
xlim([datenum(D(w+1, :)) datenum(D(r, :))]);
title('No Short Selling and Upper Bound - MS')
tick3 = datenum(2016:2022, 1, 1);
set(gca, 'xtick', tick3);
datetick('x', 'yyyy', 'keepticks', 'keeplimits')
subplot(3, 2, 3)
a = area(datenum(D(w+1:r, :)), squeeze(PortW_S8(1, :, :)));
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';  
datetick('x', 'mmm-yy')
xlim([datenum(D(w+1, :)) datenum(D(r, :))]);
title('Short selling and group constraints - GMV')
tick3 = datenum(2016:2022, 1, 1);
set(gca, 'xtick', tick3);
datetick('x', 'yyyy', 'keepticks', 'keeplimits')
subplot(3, 2, 4)
a = area(datenum(D(w+1:r, :)), squeeze(PortW_S8(2, :, :)));
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';  
xlim([datenum(D(w+1, :)) datenum(D(r, :))]);
title('Short selling and group constraints - MS')
tick3 = datenum(2016:2022, 1, 1);
set(gca, 'xtick', tick3);
datetick('x', 'yyyy', 'keepticks', 'keeplimits')
subplot(3, 2, 5)
a = area(datenum(D(w+1:r, :)), squeeze(PortW_S9(1, :, :)));
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';  
xlim([datenum(D(w+1, :)) datenum(D(r, :))]);
title('Turnover constraints - GMV')
tick3 = datenum(2016:2022, 1, 1);
set(gca, 'xtick', tick3);
datetick('x', 'yyyy', 'keepticks', 'keeplimits')
subplot(3, 2, 6)
a = area(datenum(D(w+1:r, :)), squeeze(PortW_S9(2, :, :)));
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';  
xlim([datenum(D(w+1, :)) datenum(D(r, :))]);
title('Turnover constraints - MS')
tick3 = datenum(2016:2022, 1, 1);
set(gca, 'xtick', tick3);
datetick('x', 'yyyy', 'keepticks', 'keeplimits')
legend(lab(2:16))
lgnd = legend('show');
lgnd.Position(1) = 0.03;
lgnd.Position(2) = 0.385;
movegui('center')
sgtitle('SAMPLE MOMENTS') 

% compute cumulated returns of NO SHORT SELLING and UPPER BOUND strategies
CRS4GMV = cumprod(PortRet_S4(:, 1)/100+1)-1;
CRS7GMV = cumprod(PortRet_S7(:, 1)/100+1)-1;
CRS4MS = cumprod(PortRet_S4(:, 2)/100+1)-1;
CRS7MS = cumprod(PortRet_S7(:, 2)/100+1)-1;

%% plot cumulated returns of NO SHORT SELLING and UPPER BOUND strategies
figure
set(gcf, 'Position', [100, 100, 1300, 600])
plot(datenum(D(w+1:r, :)), CRS4GMV, Color = 'r')
hold on
plot(datenum(D(w+1:r, :)), CRS7GMV, Color = 'g')
hold on
plot(datenum(D(w+1:r, :)), CRS4MS, Color = 'b')
hold on
plot(datenum(D(w+1:r, :)), CRS7MS, Color = 'm')
hold on
legend('GMV Sample Moments no short selling and upper bound', 'GMV EWMA no short selling and upper bound', 'MS Sample Moments no short selling and upper bound', 'MS EWMA no short selling and upper bound', 'Location', 'NorthWest');
tick3 = datenum(2016:2022, 1, 1);
set(gca, 'xtick', tick3);
datetick('x', 'yyyy', 'keepticks', 'keeplimits')
title('Cumulated Returns no short selling and upper bound strategies')

% compute cumulated returns of SHORT SELLING and GROUP CONSTRAINT strategies
CRS5GMV = cumprod(PortRet_S5(:, 1)/100+1)-1;
CRS8GMV = cumprod(PortRet_S8(:, 1)/100+1)-1;
CRS5MS = cumprod(PortRet_S5(:, 2)/100+1)-1;
CRS8MS = cumprod(PortRet_S8(:, 2)/100+1)-1;

%% plot cumulated returns of SHORT SELLING and GROUP CONSTRAINT strategies
figure
set(gcf, 'Position', [100, 100, 1300, 600])
plot(datenum(D(w+1:r, :)), CRS5GMV, Color = 'r')
hold on
plot(datenum(D(w+1:r, :)), CRS8GMV, Color = 'g')
hold on
plot(datenum(D(w+1:r, :)), CRS5MS, Color = 'b')
hold on
plot(datenum(D(w+1:r, :)), CRS8MS, Color = 'm')
hold on
legend('GMV Sample Moments short selling and group constrained', 'GMV EWMA short selling and group constrained', 'MS Sample Moments short selling and group constrained', 'MS EWMA short selling and group constrained', 'Location', 'NorthWest');
tick3 = datenum(2016:2022, 1, 1);
set(gca, 'xtick', tick3);
datetick('x', 'yyyy', 'keepticks', 'keeplimits')
title('Cumulated Returns short selling and group constraint strategies')

% compute cumulated returns of TURNOVER strategies
CRS6GMV = cumprod(PortRet_S6(:, 1)/100+1)-1;
CRS9GMV = cumprod(PortRet_S9(:, 1)/100+1)-1;
CRS6MS = cumprod(PortRet_S6(:, 2)/100+1)-1;
CRS9MS = cumprod(PortRet_S9(:, 2)/100+1)-1;

%% plot cumulated returns of TURNOVER strategies
figure
set(gcf, 'Position', [100, 100, 1300, 600])
plot(datenum(D(w+1:r, :)), CRS6GMV, Color = 'r')
hold on
plot(datenum(D(w+1:r, :)), CRS9GMV, Color = 'g')
hold on
plot(datenum(D(w+1:r, :)), CRS6MS, Color = 'b')
hold on
plot(datenum(D(w+1:r, :)), CRS9MS, Color = 'm')
hold on
legend('GMV Sample Moments Turnover 5% ', 'GMV EWMA Turnover 5%', 'MS Sample Moments Turnover 5% ', 'MS EWMA Turnover 5%', 'Location', 'NorthWest');
tick3 = datenum(2016:2022, 1, 1);
set(gca, 'xtick', tick3);
datetick('x', 'yyyy', 'keepticks', 'keeplimits')
title('Cumulated Returns Turnover strategies')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
save TacticalChoices.mat                                                    % save workspace