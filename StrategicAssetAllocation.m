clear;
clc;
load 'PreliminaryAnalyses.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project part II - Strategic Asset Allocation w/MONTHLY DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute MEAN and VOLATILITY of the BENCHMARK over the last 5 years  
rB = rm(:, 1);
MMrB = mean(rB(72:end, :));                                                 % from 30/11/2017 to 30/11/2022 
MVrB = var(rB(72:end, :));

% compute MEAN and VOLATILITY of the 15 indexes over the last 5 years
MMrm5 = mean(rm(72:end, 2:16));
MVrm5 = cov(rm(72:end, 2:16));                                                

%% use SAMPLE MOMENTS

% set the n° of portfolios for the EF
nport = 100;

% create an empty portfolio object assigning assets' names and moments
p1 = Portfolio('AssetList', lab(2:16), 'AssetMean', MMrm5, 'AssetCovar', MVrm5);

% set the budget UPPER and LOWER BOUNDS
p1 = setBudget(p1, 1, 1);                                                       

% set a limit to SHORT SELLING, in order to use the formula, in a way that does not change the EF weights' choice
% p1 is the UNCONSTRAINED PORTFOLIO 
p1 = setBounds(p1, -1.2);                                                      

% p2 is p1 subject to the NO SHORT SELLING CONSTRAINT
p2 = p1.setDefaultConstraints;

% compute the Global Minimum Variance portfolio
pwgtGMV = p1.estimateFrontierLimits('Min'); 
[prskGMV, pretGMV] = estimatePortMoments(p1,pwgtGMV);                       % NOT CONSTRAINED
pwgtGMV2 = p2.estimateFrontierLimits('Min'); 
[prskGMV2, pretGMV2] = estimatePortMoments(p2, pwgtGMV2);                   % CONSTRAINED

% compute the Max Sharpe portfolio
pwgtMS = estimateMaxSharpeRatio(p1); 
[prskMS, pretMS] = estimatePortMoments(p1, pwgtMS);                         % NOT CONSTRAINED
pwgtMS2 = estimateMaxSharpeRatio(p2); 
[prskMS2, pretMS2] = estimatePortMoments(p2, pwgtMS2);                      % CONSTRAINED

% pEW is the EQUALLY WEIGHTED portfolio
pEW = setInitPort(p1, 1/p1.NumAssets);
[prskEW, pretEW] = estimatePortMoments(pEW, pEW.InitPort);                      

% store weights of the EFFICIENT portfolios
pwgtNC = p1.estimateFrontier(nport);                                        % NOT CONSTRAINED                                      
pwgtC = p2.estimateFrontier(nport);                                         % CONSTRAINED

% plot the EFs and store RISKS and RETURNS of the EFFICIENT portfolios
figure
grid on
set(gcf, 'Position', [100, 100, 1300, 600]);
[priskNC, pretNC] = p1.plotFrontier(nport); 
hold on 
[priskC, pretC] = p2.plotFrontier(nport);
hold on
scatter(sqrt(diag(MVrm5)), MMrm5, 'filled', 'r'); 
text(sqrt(diag(MVrm5))+0.05, MMrm5+0.02, lab(2:16,1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12);
legend('UNCONSTRAINED', 'NO SHORT SELLING');
scatter(sqrt(MVrB), MMrB, 'filled', 'm', 'DisplayName', 'EU - BENCHMARK');
scatter(prskEW, pretEW, 'filled', 'c', 'DisplayName', 'EW');
% add GMVs and MSs to the figure 
scatter(prskGMV, pretGMV, 'filled', 'b', 'DisplayName', 'GMV-NC');
scatter(prskMS, pretMS, 'filled', 'g', 'DisplayName', 'MS-NC');
scatter(prskGMV2, pretGMV2, 'filled', 'y', 'DisplayName', 'GMV-C');
scatter(prskMS2, pretMS2, 'filled', 'k', 'DisplayName', 'MS-C');
title('EFFICIENT FRONTIERS')
xlim([0 15])
ylim([-0.5 5])
movegui('center')  

% draw weights of the two EFs w.r.t. RETURNS
figure
set(gcf, 'Position', [100, 100, 1300, 600]);
subplot(1, 2, 1);
a = area(pretNC, pwgtNC'.*(pwgtNC'>0));                      
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';   
hold on
a = area(pretNC, pwgtNC'.*(pwgtNC'<0));
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';                                                                                                          
title('UNCONSTRAINED')
xlabel('Returns');
ylabel('Weights');
axis tight
subplot(1, 2, 2);
a=area(pretC, pwgtC');
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';
title('NO SHORT SELLING')
xlabel('Returns');
ylabel('Weights');
axis tight
legend(lab(2:16));
lgnd = legend('show');
lgnd.Position(1) = 0.03;
lgnd.Position(2) = 0.385;
movegui('center') 

% draw weights of the EFs w.r.t. RISK 
figure
set(gcf, 'Position', [100, 100, 1300, 600]);
subplot(1, 2, 1);                                   
a = area(priskNC,pwgtNC'.*(pwgtNC'>0));
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';   
hold on
a = area(priskNC, pwgtNC'.*(pwgtNC'<0));
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r'; 
title('UNCONSTRAINED');
xlabel('Volatility');
ylabel('Weights');
axis tight
subplot(1, 2, 2);
a = area(priskC, pwgtC');
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';
ylim([0 1]);                                                          
title('NO SHORT SELLING')
xlabel('Volatility');
ylabel('Weights');
axis tight
legend(lab(2:16));
lgnd = legend('show');
lgnd.Position(1) = 0.03;
lgnd.Position(2) = 0.385;
movegui('center')  

%% use EQUILIBRIUM MOMENTS

% compute the EQUILIBRIUM RETURN and VARIANCE of each index over the last 5 years
rm5 = rm(72:end, 2:16);
% to find the EQUILIBRIUM RISKS and RETURNS of the BENCHMARK, we use the FULL MONTHLY BENCHMARK series
rBmeq = mean(rBm);
vBmeq = var(rBm);

% estimate the RETURNS using the CAPM and store the relevant quantities 

alpha = zeros(size(rm5, 2), 1);                                    
beta = zeros(size(rm5, 2), 1);                                       
rm2 = zeros(size(rm5, 2), 1);                                                                         
eqret = zeros(size(rm5, 2), 1);                    
resid = zeros(size(rm5, 1), size(rm5, 2));                                       

for i = 1:size(rm5, 2)                                                         
    out = regstats(rm5(:, i), rBm(575:end), 'linear', {'beta', 'r', 'rsquare', 'tstat'});  
    alpha(i, 1) = out.beta(1);                                                 
    beta(i, 1) = out.beta(2);                                                  
    rm2(i, 1) = out.rsquare;                                                   
    resid(:, i) = out.r;
    eqret(i, 1) = out.beta(2)*(rBmeq);
end

% compute the OPTIMAL portfolio by means of EQUILIBRIUM RETURNS
MMe = eqret'; 
MVe = beta*(beta')*vBmeq+diag(diag(cov(resid)));                              

% generate an UNCONSTRAINED empty portfolio and assign It the EQUILIBRIUM MOMENTS recovered from the BENCHMARK
p3 = Portfolio;
p3 = p3.setAssetMoments(MMe, MVe);
p3 = setBudget(p3, 1, 1);
p3 = setBounds(p3, -1.2);                                                     

% compute the Global Minimum Variance portfolio
pwgtGMV3 = p3.estimateFrontierLimits('Min');
[prskGMV3, pretGMV3] = estimatePortMoments(p3, pwgtGMV3);                   % UNCONSTRAINED

% compute Max Sharpe ratio
pwgtMS3 = estimateMaxSharpeRatio(p3);                                       
[prskMS3, pretMS3] = estimatePortMoments(p3, pwgtMS3);                      % UNCONSTRAINED         

% plot the UNCONSTRAINED EFs with EQUILIBRIUM MOMENTS
figure
set(gcf, 'Position', [100, 100, 1300, 600]);
[priskNC3, pretNC3] = p3.plotFrontier(nport);
hold on
p1.plotFrontier(nport);
hold on
scatter(sqrt(diag(MVrm5)), MMrm5, 'filled', 'r');
text(sqrt(diag(MVrm5))+0.05, MMrm5+0.02, lab(2:16, 1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12);
legend('CAPM EU - BENCHMARK', 'Sample moments');
scatter([prskGMV prskGMV3],[pretGMV pretGMV3],'filled','b','DisplayName', 'GMV');
scatter([prskMS prskMS3], [pretMS pretMS3], 'filled', 'g', 'DisplayName', 'MS');
title('UNCONSTRAINED EFFICIENT FRONTIERS')
hold on
scatter(prskEW, pretEW, 'filled', 'c', 'DisplayName', 'EW');
scatter(sqrt(MVrB), MMrB, 'filled', 'm', 'DisplayName', 'EU - BENCHMARK');
xlim([0 15])
ylim([-0.5 5])
movegui('center')  

% generate NO SHORT SELLING CONSTRAINED portfolio with EQUILIBRIUM MOMENTS
p4 = p3.setDefaultConstraints;                                                

% compute the Global Minimum Variance Portfolio
pwgtGMV4 = p4.estimateFrontierLimits('Min');
[prskGMV4, pretGMV4] = estimatePortMoments(p4, pwgtGMV4);                   % CONSTRAINED

% compute Max Sharpe ratio
pwgtMS4 = estimateMaxSharpeRatio(p4); 
[prskMS4, pretMS4] = estimatePortMoments(p4, pwgtMS4);                      % CONSTRAINED

% plot the NO SHORT SELLING CONSTRAINED EFs with EQUILIBRIUM MOMENTS 
figure
set(gcf, 'Position', [100, 100, 1300, 600]);
[priskC4, pretC4] = p4.plotFrontier(nport);
hold on
p2.plotFrontier(nport);
hold on
scatter(sqrt(diag(MVrm5)), MMrm5, 'filled', 'r');
text(sqrt(diag(MVrm5))+0.05, MMrm5+0.02, lab(2:16, 1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12);
legend('CAPM EU - BENCHMARK CONSTRAINED', 'Sample Moments CONSTRAINED', 'Location', 'best');
hold on
scatter([prskGMV2 prskGMV4], [pretGMV2 pretGMV4], 'filled', 'b', 'DisplayName', 'GMV');
scatter([prskMS2 prskMS4], [pretMS2 pretMS4], 'filled', 'g', 'DisplayName', 'MS');
scatter(prskEW, pretEW, 'filled', 'c', 'DisplayName', 'EW');
scatter(sqrt(MVrB), MMrB, 'filled', 'm', 'DisplayName', 'EU - BENCHMARK');
title('NO SHORT SELLING EFFICENT FRONTIERS')
xlim([0 15])
ylim([-0.5 2.5])
movegui('center') 

% store weights of the EFFICIENT portfolios 
pwgtNC3 = p3.estimateFrontier(nport);                                          
pwgtC4 = p4.estimateFrontier(nport);

% draw weights of the four EFs, the UNCONSTRAINED ones and the NO SHORT SELLING ones w.r.t. RETURNS
figure
set(gcf, 'Position', [100, 100, 1300, 600]);
subplot(1, 2, 1);
a = area(pretNC3, pwgtNC3'.*(pwgtNC3'>0));                                     
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';   
hold on
a = area(pretNC3, pwgtNC3'.*(pwgtNC3'<0));
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';                                                                                         
title('UNCONSTRAINED EU - BENCHMARK')
xlabel('Returns');
ylabel('Weights');
axis tight                                                                  % to limit the figures' axes                        
set(gcf, 'Position', [100, 100, 1300, 600]);
subplot(1, 2, 2);
a = area(pretC4, pwgtC4');
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';
title('NO SHORT SELLING EU - BENCHMARK');
xlabel('Returns');
ylabel('Weights');
axis tight
legend(lab(2:16));
lgnd=legend('show');
lgnd.Position(1) = 0.03;
lgnd.Position(2) = 0.385;
movegui('center')  

% draw weights of the four EFs, the UNCONSTRAINED ones and the NO SHORT SELLING ones w.r.t RISK
figure
set(gcf, 'Position', [100, 100, 1300, 600]);
subplot(1, 2, 1);                                   
a = area(priskNC3, pwgtNC3'.*(pwgtNC3'>0));
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';   
hold on
a = area(priskNC3,pwgtNC3'.*(pwgtNC3'<0));
a(8).FaceColor = 'y';                                                         
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r'; 
title('UNCONSTRAINED EU - BENCHMARK');
xlabel('Volatility');
ylabel('Weights');
axis tight
set(gcf, 'Position', [100, 100, 1300, 600]);
subplot(1, 2, 2);
a = area(priskC4, pwgtC4');
a(8).FaceColor = 'y';
a(9).FaceColor = 'b';
a(10).FaceColor = 'k';
a(11).FaceColor = 'r';
ylim([0 1])                                                       
title('NO SHORT SELLING EU - BENCHMARK')
xlabel('Volatility');
ylabel('Weights');
axis tight
legend(lab(2:16));
lgnd=legend('show');
lgnd.Position(1) = 0.03;
lgnd.Position(2) = 0.385;
movegui('center')

%% to find if the results are stable over time we redo the analysis on another 5-year window

% we select the period which goes from 31/10/2012 to 31/10/2017 (about 5 years)

% compute MEAN and VOLATILITY of the BENCHMARK over the new 5-year window
MMrB2 = mean(rB(11:71, :));                                                 % from 31/10/2012 to 31/10/2017                                            
MVrB2 = var(rB(11:71, :));

% compute MEAN and VOLATILITY of the assets over the new 5-year window
MMrm52 = mean(rm(11:71, 2:16));                                             % from 31/10/2012 to 31/10/2017 
MVrm52 = cov(rm(11:71, 2:16));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use SAMPLE MOMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create an empty portfolio object assigning assets' names and moments
p5 = Portfolio('AssetList', lab(2:16), 'AssetMean', MMrm52, 'AssetCovar', MVrm52);
% set the budget UPPER and LOWER BOUNDS
p5 = setBudget(p5, 1, 1);
% set a limit to SHORT SELLING,
% p5 is the UNCONSTRAINED portfolio 
p5 = setBounds(p5, -1.2);
% generate a new portfolio by imposing to p5 the NO SHORT SELLING CONSTRAINT 
p6 = p5.setDefaultConstraints;

% compute the Global Minimum Variance portfolio
pwgtGMV5 = p5.estimateFrontierLimits('Min'); 
[prskGMV5, pretGMV5] = estimatePortMoments(p5, pwgtGMV5);                   % NOT CONSTRAINED
pwgtGMV6 = p6.estimateFrontierLimits('Min'); 
[prskGMV6, pretGMV6] = estimatePortMoments(p6, pwgtGMV6);                   % CONSTRAINED

% compute the Max Sharpe
pwgtMS5 = estimateMaxSharpeRatio(p5); 
[prskMS5, pretMS5] = estimatePortMoments(p5, pwgtMS5);                      % NOT CONSTRAINED
pwgtMS6 = estimateMaxSharpeRatio(p6); 
[prskMS6, pretMS6] = estimatePortMoments(p6, pwgtMS6);                      % CONSTRAINED

% generate the EQUALLY WEIGHTED porfolio of this time window
pEW2 = setInitPort(p5, 1/p5.NumAssets);
[prskEW2, pretEW2] = estimatePortMoments(pEW2, pEW2.InitPort);

% store weights of the EFFICIENT portfolios
pwgtNC2 = p5.estimateFrontier(nport);                                          
pwgtC2 = p6.estimateFrontier(nport);

% plot the EFs and store RISKS and RETURNS of the EFFICIENT portfolios
figure
grid on 
set(gcf, 'Position', [100, 100, 1300, 600]);
[priskNC2, pretNC2] = p5.plotFrontier(nport); 
hold on 
[priskC2, pretC2] = p6.plotFrontier(nport);
hold on
scatter(sqrt(diag(MVrm52)), MMrm52, 'filled', 'r');
text(sqrt(diag(MVrm52))+0.05, MMrm52+0.02, lab(2:16,1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12);
legend('UNCONSTRAINED', 'NO SHORT SELLING');
scatter(sqrt(MVrB2), MMrB2, 'filled', 'm', 'DisplayName', 'EU - BENCHMARK');
scatter(prskEW2, pretEW2, 'filled', 'c', 'DisplayName', 'EW');
% add GMVs and TANs to the figure 
scatter(prskGMV5, pretGMV5, 'filled', 'b', 'DisplayName', 'GMV-NC');
scatter(prskMS5, pretMS5, 'filled', 'g','DisplayName', 'MS-NC');
scatter(prskGMV6, pretGMV6, 'filled', 'y', 'DisplayName', 'GMV-C');
scatter(prskMS6, pretMS6, 'filled', 'k', 'DisplayName', 'MS-C');
title('EFFICIENT FRONTIERS - Different time window')
xlim([0 15])
ylim([-0.5 5])
movegui('center') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use EQUILIBRIUM MOMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% estimate returns using the CAPM and store the relevant quantities 
rI2 = rm(11:71, 2:16);

alpha3 = zeros(size(rI2, 2), 1);                     
beta3 = zeros(size(rI2, 2), 1);                      
r2_3 = zeros(size(rI2, 2), 1);                        
eqret3 = zeros(size(rI2, 2), 1);                    
resid3 = zeros(size(rI2, 1), size(rI2, 2)); 

for i = 1:size(rI2, 2)
    out3 = regstats(rI2(:, i), rB(37:97), 'linear', {'beta', 'r', 'rsquare', 'tstat'}); 
    alpha3(i, 1) = out3.beta(1); 
    beta3(i, 1) = out3.beta(2);
    r2_3(i, 1) = out3.rsquare;
    resid3(:, i) = out3.r;
    eqret3(i, 1) = out3.beta(2)*(rBmeq);
end

% compute OPTIMAL portfolios by means of EQUILIBRIUM RETURNS  
MMe3 = eqret3'; 
MVe3 = beta3*(beta3')*vBmeq+diag(diag(cov(resid3)));                          

% generate 1 UNCONSTRAINED empty portfolio and assign It the EQUILIBRIUM MOMENTS recovered from the BENCHMARK
p7 = Portfolio;
p7 = p7.setAssetMoments(MMe3,MVe3);
p7 = setBudget(p7, 1,1);
p7 = setBounds(p7, -1.2);

% compute the Global Minimum Variance Portfolio
pwgtGMV7 = p7.estimateFrontierLimits('Min');
[prskGMV7, pretGMV7] = estimatePortMoments(p7, pwgtGMV7);

% compute Max Sharpe ratio
pwgtMS7 = estimateMaxSharpeRatio(p7); 
[prskMS7, pretMS7] = estimatePortMoments(p7, pwgtMS7); 

% plot the UNCONSTRAINED EFs with EQUILIBRIUM MOMENTS
figure
set(gcf, 'Position', [100, 100, 1300, 600]);
[priskNC4, pretNC4] = p7.plotFrontier(nport);
hold on
p5.plotFrontier(nport);
hold on
scatter(sqrt(diag(MVrm52)), MMrm52, 'filled', 'r');
text(sqrt(diag(MVrm52))+0.05, MMrm52+0.02, lab(2:16,1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12);
legend('CAPM EU - BENCHMARK', 'Sample moments');
scatter([prskGMV5 prskGMV7], [pretGMV5 pretGMV7], 'filled', 'b', 'DisplayName', 'GMV');
scatter([prskMS5 prskMS7], [pretMS5 pretMS7], 'filled', 'g', 'DisplayName', 'MS');
scatter(prskEW2, pretEW2, 'filled', 'c', 'DisplayName', 'EW');
scatter(sqrt(MVrB2), MMrB2, 'filled', 'm', 'DisplayName', 'EU - BENCHMARK');
title('UNCONSTRAINED EFFICIENT FRONTIERS - Different time window')
xlim([0 20]);
ylim([-0.5 5])
movegui('center')  

% generate NO SHORT SELLING CONSTRAINED portfolios with EQUILIBRIUM MOMENTS
p8 = p7.setDefaultConstraints;                                              % BENCHMARK values

% compute the Global Minimum Variance Portfolio
pwgtGMV8 = p8.estimateFrontierLimits('Min');
[prskGMV8, pretGMV8] = estimatePortMoments(p8, pwgtGMV8);

% compute Max Sharpe ratio
pwgtMS8 = estimateMaxSharpeRatio(p8); 
[prskMS8, pretMS8] = estimatePortMoments(p8, pwgtMS8);  

% plot the NO SHORT SELLING CONSTRAINED EFs with EQUILIBRIUM MOMENTS
figure
set(gcf, 'Position', [100, 100, 1300, 600]);
[priskC5, pretC5] = p8.plotFrontier(nport);
hold on
p6.plotFrontier(nport);
hold on
scatter(sqrt(diag(MVrm52)), MMrm52, 'filled', 'r');
text(sqrt(diag(MVrm52))+0.05, MMrm52+0.02, lab(2:16,1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12);
legend('CAPM EU - BENCHMARK CONSTRAINED', 'Sample Moments CONSTRAINED');
scatter([prskMS6 prskMS8], [pretMS6 pretMS8], 'filled', 'g', 'DisplayName', 'MS');
scatter([prskGMV6 prskGMV8], [pretGMV6 pretGMV8], 'filled', 'b', 'DisplayName', 'GMV');
hold on
scatter(prskEW2,pretEW2, 'filled', 'c', 'DisplayName', 'EW');
scatter(sqrt(MVrB2),MMrB2,'filled','m','DisplayName', 'EU - BENCHMARK');
title('NO SHORT SELLING EFFICIENT FRONTIERS - Different time window');
xlim([0 12])
ylim([-0.5 1.5])
movegui('center')   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
save StrategicAssetAllocation.mat;                                          % save workspace
