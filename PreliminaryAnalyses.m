clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Project part I - Preliminary analyses w/MONTHLY DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load MONTHLY data
[M, textM] = xlsread('MATLAB.xlsx', 'MONTHLY');                             % from 30/11/2011 to 30/11/2022

% store dates as string
D = datenum(textM(7:size(textM, 1), 1), 'dd/mm/yyyy');                           
% transform in dd/MM/yyyy format
D = datetime(D, 'ConvertFrom', 'datenum', 'Format', 'dd/MM/yyyy');

% store column labels into a string array
lab1 = textM(4, 2:size(textM, 2));                                          
lab = cell(16, 1);                                                             
% to access the contents of a cell, enclose indexes in curly brackets                
for j = 1:16
    lab{j, 1} = char(lab1(j));
    lab{j, 1} = lab{j, 1}(1:(size(lab{j, 1}, 2)-14));                       % put inside the cell the name less the last 14 digits
end

% restrict sample to a common one 
n = sum(isnan(M), 2);
p = M(n==0, :);                                                             % price indexes
d = D(n==0, :);                                                             % dates

% names
lab = {'EU - BENCHMARK', 'AT', 'BE', 'DK', 'FI', 'FR', 'DE', 'IE', 'IT', 'NL', 'NO', 'PT', 'ES', 'SE', 'CH', 'GB'}';

% plot PRICES of the BENCHMARK + 15 indexes 
figure
set(gcf, 'Position', [100, 100, 1300, 600]);
for j = 1:16
    subplot(4, 4, j);
    plot(datenum(D), p(:, j));                                                 
    tick3 = datenum(2011:3:2022, 1, 1);
    set(gca, 'xtick', tick3);
    datetick('x', 'yyyy', 'keepticks', 'keeplimits');                       % rearrange the dates
    title(lab{j, 1})
end
sgtitle('PRICES')
movegui('center')                                                           % to center the plot

% compute RETURNS                                                                                                    
rm = ((p(2:(size(p, 1)), :)./p(1:(size(p, 1)-1), :))-1)*100;

% plot RETURNS of the BENCHMARK + 15 indexes
figure                                                            
set(gcf,'Position', [100, 100, 1300, 600]);
for j = 1:16
    subplot(4, 4, j);
    plot(datenum(d(2:size(d, 1), :)), rm(:, j));
    tick3 = datenum(2011:3:2022, 1, 1);
    set(gca, 'xtick', tick3);
    datetick('x', 'yyyy', 'keepticks', 'keeplimits');
    title(lab{j, 1})
end
sgtitle('RETURNS')
movegui('center')

%% table for descriptive analyses
Mean = mean(rm)';
Median = median(rm)';
StDev = sqrt(var(rm))';
Min = min(rm)';
Max = max(rm)';
Skew = skewness(rm)';
Kurt = kurtosis(rm)';
TM = table(Mean, Median, StDev, Min, Max, Skew, Kurt, 'RowNames', lab);
format bank
TM

% store tables in Excel;
writetable(TM, 'TablesDescriptive.xlsx', 'Sheet', 'Monthly', 'WriteRowNames', 1);

%% evolution over time in relation to specific events
% MEAN                                                                      % considering the BENCHMARK 
MM = mean(rm(:, 1:16));                                                       
% plot the MEAN of each equity index (BENCHMARK + 15 indexes)               % considering the BENCHMARK 
figure
bar(1:16, MM);
set(gca, 'Xtick', 1:16, 'XTickLabel', lab);
title('TOTAL RETURNS')
movegui('center')

% VAR-COV MATRIX of each equity index (BENCHMARK + 15 indexes) 
MV = cov(rm);
% plot the STANDARD DEVIATION of each equity index (BENCHMARK + 15 indexes)
figure
bar(1:16, sqrt(diag(MV)));
set(gca, 'Xtick', 1:16, 'XtickLabel', lab);
title('TOTAL STANDARD DEVIATION')
movegui('center')

% split the sample in four parts and compare MEAN and VARIANCE moments,   
% the division was made to maintain about the same number of observations for each group 

% compute the MEAN for the 4 groups
MM1 = mean(rm(1:33, :));
MM2 = mean(rm(34:66, :));
MM3 = mean(rm(67:99, :));
MM4 = mean(rm(100:132, :));

% compute the VAR-COV MATRIX for the 4 groups
MV1 = cov(rm(1:33, :));
MV2 = cov(rm(34:66, :));
MV3 = cov(rm(67:99, :));
MV4 = cov(rm(100:132, :));

% plot the MEAN of each index for each group
figure
set(gcf, 'Position', [100, 100, 1300, 600]);                                % position of the figure window
subplot(2, 2, 1);
bar(1:16, MM1);
set(gca, 'Xtick', 1:16, 'XTickLabel', lab);
title('Returns 2011-2014')                                                
subplot(2, 2, 2);
bar(1:16, MM2);
set(gca, 'Xtick', 1:16, 'XTickLabel', lab);
title('Returns 2014-2017')                                                 
subplot(2, 2, 3);
bar(1:16, MM3);
set(gca, 'Xtick', 1:16, 'XTickLabel', lab);
title('Returns 2017-2020')                                                 
subplot(2, 2, 4);
bar(1:16, MM4);
set(gca, 'Xtick', 1:16, 'XTickLabel', lab);
title('Returns 2020-2022')                                                 
movegui('center')

% plot the STANDARD DEVIATION of each equity index (BENCHMARK + 15 assets) for each group
figure
set(gcf, 'Position', [100, 100, 1300, 600]);                                % position of the figure window
subplot(2, 2, 1);
bar(1:16, sqrt(diag(MV1)));
set(gca, 'Xtick', 1:16, 'XTickLabel', lab);
title('Standard Deviation 2011-2014')
subplot(2, 2, 2);
bar(1:16, sqrt(diag(MV2)));
set(gca, 'Xtick', 1:16, 'XTickLabel', lab);
title('Standard Deviation 2014-2017')
subplot(2, 2, 3);
bar(1:16, sqrt(diag(MV3)));
set(gca, 'Xtick', 1:16, 'XTickLabel', lab);
title('Standard Deviation 2017-2020')
subplot(2, 2, 4);
bar(1:16, sqrt(diag(MV4)));
set(gca, 'Xtick', 1:16, 'XTickLabel', lab);
title('Standard Deviation 2020-2022')
movegui('center')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the FULL MONTHLY BENCHMARK series
[BM, textBM] = xlsread('MATLAB.xlsx', 'FB MONTHLY');

% store dates as string
DB = datenum(textBM(7:size(textBM, 1), 1),'dd/mm/yyyy');                       
% transform dates in dd/MM/yyyy format
DB = datetime(DB, 'ConvertFrom', 'datenum', 'Format', 'dd/MM/yyyy');

% restrict sample to a common one  
nb = sum(isnan(BM), 2);
pb = BM(nb==0, :);                                                          % price indexes
db = DB(nb==0, :);                                                          % dates

% compute returns of the FULL MONTHLY BENCHMARK series                                                                                                  
rBm = ((pb(2:(size(pb, 1)), :)./pb(1:(size(pb, 1)-1), :))-1)*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
save PreliminaryAnalyses.mat;                                               % save workspace
