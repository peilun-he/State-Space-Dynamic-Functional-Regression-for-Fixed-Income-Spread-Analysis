clear;

% Add search path
addpath(genpath(pwd));

UK = readtable("UK_12Con.csv");
date = table2array(UK(:, 1));
date = datetime(date, 'InputFormat', 'dd/MM/yyyy');
UK = table2array(UK(:, 2: end));
UK_out = readtable("UK_12Con_Out.csv"); 
date_out = table2array(UK_out(:, 1));
date_out = datetime(date_out, 'InputFormat', 'dd/MM/yyyy');
UK_out = table2array(UK_out(:, 2: end));

US_original = readtable("US_12Con.csv"); 
US_original = table2array(US_original(:, 2: end));
US_original_out = readtable("US_12Con_Out.csv");
US_original_out = table2array(US_original_out(:, 2: end));
US_original_st1_4 = readtable("US_12Con_st1_4.csv"); 
US_original_st1_4 = table2array(US_original_st1_4(:, 2: end));
US_original_st2_4 = readtable("US_12Con_st2_4.csv"); 
US_original_st2_4 = table2array(US_original_st2_4(:, 2: end));

FR = readtable("FR_12Con.csv"); 
FR = table2array(FR(:, 2: end));
FR_out = readtable("FR_12Con_Out.csv");
FR_out = table2array(FR_out(:, 2: end));

IT = readtable("IT_12Con.csv"); 
IT = table2array(IT(:, 2: end));
IT_out = readtable("IT_12Con_Out.csv");
IT_out = table2array(IT_out(:, 2: end));

DE = readtable("DE_12Con.csv"); 
DE = table2array(DE(:, 2: end));
DE_out = readtable("DE_12Con_Out.csv");
DE_out = table2array(DE_out(:, 2: end));

JP = readtable("JP_12Con.csv"); 
JP = table2array(JP(:, 2: end));
JP_out = readtable("JP_12Con_Out.csv");
JP_out = table2array(JP_out(:, 2: end));

AU = readtable("AU_12Con.csv"); 
AU = table2array(AU(:, 2: end));
AU_out = readtable("AU_12Con_Out.csv");
AU_out = table2array(AU_out(:, 2: end));

EU = readtable("EU_12Con.csv"); 
EU = table2array(EU(:, 2: end));
EU_out = readtable("EU_12Con_Out.csv");
EU_out = table2array(EU_out(:, 2: end));

US = readtable("US_3Factors.csv");
US_st1_1 = readtable("US_3Factors_st1_1.csv"); % stress testing 1, case 1
US_st1_2 = readtable("US_3Factors_st1_2.csv"); % stress testing 1, case 2
US_st1_3 = readtable("US_3Factors_st1_3.csv"); % stress testing 1, case 3
US_st1_4 = readtable("US_3Factors_st1_4.csv"); % stress testing 1, case 4
US_st2_1 = readtable("US_3Factors_st2_1.csv"); % stress testing 2, case 1
US_st2_2 = readtable("US_3Factors_st2_2.csv"); % stress testing 2, case 2
US_st2_3 = readtable("US_3Factors_st2_3.csv"); % stress testing 2, case 3
US_st2_4 = readtable("US_3Factors_st2_4.csv"); % stress testing 2, case 4
US = table2array(US);
US_st1_1 = table2array(US_st1_1);
US_st1_2 = table2array(US_st1_2); 
US_st1_3 = table2array(US_st1_3); 
US_st1_4 = table2array(US_st1_4); 
US_st2_1 = table2array(US_st2_1); 
US_st2_2 = table2array(US_st2_2); 
US_st2_3 = table2array(US_st2_3); 
US_st2_4 = table2array(US_st2_4); 

US_factor_fore = readtable("US_3Factors_Reconstructed.csv");
US_factor_fore = table2array(US_factor_fore);
US_factor_fore_st1_4 = readtable("US_3Factors_Reconstructed_st1_4.csv");
US_factor_fore_st1_4 = table2array(US_factor_fore_st1_4);
US_factor_fore_st2_4 = readtable("US_3Factors_Reconstructed_st2_4.csv");
US_factor_fore_st2_4 = table2array(US_factor_fore_st2_4);

if size(UK, 2) == 3
    maturity = [12, 120, 360]';
elseif size(UK, 2) == 12
    maturity = [1, 3, 6, 9, 12, 24, 36, 60, 84, 120, 240, 360]';
end

if_plot = true;

st1_start = datetime("01/01/2015", 'InputFormat', 'dd/MM/yyyy'); % start date of stress test 1
st1_end = datetime("31/12/2015", 'InputFormat', 'dd/MM/yyyy'); % end date of stress test 1
st2_start = datetime("01/01/2015", 'InputFormat', 'dd/MM/yyyy'); % start date of stress test 2

if if_plot
    for i = 1: size(US, 2)
        figure;
        hold on;
        plot(date, US(:, i), "k"); 
        plot(date, US_st1_1(:, i), "r");
        plot(date, US_st1_2(:, i), "b");
        plot(date, US_st1_3(:, i), "g");
        plot(date, US_st1_4(:, i), "y");
        yLimits = ylim;
        xline(st1_start, "k", "LineWidth", 1);
        text(st1_start, yLimits(1), "Shock starts", ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left');
        xline(st1_end, "k", "LineWidth", 1);
        text(st1_end, yLimits(1), "Shock ends", ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left');
        legend(["Factor from original data", ...
            "Factor for stress testing 1, case 1", ...
            "Factor for stress testing 1, case 2", ...
            "Factor for stress testing 1, case 3", ...
            "Factor for stress testing 1, case 4"]); 

        figure;
        hold on;
        plot(date, US(:, i), "k"); 
        plot(date, US_st2_1(:, i), "r");
        plot(date, US_st2_2(:, i), "b");
        plot(date, US_st2_3(:, i), "g");
        plot(date, US_st2_4(:, i), "y");
        yLimits = ylim;
        xline(st2_start, "k", "LineWidth", 1);
        text(st2_start, yLimits(1), "Shock starts", ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left');
        legend(["Factor from original data", ...
            "Factor for stress testing 2, case 1", ...
            "Factor for stress testing 2, case 2", ...
            "Factor for stress testing 2, case 3", ...
            "Factor for stress testing 2, case 4"]); 
    end   
end

%% DNS model
yield = UK; % in-sample yields 
yield_out = UK_out; % out-of-sample yields
n_period = 12; % number of forecasting period
n_contract = size(yield, 2); % number of contracts
n_factor = 3; % number of NS factors
n_obs = size(yield, 1); % number of observations 
mdl1 = ssm(@(params) DNS(params, yield, maturity));

A0 = [3; 3; 3];   
B0 = [0.1; 0.1; 0.1];   
mu0 = [1; 1; 1];
%D0 = repelem(0.1, n_contract)'; % Diagonal covariance matrix for measurement equation
D0 = [repelem(0.1, n_contract)'; 0.1]; % Band covariance matrix
%D0 = [0.1; 0.1]; % Full covariance matrix 
param0 = [A0; B0; mu0; D0]; % initial parameters 

options = optimoptions('fminunc','MaxFunEvals',25000,'algorithm','quasi-newton', ...
    'TolFun',1e-8,'TolX',1e-8,'MaxIter',1000,'Display','off');

[est_mdl1, params1] = estimate(mdl1, yield, param0, 'Display', 'off', ...
    'options', options, 'Univariate', true);

mu = params1(2*n_factor+1: 3*n_factor)';
[~, ~, ~, ~, ~, ~, ~, deflated_yields1] = DNS(params1, yield, maturity);
deflated_states1 = filter(est_mdl1, deflated_yields1);
estimated_states1 = deflated_states1 + mu;
estimated_yields1 = estimated_states1 * est_mdl1.C';
residual1 = yield - estimated_yields1;

rmse_DNS = sqrt( mean(residual1 .^ 2) ); 

% h-step ahead forecasting
intercept = est_mdl1.C * mu';
[fore_deflated_yield1, YMSE] = forecast(est_mdl1, n_period, deflated_yields1);
fore1 = fore_deflated_yield1 + intercept'; % predicted yields 
rmse_DNS_fore = sqrt( mean( (yield_out - fore1).^2 ) );

% Save US prediction
%writematrix(fore1, "US_DNS_Prediction_st2_4.csv");

% Plots
if if_plot
    % % Estimated Nelson-Siegel factors
    figure;
    hold on;
    plot(date, estimated_states1(:, 1), "k"); 
    plot(date, estimated_states1(:, 2), "r");
    plot(date, estimated_states1(:, 3), "b");
    legend(["Estimated 1st NS factor", "Estimated 2nd NS factor", "Estimated 3rd NS factor"]);

    % Prediction 
    figure;
    hold on;
    plot(date(n_obs-12: end), yield(n_obs-12:end, 1), "k");
    plot(date_out, fore1(:, 1), "b");
    plot(date_out, fore1(:, 1) + 1.96*sqrt(YMSE(:, 1)), "r--"); 
    plot(date_out, fore1(:, 1) - 1.96*sqrt(YMSE(:, 1)), "r--");
    legend(["In-sample yields", "Predicted yields", "Upper 95% confidence bound", "Lower 95% confidence bound"]);
end

%% DNS_FR model
yield = EU; % in-sample yields
factor = US_st2_4; 
yield_out = UK_out; % out-of-sample yields
factor_out = US_factor_fore(121: end, :); % out-of-sample US factors

n_factor = 3; % number of NS factors
n_contract = size(yield, 2); % number of contracts
n_exp = size(factor, 2); % number of explanatory variables 
n_obs = size(yield, 1); % number of observations 

A0 = [3; 3; 3];    
B0 = [0.1; 0.1; 0.1];   
mu0 = [1; 1; 1];
coe0 = repelem(1, n_contract*n_exp)'; 
%D0 = repelem(0.1, n_contract)'; % Diagonal covariance matrix for measurement equation
D0 = [repelem(0.1, n_contract)'; 0.1]; % Band covariance matrix
%D0 = [0.1; 0.1]; % Full covariance matrix
param0 = [A0; B0; mu0; coe0; D0]; % initial parameters

mdl2 = ssm(@(params) DNS_FR(params, yield, factor, maturity));    

options = optimoptions('fminunc','MaxFunEvals',50000,'algorithm','quasi-newton', ...
    'TolFun',1e-8,'TolX',1e-8,'MaxIter',1000,'Display','off');

[est_mdl2, params2] = estimate(mdl2, [yield, factor], param0, 'Display', 'off', ...
    'options', options, 'Univariate', true);

mu = params2(2*n_factor+1: 3*n_factor)';
coe = zeros(n_contract, n_exp);
coe(:) = params2(3*n_factor+1: 3*n_factor+n_contract*n_exp);
[~, ~, ~, ~, ~, ~, ~, deflated_yields2] = DNS_FR(params2, yield, factor, maturity);
[deflated_states2, ~, output2] = filter(est_mdl2, deflated_yields2);
estimated_states2 = deflated_states2 + mu;
estimated_yields2 = estimated_states2 * est_mdl2.C' + factor*coe';
residual2 = yield - estimated_yields2;

rmse2 = sqrt( mean(residual2 .^ 2) ); 

% Covariance of filtered observations
cov_state = cat(3, output2.FilteredStatesCov);
cov_Y = zeros(n_contract, n_contract, n_obs);
for i = 1: n_obs
    cov_Y(:, :, i) = est_mdl2.C * cov_state(:, :, i) * est_mdl2.C' + est_mdl2.D * est_mdl2.D';
end

% Samples from state variables 
rng(1234);
n_sample = 1000;
sample_yields = zeros(n_obs, n_contract, n_sample);
for i = 61: n_obs
    %try
    %    if det(cov_state(:, :, i)) < 0
    %        if det(round(cov_state(:, :, i), 8)) < 0
    %            sample_state = mvnrnd(estimated_states2(i, :), round(cov_state(:, :, i), 4), n_sample);
    %        else
    %            sample_state = mvnrnd(estimated_states2(i, :), round(cov_state(:, :, i), 8), n_sample);
    %        end
    %    else
    %        sample_state = mvnrnd(estimated_states2(i, :), cov_state(:, :, i), n_sample);
    %    end
    %catch
    %    cov_state(:, :, i) = cov_state(:, :, i-1);
    %end
    if det(cov_state(:, :, i)) < 0 
        [e_vector, e_value] = eig(cov_state(:, :, i));
        e_value(e_value<0) = 1e-08;
        cov_state(:, :, i) = e_vector * e_value * e_vector';
    end
    sample_state = mvnrnd(estimated_states2(i, :), cov_state(:, :, i), n_sample);
    sample_yields(i, :, :) = (sample_state * est_mdl2.C' + factor(i, :)*coe')';
end

% h-step ahead forecasting
n_period = 12; % number of forecasting period
intercept = est_mdl2.C * mu';
[fore_deflated_yield2, YMSE2] = forecast(est_mdl2, n_period, deflated_yields2);
fore2 = fore_deflated_yield2 + intercept' + factor_out*coe'; % predicted yields 
resid2 = (yield_out - fore2).^2;
rmse2_fore = sqrt( mean( (yield_out - fore2).^2 ) );

% Save data
%writematrix(coe, "data_r/Coe/coe_UK.csv");

% Spread under different stress
%original = estimated_yields2; % if original US factors
%cov_Y_original = cov_Y;
%sample_original = sample_yields;

stress24 = estimated_yields2; % if US factors under stress case 1_3. Change the number "13" to the correct stress testing case. 
cov_Y_24 = cov_Y;
sample_stress_24 = sample_yields;

if if_plot
    % Estimated Nelson-Siegel factors
    figure;
    hold on;
    plot(date, estimated_states2(:, 1), "k"); 
    plot(date, estimated_states2(:, 2), "r");
    plot(date, estimated_states2(:, 3), "b");
    xlabel("Date");
    ylabel("Value");
    legend(["Estimated 1st NS factor", "Estimated 2nd NS factor", "Estimated 3rd NS factor"]);
end

%% Forecasting plot
if if_plot
    figure;
    hold on;
    plot(date, yield(:, 1), "k");
    plot(date_out, fore1(:, 1), "b");
    plot(date_out, fore1(:, 1) + 1.96*sqrt(YMSE(:, 1)), "b--"); 
    plot(date_out, fore1(:, 1) - 1.96*sqrt(YMSE(:, 1)), "b--");
    plot(date_out, fore2(:, 1), "r");
    plot(date_out, fore2(:, 1) + 1.96*sqrt(YMSE2(:, 1)), "r--"); 
    plot(date_out, fore2(:, 1) - 1.96*sqrt(YMSE2(:, 1)), "r--");
    legend(["In-sample yields", "Predicted yields - DNS model", ...
        "Upper 95% confidence bound - DNS model", "Lower 95% confidence bound - DNS model", ...
        "Predicted yields - DNS-FR model", ...
        "Upper 95% confidence bound - DNS-FR model", "Lower 95% confidence bound - DNS-FR model"]);

end

%% Spread under different stress
shock = 2;
stress = stress24;
sample_stress = sample_stress_24;
cov_Y_st = cov_Y_24;
st = "st2_4";
country = "EU";

diff = original - stress;
diff(1: 60, :) = 0;
mean_diff = zeros(size(diff, 1), 3); 
mean_diff(:, 1) = mean(diff(:, 1: 8), 2); 
mean_diff(:, 2) = mean(diff(:, 9: 10), 2);
mean_diff(:, 3) = mean(diff(:, 11: 12), 2); 

% Analytical CI
%var_original = zeros(n_obs, 3);
%var_original(:, 1) = repelem(1/8, 8) * cov_Y_original(1:8, 1:8) * repelem(1/8, 8)';
%var_original(:, 2) = repelem(1/2, 2) * cov_Y_original(9:10, 9:10) * repelem(1/2, 2)';
%var_original(:, 3) = repelem(1/2, 2) * cov_Y_original(11:12, 11:12) * repelem(1/2, 2)';
%var_original(1:60, :) = 0;

%var_st = zeros(n_obs, 3);
%var_st(:, 1) = repelem(1/8, 8) * cov_Y_st(1:8, 1:8) * repelem(1/8, 8)';
%var_st(:, 2) = repelem(1/2, 2) * cov_Y_st(9:10, 9:10) * repelem(1/2, 2)';
%var_st(:, 3) = repelem(1/2, 2) * cov_Y_st(11:12, 11:12) * repelem(1/2, 2)';
%var_st(1:60, :) = 0;

%var_diff = var_original + var_st;

% Numerical CI
%rng(1234);
%n_sample = 1000;
%sample_diff = zeros(n_obs, n_contract, n_sample);
%sample_mean_diff = zeros(n_obs, 3, n_sample); 

%for i = 61: n_obs
%    sample_original = mvnrnd(original(i, :), cov_Y_original(:, :, i), n_sample);
%    sample_st = mvnrnd(stress(i, :), cov_Y_st(:, :, i), n_sample);
%    sample_diff(i, :, :) = sample_original' - sample_st';
%    sample_mean_diff(i, 1, :) = sort(mean(sample_diff(i, 1: 8, :), 2)); 
%    sample_mean_diff(i, 2, :) = sort(mean(sample_diff(i, 9: 10, :), 2));
%    sample_mean_diff(i, 3, :) = sort(mean(sample_diff(i, 11: 12, :), 2)); 
%end

% Sampling from state variables
sample_diff2 = sample_original - sample_stress;
sample_mean_diff2 = zeros(n_obs, 3, n_sample);
sample_mean_diff2(:, 1, :) = mean(sample_diff2(:, 1: 8, :), 2);
sample_mean_diff2(:, 2, :) = mean(sample_diff2(:, 9: 10, :), 2);
sample_mean_diff2(:, 3, :) = mean(sample_diff2(:, 11: 12, :), 2);
sample_mean_diff2 = sort(sample_mean_diff2, 3); 

%for i = 1: n_obs
%    sample_mean_diff2(i, 1, :) = sort(sample_mean_diff2(i, 1, :));
%    sample_mean_diff2(i, 2, :) = sort(sample_mean_diff2(i, 2, :));
%    sample_mean_diff2(i, 3, :) = sort(sample_mean_diff2(i, 3, :));
%end

%
if if_plot
    % Mean difference of short-end, middle, and long-end maturity
    colors = [228, 26, 28;
              55, 126, 184;
              77, 175, 74];
    figure; 
    set(gcf, 'Position', [100, 100, 900, 600]);
    ax = gca;
    ax.FontSize = 22; 
    hold on;
    for i = 1: 3
        plot(date, mean_diff(:, i), "Color", colors(i, :)/255, "LineWidth", 1);
    end
    for i = 1: 3
        %plot(date, mean_diff(:, i) + 1.96*sqrt(var_diff(:, i)), "--", "Color", colors(i, :)/255, "LineWidth", 0.5);
        %plot(date, mean_diff(:, i) - 1.96*sqrt(var_diff(:, i)), "--", "Color", colors(i, :)/255, "LineWidth", 0.5);
        %plot(date, sample_mean_diff(:, i, 25), "--", "Color", colors(i, :)/255, "LineWidth", 0.5);
        %plot(date, sample_mean_diff(:, i, 975), "--", "Color", colors(i, :)/255, "LineWidth", 0.5);
        plot(date, sample_mean_diff2(:, i, 25), "--", "Color", colors(i, :)/255, "LineWidth", 0.5);
        plot(date, sample_mean_diff2(:, i, 975), "--", "Color", colors(i, :)/255, "LineWidth", 0.5);
    end
    yline(0, "LineWidth", 0.5);
    %ylim([-0.1, 0.3]);
    yLimits = ylim;
    if shock == 1
        xline(st1_start, "k", "LineWidth", 1.5);
        text(st1_start, yLimits(1), "Shock starts", ...
            'FontSize', 22, ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left');
        xline(st1_end, "k", "LineWidth", 1.5);
        text(st1_end, yLimits(1)+0.03, "Shock ends", ...
            'FontSize', 22, ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left');
    elseif shock == 2
        xline(st2_start, "k", "LineWidth", 1.5);
        text(st2_start, yLimits(1), "Shock starts", ...
            'FontSize', 22, ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left');
    end
    xlabel("Date"); 
    ylabel("Mean difference in yields");
    legend(["Short-end maturity", "Middle maturity", "Long-end maturity"], "Location", "northwest");
    
    saveas(gcf, append("Mean_diff_in_", country, "_", st, ".jpg"));

    % Difference of all curves
    colors = [166,206,227; 
              31,120,180;
              178,223,138;
              51,160,44;
              251,154,153;
              227,26,28;
              253,191,111;
              255,127,0;
              202,178,214;
              106,61,154;
              255,255,153;
              177,89,40];
    figure;
    set(gcf, 'Position', [100, 100, 900, 600]);
    hold on;
    for i = 1: size(diff, 2)
        plot(date, diff(:, i), "Color", colors(i, :)/255, "LineWidth", 1);
    end
    if shock == 1
        xline(st1_start, "k", "LineWidth", 1.5);
        text(st1_start, yLimits(1), "Shock starts", ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left');
        xline(st1_end, "k", "LineWidth", 1.5);
        text(st1_end, yLimits(1), "Shock ends", ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left');
    elseif shock == 2
        xline(st2_start, "k", "LineWidth", 1.5);
        text(st2_start, yLimits(1), "Shock starts", ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left');
    end
    xlabel("Date"); 
    ylabel("Difference in yields");
    lgd = legend(["1 month", "3 months","6 months", "9 months", "1 year", "2 years", "3 years", "5 years", ...
        "7 years", "10 years", "20 years", "30 years"], "Location", "northwest");
    title(lgd, "Time to maturity"); 
    
    % Box plot
    if shock == 1
        figure;
        boxplot(diff(61:72, :), "Labels", maturity); 
        xlabel("Time to maturity");
        ylabel("Difference in yields");

        %saveas(gcf, append("Diff_boxplot_Shock_in_", country, "_", st, ".jpg"));

        figure;
        boxplot(diff(73:end, :), "Labels", maturity); 
        xlabel("Time to maturity");
        ylabel("Difference in yields");

        %saveas(gcf, append("Diff_boxplot_AfterShock_in_", country, "_", st, ".jpg"));

    elseif shock == 2
        figure;
        boxplot(diff(61:end, :), "Labels", maturity); 
        xlabel("Time to maturity");
        ylabel("Difference in yields");

        %saveas(gcf, append("Diff_boxplot_Shock_in_", country, "_", st, ".jpg"));
    end
end

%% Moving window - DNS forecasting
country = [EU; EU_out];

n_period = 12; % number of forecasting periods
n_factor = 3; % number of NS factors

width = 60; % width of the moving window 
n_window = size(country, 1) - n_period - width + 1; % number of windows in total
rmse_in_DNS = zeros(n_window, size(country, 2)); % in-sample rmse
rmse_out_DNS = zeros(n_window, size(country, 2)); % out-of-sample rmse

for i = 1: n_window 
    i
    yield = country(i: i+width-1, :);
    yield_out = country(i + width: i+width+n_period-1, :); % out-of-sample price
    n_obs = size(yield, 1); 
    n_contract = size(yield, 2);

    mdl_mw = ssm(@(params) DNS(params, yield, maturity));

    A0 = [3; 3; 3];   
    B0 = [0.1; 0.1; 0.1];   
    mu0 = [1; 1; 1];
    %D0 = repelem(0.1, n_contract)'; % Diagonal covariance matrix for measurement equation
    D0 = [repelem(0.1, n_contract)'; 0.1]; % Band covariance matrix
    %D0 = [0.1; 0.1]; % Full covariance matrix
    param0 = [A0; B0; mu0; D0];

    options = optimoptions('fminunc','MaxFunEvals',25000,'algorithm','quasi-newton', ...
        'TolFun',1e-8,'TolX',1e-8,'MaxIter',1000,'Display','off');

    [est_mdl_mw, params_mw] = estimate(mdl_mw, yield, param0, 'Display', 'off', ...
        'options', options, 'Univariate', true);

    mu = params_mw(2*n_factor+1: 3*n_factor)';
    [~, ~, ~, ~, ~, ~, ~, deflated_yields] = DNS(params_mw, yield, maturity);
    deflated_states = filter(est_mdl_mw, deflated_yields);
    estimated_states = deflated_states + mu;
    estimated_yields = estimated_states * est_mdl_mw.C';
    residual = yield - estimated_yields;

    rmse_in_DNS(i, :) = sqrt( mean(residual .^ 2) );

    % Forecasting
    intercept = est_mdl_mw.C * mu';
    [fore_deflated_yield, YMSE] = forecast(est_mdl_mw, n_period, deflated_yields);
    fore = fore_deflated_yield + intercept';
    rmse_out_DNS(i, :) = sqrt( mean( (yield_out - fore).^2 ) );
    
    % Save data
    %writematrix(fore, append("data_r/mw/", string(i), ".csv"));
end

% Plots
if if_plot
    date_all = [date; date_out];
    
    figure;
    plot(date_all(width/2: width/2+n_window-1), rmse_in_DNS);
    xlabel("Date");
    ylabel("In-sample RMSE");

    figure;
    plot(date_all(width/2: width/2+n_window-1), mean(rmse_in_DNS, 2));
    xlabel("Date");
    ylabel("In-sample mean RMSE");

    figure;
    plot(date_all(width+n_period/2: width+n_period/2+n_window-1), rmse_out_DNS);
    xlabel("Date");
    ylabel("Out-of-sample RMSE");

    figure;
    plot(date_all(width+n_period/2: width+n_period/2+n_window-1), mean(rmse_out_DNS, 2));
    xlabel("Date");
    ylabel("Out-of-sample mean RMSE");
end

%% Moving window - DNS-FR forecasting
country = [EU; EU_out];

n_period = 12; % number of forecasting periods
n_factor = 3; % number of NS factors

width = 60; % width of the moving window 
n_window = size(country, 1) - n_period - width + 1; % number of windows in total
rmse_in_DNSFR = zeros(n_window, size(country, 2)); % in-sample rmse
rmse_out_DNSFR = zeros(n_window, size(country, 2)); % out-of-sample rmse

for i = 1: n_window
    i
    yield = country(i: i+width-1, :);
    yield_out = country(i + width: i+width+n_period-1, :); % out-of-sample price
    factor_all = readtable(append("mw_U/Recon", string(i), ".csv"));
    factor_all = table2array(factor_all);
    factor = factor_all(1: width, :);
    factor_out = factor_all(width+1: end, :); % out-of-sample US factors
    
    n_obs = size(yield, 1);
    n_contract = size(yield, 2);
    n_exp = size(factor, 2); % number of explanatory variables   
    
    A0 = [3; 3; 3];   
    B0 = [0.1; 0.1; 0.1];   
    mu0 = [1; 1; 1];
    coe0 = repelem(1, n_contract*n_exp)'; 
    %D0 = repelem(0.1, n_contract)'; % Diagonal covariance matrix for measurement equation
    D0 = [repelem(0.1, n_contract)'; 0.1]; % Band covariance matrix
    %D0 = [0.1; 0.1]; % Full covariance matrix
    param0 = [A0; B0; mu0; coe0; D0];
    
    mdl_mw = ssm(@(params) DNS_FR(params, yield, factor, maturity));    
    
    options = optimoptions('fminunc','MaxFunEvals',50000,'algorithm','quasi-newton', ...
        'TolFun',1e-8,'TolX',1e-8,'MaxIter',1000,'Display','off');
    
    [est_mdl_mw, params_mw] = estimate(mdl_mw, [yield, factor], param0, 'Display', 'off', ...
        'options', options, 'Univariate', true);
    
    mu = params_mw(2*n_factor+1: 3*n_factor)';
    [~, ~, ~, ~, ~, ~, ~, deflated_yields] = DNS_FR(params_mw, yield, factor, maturity);
    deflated_states = filter(est_mdl_mw, deflated_yields);
    estimated_states = deflated_states + mu;
    residual = deflated_yields - deflated_states * est_mdl_mw.C';
    
    rmse_in_DNSFR(i, :) = sqrt( mean(residual .^ 2) ); 
    
    coe = zeros(n_contract, n_exp);
    coe(:) = params_mw(3*n_factor+1: 3*n_factor+n_contract*n_exp);
    
    % Forecasting
    intercept = est_mdl_mw.C * mu';
    [fore_deflated_yield, YMSE] = forecast(est_mdl_mw, n_period, deflated_yields);
    fore = fore_deflated_yield + intercept' + factor_out*coe';

    rmse_out_DNSFR(i, :) = sqrt( mean( (yield_out - fore).^2 ) );
end

% Plots
if if_plot
    date_all = [date; date_out];
    
    figure;
    plot(date_all(width/2: width/2+n_window-1), rmse_in_DNSFR);
    xlabel("Date");
    ylabel("In-sample RMSE");

    figure;
    plot(date_all(width/2: width/2+n_window-1), mean(rmse_in_DNSFR, 2));
    xlabel("Date");
    ylabel("In-sample mean RMSE");

    figure;
    plot(date_all(width+n_period/2: width+n_period/2+n_window-1), rmse_out_DNSFR);
    xlabel("Date");
    ylabel("Out-of-sample RMSE");

    figure;
    plot(date_all(width+n_period/2: width+n_period/2+n_window-1), mean(rmse_out_DNSFR, 2));
    xlabel("Date");
    ylabel("Out-of-sample mean RMSE");

    colors = [5,113,176;
              202,0,32];

    figure;
    hold on;
    plot(date_all(width/2: width/2+n_window-1), mean(rmse_in_DNS, 2), "Color",colors(1,:)/255, "LineWidth", 1);
    plot(date_all(width/2: width/2+n_window-1), mean(rmse_in_DNSFR, 2), "Color",colors(2,:)/255, "LineWidth", 1);
    xlabel("Date");
    ylabel("In-sample mean RMSE");
    legend(["DNS model", "DNS-FR model"]);

    figure;
    hold on;
    plot(date_all(width+n_period/2: width+n_period/2+n_window-1), mean(rmse_out_DNS, 2), "Color",colors(1,:)/255, "LineWidth", 1);
    plot(date_all(width+n_period/2: width+n_period/2+n_window-1), mean(rmse_out_DNSFR, 2), "Color",colors(2,:)/255, "LineWidth", 1);
    xlabel("Date");
    ylabel("Out-of-sample mean RMSE");
    legend(["DNS model", "DNS-FR model"]);
end

%% Bond ladder
ex_rate = [1.28231, 1.24086, 1.25907, 1.23455, 1.23992, 1.3077,...
    1.3369, 1.29115, 1.29457, 1.33173, 1.36561, 1.36893]'; % exchange rate
ex_rate_initial = 1.32018; 
%ex_rate = repelem(1, 12); 
%ex_rate_initial = 1;
EFFR = [1.59, 1.59, 0.08, 0.05, 0.05, 0.08, 0.1, 0.09, 0.09, 0.09, 0.09, 0.09]'/100;

%EFFR = readtable("EFFR.xlsx");
%date_EFFR = table2array(EFFR(:, 1));
%date_EFFR = datetime(date_EFFR, 'InputFormat', 'MM/dd/yyyy');
%EFFR = table2array(EFFR(:, 3));
%EFFR = EFFR(ismember(date_EFFR, date_out));

mat_bl = 360; % matuirty of bond ladder portfolio
new_maturity = (1: mat_bl)';

yield = UK;
factor = US_st2_4;
yield_out = UK_out; % out-of-sample price
factor_out = US_factor_fore_st2_4(121: end, :); % out-of-sample US factors

n_factor = 3; % number of NS factors
n_contract = size(yield, 2);
n_exp = size(factor, 2); % number of explanatory variables 
n_obs = size(yield, 1);

A0 = [3; 3; 3];    
B0 = [0.1; 0.1; 0.1];   
mu0 = [1; 1; 1];
coe0 = repelem(1, n_contract*n_exp)'; 
%D0 = repelem(0.1, n_contract)';
D0 = [repelem(0.1, n_contract)'; 0.1];
%D0 = [0.1; 0.1];
param0 = [A0; B0; mu0; coe0; D0];

mdl_bl = ssm(@(params) DNS_FR(params, yield, factor, maturity));    

options = optimoptions('fminunc','MaxFunEvals',50000,'algorithm','quasi-newton', ...
    'TolFun',1e-8,'TolX',1e-8,'MaxIter',1000,'Display','off');

[est_mdl_bl, params_bl] = estimate(mdl_bl, [yield, factor], param0, 'Display', 'off', ...
    'options', options, 'Univariate', true);

mu = params_bl(2*n_factor+1: 3*n_factor)';
[~, ~, ~, ~, ~, ~, ~, deflated_yields_bl] = DNS_FR(params_bl, yield, factor, maturity);
[deflated_states_bl, ~, output_bl] = filter(est_mdl_bl, deflated_yields_bl);
estimated_states_bl = deflated_states_bl + mu;

lambda = 0.0609;
C_new = [ones(size(new_maturity)), (1 - exp(-lambda*new_maturity))./(lambda*new_maturity), ...
    ( (1 - exp(-lambda*new_maturity))./(lambda*new_maturity) - exp(-lambda*new_maturity) ) ];
deflated_yields_est = deflated_states_bl * C_new';

coe = zeros(n_contract, n_exp);
coe(:) = params_bl(3*n_factor+1: 3*n_factor+n_contract*n_exp);

% h-step ahead forecasting
n_period = 12; % number of forecasting period
intercept = est_mdl_bl.C * mu';
[forecasted_deflated_yield_bl, YMSE_bl] = forecast(est_mdl_bl, n_period, deflated_yields_bl);
forecasted_yield_bl = forecasted_deflated_yield_bl + intercept' + factor_out*coe'; % predicted yields

forecasted_yield_interp = zeros(n_period, mat_bl); % all predicted yields, including interpolated yields
forecasted_yield_interp(:, maturity(maturity <= mat_bl)) = forecasted_yield_bl(:, maturity <= mat_bl);
interp_maturity = setdiff(new_maturity, maturity); % maturity of interpolated data
linear_interp = interp1(maturity(maturity <= mat_bl), forecasted_yield_bl(:, maturity <= mat_bl)', interp_maturity); % linear interpolation
forecasted_yield_interp(:, interp_maturity) = linear_interp'; 

%n_period = 12;
%Phi = est_mdl_bl.A;
%Lambda = C_new;
%a_N = deflated_states_bl(end, :);

%forecasted_deflated_yield_bl = zeros(n_period, mat_bl);

%for i = 1: n_period
%    forecasted_deflated_yield_bl(i, :) = Lambda * Phi^i * a_N';
%end

%coe_bl = repmat(coe(maturity == mat_bl, :), mat_bl, 1); 
%coe_bl = coe;
%forecasted_yield_bl = forecasted_deflated_yield_bl + mu*Lambda' + factor_out*coe_bl';

% Portfolio price
w0 = 12000000; % initial wealth
%investment = w0/mat_bl; % investment amount per period
investment = w0 / 12;
n_inv = 13; % number of investment
fv = 100; % face value

ttm = zeros(n_period+1, n_inv); % time-to-maturity matrix, row -> dates, column -> investments
ttm(1, 1) = mat_bl; 

n_bond = zeros(n_period+1, n_inv); % number of bonds, row -> dates, column -> investments
n_bond(1, 1) = investment / ( fv*ex_rate_initial * exp(-mat_bl/12 * yield(end, maturity' == mat_bl)/100) ); % first investment, yield is known

cash = zeros(n_period+1, 1); % cash account
cash(1) = w0 - investment; % time 0

bond = zeros(n_period+1, n_inv); % bond account, row -> dates, column -> investments
bond(1, 1) = investment; % time 0

rng(1234);

for i = 1: n_period
    % Update time-to-maturity
    ttm_new = ttm(i, :) - 1;
    ttm_new(i+1) = mat_bl;
    ttm_new(ttm_new < 0) = 0;
    ttm(i+1, :) = ttm_new;

    % Update number of bonds
    n_bond(i+1, :) = n_bond(i, :);
    n_bond(i+1, i+1) = investment / ( fv*ex_rate(i) * exp(-mat_bl/12 * forecasted_yield_interp(i, mat_bl)/100) );

    % Update cash account 
    cash(i+1) = cash(i)*(1+EFFR(i)/12) - investment; 
    
    % Update bond account
    index_pos = ttm(i+1, :) > 0; % index of bonds with positive time-to-maturity
    current_yield = zeros(1, n_inv); % yields at current time point
    current_yield(index_pos) = forecasted_yield_interp(i, ttm(i+1, index_pos))/100; 
    bond(i+1, :) = n_bond(i+1, :) .* fv*ex_rate(i) .* exp(-ttm(i+1, :)/12 .* current_yield); 

    index_expiry = find(ttm(i, :) == 1); % index of expired bond
    if index_expiry
        cash(i+1) = cash(i+1) + n_bond(i+1, index_expiry) * fv*ex_rate(i); % add expired bonds into cash account
        bond(i+1, index_expiry) = 0; % delete bond
        n_bond(i+1, index_expiry) = 0; % delete bond 
    end
end

value = sum(bond,2) + cash; % portfolio value

if if_plot
    figure; 
    plot(0: 12, value); 
    xlabel("Time (in months)");
    ylabel("Portfolio value");
end

% Confidence interval
n_sample = 1000; % number of random samples 

state_cov_N = output_bl(120).FilteredStatesCov;
state_cov_fore = zeros(3, 3, n_period); % covariance matrix of the predicted state variables 
yield_cov_fore = zeros(12, 12, n_period); % covariance matrix of the predicted yields 
yield_var_fore = zeros(12, 12); % variance vector of the predicted yields 
state_cov_fore(:, :, 1) = est_mdl_bl.A * state_cov_N * est_mdl_bl.A' + est_mdl_bl.B * est_mdl_bl.B';
yield_cov_fore(:, :, 1) = est_mdl_bl.C * state_cov_fore(:, :, 1) * est_mdl_bl.C' + est_mdl_bl.D * est_mdl_bl.D';
yield_var_fore(1, :) = diag(yield_cov_fore(:, :, 1));
for i = 2: n_period
    state_cov_fore(:, :, i) = est_mdl_bl.A * state_cov_fore(:, :, i-1) * est_mdl_bl.A' + est_mdl_bl.B * est_mdl_bl.B';
    yield_cov_fore(:, :, i) = est_mdl_bl.C * state_cov_fore(:, :, i) * est_mdl_bl.C' + est_mdl_bl.D * est_mdl_bl.D';
    yield_var_fore(i, :) = diag(yield_cov_fore(:, :, i));
end

rng(1234);

sample_yields = zeros(n_period, mat_bl, n_sample);

for i = 1: n_period
    samples = mvnrnd(forecasted_deflated_yield_bl(i, :), yield_cov_fore(:, :, i), n_sample) + intercept' + factor_out(i, :)*coe';

    samples_interp = zeros(n_sample, mat_bl); % all predicted yields, including interpolated yields
    samples_interp(:, maturity(maturity <= mat_bl)) = samples(:, maturity <= mat_bl);
    interp_maturity = setdiff(new_maturity, maturity); % maturity of interpolated data
    linear_interp = interp1(maturity(maturity <= mat_bl), samples(:, maturity <= mat_bl)', interp_maturity); % linear interpolation
    samples_interp(:, interp_maturity) = linear_interp'; 

    sample_yields(i, :, :) = samples_interp';
end

value_CI = zeros(n_sample, n_period+1);

for j = 1: n_sample
    ttm = zeros(n_period+1, n_inv); % time-to-maturity matrix, row -> dates, column -> investments
    ttm(1, 1) = mat_bl; 
    
    n_bond = zeros(n_period+1, n_inv); % number of bonds, row -> dates, column -> investments
    n_bond(1, 1) = investment / ( fv*ex_rate_initial * exp(-mat_bl/12 * yield(end, maturity' == mat_bl)/100) ); % first investment, yield is known
    
    cash = zeros(n_period+1, 1); % cash account
    cash(1) = w0 - investment; % time 0
    
    bond = zeros(n_period+1, n_inv); % bond account, row -> dates, column -> investments
    bond(1, 1) = investment; % time 0
    
    for i = 1: n_period
        % Update time-to-maturity
        ttm_new = ttm(i, :) - 1;
        ttm_new(i+1) = mat_bl;
        ttm_new(ttm_new < 0) = 0;
        ttm(i+1, :) = ttm_new;
    
        % Update number of bonds
        n_bond(i+1, :) = n_bond(i, :);
        n_bond(i+1, i+1) = investment / ( fv*ex_rate(i) * exp(-mat_bl/12 * sample_yields(i, mat_bl, j)/100) );
    
        % Update cash account 
        cash(i+1) = cash(i)*(1+EFFR(i)/12) - investment; 
        
        % Update bond account
        index_pos = ttm(i+1, :) > 0; % index of bonds with positive time-to-maturity
        current_yield = zeros(1, n_inv); % yields at current time point
        current_yield(index_pos) = sample_yields(i, ttm(i+1, index_pos), j)/100; 
        bond(i+1, :) = n_bond(i+1, :) .* fv*ex_rate(i) .* exp(-ttm(i+1, :)/12 .* current_yield); 
    
        index_expiry = find(ttm(i, :) == 1); % index of expired bond
        if index_expiry
            cash(i+1) = cash(i+1) + n_bond(i+1, index_expiry) * fv*ex_rate(i); % add expired bonds into cash account
            bond(i+1, index_expiry) = 0; % delete bond
            n_bond(i+1, index_expiry) = 0; % delete bond 
        end
    end
    
    value_CI(j, :) = sum(bond,2) + cash; % portfolio value
end

if if_plot
    value_CI_sorted = sort(value_CI); 

    figure;
    hold on;
    plot(0: 12, value/1e06, "b"); 
    plot(0:12,value_CI_sorted(25, :)/1e06, "b--");
    plot(0:12,value_CI_sorted(975, :)/1e06, "b--");
    xlabel("Time (in months)");
    ylabel("Portfolio value (in million dollars)");
    legend(["Estimated portfolio values", "95% CI", "95% CI"])
end

%% Different maturities
if if_plot
    colors = [27,158,119;
              217,95,2;
              117,112,179];

    figure;
    hold on;
    plot(0: 12, value1/1e06, "Color", colors(1, :)/255, "LineWidth", 1);
    plot(0: 12, value2/1e06, "Color", colors(2, :)/255, "LineWidth", 1);
    plot(0: 12, value3/1e06, "Color", colors(3, :)/255, "LineWidth", 1);
    plot(0:12, value_CI_sorted1(25, :)/1e06, "--", "Color", colors(1, :)/255, "LineWidth", 1);
    plot(0:12, value_CI_sorted1(975, :)/1e06, "--", "Color", colors(1, :)/255, "LineWidth", 1);
    plot(0:12, value_CI_sorted2(25, :)/1e06, "--", "Color", colors(2, :)/255, "LineWidth", 1);
    plot(0:12, value_CI_sorted2(975, :)/1e06, "--", "Color", colors(2, :)/255, "LineWidth", 1);
    plot(0:12, value_CI_sorted3(25, :)/1e06, "--", "Color", colors(3, :)/255, "LineWidth", 1);
    plot(0:12, value_CI_sorted3(975, :)/1e06, "--", "Color", colors(3, :)/255, "LineWidth", 1);
    xlabel("Time (in months)");
    ylabel("Portfolio value (in million dollars)");
    legend(["6 months maturity", "1 year maturity", "30 years maturity"], "Location", "northwest");
    %legend(["Original data", "Temporary shock", "Permanent shock"], "Location", "northwest");
end

%% Different shocks 
if if_plot
    %diff_CI_sorted2 = sort(value_CI2 - value_CI1);
    %diff_CI_sorted3 = sort(value_CI3 - value_CI1);

    colors = [202,0,32;
              64,64,64];

    figure;
    ax = gca;
    ax.FontSize = 14; 
    hold on;
    plot(0:12, value2 - value1, "Color", colors(1, :)/255, "LineWidth", 1);
    plot(0:12, value3 - value1, "Color", colors(2, :)/255, "LineWidth", 1);
    %plot(0:12, diff_CI_sorted2(25, :), "--", "Color", colors(1, :)/255, "LineWidth", 1);
    %plot(0:12, diff_CI_sorted2(975, :), "--", "Color", colors(1, :)/255, "LineWidth", 1);
    %plot(0:12, diff_CI_sorted3(25, :), "--", "Color", colors(2, :)/255, "LineWidth", 1);
    %plot(0:12, diff_CI_sorted3(975, :), "--", "Color", colors(2, :)/255, "LineWidth", 1);
    xlabel("Time (in months)");
    ylabel("Difference in portfolio value (in dollars)");
    legend(["Temporary shock", "Permanent shock"]);
end

%% Value at Risk
if if_plot
    colors = [27,158,119;
              217,95,2;
              117,112,179];

    figure;
    hold on;
    plot(0:12, value_CI_sorted1(50, :)/1e06, "Color", colors(1, :)/255, "LineWidth", 1);
    plot(0:12, value_CI_sorted2(50, :)/1e06, "Color", colors(2, :)/255, "LineWidth", 1);
    plot(0:12, value_CI_sorted3(50, :)/1e06, "Color", colors(3, :)/255, "LineWidth", 1);
    xlabel("Time (in months)");
    ylabel("5% value at risk (in million dollars)");
    legend(["6 months maturity", "1 year maturity", "30 years maturity"], "Location", "northwest");
end

%% Difference of VaR
if if_plot
    colors = [202,0,32;
              64,64,64];

    figure;
    ax = gca;
    ax.FontSize = 14; 
    hold on;
    plot(0:12, value_CI_sorted2(50, :) - value_CI_sorted1(50, :), "Color", colors(1, :)/255, "LineWidth", 1);
    plot(0:12, value_CI_sorted3(50, :) - value_CI_sorted1(50, :), "Color", colors(2, :)/255, "LineWidth", 1);
    xlabel("Time (in months)");
    ylabel("Difference in 5% value at risk (in dollars)");
    legend(["Temporary shock", "Permanent shock"]);
end

%% DNS factor loadings
tau = 1: 400;
lambda = 0.0609;
beta1 = repelem(1, length(tau));
beta2 = (1-exp(-lambda*tau)) ./ (lambda*tau);
beta3 = (1-exp(-lambda*tau)) ./ (lambda*tau) - exp(-lambda*tau); 

figure;
hold on;
plot(tau, beta1, "k", "LineWidth", 1);
plot(tau, beta2, "r", "LineWidth", 1);
plot(tau, beta3, "b", "LineWidth", 1);
ylim([0, 1.1]);
xlabel("Maturity in months");
ylabel("Loadings");
legend(["Level", "Slope", "Curvature"]);

%% Data
[X1, Y1] = meshgrid([date;date_out], maturity); 
Z1 = [UK; UK_out]'; 
figure;
surf(X1, Y1, Z1);
colorbar;
xlabel("Date");
ylabel("Maturity in months");
zlabel("Yield");

index = [1,2,3,5,6,7,8,10,11,12];
[X2, Y2] = meshgrid([date;date_out], maturity(index)); 
Z2 = [UK(:, index); UK_out(:, index)]'; 
figure;
surf(X2, Y2, Z2);
colorbar;
xlabel("Date");
ylabel("Maturity in months");
zlabel("Yield");

%% Yield curve
price = UK;
index1 = find(date == datetime("2011-01-31"));
index2 = find(date == datetime("2015-08-02"));
index3 = find(date == datetime("2019-09-01"));

colors = [178,223,138;
          166,206,227;
          31,120,180];

figure;
%set(gcf, 'Position', [100, 100, 900, 600]);
ax = gca;
ax.FontSize = 14; 
hold on;
plot(maturity, price(index1, :), "Color", colors(1, :)/255, "LineWidth", 1);
plot(maturity, estimated_yields1(index1, :), "Color", colors(2, :)/255, "LineWidth", 1);
plot(maturity, estimated_yields2(index1, :), "Color", colors(3, :)/255, "LineWidth", 1);
xlabel("Maturity in months");
ylabel("Yield");
legend(["Original data", "DNS estimation", "DNS-FR estimation"]);

figure;
ax = gca;
ax.FontSize = 14; 
hold on;
plot(maturity, price(index2, :), "Color", colors(1, :)/255, "LineWidth", 1);
plot(maturity, estimated_yields1(index2, :), "Color", colors(2, :)/255, "LineWidth", 1);
plot(maturity, estimated_yields2(index2, :), "Color", colors(3, :)/255, "LineWidth", 1);
xlabel("Maturity in months");
ylabel("Yield");
legend(["Original data", "DNS estimation", "DNS-FR estimation"]);
axes("Position", [0.7, 0.2, 0.2, 0.2]);
box on;
hold on;
plot(maturity(1:6), price(index2, 1:6), "Color", colors(1, :)/255, "LineWidth", 1);
plot(maturity(1:6), estimated_yields1(index2, 1:6), "Color", colors(2, :)/255, "LineWidth", 1);
plot(maturity(1:6), estimated_yields2(index2, 1:6), "Color", colors(3, :)/255, "LineWidth", 1);

figure;
ax = gca;
ax.FontSize = 14; 
hold on;
plot(maturity, price(index3, :), "Color", colors(1, :)/255, "LineWidth", 1);
plot(maturity, estimated_yields1(index3, :), "Color", colors(2, :)/255, "LineWidth", 1);
plot(maturity, estimated_yields2(index3, :), "Color", colors(3, :)/255, "LineWidth", 1);
xlabel("Maturity in months");
ylabel("Yield");
legend(["Original data", "DNS estimation", "DNS-FR estimation"]);

%% Yield curve under different stress testings
price = UK;
index1 = find(date == datetime("2011-01-31"));
index2 = find(date == datetime("2015-08-02"));
index3 = find(date == datetime("2019-09-01"));

figure;
hold on;
plot(maturity, price(index1, :)); 
plot(maturity, original(index1, :)); 
plot(maturity, stress21(index1, :)); 
plot(maturity, stress22(index1, :)); 
plot(maturity, stress23(index1, :)); 
plot(maturity, stress24(index1, :)); 
xlabel("Maturity in months");
ylabel("Yield");
legend(["Real yields", "Estimation using original US factors", "Stress 1.1", "Stress 1.2", "Stress 1.3", "Stress 1.4"]);

figure;
hold on;
plot(maturity, price(index2, :)); 
plot(maturity, original(index2, :)); 
plot(maturity, stress21(index2, :)); 
plot(maturity, stress22(index2, :)); 
plot(maturity, stress23(index2, :)); 
plot(maturity, stress24(index2, :)); 
xlabel("Maturity in months");
ylabel("Yield");
legend(["Real yields", "Estimation using original US factors", "Stress 1.1", "Stress 1.2", "Stress 1.3", "Stress 1.4"]);

figure;
hold on;
plot(maturity, price(index3, :)); 
plot(maturity, original(index3, :)); 
plot(maturity, stress21(index3, :)); 
plot(maturity, stress22(index3, :)); 
plot(maturity, stress23(index3, :)); 
plot(maturity, stress24(index3, :)); 
xlabel("Maturity in months");
ylabel("Yield");
legend(["Real yields", "Estimation using original US factors", "Stress 1.1", "Stress 1.2", "Stress 1.3", "Stress 1.4"]);
