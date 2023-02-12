% Off-grid Microgrid Supply-Demand Model V2
% Description: Supply-demand simulation for off-grid rural community
% projects following a bottom-up approach.
% Author: David Sanchinelli
% Last update: 10/02/23

%==========================================================================
% INPUT SIMULATION DATA (MODIFIABLE / CHECK MODEL MANUAL FOR DETAILS)
%==========================================================================
% Topology scenario to evaluate
topology = 1; % 0: Stand-alone, 1: Centralized, 2: Distributed
% Community to evaluate
community_type = 1;  % 1: San Pedro Carcha, 2: San Bartolome Jocotenango, 3: San Gaspar Ixchil, 4: El Chal
% Type of data set available
dataset_type = 0; % 0: Averaged (monthly average data), 1: Historic-Averaged(hourly data - few years), 2: Historic(hourly data - all years)
data_path = 'C:\Users\david\Desktop'; % Location of data files
% Run supply optimization program or integrated analysis simulation
mix_test = false; % true = run the supply optimization program (Runtime will be set to 365 days, no integrated analysis)
R_factor = 1.5; % Reliability factor to adjust optimization results (default = 1.5)
% Technology mix to evaluate (if integrated analysis simulation is selected)
solar_n = 44; % If mix_test = false, set number of solar panels
wind_n = 9; % If mix_test = false, set number of wind turbines
batteries_n = 9; % If mix_test = false, set number of batteries
% Total time periods for integrated analysis (set time in days)
short_time = 365*2;
mid_time = 365*10;
long_time = 365*15;
% Randomizing parameters 'Rh' and 'Rw'
Rh = 0.1; % Max. percentage of 'h' to add/remove
Rw = 0.5; % Max. percentage of 1 hour to add/remove
% Initial standard deviation for normal distribution of peak appliances (peak-value correction)
sigma = 2; % Initial value
% Plot graphs for resulting scenario characteristics
plot_results = true; % Choose to plot results or not
plot_path = 'Final Results\Scenario 1\Centralized'; % Folder to store results and graphs
% Profiles to generate per day for convergence
convergence_test = true; % Set convergence test or set profiles to generate
convergence_percent = 0.25; % Average value and standard deviation of convergence test (default = 0.25)
profiles_per_day = 3; % If test is false, set specific number of profiles (min. set to 3)
%==========================================================================

%==========================================================================
% INPUT COMMUNITY DATA FROM EXCEL FILES (MODIFIABLE)
%==========================================================================
switch community_type
    case 1
        pop_data = readmatrix(data_path+"\COMMUNITY_DATA.xlsx",'Sheet','A. General Community','Range','C10:O17');
        app_data = readmatrix(data_path+"\COMMUNITY_DATA.xlsx",'Sheet','1. San Pedro Carcha','Range','A2:M37');
        if dataset_type == 0
            weather_data = readmatrix(data_path+"\COMMUNITY_DATA.xlsx",'Sheet','A. General Community','Range','C2:N9');
        else
            if mix_test == true
                weather_data = readmatrix(data_path+"\COMMUNITY_DATA_HISTORIC3Y.xlsx",'Sheet','1. San Pedro Carcha','Range','B11');
            else
                weather_data = readmatrix(data_path+"\COMMUNITY_DATA_HISTORIC4Y.xlsx",'Sheet','1. San Pedro Carcha','Range','B11');
            end
        end
    case 2
        pop_data = readmatrix(data_path+"\COMMUNITY_DATA.xlsx",'Sheet','A. General Community','Range','C27:O34');
        app_data = readmatrix(data_path+"\COMMUNITY_DATA.xlsx",'Sheet','2. San Bartolome Jocotenango','Range','A2:M37');
        if dataset_type == 0
            weather_data = readmatrix(data_path+"\COMMUNITY_DATA.xlsx",'Sheet','A. General Community','Range','C19:N26');
        else
            if mix_test == true
                weather_data = readmatrix(data_path+"\COMMUNITY_DATA_HISTORIC3Y.xlsx",'Sheet','2. San Bartolome Jocotenango','Range','B11');
            else
                weather_data = readmatrix(data_path+"\COMMUNITY_DATA_HISTORIC4Y.xlsx",'Sheet','2. San Bartolome Jocotenango','Range','B11');
            end
        end
    case 3
        pop_data = readmatrix(data_path+"\COMMUNITY_DATA.xlsx",'Sheet','A. General Community','Range','C45:O52');
        app_data = readmatrix(data_path+"\COMMUNITY_DATA.xlsx",'Sheet','3. San Gaspar Ixchil','Range','A2:M37');
        if dataset_type == 0
            weather_data = readmatrix(data_path+"\COMMUNITY_DATA.xlsx",'Sheet','A. General Community','Range','C37:N44');
        else
            if mix_test == true
                weather_data = readmatrix(data_path+"\COMMUNITY_DATA_HISTORIC3Y.xlsx",'Sheet','3. San Gaspar Ixchil','Range','B11');
            else
                weather_data = readmatrix(data_path+"\COMMUNITY_DATA_HISTORIC4Y.xlsx",'Sheet','3. San Gaspar Ixchil','Range','B11');
            end
        end
    case 4
        pop_data = readmatrix(data_path+"\COMMUNITY_DATA.xlsx",'Sheet','A. General Community','Range','C63:O70');
        app_data = readmatrix(data_path+"\COMMUNITY_DATA.xlsx",'Sheet','4. El Chal','Range','A2:M37');
        if dataset_type == 0
            weather_data = readmatrix(data_path+"\COMMUNITY_DATA.xlsx",'Sheet','A. General Community','Range','C55:N62');
        else
            if mix_test == true
                weather_data = readmatrix(data_path+"\COMMUNITY_DATA_HISTORIC3Y.xlsx",'Sheet','4. El Chal','Range','B11');
            else
                weather_data = readmatrix(data_path+"\COMMUNITY_DATA_HISTORIC4Y.xlsx",'Sheet','4. El Chal','Range','B11');
            end
        end
end
topology_data = readmatrix(data_path+"\COMMUNITY_DATA.xlsx",'Sheet','B. Supply Data','Range','B2:H4');
price_data = readmatrix(data_path+"\COMMUNITY_DATA.xlsx",'Sheet','B. Supply Data','Range','B7:F9');
%==========================================================================

%==========================================================================
% INITIALIZE SUPPLY VARIABLES
%==========================================================================
% Set time for supply optimization program (if selected)
if mix_test == true
    short_time = 50;
    mid_time = 250;
    long_time = 365;
end
% Indicators
total_time = long_time;
time_step = 60; % Calculations are done per minute
reliability_counter = 0;
reliability_hours = 0;
avg_households_affected_rel = zeros(total_time*1440/time_step,1);
prev_down_rel = 0;
resilience_counter = 0;
resilience_hours = 0;
avg_households_affected_res = zeros(total_time*1440/time_step,1);
prev_down_res = 0;
results_short = zeros(3,3);
results_mid = zeros(3,3);
results_long = zeros(3,3);
if dataset_type == 1
    dataset_size = floor(size(weather_data,1)/8760); % How many years in dataset
    if dataset_size > 1
        mean_weather = zeros(8760,5);
        total_weather = zeros(dataset_size,1);
        for i = 1:8760
            for j = 1:5
                for k = 1:dataset_size
                    total_weather(k) = weather_data(((k-1)*8760)+i,j);
                end
                mean_weather(i,j) = mean(total_weather,'all');
            end
        end
        weather_data = mean_weather;
    end
end
% Solar Variables
solar_failure_rate = zeros(solar_n,1);
solar_failure_rate = solar_failure_rate + 0.0005;
solar_dirt_r = 0; % Efficiency drop from dirty panels
solar_down = zeros(solar_n,1); % Status of each panel (up or down)
used_solar_n = solar_n; % Total panels online
if dataset_type == 0
    solar_data = weather_data(5:6,:);
    solar_curve = normpdf(1:1440/time_step,1440/time_step/2,(360/time_step)/3);
    solar_curve(1:300/time_step) = 0;
    solar_curve(1140/time_step:1440/time_step) = 0;
    total_solarscale = [];
else
    solar_data = weather_data(:,2);
    solar_data = solar_data/1000;
end
% Wind Variables
wind_failure_rate = zeros(wind_n,1);
wind_failure_rate = wind_failure_rate + 0.0054;
wind_down = zeros(wind_n,1); % Status of each turbine (up or down)
used_wind_n = wind_n; % Total turbines online
if dataset_type == 0
    wind_data = weather_data(7:8,:);
    wind_curve = gampdf(1:1440/60,5,3); % ONLY SUITABLE FOR 1H TIMESTEPS
    wind_to_power = zeros(size(wind_curve,2),1);
    total_windscale = [];
else
    wind_data = weather_data(:,4);
    wind_to_power = zeros(24,1);
end
% Total Supply/Battery Variables
total_panelsperday = zeros(total_time,1);
total_turbinesperday = zeros(total_time,1);
total_supply = zeros(1440/time_step,1);
total_supply_oandm = 0;
total_battery = zeros(1440/time_step,1);
total_balance = zeros(1440/time_step,1);
solarrate_history = zeros(total_time,solar_n);
windrate_history = zeros(total_time,wind_n);
ttotal_supply = zeros(total_time,1440/time_step);
ttotal_battery = zeros(total_time,1440/time_step);
ttotal_balance = [];
% Weather Variables
if dataset_type == 0
    rain_data = weather_data(4,:);
else
    rain_data = weather_data(:,3);
end
% Topology Variables
if topology == 0
    initial_solarP = topology_data(1,1)/1000;
    solar_panel_P = initial_solarP;
    panelsize = topology_data(1,7);
    initial_windP = topology_data(1,2)/1000;
    wind_P = initial_windP;
    wind_cutin = topology_data(1,5);
    wind_rated = topology_data(1,6);
    hydro_P(1:1440/time_step) = topology_data(1,3)/1000;
    battery_size = topology_data(1,4);
    battery_charge = battery_size * batteries_n;
    max_battery = battery_charge * 0.8; % Include round-trip efficiency
    battery_charge_h = max_battery;
    prev_down_res_h = 0;
    prev_down_rel_h = 0;
    investment_price = price_data(1,:);
elseif topology == 1
    initial_solarP = topology_data(2,1)/1000;
    solar_panel_P = initial_solarP;
    panelsize = topology_data(2,7);
    initial_windP = topology_data(2,2)/1000;
    wind_P = initial_windP;
    hydro_P(1:1440/time_step) = topology_data(2,3)/1000;
    wind_cutin = topology_data(2,5);
    wind_rated = topology_data(2,6);
    battery_size = topology_data(2,4);
    battery_charge = battery_size * batteries_n;
    max_battery = battery_charge * 0.8; % Include round-trip efficiency
    investment_price = price_data(2,:);
elseif topology == 2
    initial_solarP = topology_data(3,1)/1000;
    solar_panel_P = initial_solarP;
    panelsize = topology_data(3,7);
    initial_windP = topology_data(3,2)/1000;
    wind_P = initial_windP;
    hydro_P(1:1440/time_step) = topology_data(3,3)/1000;
    wind_cutin = topology_data(3,5);
    wind_rated = topology_data(3,6);
    battery_size = topology_data(3,4);
    battery_charge = battery_size * batteries_n;
    max_battery = battery_charge * 0.8; % Include round-trip efficiency
    investment_price = price_data(3,:);
    topology = 1;
end
%==========================================================================

%==========================================================================
% INITIALIZE DEMAND VARIABLES
%==========================================================================
% Obtain holiday list
holiday_tt = pop_data(1,1:13);
% Obtain total population of community
households = pop_data(3,1);
people_perh = pop_data(4,1);
tpop = households * people_perh;
% Obtain population growth rate
pop_growth = pop_data(2,1:6);
pop_growth_r = zeros(5,1);
for i = 1:5
    pop_growth_r(i) = (pop_growth(i+1)-pop_growth(i))/pop_growth(i);
end
pop_growth_r = mean(pop_growth_r);
% Obtain household power demand increase rate when there is rainy weather
pop_r_data = pop_data(5:8,1);
rain_growth_r = ((pop_r_data(4)*pop_r_data(2)) + pop_r_data(1) + pop_r_data(3))/((1-pop_r_data(4))*pop_r_data(2));
% Obtain total number of user classes, users per class, and appliances per class. 
N_tt = app_data(:,2);
N_nan = ~isnan(N_tt);
j_tt = sum(N_nan); % Total number of user classes
N_nan_end = [N_nan; true];
i_tt = zeros(j_tt,1);
current_j = 1;
redo = true;
current_state = true;
i_count = 0;
for i = 1:size(N_tt,1)+1
    if N_nan_end(i) == current_state
        current_state = false;
        i_count = i_count + 1;
    else
        i_tt(current_j) = i_count;
        i_count = 1;
        current_j = current_j + 1;
    end
end
N_tt(isnan(N_tt)) = []; % Total users per class
N_tt_init = sum(N_tt(1:3));
i_tt(:,2) = i_tt(:,1);
i_tt(:,1) = []; % Total appliances per class

% Initialization of other variables
sigma_init = sigma;
monthday_c = 0;
holiday = 1;
total_load = zeros(1440/time_step,j_tt);
total_load_max = zeros(1440/time_step,j_tt);
ttotal_load = zeros(total_time,1440/time_step);
ttotal_load_houses = zeros(total_time,1440/time_step);
ttotal_load_perclass = zeros(total_time,1440/time_step,j_tt);
total_load_conv = zeros(1440/time_step,j_tt,2);
rain_year  = zeros(365,1);
rain_month = zeros(30,1);
%==========================================================================

%==========================================================================
% SET CURRENT DAY'S PARAMETERS FOR DAILY LOAD PROFILE
%==========================================================================
for day = 1:total_time % day: current day being evaluated
    % Day and month counter
    monthday_c = monthday_c + 1;
    month = ceil(monthday_c/31);
    % Convergence test variables
    current_conv = 1;
    mean_test = 1;
    std_test = 1;
    total_load_sum_conv = zeros(1440/time_step,2);
    total_load_sum_houses_conv = zeros(1440/time_step,2);
    mean_test_total_prev = zeros(1440/time_step,1);
    std_test_total_prev = zeros(1440/time_step,1);
    mean_test_total = zeros(1440/time_step,1);
    mean_test_percentage = zeros(1440/time_step,1);
    std_test_total = zeros(1440/time_step,1);
    std_test_percentage = zeros(1440/time_step,1);
    success_mean_k = 0;
    success_std_k = 0;
    % Update population and supply data once a year
    if monthday_c > 365
        monthday_c = 1;
        % Update total population and number of households
        N_tt(3) = N_tt(3) + ceil(pop_growth_r*tpop/people_perh); % Add new households (low-consumption tier)
        tpop = ceil(tpop + (pop_growth_r*tpop));
        % Update solar panel, wind turbine, and battery efficiency
        solar_panel_P = solar_panel_P - (initial_solarP*0.005);
        max_battery = max_battery - (battery_size*batteries_n*(1-0.999772));
        wind_P = wind_P - (initial_windP*(1-0.984));
    end
    % Set electricity increase according to holiday
    if any(ismember(monthday_c,holiday_tt))
        holiday = 0.7112; % 0.2888 reduction in consumption
    else
        holiday = 1;
    end
    % Set electricity increase according to rain
    if dataset_type == 0
        if monthday_c == 1
            for p = 0:11 % Randomize rainy days per month
                precip = ceil(rain_data(p+1));
                rain_month = zeros(30,1);
                rain_month(randperm(30, precip)) = 1;
                rain_year((p*30)+1:(p*30)+30) = rain_month;
            end
        end
        % Adjust daily rain rate to increase the functioning window and time
        if rain_year(monthday_c) == 1
            R_rain = 1+((rain_growth_r-1)*rand(1,1)); % Based on pop. not at home now staying home
            rain_day = 0.5+((1-0.5)*rand(1,1)); % Select rain intensity randomly
            R_rain = ((R_rain-1)*rain_day)+1; % Set rain intensity factor
        else
            R_rain = 1;
        end
    else
        if monthday_c == 1
            for p = 1:8760 % Determine rain intensity per hour
                if rain_data(p) < 2.0
                    rain_year(p) = 0;
                elseif rain_data(p) < 2.5
                    rain_year(p) = 0.25;
                elseif rain_data(p) < 7.5
                    rain_year(p) = 0.5;
                elseif rain_data(p) < 50
                    rain_year(p) = 0.75;
                else
                    rain_year(p) = 1;
                end
            end
        end
        rain_day = max(rain_year((((monthday_c-1)*24)+1):(((monthday_c-1)*24)+24))); % Calculate type of rain for the day
        if rain_day > 0
            R_rain = 1+((rain_growth_r-1)*rand(1,1)); % Based on pop. not at home now staying home
            R_rain = ((R_rain-1)*rain_day)+1; % Set rain intensity factor
        else
            R_rain = 1;
        end
    end
    % Set electricity increase according to daily temperature
    if dataset_type == 0
        temp = weather_data(1:3,month);
        day_temp = round(normrnd(temp(2),temp(2)/3,[1 1]));
        if day_temp > temp(1)
            day_temp = temp(1);
        elseif day_temp < temp(3)
            day_temp = temp(3);
        end
        day_temp_c = 1 + (0.015*abs(day_temp - temp(2))); % 0.015/C increase in functioning time
    else
        temp = weather_data(:,1);
        day_temp = mean(temp((((monthday_c-1)*24)+1):(((monthday_c-1)*24)+24)));
        max_day_temp = max(temp((((monthday_c-1)*24)+1):(((monthday_c-1)*24)+24))) - day_temp;
        min_day_temp = day_temp - min(temp((((monthday_c-1)*24)+1):(((monthday_c-1)*24)+24)));
        % 0.015/C increase/decrease in functioning time
        if max_day_temp > min_day_temp
            day_temp_c = 1 + (0.015*round(max_day_temp));
        else
            day_temp_c = 1 + (0.015*round(min_day_temp));
        end
    end
%==========================================================================

%==========================================================================
% SET CONVERGENCE RUNS FOR EACH DAILY LOAD PROFILE
%==========================================================================
    disp("Day "+day+":");
    fprintf('Convergence runs: 000');
    %for conv = 1:profiles_per_day
    while mean_test == 1 && std_test == 1
        if current_conv < 10
            fprintf('\b%u',current_conv);
        elseif current_conv < 100
            fprintf('\b\b%u',current_conv);
        elseif current_conv < 1000
            fprintf('\b\b\b%u',current_conv);
        end
        init_i = 1;
        final_i = 0;
        % Load Profile cycle for each user class (current_j: current class evaluated)
        for z = 1:j_tt
            current_j = z;
            % Initialize remaining data for current_j
            i_t = i_tt(current_j);
            N = N_tt(current_j);
            final_i = final_i + i_t;
            n = app_data(init_i:final_i,4);
            P = app_data(init_i:final_i,5);
            d = app_data(init_i:final_i,6);
            h = app_data(init_i:final_i,7);
            wf1 = app_data(init_i:final_i,8:9)*60;
            wf2 = app_data(init_i:final_i,10:11)*60;
            wf3 = app_data(init_i:final_i,12:13)*60;
            wf1(wf1==0) = 1;
            wf2(wf2==0) = 1;
            wf3(wf3==0) = 1;
            init_i = init_i + i_t;
%==========================================================================

%==========================================================================
% CALCULATE TOTAL REQUIRED DAILY ELECTRICITY FOR EACH USER CLASS
%==========================================================================
            Ec = 0;
            for j = 1:i_t
                Ec = Ec + (n(j)*P(j)*h(j)*holiday*R_rain*day_temp_c/60);
            end
            Ec = Ec*N; % daily energy consumption for all users in class
%==========================================================================

%==========================================================================
% FUNCTIONING TIME AND WINDOWS RANDOMIZATION
%==========================================================================
            nt = zeros(3,1);
            for i = 1:i_t
                % Randomize functioning time 'h' (as 'nt' cycles)
                nt(i) = round(((h(i)*holiday*R_rain*day_temp_c) + randi([round(-Rh*h(i)) round(Rh*h(i))],1,1)) / d(i));
                % Randomize functioning windows 'wf'
                if sum(isnan(wf1(i,:))) == 0
                    wf1(i,:) = wf1(i,:) + randi([round(-60*Rw*R_rain) round(60*Rw*R_rain)],1,1);
                    for m = 1:2
                        if wf1(i,m) < 1
                            wf1(i,m) = 1;
                        elseif wf1(i,m) > 1440
                            wf1(i,m) = 1440;
                        end
                    end
                end
                if sum(isnan(wf2(i,:))) == 0
                    wf2(i,:) = wf2(i,:) + randi([round(-60*Rw*R_rain) round(60*Rw*R_rain)],1,1);
                    for m = 1:2
                        if wf1(i,m) < 1
                            wf1(i,m) = 1;
                        elseif wf1(i,m) > 1440
                            wf1(i,m) = 1440;
                        end
                    end
                end
                if sum(isnan(wf3(i,:))) == 0
                    wf2(i,:) = wf2(i,:) + randi([round(-60*Rw*R_rain) round(60*Rw*R_rain)],1,1);
                    for m = 1:2
                        if wf1(i,m) < 1
                            wf1(i,m) = 1;
                        elseif wf1(i,m) > 1440
                            wf1(i,m) = 1440;
                        end
                    end
                end
            end
            wf1_nan = isnan(wf1);
            wf2_nan = isnan(wf2);
            wf3_nan = isnan(wf3);
            wf1(wf1_nan)=-1;
            wf2(wf2_nan)=-1;
            wf3(wf3_nan)=-1;
%==========================================================================

%==========================================================================
% CALCULATE POWER PEAK WINDOW
%==========================================================================
            wf_t = [wf1 ; wf2 ; wf3];
            wf=sortrows(wf_t);
            c = 1;
            ct = 1;
            while ct <= size(wf_t,1)
                if wf(c,1) == -1
                    wf(c,:) = [];
                    c = c - 1;
                end
                c = c + 1;
                ct = ct + 1;
            end
            % Interval intersection algorithm
            S = sort(unique(wf(:)));
            count = zeros(1,numel(S)-1);
            for i = 1:numel(S)-1
              I1 = S(i);
              I2 = S(i+1);
              for j = 1:size(wf,1)
                if I1 >= wf(j,1) && I2 <= wf(j,2) 
                  count(i) = count(i) + 1;
                end
              end  
            end
            k = max(count);
            f = find(diff([false,count==k,false])~=0);
            [~,p] = max(S(f(2:2:end))-S(f(1:2:end-1)));
            wf_sorted = S(f(2*p-1):f(2*p))'; % Power peak window
%==========================================================================

%==========================================================================
% CALCULATE DAILY POWER PEAK VALUE AND OCURRENCE TIME
%==========================================================================
            peak_power = 0;
            wf_start = [wf1(:,1) , wf2(:,1) , wf3(:,1)];
            wf_end = [wf1(:,2) , wf2(:,2) , wf3(:,2)];
            app_amount = size(wf1,1);
            wf_peak = zeros(app_amount,1);
            peak_apps = zeros(i_t,1);
            for i = 1:app_amount
                if wf_start(i,1) ~= -1 
                    if wf_start(i,1) < wf_sorted(2) && wf_end(i,1) > wf_sorted(1)
                        if peak_apps(i) == 0
                            peak_power = peak_power + (n(i,1)*P(i,1));
                            peak_apps(i) = 1;
                        end
                    elseif wf_start(i,2) < wf_sorted(2) && wf_end(i,2) > wf_sorted(1)
                        if peak_apps(i) == 0
                            peak_power = peak_power + (n(i)*P(i));
                            peak_apps(i) = 1;
                        end
                    elseif wf_start(i,3) < wf_sorted(2) && wf_end(i,3) > wf_sorted(1)
                        if peak_apps(i) == 0
                            peak_power = peak_power + (n(i,1)*P(i,1));
                            peak_apps(i) = 1;
                        end
                    end
                end
            end
            peak_power = peak_power*N; % theoretical peak power for all users in class
            tp = randi([wf_sorted(1) wf_sorted(2)],1,1); % Peak time randomly chosen with uniform prob.
%==========================================================================

%==========================================================================
% CALCULATE COINCIDENCE FACTOR, LOAD FACTOR, AND ACTUAL POWER PEAK
%==========================================================================
            F = @(x) [(((1/(0.187 + 0.813*exp((-4)*(((1-x(2))^2)+((1-x(2))^16))))) * (1-((1-(0.187 + 0.813*exp((-4)*(((1-x(2))^2)+((1-x(2))^16)))))^(1/x(2))))) * x(2) + (1-((1/(0.187 + 0.813*exp((-4)*(((1-x(2))^2)+((1-x(2))^16))))) * (1-((1-(0.187 + 0.813*exp((-4)*(((1-x(2))^2)+((1-x(2))^16)))))^(1/x(2)))))*x(2))*N^(-1/((1/(0.187 + 0.813*exp((-4)*(((1-x(2))^2)+((1-x(2))^16))))) * (1-((1-(0.187 + 0.813*exp((-4)*(((1-x(2))^2)+((1-x(2))^16)))))^(1/x(2)))))))-x(1);
                (24*x(2)*x(1)*peak_power)-Ec];
            x_init = [0.5;0.5];
            options = optimoptions('fsolve','Display','none');
            x = fsolve(F,x_init,options); % x(1) = fc; x(2) = fl
            pl = peak_power*x(1); % Actual peak power for all users in class
%==========================================================================

%==========================================================================
% LOAD PROFILE COMPUTATION BY SAMPLING NON-PEAK AND PEAK APPLIANCES
%==========================================================================
            while redo == true
                t = zeros(i_t,max(nt, [],"all")*max(n)*max(N_tt));
                wf1(wf1_nan)=0;
                wf2(wf2_nan)=0;
                wf3(wf3_nan)=0;
                for i = 1:i_t
                    range1 = wf1(i,2)-wf1(i,1);
                    range2 = wf2(i,2)-wf2(i,1);
                    range3 = wf3(i,2)-wf3(i,1);
                    wf_i = range1 + range2 + range3;
                    if peak_apps(i) == 0
                        t_temp = randi(wf_i,nt(i)*n(i)*N,1)'; % Non-peak appliances uniform dist.
                    else
                        if tp < wf1(i,2)
                            norm_t = tp-wf1(i,1);
                        elseif tp < wf2(i,2)
                            norm_t = range1 + (tp-wf2(i,1));
                        elseif tp < wf3(i,2)
                            norm_t = range1 + range2 + (tp-wf3(i,1));
                        end
                        t_temp = round(normrnd(norm_t,norm_t/sigma,[(nt(i)*n(i)*N) 1])'); % Peak appliances normal dist.
                    end
                    t_temp(isnan(t_temp))=0;
                    for o = 1:nt(i)*n(i)*N
                        if t_temp(o) <= range1
                            if t_temp(o) < 1
                                t_temp(o) = 1;
                            end
                            t_temp(o) = t_temp(o) + wf1(i,1);
                        elseif t_temp(o) <= (range1 + range2)
                            t_temp(o) = t_temp(o) - range1 + wf2(i,1);
                        else
                            if t_temp(o) > 1440
                                t_temp(o)= 1440;
                            end
                            t_temp(o) = t_temp(o) - range1 - range2 + wf3(i,1);
                        end
                    end
                    for m = 1:size(t_temp,2)
                        t(i,m) = t_temp(m);
                    end
                end
                load_profile = zeros(i_t,1440);
                load_p_temp = zeros(i_t,1440);
                for i = 1:i_t
                    for m = 1:(max(nt, [],"all")*max(n)*max(N_tt))
                        load_start = t(i,m);
                        if load_start ~= 0
                            load_end = t(i,m)+d(i)-1;
                            if load_end > 1440
                                load_end = 1440;
                            end
                            load_p_temp(i,load_start:load_end) = load_p_temp(i,load_start:load_end) + P(i);
                        end
                    end
                    for o = 1:1440
                        if load_p_temp(i,o) > n(i)*P(i)*N
                            load_p_temp(i,o) = n(i)*P(i)*N;
                        end
                    end
                    load_profile(i,:) = load_p_temp(i,:);
                end
                % Final daily load profile for current_j adjusted for selected timestep
                i_k = 0;
                load_sum = zeros(1440,1);
                for i = 1:1440
                    load_sum(i) = load_sum(i) + sum(load_profile(:,i));
                end
                for i = 1:1440/time_step
                    total_load(i,current_j) = mean(load_sum(((time_step*i_k)+1):(time_step*i)));
                    total_load_max(i,current_j) = max(load_sum(((time_step*i_k)+1):(time_step*i)));
                    i_k = i_k + 1;
                end
%==========================================================================

%==========================================================================
% PEAK VALUE CROSSCHECKING AND STANDARD DEVIATION CORRECTION (SHIFTING)
%==========================================================================
                % Crosscheck generated peak value with real value at 25% error margin
                if max(total_load_max(:,current_j)) < pl*0.75
                    sigma = sigma + 0.5;
                    if sigma > 50
                        sigma = sigma_init;
                        redo = false;
                    end
                elseif max(total_load_max(:,current_j)) > pl*1.25
                    sigma = sigma - 0.1;
                    if sigma < 0.1
                        sigma = sigma_init;
                        redo = false;
                    end
                else
                    %disp("Sigma:"+sigma);
                    redo = false;
                    sigma = sigma_init;
                end
            end
            redo = true;            
        end
%==========================================================================

%==========================================================================
% DAILY LOAD PROFILES AGGREGATION OF ALL USER CLASSES
%==========================================================================
        total_load_sum = zeros(1440/time_step,1);
        total_load_sum_houses = zeros(1440/time_step,1);
        for i = 1:1440/time_step
            total_load_sum(i) = sum(total_load(i,:));
            total_load_sum_houses(i) = sum(total_load(i,1:3));
        end
        total_load_sum = total_load_sum/1000; % W → kW for total aggregation
        total_load_sum_houses = total_load_sum_houses/1000; % W → kW for total aggregation
        
        total_load_sum_conv(:,current_conv) = total_load_sum;
        total_load_sum_houses_conv(:,current_conv) = total_load_sum_houses;
        for total_load_j = 1:j_tt
            for total_load_i = 1:1440/time_step
                total_load_conv(total_load_i,total_load_j,current_conv) = total_load(total_load_i,total_load_j);
            end
        end
        switch current_conv
            case 2
                for k = 1:(1440/time_step)
                    mean_test_total_prev(k) = mean(total_load_sum_conv(k,current_conv));
                    std_test_total_prev(k) = std(total_load_sum_conv(k,current_conv));
                end
            otherwise
                for k = 1:(1440/time_step)
                    mean_test_total(k) = mean(total_load_sum_conv(k,1:current_conv));
                    mean_test_percentage(k) = abs((mean_test_total(k) - mean_test_total_prev(k)) / mean_test_total(k));
                    if mean_test_percentage(k) < (convergence_percent/100)
                        success_mean_k = success_mean_k + 1; 
                    end
                    std_test_total(k) = std(total_load_sum_conv(k,1:current_conv));
                    std_test_percentage(k) = abs((std_test_total(k) - std_test_total_prev(k)) / std_test_total(k));
                    if mean_test_percentage(k) < (convergence_percent/100)
                        success_std_k = success_std_k + 1; 
                    end
                end
                %disp("mean test: " + success_mean_k + " & std. test: " + success_std_k);
                if convergence_test == true
                    if success_mean_k >= ((1440/time_step)*0.95) && success_std_k >= ((1440/time_step)*0.95)
                        mean_test = 0;
                        std_test = 0;
                    else
                        success_mean_k = 0; 
                        success_std_k = 0;
                    end
                    mean_test_total_prev = mean_test_total;
                    std_test_total_prev = std_test_total;
                else
                    if current_conv == profiles_per_day
                        mean_test = 0;
                        std_test = 0;
                    end
                end
        end
        current_conv = current_conv + 1;
    end

    total_load_sum_conv(:,1) = [];
    mean_total_sum = mean(total_load_sum_conv,2);

    total_load_sum_houses_conv(:,1) = [];
    mean_total_sum_houses = mean(total_load_sum_houses_conv,2);
    
    mean_total = zeros(1440/time_step,j_tt);
    for total_load_j = 1:j_tt
        for total_load_i = 1:1440/time_step
            mean_total(total_load_i,total_load_j) = mean(total_load_conv(total_load_i,total_load_j,:));
        end
    end
    ttotal_load_perclass(day,:,:) = mean_total/1000;
    ttotal_load(day,:) = mean_total_sum;
    ttotal_load_houses(day,:) = mean_total_sum_houses;
%==========================================================================

%==========================================================================
% CALCULATE SUPPLY-SIDE DAILY PARAMETERS (SOLAR & WIND ONLY)
%==========================================================================
% SOLAR PANELS
    solar_failure_rate = solar_failure_rate + 0.000193656; % Adjust failure rate per day
    if dataset_type == 0
        solar_data_month = solar_data(:,month); % Set solar incidence of day
        solar_data_day = (((solar_data_month(2)-solar_data_month(1))/31)*mod(day,31))+solar_data_month(1);
        solar_curve_scale = rescale(solar_curve,0,solar_data_day);
    else
        solar_curve_scale = solar_data((((monthday_c-1)*24)+1):(((monthday_c-1)*24)+24)); % Set solar incidence of day w/historic data
    end
    for solar_i = 1:solar_n
        r = rand(1,1);
        if solar_failure_rate(solar_i) >= r && solar_down(solar_i) == 0 % Check if panel suffers malfunction
            mu = 2;
            rnexp = exprnd(mu, 1, 50);
            fixlim = fix(rnexp(rnexp>=1 & rnexp<=7)); % How long to repair
            solar_down(solar_i) = fixlim(randi(length(fixlim)));
            used_solar_n = used_solar_n - 1;
            solar_failure_rate(solar_i) = (0.0005*ceil(day/365)) + (((solar_failure_rate(solar_i))-(0.0005*ceil(day/365)))*rand(1,1));
        end
    end
    if R_rain == 1 % Check if sunny day
        solar_dirt_r = solar_dirt_r + 0.00235; % Dirt accumulation until rain
        solar_curve_c = solar_curve_scale .* (solar_panel_P * used_solar_n * panelsize * (1-solar_dirt_r) * 0.8);
    else % Check if rainy day
        if solar_dirt_r < 0.01
            solar_dirt_r = 0;
        else
            solar_dirt_r = 0.01;
        end
        random_raineffect = 0.3 + (0.2)*rand(1,1);
        solar_curve_c = solar_curve_scale .* (solar_panel_P * used_solar_n * panelsize * (1-solar_dirt_r) * 0.8 * random_raineffect);
    end
    if sum(solar_down) > 0 % Reduce days for malfunctioning panel
        solar_down(~solar_down)=1;
        solar_down = solar_down - 1;
        used_solar_n = numel(solar_down(~solar_down));
    end
% WIND TURBINES
    wind_failure_rate = wind_failure_rate + 0.00274; % Adjust failure rate per day
    if dataset_type == 0
        wind_data_month = wind_data(:,month); % Set max wind speed of day
        wind_data_day = (((wind_data_month(2)-wind_data_month(1))/31)*mod(day,31))+wind_data_month(1);
        wind_curve_scale = rescale(wind_curve,1.2,wind_data_day);
    else
        wind_curve_scale = wind_data((((monthday_c-1)*24)+1):(((monthday_c-1)*24)+24)); % Set wind speed of day w/historic data
    end
    for wind_i = 1:size(wind_curve_scale,2)
        if wind_curve_scale(wind_i) > wind_cutin
            wind_to_power(wind_i) = (((wind_P)/(wind_rated-wind_cutin))*(wind_curve_scale(wind_i)-wind_cutin)); % Wind to power eq.
        else
            wind_to_power(wind_i) = 0;
        end
    end
    for wind_i = 1:wind_n
        r = rand(1,1);
        if wind_failure_rate(wind_i) >= r && wind_down(wind_i) == 0 % Check if turbine suffers malfunction
            mu = 2;
            rnexp = exprnd(mu, 1, 50);
            fixlim = fix(rnexp(rnexp>=1 & rnexp<=7)); % How long to repair
            wind_down(wind_i) = fixlim(randi(length(fixlim)));
            used_wind_n = used_wind_n - 1;
            wind_failure_rate(wind_i) = (0.0054*ceil(day/365)) + (((wind_failure_rate(wind_i))-(0.0054*ceil(day/365)))*rand(1,1));
        end
    end
    if R_rain == 1 % Check if sunny day
        wind_curve_c = wind_to_power .* (used_wind_n);
    else % Check if rainy day
        random_raineffect = 0.7 + (0.1)*rand(1,1);
        %wind_curve_c = wind_to_power .* (used_wind_n * random_raineffect);
        wind_curve_c = wind_to_power .* (used_wind_n); % Effect isn't considered
    end
    if sum(wind_down) > 0 % Reduce days for malfunctioning turbine
        wind_down(~wind_down)=1;
        wind_down = wind_down - 1;
        used_wind_n = numel(wind_down(~wind_down));
    end
%==========================================================================

%==========================================================================
% CENTRALIZED & DISTRIBUTED TOPOLOGY BALANCING AND INDICATORS CALCULATION
%==========================================================================
    if topology == 1
        for k = 1:1440/time_step
            total_supply(k) = solar_curve_c(k) + wind_curve_c(k) + hydro_P(k);
            if total_supply(k) >= mean_total_sum(k)
                total_balance(k) = mean_total_sum(k); % Total covered demand
                total_supply(k) = total_supply(k) - mean_total_sum(k);
                battery_charge = battery_charge + (total_supply(k)*(time_step/60));
                if battery_charge > max_battery
                    battery_charge = max_battery;
                end
                total_battery(k) = battery_charge/(time_step/60);
                prev_down_res = 0;
                prev_down_rel = 0;
            elseif (total_supply(k)+(battery_charge/(time_step/60))) >= mean_total_sum(k)
                total_balance(k) = mean_total_sum(k); % Total covered demand
                battery_charge = battery_charge - ((mean_total_sum(k)-total_supply(k))/(time_step/60));
                total_supply(k) = 0;
                total_battery(k) = battery_charge/(time_step/60);
                prev_down_res = 0;
                prev_down_rel = 0;
            else
                % Resilience/Reliability Error
                total_supply(k) = total_supply(k) + (battery_charge/(time_step/60)) - mean_total_sum(k); % Include round-trip efficiency
                total_balance(k) = total_supply(k); % Total covered demand
                battery_charge = 0;
                total_battery(k) = battery_charge;
                affected_households = abs(total_supply(k))/(mean_total_sum_houses(k)/sum(N_tt(1:3)));
                if affected_households > sum(N_tt(1:3))
                    affected_households = sum(N_tt(1:3));
                end
                if R_rain == 1 && sum(solar_down) == 0 && sum(wind_down) == 0
                    if prev_down_rel == 0
                        reliability_counter = reliability_counter + 1;
                        prev_down_rel = 1;
                    end
                    reliability_hours = reliability_hours + 1;
                    avg_households_affected_rel(k*day) = affected_households;
                else
                    if prev_down_res == 0
                        resilience_counter = resilience_counter + 1;
                        prev_down_res = 1;
                    end
                    resilience_hours = resilience_hours + 1;
                    avg_households_affected_res(k*day) = affected_households;
                end
            end
        end
    end
%==========================================================================

%==========================================================================
% STAND-ALONE TOPOLOGY BALANCING AND INDICATORS CALCULATION
%==========================================================================
    if topology == 0
        total_supply_h = zeros(1440/time_step,1);
        for k = 1:1440/time_step
            total_supply_h(k) = solar_curve_c(k) + wind_curve_c(k);
            % FOR HOUSES
            if total_supply_h(k) >= mean_total_sum_houses(k)
                total_balance(k) = mean_total_sum_houses(k); % Total covered demand
                total_supply_h(k) = total_supply_h(k) - mean_total_sum_houses(k);
                total_battery(k) = battery_charge_h/(time_step/60);
                battery_charge_h = battery_charge_h + (total_supply_h(k)*(time_step/60));
                if battery_charge_h > max_battery
                    battery_charge_h = max_battery;
                end
                prev_down_res_h = 0;
                prev_down_rel_h = 0;
            elseif (total_supply_h(k)+(battery_charge_h/(time_step/60))) >= mean_total_sum_houses(k)
                total_balance(k) = mean_total_sum_houses(k); % Total covered demand
                battery_charge_h = battery_charge_h - ((mean_total_sum_houses(k)-total_supply_h(k))/(time_step/60)); % Include round-trip efficiency
                total_supply_h(k) = 0;
                total_battery(k) = battery_charge_h/(time_step/60);
                prev_down_res_h = 0;
                prev_down_rel_h = 0;
            else
                % Resilience/Reliability Error
                total_supply_h(k) = total_supply_h(k) + (battery_charge_h/(time_step/60)) - mean_total_sum_houses(k);
                total_balance(k) = total_supply_h(k); % Total covered demand
                battery_charge_h = 0;
                total_battery(k) = battery_charge_h;
                affected_households = abs(total_supply_h(k))/(mean_total_sum_houses(k)/sum(N_tt(1:3)));
                if affected_households > sum(N_tt(1:3))
                    affected_households = sum(N_tt(1:3));
                end
                if R_rain == 1 && sum(solar_down) == 0 && sum(wind_down) == 0
                    if prev_down_rel_h == 0
                        reliability_counter = reliability_counter + 1;
                        prev_down_rel_h = 1;
                    end
                    reliability_hours = reliability_hours + 1;
                    avg_households_affected_rel(k*day) = affected_households;
                else
                    if prev_down_res_h == 0
                        resilience_counter = resilience_counter + 1;
                        prev_down_res_h = 1;
                    end
                    resilience_hours = resilience_hours + 1;
                    avg_households_affected_res(k*day) = affected_households;
                end
            end
        end
        total_supply = total_supply_h;
    end
%==========================================================================

%==========================================================================
% DAILY SUPPLY BALANCING AND INDICATORS AGGREGATION FOR ALL USER CLASSES
%==========================================================================
    total_panelsperday(day) = used_solar_n;
    total_turbinesperday(day) = used_wind_n;
    if dataset_type == 0
        total_windscale = [total_windscale; wind_curve_scale'];
        total_solarscale = [total_solarscale; solar_curve_scale'];
    end
    for i = 1:1440/time_step
        if total_balance(i) > 0
            total_supply_oandm =  total_supply_oandm + sum(total_balance(i));
        end
    end
    ttotal_balance = [ttotal_balance, total_balance'];
    if day == short_time
        total_supply_oandm_short = total_supply_oandm;
        resilience_counter_short = resilience_counter;
        resilience_hours_short = resilience_hours;
        reliability_counter_short = reliability_counter;
        reliability_hours_short = reliability_hours;
        avg_households_affected_res_short = avg_households_affected_res;
        avg_households_affected_rel_short = avg_households_affected_rel;
        total_panelsperday_short = total_panelsperday(1:short_time);
        total_turbinesperday_short = total_turbinesperday(1:short_time);
        ttotal_balance_short = ttotal_balance;
    elseif day == mid_time
        total_supply_oandm_mid = total_supply_oandm;
        resilience_counter_mid = resilience_counter;
        resilience_hours_mid = resilience_hours;
        reliability_counter_mid = reliability_counter;
        reliability_hours_mid = reliability_hours;
        avg_households_affected_res_mid = avg_households_affected_res;
        avg_households_affected_rel_mid = avg_households_affected_rel;
        total_panelsperday_mid = total_panelsperday(1:mid_time);
        total_turbinesperday_mid = total_turbinesperday(1:mid_time);
        ttotal_balance_mid = ttotal_balance;
    end
    solarrate_history(day,:) = solar_failure_rate*100;
    windrate_history(day,:) = wind_failure_rate*100;
    ttotal_supply(day,:) = total_supply;
    ttotal_battery(day,:) = total_battery/(time_step/60);
    fprintf('\n');
    disp("Day "+day+" completed");
    fprintf('\n');
end
%==========================================================================

%==========================================================================
% TOTAL SUPPLY BALANCING AND INDICATORS AGGREGATION FOR ALL USER CLASSES
%==========================================================================
ttotal_supply_mean = zeros(1440/time_step,1);
ttotal_battery_mean = zeros(1440/time_step,1);
for i = 1:1440/time_step
    ttotal_supply_mean(i) = mean(ttotal_supply(:,i));
    ttotal_battery_mean(i) = mean(ttotal_battery(:,i));
end

% Affected households reliability calculation
avg_households_affected_rel_short(~avg_households_affected_rel_short) = [];
avg_households_affected_rel_short = mean(avg_households_affected_rel_short);
avg_households_affected_rel_mid(~avg_households_affected_rel_mid) = [];
avg_households_affected_rel_mid = mean(avg_households_affected_rel_mid);
avg_households_affected_rel(~avg_households_affected_rel) = [];
avg_households_affected_rel = mean(avg_households_affected_rel);
if isnan(avg_households_affected_rel_short)
    avg_households_affected_rel_short = 0;
else
    avg_households_affected_rel_short = (avg_households_affected_rel_short/sum(N_tt(1:3)))*100;
end
if isnan(avg_households_affected_rel_mid)
    avg_households_affected_rel_mid = 0;
else
    avg_households_affected_rel_mid = (avg_households_affected_rel_mid/sum(N_tt(1:3)))*100;
end
if isnan(avg_households_affected_rel)
    avg_households_affected_rel = 0;
else
    avg_households_affected_rel = (avg_households_affected_rel/sum(N_tt(1:3)))*100;
end
% Affected households resilience calculation
avg_households_affected_res_short(~avg_households_affected_res_short) = [];
avg_households_affected_res_short = mean(avg_households_affected_res_short);
avg_households_affected_res_mid(~avg_households_affected_res_mid) = [];
avg_households_affected_res_mid = mean(avg_households_affected_res_mid);
avg_households_affected_res(~avg_households_affected_res) = [];
avg_households_affected_res = mean(avg_households_affected_res);
if isnan(avg_households_affected_res_short)
    avg_households_affected_res_short = 0;
else
    avg_households_affected_res_short = (avg_households_affected_res_short/sum(N_tt(1:3)))*100;
end
if isnan(avg_households_affected_res_mid)
    avg_households_affected_res_mid = 0;
else
    avg_households_affected_res_mid = (avg_households_affected_res_mid/sum(N_tt(1:3)))*100;
end
if isnan(avg_households_affected_res)
    avg_households_affected_res = 0;
else
    avg_households_affected_res = (avg_households_affected_res/sum(N_tt(1:3)))*100;
end
% Calculate initial investment and LCOE
if topology == 0
    initial_investment = (investment_price(1)*solar_n)+(investment_price(2)*wind_n)+(investment_price(3)*households)+(investment_price(4)*batteries_n)+investment_price(5);
else
    initial_investment = (investment_price(1)*solar_n)+(investment_price(2)*wind_n)+investment_price(3)+(investment_price(4)*batteries_n)+investment_price(5);
end
initial_investment_household = initial_investment/N_tt_init;
total_oandm_short = (17*initial_solarP*solar_n*(short_time/365))+(26.2*initial_windP*wind_n*(short_time/365));
total_oandm_mid = (17*initial_solarP*solar_n*(mid_time/365))+(26.2*initial_windP*wind_n*(mid_time/365));
total_oandm = (17*initial_solarP*solar_n*(long_time/365))+(26.2*initial_windP*wind_n*(long_time/365));
lcoe_short = (initial_investment + total_oandm_short) / total_supply_oandm_short;
lcoe_mid = (initial_investment + total_oandm_mid) / total_supply_oandm_mid;
lcoe = (initial_investment + total_oandm) / total_supply_oandm;
% Print techno-economic indicators
if mix_test == false
    disp("------------------------");
    disp("SHORT TERM:");
    disp("------------------------");
    disp("ECONOMIC INDICATORS:");
    disp("------------------------");
    disp("Initial Investment: $"+initial_investment+" or $"+initial_investment_household+"/household");
    disp("LCOE: $"+lcoe_short+"/kWh");
    disp("RESILIENCE INDICATORS:");
    disp("Blackouts: "+resilience_counter_short);
    disp("Total off time: "+resilience_hours_short);
    disp("Avg. impacted households: "+ceil(avg_households_affected_res_short)+"%");
    disp("------------------------");
    disp("RELIABILITY INDICATORS:");
    disp("------------------------");
    disp("Blackouts: "+reliability_counter_short);
    disp("Total off time: "+reliability_hours_short);
    disp("Avg. impacted households: "+ceil(avg_households_affected_rel_short)+"%");
    disp("------------------------");
    results_short(1,:) = [initial_investment initial_investment_household lcoe_short];
    results_short(2,:) = [resilience_counter_short resilience_hours_short ceil(avg_households_affected_res_short)];
    results_short(3,:) = [reliability_counter_short reliability_hours_short ceil(avg_households_affected_rel_short)];
    disp("------------------------");
    disp("MID TERM:");
    disp("------------------------");
    disp("ECONOMIC INDICATORS:");
    disp("------------------------");
    disp("Initial Investment: $"+initial_investment+" or $"+initial_investment_household+"/household");
    disp("LCOE: $"+lcoe_mid+"/kWh");
    disp("RESILIENCE INDICATORS:");
    disp("Blackouts: "+resilience_counter_mid);
    disp("Total off time: "+resilience_hours_mid);
    disp("Avg. impacted households: "+ceil(avg_households_affected_res_mid)+"%");
    disp("------------------------");
    disp("RELIABILITY INDICATORS:");
    disp("------------------------");
    disp("Blackouts: "+reliability_counter_mid);
    disp("Total off time: "+reliability_hours_mid);
    disp("Avg. impacted households: "+ceil(avg_households_affected_rel_mid)+"%");
    disp("------------------------");
    results_mid(1,:) = [initial_investment initial_investment_household lcoe_mid];
    results_mid(2,:) = [resilience_counter_mid resilience_hours_mid ceil(avg_households_affected_res_mid)];
    results_mid(3,:) = [reliability_counter_mid reliability_hours_mid ceil(avg_households_affected_rel_mid)];
    disp("------------------------");
    disp("LONG TERM:");
    disp("------------------------");
    disp("ECONOMIC INDICATORS:");
    disp("------------------------");
    disp("Initial Investment: $"+initial_investment+" or $"+initial_investment_household+"/household");
    disp("LCOE: $"+lcoe+"/kWh");
    disp("RESILIENCE INDICATORS:");
    disp("Blackouts: "+resilience_counter);
    disp("Total off time: "+resilience_hours);
    disp("Avg. impacted households: "+ceil(avg_households_affected_res)+"%");
    disp("------------------------");
    disp("RELIABILITY INDICATORS:");
    disp("------------------------");
    disp("Blackouts: "+reliability_counter);
    disp("Total off time: "+reliability_hours);
    disp("Avg. impacted households: "+ceil(avg_households_affected_rel)+"%");
    disp("------------------------");
    results_long(1,:) = [initial_investment initial_investment_household lcoe];
    results_long(2,:) = [resilience_counter resilience_hours ceil(avg_households_affected_res)];
    results_long(3,:) = [reliability_counter reliability_hours ceil(avg_households_affected_rel)];
    writematrix(results_short,fullfile(plot_path,'results_2Y.txt'));
    writematrix(results_mid,fullfile(plot_path,'results_10Y.txt'));
    writematrix(results_long,fullfile(plot_path,'results_15Y.txt'));
end
%==========================================================================

%==========================================================================
% TOTAL LOAD PROFILES AGGREGATION OF ALL USER CLASSES
%==========================================================================
ttotal_load_mean_short = zeros(1440/time_step,1);
ttotal_load_max_short = zeros(1440/time_step,1);
ttotal_load_min_short = zeros(1440/time_step,1);
ttotal_load_mean_mid = zeros(1440/time_step,1);
ttotal_load_max_mid = zeros(1440/time_step,1);
ttotal_load_min_mid = zeros(1440/time_step,1);
ttotal_load_mean = zeros(1440/time_step,1);
ttotal_load_max = zeros(1440/time_step,1);
ttotal_load_min = zeros(1440/time_step,1);
for i = 1:1440/time_step
    if topology == 0
        ttotal_load_mean_short(i) = mean(ttotal_load_houses(1:short_time,i));
        ttotal_load_max_short(i) = max(ttotal_load_houses(1:short_time,i));
        ttotal_load_min_short(i) = min(ttotal_load_houses(1:short_time,i));

        ttotal_load_mean_mid(i) = mean(ttotal_load_houses(1:mid_time,i));
        ttotal_load_max_mid(i) = max(ttotal_load_houses(1:mid_time,i));
        ttotal_load_min_mid(i) = min(ttotal_load_houses(1:mid_time,i));

        ttotal_load_mean(i) = mean(ttotal_load_houses(:,i));
        ttotal_load_max(i) = max(ttotal_load_houses(:,i));
        ttotal_load_min(i) = min(ttotal_load_houses(:,i));
    else
        ttotal_load_mean_short(i) = mean(ttotal_load(1:short_time,i));
        ttotal_load_max_short(i) = max(ttotal_load(1:short_time,i));
        ttotal_load_min_short(i) = min(ttotal_load(1:short_time,i));

        ttotal_load_mean_mid(i) = mean(ttotal_load(1:mid_time,i));
        ttotal_load_max_mid(i) = max(ttotal_load(1:mid_time,i));
        ttotal_load_min_mid(i) = min(ttotal_load(1:mid_time,i));

        ttotal_load_mean(i) = mean(ttotal_load(:,i));
        ttotal_load_max(i) = max(ttotal_load(:,i));
        ttotal_load_min(i) = min(ttotal_load(:,i));
    end
end

ttotal_load_mean_perclass = zeros(1440/time_step,j_tt);
for j = 1:j_tt
    for i = 1:1440/time_step
    ttotal_load_mean_perclass(i,j) = mean(ttotal_load_perclass(:,i,j));
    end
end

load_household = zeros(1440/time_step,1);
load_community = zeros(1440/time_step,1);
for i = 1:1440/time_step
    load_household(i) = sum(ttotal_load_mean_perclass(i,1:3));
    load_community(i) = sum(ttotal_load_mean_perclass(i,4:7));
end
%==========================================================================

%==========================================================================
% SET TECHNOLOGY MIX
%==========================================================================
if mix_test == true
    ttotal_wind_mean = zeros(1440/time_step,1);
    ttotal_solar_mean = zeros(1440/time_step,1);
    hours_s = zeros(1440/time_step,1);
    hours_n = zeros(1440/time_step,1);
    meansolarPperunit = zeros(1440/time_step,1); 
    meanwindPperunit = zeros(1440/time_step,1);
    if dataset_type == 0
        wind_data = total_windscale;
        solar_data = total_solarscale;
    end
    for i = 1:1440/time_step
        total_wind_hour = [];
        total_solar_hour = [];
        for j = 0:364
            total_wind_hour = [total_wind_hour wind_data(i+(24*j))];
            total_solar_hour = [total_solar_hour solar_data(i+(24*j))];
        end
        ttotal_wind_mean(i) = mean(total_wind_hour);
        ttotal_solar_mean(i) = mean(total_solar_hour);
    end
    for i = 1:1440/time_step
        meanwindPperunit(i) = ((initial_windP)/(wind_rated-wind_cutin))*(ttotal_wind_mean(i)-wind_cutin);
        meansolarPperunit(i) = ttotal_solar_mean(i)*initial_solarP*panelsize*0.8;
        if ttotal_solar_mean(i) < 0.05
            hours_n(i) = 1;
        end
        if meanwindPperunit(i) < 0
            meanwindPperunit(i) = 0;
        end
    end
    load_down_kwh = sum(ttotal_load_max.*hours_n);
    [max_load_down,max_down_h] = max(ttotal_load_max.*hours_n);
    hours_s(:) = ~hours_n;
    load_up_kwh = sum(ttotal_load_max.*hours_s);
    [max_load_up,max_up_h] = max(ttotal_load_max.*hours_s);
    total_wind_kwh_up = sum(meanwindPperunit.*hours_s);
    total_wind_kwh_down = sum(meanwindPperunit.*hours_n);
    total_hydro_kwh_up = sum(hydro_P.*hours_s');
    total_hydro_kwh_down = sum(hydro_P.*hours_n');
    total_solar_kwh = sum(meansolarPperunit.*hours_s);
    % Optimization problem
    if topology == 0
        total_hydro_kwh_up = 0;
        total_hydro_kwh_down = 0;
    end
    s = optimvar('s',1,'LowerBound',1);
    w = optimvar('w',1,'LowerBound',1);
    b = optimvar('b',1,'LowerBound',1);
    prob = optimproblem('ObjectiveSense','min');
    prob.Objective = (investment_price(1)*s)+(investment_price(2)*w)+(investment_price(4)*b);
    prob.Constraints.cons1 = (total_wind_kwh_up*w) + (total_solar_kwh*s) - (battery_size*b) >= (load_up_kwh - total_hydro_kwh_up);
    prob.Constraints.cons2 = (total_wind_kwh_down*w) + (battery_size*b) >= (load_down_kwh - total_hydro_kwh_down);
    prob.Constraints.cons3 = (meanwindPperunit(max_up_h)*w) + (meansolarPperunit(max_up_h)*s) >= (max_load_up - hydro_P);
    prob.Constraints.cons4 = (meanwindPperunit(max_down_h)*w) + ((battery_size/sum(hours_n))*b) >= (max_load_down - hydro_P);
    sol = solve(prob);
    if sol.s == 1 && total_solar_kwh > 0.1
        auto_result = [ceil(sol.w*R_factor/6) ceil(sol.w*R_factor/1.2) ceil(sol.b*R_factor)];
    elseif sol.w == 1 && (total_wind_kwh_down+total_wind_kwh_up) > 0.1
        auto_result = [ceil(sol.s*R_factor/1.2) ceil(sol.s*R_factor/6) ceil(sol.b*R_factor)];
    else
        auto_result = [ceil(sol.s*R_factor/1.2) ceil(sol.w*R_factor/1.2) ceil(sol.b*R_factor)];
    end
    disp("----------------------");
    disp("SUGGESTED TECH. MIX:");
    disp("----------------------");
    disp("Solar panels: "+auto_result(1));
    disp("Wind turbines: "+auto_result(2));
    disp("Batteries: "+auto_result(3));
    disp("----------------------");
end
%==========================================================================

%==========================================================================
% PLOT CURVES FOR DIFFERENT CHARACTERISTICS
%==========================================================================
% SHORT TERM PLOTS
beep
if plot_results == true && mix_test == false
    X = 1:1440/time_step;
% Plot curve for class sum
    figure(1)
    fitobject = fit(X',ttotal_load_max_short,'cubicspline');
    fitobject3 = fit(X',ttotal_load_min_short,'cubicspline');
    fitobject2 = fit(X',ttotal_load_mean_short,'cubicspline');
    p1 = plot(fitobject);
    hold on
    p1.Color = "#939393";
    p2 = plot(fitobject2);
    p2.Color = "black";
    p3 = plot(fitobject3);
    p3.Color = "#474747";
    title('Average Total Load Forecast Profile - Short-term');
    xlabel('Hour');
    ylabel('Power Demand (kW)');
    legend({'Max Load','Mean Load','Min Load'},'Location','northwest');
    xticks(linspace(1,1440/time_step,24));
    xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'});
    set(gcf, 'Position',  [800, 500, 1000, 500]);
    xlim([1 1440/time_step]);
    ylim([0 max(ttotal_load_max_short)+1]);
    ax = gca;
    ax.YGrid = 'on';
    filename = 'Total daily load forecast 2Y';
    saveas(figure(1),fullfile(plot_path, filename));
    saveas(figure(1),fullfile(plot_path, filename),'bmp');
% Plot curve for solar panel status
    figure(2)
    fig_bins = max(total_panelsperday_short)-min(total_panelsperday_short)+1;
    histogram(total_panelsperday_short,fig_bins,'BinLimits',[min(total_panelsperday_short)-0.5 max(total_panelsperday_short)+0.5])
    title('Frequency of Solar Panels in Operation - Short-term');
    xlabel('Operating Units');
    ylabel('Total Days');
    set(gcf, 'Position',  [800, 500, 1000, 500]);
    xlim([min(total_panelsperday_short)-1 solar_n+1]);
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'on';
    filename = 'Solar panel frequency 2Y';
    saveas(figure(2),fullfile(plot_path, filename));
    saveas(figure(2),fullfile(plot_path, filename),'bmp');
% Plot curve for wind turbine status
    figure(3)
    fig_bins = max(total_turbinesperday_short)-min(total_turbinesperday_short)+1;
    histogram(total_turbinesperday_short,fig_bins,'BinLimits',[min(total_turbinesperday_short)-0.5 max(total_turbinesperday_short)+0.5])
    title('Frequency of Wind Turbines in Operation - Short-term');
    xlabel('Operating Units');
    ylabel('Total Days');
    set(gcf, 'Position',  [800, 500, 1000, 500]);
    xlim([min(total_turbinesperday_short)-1 wind_n+1]);
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'on';
    filename = 'Wind turbine frequency 2Y';
    saveas(figure(3),fullfile(plot_path, filename));
    saveas(figure(3),fullfile(plot_path, filename),'bmp');
% Plot curve for supply-demand balance history
    figure(4)
    fig_bins = ceil(max(ttotal_balance_short))-floor(min(ttotal_balance_short));
    histogram(ttotal_balance_short, fig_bins,'BinLimits',[floor(min(ttotal_balance_short)) ceil(max(ttotal_balance_short))],'Orientation','horizontal');
    title('Supply-Demand Balance - Short-term');
    xlabel('Total Hours');
    ylabel('Electricity Demand (kWh)');
    set(gcf, 'Position',  [800, 500, 1000, 500]);
    set(gca, 'XScale', 'log');
    xlim([0.9 inf]);
    ax = gca;
    ax.YGrid = 'on';
    ax.XGrid = 'on';
    ax.XMinorGrid = 'on';
    filename = 'Supply demand balance 2Y';
    saveas(figure(4),fullfile(plot_path, filename));
    saveas(figure(4),fullfile(plot_path, filename),'bmp');
% MID TERM PLOTS
% Plot curve for class sum
    figure(5)
    fitobject = fit(X',ttotal_load_max_mid,'cubicspline');
    fitobject3 = fit(X',ttotal_load_min_mid,'cubicspline');
    fitobject2 = fit(X',ttotal_load_mean_mid,'cubicspline');
    p1 = plot(fitobject);
    hold on
    p1.Color = "#939393";
    p2 = plot(fitobject2);
    p2.Color = "black";
    p3 = plot(fitobject3);
    p3.Color = "#474747";
    title('Average Total Load Forecast Profile - Mid-term');
    xlabel('Hour');
    ylabel('Power Demand (kW)');
    legend({'Max Load','Mean Load','Min Load'},'Location','northwest');
    xticks(linspace(1,1440/time_step,24));
    xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'});
    set(gcf, 'Position',  [800, 500, 1000, 500]);
    xlim([1 1440/time_step]);
    ylim([0 max(ttotal_load_max_mid)+1]);
    ax = gca;
    ax.YGrid = 'on';
    filename = 'Total daily load forecast 10Y';
    saveas(figure(5),fullfile(plot_path, filename));
    saveas(figure(5),fullfile(plot_path, filename),'bmp');
% Plot curve for solar panel status
    figure(6)
    fig_bins = max(total_panelsperday_mid)-min(total_panelsperday_mid)+1;
    histogram(total_panelsperday_mid,fig_bins,'BinLimits',[min(total_panelsperday_mid)-0.5 max(total_panelsperday_mid)+0.5])
    title('Frequency of Solar Panels in Operation - Mid-term');
    xlabel('Operating Units');
    ylabel('Total Days');
    set(gcf, 'Position',  [800, 500, 1000, 500]);
    xlim([min(total_panelsperday_mid)-1 solar_n+1]);
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'on';
    filename = 'Solar panel frequency 10Y';
    saveas(figure(6),fullfile(plot_path, filename));
    saveas(figure(6),fullfile(plot_path, filename),'bmp');
% Plot curve for wind turbine status
    figure(7)
    fig_bins = max(total_turbinesperday_mid)-min(total_turbinesperday_mid)+1;
    histogram(total_turbinesperday_mid,fig_bins,'BinLimits',[min(total_turbinesperday_mid)-0.5 max(total_turbinesperday_mid)+0.5])
    title('Frequency of Wind Turbines in Operation - Mid-term');
    xlabel('Operating Units');
    ylabel('Total Days');
    set(gcf, 'Position',  [800, 500, 1000, 500]);
    xlim([min(total_turbinesperday_mid)-1 wind_n+1]);
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'on';
    filename = 'Wind turbine frequency 10Y';
    saveas(figure(7),fullfile(plot_path, filename));
    saveas(figure(7),fullfile(plot_path, filename),'bmp');
% Plot curve for supply-demand balance history
    figure(8)
    fig_bins = ceil(max(ttotal_balance_mid))-floor(min(ttotal_balance_mid));
    histogram(ttotal_balance_mid, fig_bins,'BinLimits',[floor(min(ttotal_balance_mid)) ceil(max(ttotal_balance_mid))],'Orientation','horizontal');
    title('Supply-Demand Balance - Mid-term');
    xlabel('Total Hours');
    ylabel('Electricity Demand (kWh)');
    set(gcf, 'Position',  [800, 500, 1000, 500]);
    set(gca, 'XScale', 'log');
    xlim([0.9 inf]);
    ax = gca;
    ax.YGrid = 'on';
    ax.XGrid = 'on';
    ax.XMinorGrid = 'on';
    filename = 'Supply demand balance 10Y';
    saveas(figure(8),fullfile(plot_path, filename));
    saveas(figure(8),fullfile(plot_path, filename),'bmp');
% LONG TERM PLOTS
    % Plot curve for each class
    figure(9)
    plot(X,load_household);
    hold on
    plot(X,load_community);
    legend({'Household Load','Community Load'},'Location','northwest');
    title('Daily Load Forecast Profile per Class - Long-term');
    xlabel('Hours');
    ylabel('Power Demand (kW)');
    xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'});
    set(gcf, 'Position',  [800, 500, 1000, 500]);
    xlim([1 1440/time_step]);
    ax = gca;
    ax.YGrid = 'on';
    filename = 'Daily load forecast per class 15Y';
    saveas(figure(9),fullfile(plot_path, filename));
    saveas(figure(9),fullfile(plot_path, filename),'bmp');
% Plot curve for class sum
    figure(10)
    fitobject = fit(X',ttotal_load_max,'cubicspline');
    fitobject3 = fit(X',ttotal_load_min,'cubicspline');
    fitobject2 = fit(X',ttotal_load_mean,'cubicspline');
    p1 = plot(fitobject);
    hold on
    p1.Color = "#939393";
    p2 = plot(fitobject2);
    p2.Color = "black";
    p3 = plot(fitobject3);
    p3.Color = "#474747";
    title('Average Total Load Forecast Profile - Long-term');
    xlabel('Hour');
    ylabel('Power Demand (kW)');
    legend({'Max Load','Mean Load','Min Load'},'Location','northwest');
    xticks(linspace(1,1440/time_step,24));
    xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'});
    set(gcf, 'Position',  [800, 500, 1000, 500]);
    xlim([1 1440/time_step]);
    ylim([0 max(ttotal_load_max)+1]);
    ax = gca;
    ax.YGrid = 'on';
    filename = 'Total daily load forecast 15Y';
    saveas(figure(10),fullfile(plot_path, filename));
    saveas(figure(10),fullfile(plot_path, filename),'bmp');
% Plot curve for solar panel status
    figure(11)
    fig_bins = max(total_panelsperday)-min(total_panelsperday)+1;
    histogram(total_panelsperday,fig_bins,'BinLimits',[min(total_panelsperday)-0.5 max(total_panelsperday)+0.5])
    title('Frequency of Solar Panels in Operation - Long-term');
    xlabel('Operating Units');
    ylabel('Total Days');
    set(gcf, 'Position',  [800, 500, 1000, 500]);
    xlim([min(total_panelsperday)-1 solar_n+1]);
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'on';
    filename = 'Solar panel frequency 15Y';
    saveas(figure(11),fullfile(plot_path, filename));
    saveas(figure(11),fullfile(plot_path, filename),'bmp');
% Plot curve for wind turbine status
    figure(12)
    fig_bins = max(total_turbinesperday)-min(total_turbinesperday)+1;
    histogram(total_turbinesperday,fig_bins,'BinLimits',[min(total_turbinesperday)-0.5 max(total_turbinesperday)+0.5])
    title('Frequency of Wind Turbines in Operation - Long-term');
    xlabel('Operating Units');
    ylabel('Total Days');
    set(gcf, 'Position',  [800, 500, 1000, 500]);
    xlim([min(total_turbinesperday)-1 wind_n+1]);
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'on';
    filename = 'Wind turbine frequency 15Y';
    saveas(figure(12),fullfile(plot_path, filename));
    saveas(figure(12),fullfile(plot_path, filename),'bmp');
% Plot curve for supply-demand balance history
    figure(13)
    fig_bins = ceil(max(ttotal_balance))-floor(min(ttotal_balance));
    histogram(ttotal_balance, fig_bins,'BinLimits',[floor(min(ttotal_balance)) ceil(max(ttotal_balance))],'Orientation','horizontal');
    title('Supply-Demand Balance - Long-term');
    xlabel('Total Hours');
    ylabel('Electricity Demand (kWh)');
    set(gcf, 'Position',  [800, 500, 1000, 500]);
    set(gca, 'XScale', 'log');
    xlim([0.9 inf]);
    ax = gca;
    ax.YGrid = 'on';
    ax.XGrid = 'on';
    ax.XMinorGrid = 'on';
    filename = 'Supply demand balance 15Y';
    saveas(figure(13),fullfile(plot_path, filename));
    saveas(figure(13),fullfile(plot_path, filename),'bmp');
end
%==========================================================================