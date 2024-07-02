classdef Model_NMC_Gr
    %MODEL_NMC_Gr Kokam NMC-Gr dynamic aging model
    % Fit to accelerated aging data conducted at NREL for Sunpower
    % See data reference at https://ieeexplore.ieee.org/iel7/7951530/7962914/07963578.pdf
    % Calendar aging with T, SOC dependence
    % Cycle aging with T, estimated C-rate and DOD dependence
    % Capacity degradation model is a competitive model of the form:
    %   q = min(q_LLI, q_LAM)
    %   q_LLI = 1 - qLoss_calendar(t, T, SOC) - qLoss_cycling(N, T, DOD)
    %   q_LAM = 1.01 - qLoss_LAM(N, T, DOD)
    % 3/11/2024 - NOTE
    % Resistance model is only for total DC resistance. 

    properties
        % Capacity degradation model parameters and functions
        params_q
        submodels_q
        % Resistance degradation model parameters and functions
        params_r
        submodels_r
        % Degradation model outputs and states
        Outputs
        States
        % Nominal values for capacity and resistnace
        Battery_Properties
    end
    
    methods
        function obj = Model_NMC_Gr()
            %MODEL_NMC_Gr Construct an instance of this class
            %%% Nominal battery properties %%%
            % these aren't really used, BatteryNMC overrides
            obj.Battery_Properties.Q_Ah = 75;
            obj.Battery_Properties.R_Ohms = 1e-3;

            %%% Load test data, just need variable names
            load('nmc_gr\NREL_Kokam75Ah_LifeTest_04Apr2016_add.mat', 'kokam') % data, need data variables to export function handle
            Data = dataStructToTable(kokam); 
            clearvars kokam
            Data.EFC = Data.Ahtp ./ 75;
            Data.soc = Data.socavg;
            x_a_eq = @(SOC) 8.5e-3 + SOC.*(7.8e-1 - 8.5e-3);
            Ua_eq = @(x_a) 0.6379 + 0.5416.*exp(-305.5309.*x_a) + 0.044.*tanh(-1.*(x_a-0.1958)./0.1088) - 0.1978.*tanh((x_a-1.0571)./0.0854) - 0.6875.*tanh((x_a+0.0117)./0.0529) - 0.0175.*tanh((x_a-0.5692)./0.0875);
            Data.Ua = Ua_eq(x_a_eq(Data.soc));
            Data.TdegKN = Data.TdegK ./ (273.15+35); % normalize at 35C (mean temp of aging)
            Data.UaN = Data.Ua ./ 0.123; % normalize at 50% soc
            data_vars = Data.Properties.VariableNames;
            data_vars = data_vars(3:end);
            clearvars Data
            
            %%% Parameters and funtions %%%
            % Load ReducedOrderModel objects for parameter values and
            % exported function handles for rate equations
            
            % Capacity fade parameters and functions
            load('nmc_gr\Model_q.mat') % model object
            % degradation rate equation function handles
            fh_q1 = exportFuncHandle(Mdl_q.SubModels(1).SubModels(1), data_vars);
            fh_q3 = exportFuncHandle(Mdl_q.SubModels(1).SubModels(2), data_vars);
            fh_q5 = exportFuncHandle(Mdl_q.SubModels(2).SubModels(1), data_vars);
            % Load parameter names and values
            pvars  = getAllParamCellstr(Mdl_q);
            [p, ~] = getParamValues(Mdl_q, pvars);
            p      = array2table(p, 'VariableNames', pvars);
            % Store parameters and function handles
            obj.params_q = p;
            obj.submodels_q = {fh_q1, fh_q3, fh_q5};

            % Resistance growth parameters
            load('nmc_gr\Model_r.mat', 'Mdl_r') % model object
            pvars  = getAllParamCellstr(Mdl_r);
            [p, ~] = getParamValues(Mdl_r, pvars);
            p      = array2table(p, 'VariableNames', pvars);
            obj.params_r = p;
            obj.submodels_r = []; % no unique submodels; uses capacity fade submodels

            %%% Outputs: initial values %%%
            outputs = {'q','r'};
            intial_outputs = ones(1, length(outputs));
            obj.Outputs = array2table(intial_outputs, 'VariableNames', outputs);
            %%% States: initial values %%%
            states = {'q_loss_LLI_cal', 'q_loss_LLI_cyc', 'q_loss_LAM', 'r_gain_LLI_cal', 'r_gain_LLI_cyc', 'r_gain_LAM'};
            initial_states = zeros(1, length(states));
            obj.States = array2table(initial_states, 'VariableNames', states);
        end

        function obj = solve_aging(obj, ElectroThermalResponse)
            %SOLVE_AGING Solve aging model
            %   Extract battery stress from its ElectroThermalResponse from
            %   a simulation and update battery health states like capacity
            %   and resistance.
            arguments
                obj Model_NMC_Gr
                ElectroThermalResponse table
            end

            % Extract stressors
            Stressors = obj.extract_stressors(ElectroThermalResponse);
            % Prior q_LAM value for resistance model
            q_LAM = 1.01 - obj.States.q_loss_LAM;
            % Update capacity
            obj = degrade_capacity(obj, Stressors);
            % Update resistance
            obj = degrade_resistance(obj, Stressors, q_LAM);
        end

        function obj = degrade_capacity(obj, Stressors)
            % States = obj.States;
            Params = obj.params_q;
            SubModels = obj.submodels_q;

            % Loss of lithium inventory, calendar aging: qLossCal
            % Rate and order parameters
            k = SubModels{1}; % function_handle
            p = Params.q2;
            % Prior value of state and independent variable change
            y = obj.States.q_loss_LLI_cal;
            dx = Stressors.dt;
            % Calculate state change
            dy = dynamicPower(y, dx, [k(Stressors, Params), p]);
            % Update state value
            obj.States.q_loss_LLI_cal = y + dy;

            % Loss of lithium inventory, cycle aging: qLossCyc
            % Rate and order parameters
            k = SubModels{2}; % function_handle
            p = Params.q4;
            % Prior value of state and independent variable change
            y = obj.States.q_loss_LLI_cyc;
            dx = Stressors.dEFC;
            % Calculate state change
            dy = dynamicPower(y, dx, [abs(k(Stressors, Params)), p]);
            % Update state value
            obj.States.q_loss_LLI_cyc = y + dy;

            % Loss of active material (accelerating fade): LAM
            % Sigmoid parameters:
            %   yInf: value of y at x=Inf
            %   tau: 1/x_halfMax
            %   p: order
            yInf = 1;
            tau0 = 8e3;
            tau = SubModels{3}; % function_handle
            tau = tau(Stressors, Params);
            tau = max(tau, tau0); % prevent tau from becoming unrealistically large at low DOD or high temp
            tau = 1 ./ tau;
            p = Params.pLAM;
            % Prior value of state and independent variable change
            y = obj.States.q_loss_LAM;
            dx = Stressors.dEFC;
            % Calculate state change
            dy = dynamicSigmoid(y, dx, [yInf, tau, p]);
            % Update state value
            obj.States.q_loss_LAM = y + dy;

            % Capacity loss
            q_LLI = 1 - obj.States.q_loss_LLI_cal - obj.States.q_loss_LLI_cyc;
            q_LAM = 1.01 - obj.States.q_loss_LAM;
            q = min([q_LLI, q_LAM], [], 2);
            obj.Outputs.q = q;
        end

        function obj = degrade_resistance(obj, Stressors, q_LAM)
            % Calculate change in resistance and update state and output
            % values

            % all submodels are the same as the capacity fade models
            SubModels = obj.submodels_q;
            % need both resistance and capacity fade model parameters
            Params = [obj.params_q, obj.params_r];

            % Loss of lithium inventory, calendar aging: r_gain_LLI_cal
            % Rate and order parameters
            k = SubModels{1}; % function_handle
            p = Params.r2;
            % Prior value of state and independent variable change
            y = obj.States.r_gain_LLI_cal;
            dx = Stressors.dt;
            % Calculate state change
            dy = dynamicPower(y, dx, [Params.r1 .* k(Stressors, Params), p]);
            % Update state value
            obj.States.r_gain_LLI_cal = y + dy;

            % Loss of lithium inventory, cycle aging: r_gain_LLI_cyc
            % Rate and order parameters
            k = SubModels{2}; % function_handle
            p = Params.r4;
            % Prior value of state and independent variable change
            y = obj.States.r_gain_LLI_cyc;
            dx = Stressors.dEFC;
            % Calculate state change
            dy = dynamicPower(y, dx, [Params.r3 .* abs(k(Stressors, Params)), p]);
            % Update state value
            obj.States.r_gain_LLI_cyc = y + dy;

            % Loss of active material (accelerating fade): r_LAM
            y = obj.States.r_gain_LAM;
            % Simply linear with LAM
            x0 = 1.0123 - q_LAM;
            x = 1.0123 - (1.01 - obj.States.q_loss_LAM);
            dx = x - x0;
            % Calculate state change
            dy = Params.r6 * dx;
            % Update state value
            obj.States.r_gain_LAM = y + dy;

            % Resistance growth
            r_LLI = 1 + obj.States.r_gain_LLI_cal + obj.States.r_gain_LLI_cyc;
            r_LAM = 1 + obj.States.r_gain_LAM;
            r = max([r_LLI, r_LAM], [], 2);
            obj.Outputs.r = r;
        end
    end

    methods (Static)
        function Stressors = extract_stressors(ElectroThermalResponse, nominal_capacity)
        % EXTRACT_STRESSORS Output battery stress values from battery measurements
        %{
        Input timeseries:
            cycle.t:     time in days
            cycle.soc:   soc in normal units, 0 to 1
            cycle.TdegC: temperature in celsius
        Outputs:
            cycle.dt:    Time difference in days
            cycle.dEFC:  Equvialent full cycles during this period
            cycle.TdegK: Time weighted average of temperature
            cycle.soc:   Time weighted average of soc
            cycle.Ua:    Time weighted average of Ua
            cycle.dod:   Difference b/w max and min SOC
            cycle.Crate: Time weighted average of Crate in non-resting periods
        
        AuDOD variables
        t, N, DOD, TdegK, soc, Ua
        %}
        arguments
            ElectroThermalResponse table
            nominal_capacity = 1
        end
        % Time
        dt_s = ElectroThermalResponse.Time_s(end) - ElectroThermalResponse.Time_s(1);
        dt = dt_s ./ (24*3600); % Time in seconds to time in days
        
        % Temperature
        TdegC = trapz(ElectroThermalResponse.Time_s, ElectroThermalResponse.Temperature_C) / dt_s;
        TdegK = TdegC + 273.15;
        TdegKN = TdegK ./ (273.15+35); % Normalized temperature
        
        % SOC and voltage
        soc = trapz(ElectroThermalResponse.Time_s, ElectroThermalResponse.SOC) / dt_s;
        x_a_eq = @(SOC) 8.5e-3 + SOC.*(7.8e-1 - 8.5e-3);
        Ua_eq = @(x_a) 0.6379 + 0.5416.*exp(-305.5309.*x_a) + 0.044.*tanh(-1.*(x_a-0.1958)./0.1088) - 0.1978.*tanh((x_a-1.0571)./0.0854) - 0.6875.*tanh((x_a+0.0117)./0.0529) - 0.0175.*tanh((x_a-0.5692)./0.0875);
        Ua = Ua_eq(x_a_eq(soc));
        UaN = Ua ./ 0.123;
        
        % Cycling stressors:
        % dod
        dod = max(ElectroThermalResponse.SOC) - min(ElectroThermalResponse.SOC);
        % charge throughput
        amp_seconds = trapz(ElectroThermalResponse.Time_s, abs(ElectroThermalResponse.Current_A));
        amp_hours = amp_seconds / 3600;
        dEFC = amp_hours / (2*nominal_capacity);
        
        Stressors = table(dt, dEFC, TdegC, TdegK, TdegKN, soc, Ua, UaN, dod);
        end
    end
end

