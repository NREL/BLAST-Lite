function xx = func_LifeMdl_UpdateStates(lifeMdl,cycle,xx0)
% Calculates capacity fade & resistance growth rates & integrates rates
%   forward in time to update states for a complex cycling profile
%   at a given state-of-life.
%
%   - Author: Kandler Smith, NREL, Copyright 2010
%   - Based on paper by EVS-24 paper by Smith,Markel,Pesaran, 2009
%   - Modified for new NCA model, KAS 6/2010
%   - Added FeP model based on A123 2.3Ah 26650 cell (KAS, 4/28/15)
%   - Adapted FeP model for Kokam 75Ah NMC cell (KAS,8/1/16)
%   - Adapted for Shell Battery Life Models (Paul Gasper, 2021)
%
% Inputs:
%   perfMdl - structure containing battery performance model parameters
%
%   lifeMdl - structure containing battery life model parameters
%           * see func_LoadLifeParameters.m for definition of life-parameters
%
%   cycle   - structure containing duty-cycle & temperature variables, e.g.
%           cycle.tsec  - time stamp array
%           cycle.soc   - state-of-charge array corresponding to time stamps (0=empty, 1=full charge)
%           cycle.TdegC - temperature array corresponding to time stamps (or use a
%                         single value to represent constant temperature) (degrees Celcius)
%           cycle.dsoc_i- array of delta-state-of-charge swings encountered by
%                         the battery during cycling calulated by Rainflow algorithm
%           cycle.ncycle_i-array of number of cycles corresponding to each
%                         delta-state-of-charge calculated by Rainflow algorithm
%
%   dtdays_life - timestep between previous and present battery
%            state-of-life (days)
%
%   xx0  - vector of life model states at previous state-of-life
%
% Outputs:
%   xx    - vector of life model states integrated to next state-of-life
%
%   dxxdt - vector of life model states time rate of change

% Rate equations are derived from trajectory equations in this manner:
%   Traj. eq: y = f(x)
%   Invert to solve for x: x_hat = g(y)
%   Rate eq: dy/dx = df(x_hat)/dx
%   The change of the state in each timestep is then dy/dx * dx (x = time
%       or energy throughput or whatever)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalized stressor values
cycle.TdegKN = cycle.TdegK ./ (273.15 + 35);
cycle.TdegC = cycle.TdegK - 273.15;
cycle.UaN = cycle.Ua ./ 0.123;

switch lifeMdl.model
    case 'NMC|Gr'
        % Calculate degradation rates
        q1 = lifeMdl.q1_0 .* exp(lifeMdl.q1_1 .* (1 ./ cycle.TdegKN)) .* exp(lifeMdl.q1_2 .* (cycle.UaN ./ cycle.TdegKN));
        q3 = lifeMdl.q3_0 .* exp(lifeMdl.q3_1 .* (1 ./ cycle.TdegKN)) .* exp(lifeMdl.q3_2 .* exp(cycle.dod.^2));
        q5 = lifeMdl.q5_0 + lifeMdl.q5_1 .* (cycle.TdegC - 55) .* cycle.dod;

        % Calculate incremental state changes
        dq_LLI_t = dynamicPowerLaw(xx0(1), cycle.dt, 2*q1, lifeMdl.q2);
        dq_LLI_EFC = dynamicPowerLaw(xx0(2), cycle.dEFC, q3, lifeMdl.q4);
        dq_LAM = dynamicSigmoid(xx0(3), cycle.dEFC, 1, 1/q5, lifeMdl.p_LAM);
        dr_LLI_t = dynamicPowerLaw(xx0(4), cycle.dt, lifeMdl.r1 * q1, lifeMdl.r2);
        dr_LLI_EFC = dynamicPowerLaw(xx0(5), cycle.dEFC, lifeMdl.r3 * q3, lifeMdl.r4);

        % Pack up
        dxx = [dq_LLI_t; dq_LLI_EFC; dq_LAM; dr_LLI_t; dr_LLI_EFC; xx0(6)];
        xx = xx0 + dxx;

    case 'LFP|Gr'
        % capacity
        p = lifeMdl.p;
        
        % Degradation rate equations
        kcal = @(S,p) abs(p.p1) * exp(p.p2 / S.TdegK) * exp(p.p3 * S.Ua);
        kcyc = @(S,p) (abs(p.p4) + abs(p.p5) * S.dod + abs(p.p6) * S.Crate);
        
        % Calendar loss
        dqLossCal = dynamicPowerLaw(xx0(1), cycle.dt, kcal(cycle,p), p.pcal);
        % Cycling loss
        dqLossCyc = dynamicPowerLaw(xx0(2), cycle.dEFC, kcyc(cycle,p), p.pcyc);
        
        % resistance
        p = lifeMdl.p_rdc;
        drGainCal = dynamicPowerLaw(xx0(3), cycle.dt, kcal(cycle,p), p.pcal);
        drGainCyc = dynamicPowerLaw(xx0(4), cycle.dEFC, kcyc(cycle,p), p.pcyc);
        
        % Pack up & accumulate state vector forward in time
        dxx = [dqLossCal; dqLossCyc; drGainCal; drGainCyc];
        xx = xx0 + dxx;
end

    function dy = dynamicPowerLaw(y0, dx, k, p)
        % DYNAMICPOWERLAW calculates the change of state for states modeled by a
        % power law equation, y = k*x^p. The output is instantaneous slope of dy
        % with respect to dx.
        if y0 == 0
            if dx == 0
                dydx = 0;
            else
                y0   = k * dx^p;
                dydx = y0 / dx;
            end
        else
            if dx == 0
                dydx = 0;
            else
                dydx = k * p * (y0 / k)^((p-1) / p);
            end
        end
        dy = dydx * dx;
    end

    function dy = dynamicSigmoid(y0, dx, y_inf, k, p)
        if y0 == 0
            if dx == 0
                dydx = 0;
            else
                dy = 2 * y_inf * (1/2 - 1 / (1 + exp((k * dx) ^ p)));
                dydx = dy / dx;
            end
        else
            if dx == 0
                dydx = 0;
            else
                x_inv = (1 / k) * ((log(-(2 * y_inf/(y0-y_inf)) - 1)) ^ (1 / p) );
                z = (k * x_inv) ^ p;
                dydx = (2 * y_inf * p * exp(z) * z) / (x_inv * (exp(z) + 1) ^ 2);
            end
        end
        dy = dydx * dx;
    end
end