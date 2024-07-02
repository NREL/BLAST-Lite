function yy = func_LifeMdl_UpdateOutput(lifeMdl, xx)
% Calculates capacity fade & resistance growth at a given state-of-life
%
%   - Author: Kandler Smith, NREL, Copyright 2010
%   - Based on paper by EVS-24 paper by Smith,Markel,Pesaran, 2009
%   - Modified for new NCA model, KAS 6/2010
%   - Added FeP model based on A123 2.3Ah 26650 cell (KAS, 4/28/15)
%   - Adapted FeP model for Kokam 75Ah NMC cell (KAS,8/1/16)
%
% Inputs:
%   perfMdl - structure containing battery performance model parameters
%
%   lifeMdl - structure containing battery life model parameters
%           * see func_LoadLifeParameters.m for definition of life-parameters
%
%   xx      - vector of life model states at present state-of-life
%
% Outputs:
%   yy      - vector of life model capacity fade and resistance growth

switch lifeMdl.model
    case 'NMC|Gr'
        q_LLI = 1 - xx(1) - xx(2);
        q_LLI_t = 1 - xx(1);
        q_LLI_EFC = 1 - xx(2);
        q_LAM = 1.01 - xx(3);
        q = min([q_LLI, q_LAM]);

        % Resistance
        r_LLI = 1 + xx(4) + xx(5);
        r_LLI_t = 1 + xx(4);
        r_LLI_EFC = 1 + xx(5);
        r_LAM = lifeMdl.r5 + lifeMdl.r6 * (1 / q_LAM);
        r = max([r_LLI, r_LAM]);

        % Assemble output
        yy = [q, q_LLI, q_LLI_t, q_LLI_EFC, q_LAM, r, r_LLI, r_LLI_t, r_LLI_EFC, r_LAM];
    case 'LFP|Gr'
        qLossCal = xx(1);
        % Scale cycling loss to assume cells aren't so long lived as the
        % Sony Murata cells are (10-15 thousand EFCs).
        scalingfactor = 4;
        qLossCyc = scalingfactor*xx(2);
        q = 1 - qLossCal - qLossCyc;
        rGainCal = xx(3);
        rGainCyc = scalingfactor*xx(4);
        r = 1 + rGainCal + rGainCyc;
        yy = [q qLossCal qLossCyc r rGainCal rGainCyc];
end
end