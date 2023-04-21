function cycle = func_LifeMdl_StressStatistics(cycle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adds to "cycle" structure additional aging stress statistics.
%  Scalars:
%   cycle.Crate_rms,Crate_disrms,Crate_chgrms
%  Vectors: (using rainflow algorithm)
%   cycle.dsoc_i,ncyc_i,Crate_dis_i,Crate_chg_i,TavgK_i,socavg_i...
%
% Copyright Kandler Smith & Ying Shi, NREL 6/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%}
% dt
dt = cycle.t(end) - cycle.t(1);

% TdegK
TdegK = trapz(cycle.t, cycle.TdegC + 273.15)/dt;

% soc
soc = trapz(cycle.t, cycle.soc)/dt;

% Ua:
x_a_eq = @(SOC) 8.5e-3 + SOC.*(7.8e-1 - 8.5e-3);
Ua_eq = @(x_a) 0.6379 + 0.5416.*exp(-305.5309.*x_a) + 0.044.*tanh(-1.*(x_a-0.1958)./0.1088) - 0.1978.*tanh((x_a-1.0571)./0.0854) - 0.6875.*tanh((x_a+0.0117)./0.0529) - 0.0175.*tanh((x_a-0.5692)./0.0875);
Ua = Ua_eq(x_a_eq(soc));

% DOD
dod = max(cycle.soc) - min(cycle.soc);

% Crate
Crate = diff(cycle.soc)./diff(cycle.t .* 24);
% Crate is the average of the non-zero rates (Absolute value)
% Any rate smaller than a 1000 hour chg/dischg rate is
% considered to be negligible.
Crate = abs(Crate);
maskNonResting = Crate > 1e-3;
if any(maskNonResting)
    Crate = Crate(maskNonResting);
    difft = diff(cycle.t);
    difft = difft(maskNonResting);
    tNonResting = cumsum(difft);
    Crate = trapz(tNonResting, Crate)/(tNonResting(end)-tNonResting(1));
else
    Crate = 0;
end
% dEFC
dEFC = sum(abs(diff(cycle.soc)))/2;

% output the extracted stress statistics
cycle = struct();
cycle.dt =      dt;
cycle.dEFC =    dEFC;
cycle.TdegK =   TdegK;
cycle.soc =     soc;
cycle.Ua =      Ua;
cycle.dod =     dod;
cycle.Crate =   Crate;
end
