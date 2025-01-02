function tickVarStr = get_tick_label(tickVar)
switch tickVar
    case 'Cycle_Index'
        tickVarStr = 'Cycle';
    case 'Charge_Throughput'
        tickVarStr = 'Charge throughput (Ah)';
    case 'Energy_Throughput'
        tickVarStr = 'Energy throughput (Wh)';
    case 'Equivalent_Full_Cycles'
        tickVarStr = 'Equiv. full cycles';
    case {'Time_s', 'tsecs_start', 'tsecs_end'}
        tickVarStr = 'Time (s)';
    case 'Time_d'
        tickVarStr = 'Time (days)';
    case 'Datenum_d'
        tickVarStr = 'Test datenum (days)';
    case 'datenum_d'
        tickVarStr = 'Relative test date (days)';
    case 'tsecs_cycle'
        tickVarStr = 'Cycle duration (s)';
    case 'Q_chg'
        tickVarStr = 'Charge capacity (Ah)';
    case 'q_chg'
        tickVarStr = 'Relative charge capacity';
    case 'Q_dis'
        tickVarStr = 'Discharge capacity (Ah)';
    case 'q_dis'
        tickVarStr = 'Relative discharge capacity';
    case 'CE'
        tickVarStr = 'Coulombic efficieny';
    case 'E_chg'
        tickVarStr = 'Charge energy (Wh)';
    case 'e_chg'
        tickVarStr = 'Relative charge energy';
    case 'E_dis'
        tickVarStr = 'Discharge energy (Wh)';
    case 'e_dis'
        tickVarStr = 'Relative discharge energy';
    case 'EE'
        tickVarStr = 'Energy efficiency';
    case 'DeltaV'
        tickVarStr = 'Voltage hysteresis, \DeltaV (V)';
    case 'deltaV'
        tickVarStr = 'Relative voltage hysteresis (V)';
    case 'V_min'
        tickVarStr = 'Minimum voltage (V)';
    case 'V_max'
        tickVarStr = 'Maximum voltage (V)';
    case 'V_avg'
        tickVarStr = 'Average voltage (V)';
    case 'I_min'
        tickVarStr = 'Minimum current (A)';
    case 'I_max'
        tickVarStr = 'Maximum current (A)';
    case 'I_avg'
        tickVarStr = 'Average current (A)';
    case 'P_min'
        tickVarStr = 'Minimum power (W)';
    case 'P_max'
        tickVarStr = 'Maximum power (W)';
    case 'P_avg'
        tickVarStr = 'Average power (W)';
    case 'T_min'
        tickVarStr = 'Minimum Temperature (^{\circ}C)';
    case 'T_max'
        tickVarStr = 'Maximum Temperature (^{\circ}C)';
    case 'T_avg'
        tickVarStr = 'Average Temperature (^{\circ}C)';
    otherwise
        tickVarStr = tickVar;
end
tickVarStr = string(tickVarStr);
end