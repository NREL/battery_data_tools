function check_IVt_data_formatting(RawData)
% CHECK_IVT_DATA_FORMATTING Checks raw current-voltage-time data formatting
%   Checks the formatting of the raw data table for some common formatting
%   errors such as missing variables, nonnegativity, and monotonicity.
arguments
    RawData table % Raw current-voltage-time battery data
end

% Check for required variables
required_variables = ["Datenum_d", "Time_s", "Cycle_Index", "Step", "Current_A",...
    "Voltage_V", "Power_W", "Charge_Capacity_Ah", "Discharge_Capacity_Ah",...
    "Charge_Energy_Wh", "Discharge_Energy_Wh"];
common_variables = intersect(required_variables, RawData.Properties.VariableNames);
missing_variables = setdiff(required_variables, common_variables);
if length(missing_variables) > 0
    error("The following variables must be added to the raw data table: " + join(missing_variables, ", "))
end
% Check impedance variables formatting, if they are there
if any(contains(RawData.Properties.VariableNames, 'Frequency', 'IgnoreCase', true))
    required_eis_vars = ["Frequency_Hz", "Z_Real_Ohm", "Z_Imag_Ohm", "Z_Mag_Ohm", "Z_Phase_Degree"];
    if ~all(contains(required_eis_vars, RawData.Properties.VariableNames))
        error("Check formatting of impedance variable names, should be Frequency_Hz, Z_Real_Ohm, Z_Imag_Ohm, Z_Mag_Ohm, and Z_Phase_Degree")
    end
end

% Charge/Discharge Ah and Wh values are nonnegative
if any(RawData.Charge_Capacity_Ah < 0)
    error("Charge_Capacity_Ah must be nonnegative.")
end
if any(RawData.Charge_Energy_Wh < 0)
    error("Charge_Energy_Wh must be nonnegative.")
end
if any(RawData.Discharge_Capacity_Ah < 0)
    error("Discharge_Capacity_Ah must be nonnegative.")
end
if any(RawData.Discharge_Energy_Wh < 0)
    error("Discharge_Energy_Wh must be nonnegative.")
end

% Cycle indices are monotonically increasing, always positive
if any(RawData.Cycle_Index < 0)
    error("Cycle_Index must be nonnegative")
end
if ~issorted(RawData.Cycle_Index, 'ascend')
    error("Cycle_Index must be monotonically increasing.")
end

% Time_s is monotonically increasing, always positive
if any(RawData.Time_s < 0)
    error("Time_s must be nonnegative")
end
if ~issorted(RawData.Time_s, 'ascend')
    error("Time_s must be monotonically increasing.")
end

% Datenum is monotonically increasing, always positive
if any(RawData.Datenum_d < 0)
    error("Datenum_d must be nonnegative")
end
if ~issorted(RawData.Datenum_d, 'ascend')
    % can be identical datenums, depending on time sampling rate and
    % datetime resolution, which will flag as not sorted
    if min(diff(RawData.Datenum_d)) < 0
        error("Datenum_d must be monotonically increasing.")
    end
end

% Charge_Throughput_Ah is monotonically increasing, always positive
if any(strcmp(RawData.Properties.VariableNames, 'Charge_Throughput_Ah'))
    if any(RawData.Charge_Throughput_Ah < 0)
        error("Charge_Throughput_Ah must be nonnegative")
    end
    if ~issorted(RawData.Charge_Throughput_Ah, 'ascend')
        error("Charge_Throughput_Ah must be monotonically increasing.")
    end
end

% Energy_Throughput_Wh is monotonically increasing, always positive
if any(strcmp(RawData.Properties.VariableNames, 'Energy_Throughput_Wh'))
    if any(RawData.Energy_Throughput_Wh < 0)
        error("Energy_Throughput_Wh must be nonnegative")
    end
    if ~issorted(RawData.Energy_Throughput_Wh, 'ascend')
        error("Energy_Throughput_Wh must be monotonically increasing.")
    end
end

% Equivalent Full Cycles is monotonically increasing, always positive
if any(strcmp(RawData.Properties.VariableNames, 'Equivalent_Full_Cycles'))
    if any(RawData.Equivalent_Full_Cycles < 0)
        error("Equivalent_Full_Cycles must be nonnegative")
    end
    if ~issorted(RawData.Equivalent_Full_Cycles, 'ascend')
        error("Equivalent_Full_Cycles must be monotonically increasing.")
    end
end
end