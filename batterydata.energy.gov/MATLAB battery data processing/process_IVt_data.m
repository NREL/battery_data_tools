function [ProcessedData, SummaryData] = process_IVt_data(RawData, options)
% PROCESS_IVT_DATA Process and summarize raw current-voltage-time data
%   Processes raw current-voltage-time data and then extracts summary
%   metrics. Processing first concatenates raw data tables if more than one
%   is provided, automatically handling the concatenation of tables with
%   different meaured variables (for example, tables with AC and DC
%   variables). Formatting of the raw data tables are checked to ensure
%   proper function of the procesing and plotting methods in the 'battery
%   data tools'. Processing accumulates charge and energy throughput and
%   calculates equivalent full cycles using a nominal capacity value,
%   either extracted from the 'NominalCycle' (default first cycle) or using
%   an input nominal capacity value. Processing also calculates 
%   differential capacity and differential voltage for any charge or
%   discharge. Differential analysis on low resolution voltage/current
%   measurement or noisy data is helped by data smoothing; charge and
%   discharge steps are interpolated into 1e4 evenly spaced segments on the
%   2-D curve and then the differential values are smoothed using an
%   Savtizky-Golay filter, with the filter window width set by the input
%   'DifferentialAnalysisSmoothingWindow'. If variables 
%   'Absolute_Charge_Throughput_Ah', 'Absolute_Energy_Throughput_Wh',
%   'Charge_Throughput_Ah', 'Energy_Throughput_Wh',
%   'Equivalent_Full_Cycles', 'Differential_Capacity_Ah_V', or
%   'Differential_Voltage_V_Ah' are already in the RawData table, they will
%   not be calculated (user input overrides automatic data processing).
%   Summary metrics are extracted and output as a new table with one row
%   per cycle. This fuction assumes that each 'cycle' has only one
%   contiguous charge and one contiguous discharge; cycles that violate
%   this assumption should not cause any errors, but the output summary
%   metrics or differential capacity/voltage calculations may be
%   nonsensical for cycles that violate these assumptions.
%
%   Required input:
%       RawData (table or cell): Raw current-voltage-time battery data.
%           Data formatting requirements are automatically checked by
%           calling the 'check_IVt_data_formatting' function. If a cell
%           array is input, the data tables in the array are sorted by
%           the Datenum_d variable and concatenated before processing.
%   Optional inputs (options):
%       NominalCycleIndex (int): cycle index to treat as 'nominal' for
%           calculating relative values of capacity and energy. Discharge
%           capacity of the nominal cycle is used for calculating
%           equivalent-full-cycles.
%       NominalCapacity (double): Nominal cell capacity in Ah, overriding
%           the nominal discharge capacity meausred from the 'NominalCycle'
%           specified above. Used for calculating equivalent-full-cycles.
%       CycleLabels (cell):
%       DifferentialAnalysisSmoothingWindow (int): window width for data
%           smoothing on differential capacity and voltage analysis. The
%           window length per unit is ~mV; the raw data is interpolated at
%           1e4 evenly spaced chords in the V-Ah basis before smoothing so
%           that both voltage and capacity vectors are sampled at approx.
%           the same sampling rate to ensure best quality on both dQdV and
%           dVdQ. Default is 250, channel voltage/current resolution and
%           sampling frequency will impact the necessary smoothing window
%           size. Larger windows for noisier data.
%       CurrentResolution (double): Threshold for determining zero current
%           for detecting charge and discharging steps. Default is 0.
%   Outputs:
%       ProcessedData (table): Raw current-voltage-time data table after
%           concatenation of input data tables (if there are more than 1),
%           accumulation of charge/energy throughput and equivalent full
%           cycle count, and calculation of differential capacity and
%           differential voltage.
%       SummaryData (table): table of cycle-by-cycle performance metrics.
%
%   Description of SummaryData table variables:
%       

arguments
    RawData {mustBeA(RawData, ["table","cell"])}
    options.NominalCycleIndex (1,1) double {mustBeInteger, mustBeNonnegative} = 1
    options.NominalCapacity (1,1) double = 0
    options.CycleLabels cell = {[]}
    options.DifferentialAnalysisSmoothingWindow (1,1) double {mustBeInteger, mustBeNonnegative} = 250
    options.CurrentResolution (1,1) double {mustBeNonnegative} = 0
end
nom_cyc_index = options.NominalCycleIndex;
diff_analysis_smoothing_window = options.DifferentialAnalysisSmoothingWindow;
cycle_labels = options.CycleLabels;
cycle_label_names = [cycle_labels{:,1}]';
cycle_label_rules = cycle_labels(:,2);
segment_labels = cycle_labels(:,3);

if isa(RawData, 'table')
    % Check that the raw dats is formatted correctly for processing
    check_IVt_data_formatting(RawData)
else %cell array of tables
    % Concat function checks formatting as well.
    RawData = concat_data_tables(RawData{:});
end
% Calculate charge-throughput, energy-throughput if they don't exist
if ~any(strcmp(RawData.Properties.VariableNames, 'Absolute_Charge_Throughput_Ah'))
    Ah_Chg = [RawData.Charge_Capacity_Ah(1); diff(RawData.Charge_Capacity_Ah)];
    Ah_Dis = [RawData.Discharge_Capacity_Ah(1); diff(RawData.Discharge_Capacity_Ah)];
    % Erase negative values at variable reset between cycles
    Ah_Chg(Ah_Chg < 0) = 0; 
    Ah_Dis(Ah_Dis < 0) = 0;
    % Accumulate
    RawData.Absolute_Charge_Throughput_Ah = cumsum(Ah_Chg + Ah_Dis);
end
if ~any(strcmp(RawData.Properties.VariableNames, 'Charge_Throughput_Ah'))
    Ah_Chg = [RawData.Charge_Capacity_Ah(1); diff(RawData.Charge_Capacity_Ah)];
    Ah_Dis = [RawData.Discharge_Capacity_Ah(1); diff(RawData.Discharge_Capacity_Ah)];
    % Erase negative values at variable reset between cycles
    Ah_Chg(Ah_Chg < 0) = 0; 
    Ah_Dis(Ah_Dis < 0) = 0;
    % Make Ah_Dis negative
    Ah_Dis = -1 .* Ah_Dis;
    % Accumulate
    RawData.Charge_Throughput_Ah = cumsum(Ah_Chg + Ah_Dis);
end
if ~any(strcmp(RawData.Properties.VariableNames, 'Absolute_Energy_Throughput_Wh'))
    Wh_Chg = [RawData.Charge_Energy_Wh(1); diff(RawData.Charge_Energy_Wh)];
    Wh_Dis = [RawData.Discharge_Energy_Wh(1); diff(RawData.Discharge_Energy_Wh)];
    % Erase negative values at variable reset between cycles
    Wh_Chg(Wh_Chg < 0) = 0; 
    Wh_Dis(Wh_Dis < 0) = 0;
    % Accumulate
    RawData.Absolute_Energy_Throughput_Wh = cumsum(Wh_Chg + Wh_Dis);
end
if ~any(strcmp(RawData.Properties.VariableNames, 'Energy_Throughput_Wh'))
    Wh_Chg = [RawData.Charge_Energy_Wh(1); diff(RawData.Charge_Energy_Wh)];
    Wh_Dis = [RawData.Discharge_Energy_Wh(1); diff(RawData.Discharge_Energy_Wh)];
    % Erase negative values at variable reset between cycles
    Wh_Chg(Wh_Chg < 0) = 0; 
    Wh_Dis(Wh_Dis < 0) = 0;
    % Make Wh_Dis negative
    Wh_Dis = -1 .* Wh_Dis;
    % Accumulate
    RawData.Energy_Throughput_Wh = cumsum(Wh_Chg + Wh_Dis);
end

% Extract metrics from each cycle
% Calculate dQdV and dVdQ for each cycle
cycle = unique(RawData.Cycle_Index);
n = length(cycle);
% Variables to pull from each cycle
tsecs_start = zeros(n,1);
tsecs_end = zeros(n,1);
tsecs_cycle = zeros(n,1);
Datenum_d = zeros(n,1);
Cycle_Index = cycle; % we already know this one
Absolute_Charge_Throughput_Ah = zeros(n,1);
Absolute_Energy_Throughput_Wh = zeros(n,1);
Charge_Throughput_Ah = zeros(n,1);
Energy_Throughput_Wh = zeros(n,1);
Q_dis = zeros(n,1);
Q_chg = zeros(n,1);
CE = zeros(n,1);
E_dis = zeros(n,1);
E_chg = zeros(n,1);
EE = zeros(n,1);
T_min = zeros(n,1);
T_max = zeros(n,1);
T_avg = zeros(n,1);
P_min = zeros(n,1);
P_max = zeros(n,1);
P_avg = zeros(n,1);
V_min = zeros(n,1);
V_max = zeros(n,1);
V_avg = zeros(n,1);
V_chg_avg = zeros(n,1); % voltage at 50% SOC
V_dis_avg = zeros(n,1); % voltage at 50% SOC
DeltaV = zeros(n,1); % difference of last two variables
I_max = zeros(n,1);
I_min = zeros(n,1);
I_avg = zeros(n,1);
Cycle_Label = strings(n,1);

RawData.Cycle_Label(:) = "";
RawData.Segment_Label(:) = "";

for i = 1:length(cycle)
    % Get data from just this cycle
    this_cycle = cycle(i);
    mask_this_cycle = RawData.Cycle_Index == this_cycle;
    data_this_cycle = RawData(mask_this_cycle, :);
    
    % Charge and discharge data
    mask_chg = data_this_cycle.Current_A > options.CurrentResolution;
    data_chg = data_this_cycle(mask_chg, :);
    mask_dis = data_this_cycle.Current_A < -options.CurrentResolution;
    data_dis = data_this_cycle(mask_dis, :);

    % Extract metrics
    tsecs_start(i) = data_this_cycle.Time_s(1);
    tsecs_end(i) = data_this_cycle.Time_s(end);
    tsecs_cycle(i) = tsecs_end(i) - tsecs_start(i);
    Datenum_d(i) = data_this_cycle.Datenum_d(1);
    % Charge:
    if any(mask_chg) && length(mask_chg(mask_chg)) > 1
        Q_chg(i) = max(data_this_cycle.Charge_Capacity_Ah);
        E_chg(i) = max(data_this_cycle.Charge_Energy_Wh);
        % V average calculation
        data_chg.soc = data_chg.Charge_Capacity_Ah ./ Q_chg(i);
        % Can only interpolate b/w unique values
        [~, idx, ~] = unique(data_chg.soc, 'last');
        try
            V_chg_avg(i) = interp1(data_chg.soc(idx), data_chg.Voltage_V(idx), 0.5);
        catch
            V_chg_avg(1) = NaN;
        end
    else
        Q_chg(i) = NaN;
        E_chg(i) = NaN;
        V_chg_avg(i) = NaN;
    end
    % Discharge:
    if any(mask_dis) && length(mask_dis(mask_dis)) > 1
        Q_dis(i) = max(data_this_cycle.Discharge_Capacity_Ah);
        E_dis(i) = max(data_this_cycle.Discharge_Energy_Wh);
        % V average calculation 
        data_dis.soc = data_dis.Discharge_Capacity_Ah ./ Q_dis(i);
        [~, idx, ~] = unique(data_dis.soc, 'last');
        try
            V_dis_avg(i) = interp1(data_dis.soc(idx), data_dis.Voltage_V(idx), 0.5);
        catch
            V_dis_avg(i) = NaN;
        end
    else
        Q_dis(i) = NaN;
        E_dis(i) = NaN;
        V_dis_avg(i) = NaN;
    end

    % Efficiency metrics:
    CE(i) = Q_dis(i) / Q_chg(i);
    EE(i) = E_dis(i) / E_chg(i);
    DeltaV(i) = V_chg_avg(i) - V_dis_avg(i);

    % Other metrics that are unlikely to cause errors:
    Absolute_Charge_Throughput_Ah(i) = max(data_this_cycle.Absolute_Charge_Throughput_Ah);
    Absolute_Energy_Throughput_Wh(i) = max(data_this_cycle.Absolute_Energy_Throughput_Wh);
    Charge_Throughput_Ah(i) = data_this_cycle.Charge_Throughput_Ah(end);
    Energy_Throughput_Wh(i) = data_this_cycle.Energy_Throughput_Wh(end);
    
    % Time-based average, handles unevenly sampled data
    V_min(i) = min(data_this_cycle.Voltage_V);
    V_max(i) = max(data_this_cycle.Voltage_V);
    V_avg(i) = trapz(data_this_cycle.Time_s, data_this_cycle.Voltage_V) ./ tsecs_cycle(i);
    I_min(i) = min(data_this_cycle.Current_A);
    I_max(i) = max(data_this_cycle.Current_A);
    I_avg(i) = trapz(data_this_cycle.Time_s, data_this_cycle.Current_A) ./ tsecs_cycle(i);
    P_min(i) = min(data_this_cycle.Power_W);
    P_max(i) = max(data_this_cycle.Power_W);
    P_avg(i) = trapz(data_this_cycle.Time_s, data_this_cycle.Power_W) ./ tsecs_cycle(i);

    % If temperature exists, calculate temperature variables
    if any(strcmp(RawData.Properties.VariableNames, 'Cell_Temperature_C'))
        T_min(i) = min(data_this_cycle.Cell_Temperature_C);
        T_max(i) = max(data_this_cycle.Cell_Temperature_C);
        T_avg(i) = trapz(data_this_cycle.Time_s, data_this_cycle.Cell_Temperature_C) ./ tsecs_cycle(i);
    else
        T_min(i) = NaN;
        T_max(i) = NaN;
        T_avg(i) = NaN;
    end

    % Check cycle name
    if length(cycle_label_rules) == 1 && isempty(cycle_label_rules{1})
        Cycle_Label(i) = cycle_label_names;
    else
        isLabel = false(length(cycle_label_names), 1);
        for iLabel = 1:length(cycle_label_rules)
            check_label = cycle_label_rules{iLabel};
            isLabel(iLabel) = all(check_label(data_this_cycle));
        end
        idxLabel = find(isLabel);
        if length(idxLabel) > 1
            error("Check cycle label rules, only one can be valid at a time, the labels {'" + join(cycle_label_names(idxLabel), "', '") + "'} were all valid.")
        elseif isempty(idxLabel)
            Cycle_Label(i) = "";
        else
            Cycle_Label(i) = cycle_label_names(idxLabel);
            % Label the measurement data table
            RawData.Cycle_Label(mask_this_cycle) = cycle_label_names(idxLabel);
            % Check the segment labels for this cycle
            segment_labels_this_cycle = segment_labels{idxLabel};
            if ~isempty(segment_labels_this_cycle)
                segment_label_names = [segment_labels_this_cycle{:,1}]';
                segment_label_rules = segment_labels_this_cycle(:,2);
                for iSegmentLabel = 1:length(segment_label_names)
                    check_segment_label = segment_label_rules{iSegmentLabel};
                    mask_this_label = check_segment_label(RawData);
                    RawData.Segment_Label(mask_this_cycle & mask_this_label) = segment_label_names(iSegmentLabel);
                end
            end
        end
    end

    % Calculate dQdV and dVdQ and insert into the raw_data table:
    if any(mask_chg) && length(mask_chg(mask_chg)) > 3
        % Charge:        
        % Cut out constant-voltage phase, we can't process this data using
        % dQdV or dVdQ.
        dV = [0; diff(data_chg.Voltage_V)];
        idx_cc = find(dV, 1, "last");
        if ~isempty(idx_cc)
            data_analysis = data_chg(1:idx_cc, :);
        else
            data_analysis = data_chg;
        end
        % Interpolate the V-Ah arc at evenly spaced intervals ~100 microvolts
        V_range = range(data_analysis.Voltage_V);
        n = ceil(V_range*1e4);
        [~, idx] = unique(data_analysis.Charge_Capacity_Ah);
        if length(idx) > 2
            [xx, dxxdt] = interparc(n, data_analysis.Charge_Capacity_Ah(idx), data_analysis.Voltage_V(idx), 'linear');
            % Smooth the diffs and the use chain rule to calculate differences
            dxxdt = smoothdata(dxxdt, 'sgolay', diff_analysis_smoothing_window);
            dQdV = dxxdt(:,1) ./ dxxdt(:,2);
            dVdQ = dxxdt(:,2) ./ dxxdt(:,1);
            % interparc does interpolation over evenly spaced segment arcs.
            % Calculate the segment distances from the original vector to
            % return to the original basis.
            segmentlen = sqrt(sum(diff([data_analysis.Charge_Capacity_Ah, data_analysis.Voltage_V],[],1).^2,2));
            segmentlen = segmentlen/sum(segmentlen);
            segmentdists = [0;cumsum(segmentlen)];
            if length(dQdV) > 1 % Occasionally there's only 1 segment and this breaks
                dQdV = interp1(linspace(0,1,n), dQdV, segmentdists); dQdV(idx_cc+1:height(data_chg)) = NaN;
                dQdV = filloutliers(dQdV, 'spline', 'movmedian', ceil(length(dQdV)/10));
                dVdQ = interp1(linspace(0,1,n), dVdQ, segmentdists); dVdQ(idx_cc+1:height(data_chg)) = NaN;
                dVdQ = filloutliers(dVdQ, 'spline', 'movmedian', ceil(length(dQdV)/10));
            end
    
            % Store at original measurement points
            data_this_cycle.dQdV(mask_chg) = dQdV;
            data_this_cycle.dVdQ(mask_chg) = dVdQ;
        end
    end
    if any(mask_dis) && length(mask_dis(mask_dis)) > 3
        % Discharge:
        % Cut out constant-voltage phase, we can't process this data using
        % dQdV or dVdQ.
        dV = [0; diff(data_dis.Voltage_V)];
        idx_cc = find(dV, 1, "last");
        if ~isempty(idx_cc)
            data_analysis = data_dis(1:idx_cc, :);
        else
            data_analysis = data_dis;
        end
        % Interpolate the V-Ah arc at evenly spaced intervals ~100 microvolts
        V_range = range(data_analysis.Voltage_V);
        n = ceil(V_range*1e4);
        [~, idx] = unique(data_analysis.Discharge_Capacity_Ah);
        if length(idx) > 2
            [xx, dxxdt] = interparc(n, data_analysis.Discharge_Capacity_Ah(idx), data_analysis.Voltage_V(idx), 'linear');
            % Smooth the diffs and the use chain rule to calculate differences
            dxxdt = smoothdata(dxxdt, 'sgolay', diff_analysis_smoothing_window);
            dQdV = dxxdt(:,1) ./ dxxdt(:,2);
            dVdQ = dxxdt(:,2) ./ dxxdt(:,1);
                
            % interparc does interpolation over evenly spaced segment arcs.
            % Calculate the segment distances from the original vector to
            % return to the original basis.
            segmentlen = sqrt(sum(diff([data_analysis.Discharge_Capacity_Ah, data_analysis.Voltage_V],[],1).^2,2));
            segmentlen = segmentlen/sum(segmentlen);
            segmentdists = [0;cumsum(segmentlen)];
            if length(dQdV) > 1 % Occasionally there's only 1 segment and this breaks
                dQdV = interp1(linspace(0,1,n), dQdV, segmentdists); dQdV(idx_cc+1:height(data_dis)) = NaN;
                dQdV = filloutliers(dQdV, 'spline', 'movmedian', ceil(length(dQdV)/10));
                dVdQ = interp1(linspace(0,1,n), dVdQ, segmentdists); dVdQ(idx_cc+1:height(data_dis)) = NaN;
                dVdQ = filloutliers(dVdQ, 'spline', 'movmedian', ceil(length(dQdV)/10));
            end
    
            % Store at original measurement points
            data_this_cycle.dQdV(mask_dis) = dQdV;
            data_this_cycle.dVdQ(mask_dis) = dVdQ;
        end
    end
    % Insert into data_raw
    data_this_cycle.dQdV(~(mask_chg | mask_dis)) = NaN;
    data_this_cycle.dVdQ(~(mask_chg | mask_dis)) = NaN;
    RawData.Differential_Capacity_Ah_V(mask_this_cycle) = data_this_cycle.dQdV;
    RawData.Differential_Voltage_V_Ah(mask_this_cycle) = data_this_cycle.dVdQ; 
end

% Rename Raw to Processed
ProcessedData = RawData; clearvars RawData;

% Create SummaryData table
% EFC from nominal capacity
if options.NominalCapacity == 0
    Equivalent_Full_Cycles = Absolute_Charge_Throughput_Ah ./ (2*Q_dis(nom_cyc_index));
    Equivalent_Charge_Slippage = Charge_Throughput_Ah ./ (2*Q_dis(nom_cyc_index));
else
    Equivalent_Full_Cycles = Absolute_Charge_Throughput_Ah ./ (2*options.NominalCapacity);
    Equivalent_Charge_Slippage = Charge_Throughput_Ah ./ (2*options.NominalCapacity);
end
% Relative metrics
q_chg = Q_chg ./ Q_chg(nom_cyc_index);
q_dis = Q_dis ./ Q_dis(nom_cyc_index);
e_chg = E_chg ./ E_chg(nom_cyc_index);
e_dis = E_dis ./ E_dis(nom_cyc_index);
deltaV = DeltaV ./ DeltaV(nom_cyc_index);

% Copy tsecs_start with a nicer name and convert to days
Time_s = tsecs_start;
Time_d = Time_s ./ (24*3600);

% Make a relative datenum variable
datenum_d = Datenum_d - Datenum_d(1);

% Output table
SummaryData = table(Cycle_Index, Cycle_Label,...
    Time_s, Time_d, Datenum_d, datenum_d,...
    Absolute_Charge_Throughput_Ah, Absolute_Energy_Throughput_Wh, Equivalent_Full_Cycles,...
    Charge_Throughput_Ah, Energy_Throughput_Wh, Equivalent_Charge_Slippage,...
    tsecs_start, tsecs_end, tsecs_cycle,...
    Q_chg, q_chg, Q_dis, q_dis, CE,...
    E_chg, e_chg, E_dis, e_dis, EE,...
    DeltaV, deltaV,...
    V_min, V_max, V_avg,...
    I_min, I_max, I_avg,...
    P_min, P_max, P_avg,...
    T_min, T_max, T_avg);
end