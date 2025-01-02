function varargout = plot_measurements(ProcessedData, SummaryData, options)
% PLOT_MEASUREMENTS Plots data from various battery measurement types.
%   Plots raw current-voltage-time battery data for various types of
%   measurements, either plotting all recorded measurements or from up to 5
%   specific cycle numbers. Can plot to an input axis handle, and
%   optionally outputs the plotted axis.
%   Required inputs:
%       ProcessedData (table): table of raw current-voltage-time timeseries
%           data. Must have been processed using the 'process_IVt_data'
%           function before using plot_measurements.
%       SummaryData (table): table of summarized battery performance
%           metrics, calculated using process_IVt_data. Any SummaryData 
%           field may be used to denote each measurement as set by the
%           optional input 'ColorVariable'. Cycles to plot can be set by
%           any Cycle_Label in the SummaryData.
%       options (name-value input):
%           Measurement (text): type of measurement to plot:
%               'Voltage-capacity', 'Voltage-time', 'Charge slippage'....
%               'Polarization-time', 'Differential capacity', 'Differential voltage', ...
%               'EIS Nyquist', 'EIS Bode Real-Imaginary', 'EIS Bode Mag-Phase'
%   Optional inputs:
%       options (name-value inputs):
%           Direction (text): 'Charge', 'Discharge', or 'Both', used for
%               plotting measurements with clear charge or discharge
%               directions (default 'Both').
%           Cycle: 'All', a Cycle_Label from SummaryData, or a vector 
%               of cycle numbers (default 'All').
%           Segment: a Segment_Label to plot for each 'Cycle' being plotted
%           LineStyle (text): MATLAB linestyle, as input to plot() (default '-')
%           LineStyle_Secondary (text): Secondary line style for secondary
%               y-axis, used for EIS Bode plots (default '--').
%           LineWidth (double): Line width (default 1)
%           AxisHandle (Axes): An axis handle to plot on. If empty, plots
%               on a new figure like default MATLAB plot behavior.
%           AxisLimits (double): Axis limits, [x_min x_max y_min y_max]
%           LegendStyle (text): Sets behavior for denoting each cycle as
%               'Legend' (in a legend), 'Colorbar' (using a colorbar), or
%               'Off' (none) (default 'Colorbar')
%           LegendLocation (text): Axis location for the legend/colorbar.
%               In addition to normal MATLAB inputs, can also be
%               'tile_...', which places the legend on the specified
%               tiledlayout grid location (default 'eastoutside').
%           ColorVariable (text): A variable from SummaryData table to name
%               legend entries or put as colorbar ticks (default 'Cycle').
%           Colormap (function_handle): function_handle that generates a
%               arbitrarily long set of RGB colors (default @(n) cool(n)).
%   Outputs:
%       ax (Axes): the axis handle that was plotted on
arguments
    ProcessedData table
    SummaryData table
    options.Measurement {mustBeMember(options.Measurement, {...
        'Voltage-capacity', 'Voltage-time','Polarization-time',...
        'Charge slippage', 'Differential capacity', 'Differential voltage',...
        'EIS Nyquist', 'EIS Bode Real-Imaginary', 'EIS Bode Mag-Phase'})}
    options.Direction {mustBeMember(options.Direction, {'Charge', 'Discharge', 'Both'})} = 'Both'
    % Optional inputs
    options.Cycle = 'All'; % 'All', a Cycle_Label from SummaryData, or a vector of cycle numbers.
    options.Segment = ''; % A Segment_Label from the 'Cycle' being plotted
    % Line options
    options.LineStyle {mustBeMember(options.LineStyle, {'-','--',':','.-'})} = '-'
    options.LineStyle_Secondary {mustBeMember(options.LineStyle_Secondary, {'-','--',':','.-'})} = '--'
    options.LineWidth (1,1) double  = 1
    % Plotting decorations
    options.AxisHandle = [] % An axis handle to plot on
    options.AxisLimits = [] % Axis limits, [x_min x_max y_min y_max]
    % Legend or colorbar location to show cycle number of the plotted meausrement,
    % can also be in the format 'tile_LOC' where LOC is east, north, west, or south
    % if this axis is on a tiledlayout
    options.LegendStyle {mustBeMember(options.LegendStyle, {'Legend','Colorbar','Off'})} = 'Colorbar'
    options.LegendLocation {mustBeMember(options.LegendLocation, {'east','west','north','south',...
        'eastoutside','westoutside','northoutside','southoutside',...
        'tile_east','tile_west','tile_north','tile_south'})} = 'eastoutside' 
    options.ColorVariable {mustBeMember(options.ColorVariable, {'Cycle_Index','Absolute_Charge_Throughput_Ah','Absolute_Energy_Throughput_Wh',...
        'Equivalent_Full_Cycles','Time_s','Time_d','Datenum_d','datenum_d','tsecs_start','tsecs_end','tsecs_cycle','Q_chg','q_chg','Q_dis','q_dis','CE',...
        'E_chg','e_chg','E_dis','e_dis','EE','DeltaV','deltaV','V_min','V_max','V_avg','I_min','I_max','I_avg',...
        'P_min','P_max','P_avg','T_min','T_max','T_avg'})} = 'Cycle_Index'
    options.Colormap = @(n) cool(n) % function_handle to a valid colormap
    options.CurrentResolution (1,1) double {mustBePositive} = 1e-6
end
% Input checking for ColorVariable
if isstring(options.Cycle) || ischar(options.Cycle)
    if strcmp(options.Cycle, 'All')
        Cycles = unique(ProcessedData.Cycle_Index);
    else
        validCycleLabels = unique(SummaryData.Cycle_Label);
        if ~any(strcmp(options.Cycle, validCycleLabels))
            error("Must specify one Cycle_Label out of: " + join(validCycleLabels, ", ") + ".")
        else
            Cycles = SummaryData.Cycle_Index(strcmp(SummaryData.Cycle_Label, options.Cycle), :);
        end
    end
else
    assert(isvector(options.Cycle) && isnumeric(options.Cycle), 'The input options.Cycle must be a numeric vector.');
    Cycles = options.Cycle;
end
if iscolumn(Cycles)
    Cycles = Cycles';
end
n = length(Cycles);
gobj = gobjects(n,1);
if isempty(options.AxisHandle)
    figure; hold on; box on; grid on; ax = gca;
else
    ax = options.AxisHandle; hold on; box on; grid on;
end
% If n is small and this is an impedance bode plot, use the legend
if n < 10 && contains(options.Measurement, 'EIS Bode')
    options.LegendStyle = 'Legend';
end
% Line styles
colors = options.Colormap(n);
ls = options.LineStyle;
lw = options.LineWidth;
tick_variable = options.ColorVariable;
tickVarStr = get_tick_label(tick_variable);
% Plot lines
for i = 1:n
    this_cycle = Cycles(i);
    mask_this_cycle = ProcessedData.Cycle_Index == this_cycle;
    Data_this_cycle = ProcessedData(mask_this_cycle, :);
    if ~isempty(options.Segment)
        mask_this_segment = strcmp(Data_this_cycle.Segment_Label, options.Segment);
        idx_segments = contiguous(double(mask_this_segment));
        if size(idx_segments, 1) > 1
            idx_segments = idx_segments{2,2};
            segment_starts = idx_segments(:,1);
            segment_stops = idx_segments(:,2);
        else
            segment_starts = 1;
            segment_stops = height(Data_this_cycle);
        end
        hold(ax, 'on')
    else
        segment_starts = 1;
        segment_stops = height(Data_this_cycle);
    end
    for i_segment = 1:length(segment_starts)
        Data_this_segment = Data_this_cycle(segment_starts(i_segment):segment_stops(i_segment), :);
        Data_Discharge = Data_this_segment(Data_this_segment.Current_A < -options.CurrentResolution, :);
        Data_Charge = Data_this_segment(Data_this_segment.Current_A > options.CurrentResolution, :);
        % Line label
        if ~strcmp(options.ColorVariable, 'Cycle_Index')
            this_colorVar = SummaryData{i, tick_variable};
        else
            this_colorVar = this_cycle;
        end
        switch options.Measurement
            case 'Voltage-capacity'
                switch options.Direction
                    case 'Discharge'
                        if ~isempty(Data_Discharge)
                            gobj(i) = plot(ax, Data_Discharge.Discharge_Capacity_Ah, Data_Discharge.Voltage_V,...
                                ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar, 'LineWidth', lw);
                        end
                    case 'Charge'
                        if ~isempty(Data_Charge)
                            gobj(i) = plot(ax, Data_Charge.Charge_Capacity_Ah, Data_Charge.Voltage_V,...
                                ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar);
                        end
                    case 'Both'
                        if ~isempty(Data_Charge)
                            Q_chg = max(Data_Charge.Charge_Capacity_Ah);
                            plot(ax, Data_Charge.Charge_Capacity_Ah, Data_Charge.Voltage_V,...
                                ls, 'Color', colors(i,:), 'LineWidth', lw);
                            is_empty_charge = false;
                        else
                            is_empty_charge = true;
                        end
                        if ~isempty(Data_Discharge)
                            if ~is_empty_charge
                                gobj(i) = plot(ax, Q_chg - Data_Discharge.Discharge_Capacity_Ah, Data_Discharge.Voltage_V,...
                                    ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar, 'LineWidth', lw);
                            else
                                gobj(i) = plot(ax, Data_Discharge.Discharge_Capacity_Ah, Data_Discharge.Voltage_V,...
                                    ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar, 'LineWidth', lw);
                            end
                        end
                end
            case 'Voltage-time'
                gobj(i) = plot(ax, Data_this_segment.Time_s - Data_this_segment.Time_s(1), Data_this_segment.Voltage_V,...
                    ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar, 'LineWidth', lw);
            case 'Polarization-time'
                gobj(i) = plot(ax, Data_this_segment.Time_s - Data_this_segment.Time_s(1), Data_this_segment.Voltage_V - Data_this_segment.Voltage_V(1),...
                    ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar, 'LineWidth', lw);
            case 'Charge slippage'
                gobj(i) = plot(ax, Data_this_segment.Charge_Throughput_Ah, Data_this_segment.Voltage_V,...
                    ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar, 'LineWidth', lw);
            case 'Differential capacity'
                switch options.Direction
                    case 'Discharge'
                        if ~isempty(Data_Discharge)
                            gobj(i) = plot(ax, Data_Discharge.Voltage_V, Data_Discharge.Differential_Capacity_Ah_V,...
                                ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar, 'LineWidth', lw);
                        end
                    case 'Charge'
                        if ~isempty(Data_Charge)
                            gobj(i) = plot(ax, Data_Charge.Voltage_V, Data_Charge.Differential_Capacity_Ah_V,...
                                ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar, 'LineWidth', lw);
                        end
                    case 'Both'
                        if ~isempty(Data_Discharge)
                            gobj(i) = plot(ax, Data_Discharge.Voltage_V, Data_Discharge.Differential_Capacity_Ah_V,...
                                ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar, 'LineWidth', lw);
                        end
                        if ~isempty(Data_Charge)
                            plot(ax, Data_Charge.Voltage_V, Data_Charge.Differential_Capacity_Ah_V,...
                                ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar, 'LineWidth', lw);
                        end
                end
            case 'Differential voltage'
                switch options.Direction
                    case 'Discharge'
                        if ~isempty(Data_Discharge)
                            gobj(i) = plot(ax, Data_Discharge.Discharge_Capacity_Ah, Data_Discharge.Differential_Voltage_V_Ah,...
                                ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar, 'LineWidth', lw);
                        end
                    case 'Charge'
                        if ~isempty(Data_Charge)
                            gobj(i) = plot(ax, Data_Charge.Charge_Capacity_Ah, Data_Charge.Differential_Voltage_V_Ah,...
                                ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar, 'LineWidth', lw);
                        end
                    case 'Both'
                        if ~isempty(Data_Discharge)
                            gobj(i) = plot(ax, Data_Discharge.Discharge_Capacity_Ah, Data_Discharge.Differential_Voltage_V_Ah,...
                                ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar, 'LineWidth', lw);
                        end
                        if ~isempty(Data_Charge)
                            plot(ax, Data_Charge.Charge_Capacity_Ah, Data_Charge.Differential_Voltage_V_Ah,...
                                ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar, 'LineWidth', lw);
                        end
                end
            case 'EIS Nyquist'
                gobj(i) = plot(ax, Data_this_segment.Z_Real_Ohm, Data_this_segment.Z_Imag_Ohm, ...
                    ls, 'Color', colors(i,:), 'DisplayName', tickVarStr + " " + this_colorVar, 'LineWidth', lw);
                % Markers on frequency decades
                mask = ~isnan(Data_this_segment.Frequency_Hz); % avoid NAN rows if EIS is combined with a DC measurement like a rest
                logfreq = round(log10(Data_this_segment.Frequency_Hz(mask)));
                freq_decades = logspace(min(logfreq), max(logfreq), range(logfreq)+1);
                Z_Real_decades = interp1(Data_this_segment.Frequency_Hz(mask), Data_this_segment.Z_Real_Ohm(mask), freq_decades);
                Z_Imag_decades = interp1(Data_this_segment.Frequency_Hz(mask), Data_this_segment.Z_Imag_Ohm(mask), freq_decades);
                plot(ax, Z_Real_decades, Z_Imag_decades, 'd', 'Color', colors(i,:), 'LineWidth', lw)
                % Axis formatting for Nyquist plots
                set(ax, 'YDir', 'reverse')
                axis square            
            case 'EIS Bode Real-Imaginary'
                gobj(i) = plot(ax, Data_this_segment.Frequency_Hz, Data_this_segment.Z_Real_Ohm, ...
                    ls, 'Color', colors(i,:), 'LineWidth', lw, ...
                    'DisplayName', tickVarStr + " " + this_colorVar + " (Real: " + ls + ", Imag: " + options.LineStyle_Secondary + ")");
                set(ax, 'XScale', 'log')
                colororder({'k','k'})
                yyaxis right
                plot(Data_this_segment.Frequency_Hz, Data_this_segment.Z_Imag_Ohm, ...
                    options.LineStyle_Secondary, 'Color', colors(i,:), 'LineWidth', lw);
                set(gca, 'YDir', 'reverse')
                yyaxis left
            case 'EIS Bode Mag-Phase'
                gobj(i) = plot(ax, Data_this_segment.Frequency_Hz, Data_this_segment.Z_Mag_Ohm, ...
                    ls, 'Color', colors(i,:), 'LineWidth', lw, ...
                    'DisplayName', tickVarStr + " " + this_colorVar + " (|Z|: " + ls + ", \Theta: " + options.LineStyle_Secondary + ")");
                set(ax, 'XScale', 'log')
                colororder({'k','k'})
                yyaxis right
                plot(Data_this_segment.Frequency_Hz, Data_this_segment.Z_Phase_Degree, ...
                    options.LineStyle_Secondary, 'Color', colors(i,:), 'LineWidth', lw);
                yyaxis left
        end
    end
end
% Decorations
switch options.Measurement
    case 'Voltage-capacity'
        xlabel('Capacity (Ah)')
        ylabel('Voltage (V)')
    case 'Voltage-time'
        xlabel('Time (s)')
        ylabel('Voltage (V)')
    case 'Polarization-time'
        xlabel('Time (s)')
        ylabel('Polarization [V - V_0] (V)')
    case 'Charge slippage'
        xlabel('Cumulative charge (Ah)')
        ylabel('Voltage (V)')
    case 'Differential capacity'
        xlabel('Voltage (V)')
        ylabel('Differential capacity (Ah\cdotV^{-1})')
    case 'Differential voltage'
        xlabel('Capacity (Ah)')
        ylabel('Differential voltage (V\cdotAh^{-1})')
    case 'EIS Nyquist'
        xlabel('Z_{Real} (\Omega)')
        ylabel('Z_{Imaginary} (\Omega)')
    case 'EIS Bode Real-Imaginary'
        xlabel('Frequency (Hz)')
        yyaxis left
        ylabel('Z_{Real} (\Omega)')
        yyaxis right
        ylabel('Z_{Imaginary} (\Omega)')
    case 'EIS Bode Mag-Phase'
        xlabel('Frequency (Hz)')
        yyaxis left
        ylabel('|Z| (\Omega)')
        yyaxis right
        ylabel('\Theta (\circ)')
end
if strcmp(options.LegendStyle, 'Legend')
    if contains(options.LegendLocation, 'tile_')
        loc = split(options.LegendLocation, '_');
        loc = loc{2};
        lgd = legend(gobj);
        lgd.Layout.Tile = loc;
    else
        legend(gobj, 'Location', options.LegendLocation)
    end
elseif strcmp(options.LegendStyle, 'Colorbar')
    colormap(ax, colors)
    if contains(options.LegendLocation, 'tile_')
        loc = split(options.LegendLocation, '_');
        loc = loc{2};
        cb = colorbar();
        cb.Layout.Tile = loc;
    else
        cb = colorbar('Location', options.LegendLocation);
    end
    cb.Label.String = tickVarStr;
    if n <= 5
        cb.Ticks = ((((1:n)-1).*2)+1)./(n*2);
        if strcmp(options.ColorVariable, 'Cycle_Index')
            cb.TickLabels = compose("%d", Cycles);
        else
            colorVar = SummaryData.(tick_variable);
            mask_cycles = any(SummaryData.Cycle_Index == Cycles, 2);
            colorVar = colorVar(mask_cycles);
            cb.TickLabels = compose("%0.3g", colorVar);
        end
    else
        cb.Ticks = 0:0.2:1;
        if strcmp(options.ColorVariable, 'Cycle_Index')
            cb.TickLabels = compose("%d", round(linspace(Cycles(1), Cycles(end), 6)));
        else
            colorVar = SummaryData.(tick_variable);
            mask_cycles = any(SummaryData.Cycle_Index == Cycles, 2);
            colorVar = colorVar(mask_cycles);
            cb.TickLabels = compose("%0.3g", linspace(colorVar(1), colorVar(end), 6));
        end
    end
    cb.Direction = 'reverse';
end
if ~isempty(options.AxisLimits)
    axis(ax, options.AxisLimits);
end

% Output ax only if it's requested
if nargout > 0
    varargout = {ax};
else
    varargout = {};
end
end