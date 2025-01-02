function varargout = plot_metric(SummaryData, Xvar, Yvar, LineSpec, options)
% PLOT_METRIC Plots summary data of battery performance metrics.
%   Plots summary data like capacity, efficiency, or voltage extracted from
%   raw cycling data. Optionally outputs the plotted axis. Formatted X and
%   Y axis labels are automatically added to the plot based on the Xvar and
%   Yvar inputs.
%   Required inputs:
%       Data (table): table of summary data
%       Xvar (text): variable name to plot on the x-axis. Can be any
%           variable from the SummaryData table.
%       Yvar (text): variable name to plot on the y-axis. Can be any
%           variable from the SummaryData table.
%   Optional inputs:
%       LineSpec (char): MATLAB line specification, like 'ok' or 'b-' (as when using plot)
%       options:
%           AxisHandle (Axes): An axis handle to plot on
%           AxisLimits (double): Axis limits, [x_min x_max y_min y_max]
%           Cycle: 'All', a Cycle_Label from SummaryData, or a vector 
%               of cycle numbers (default 'All').
%           Any properties of Line to specify plotting style (as when using plot).
%   Outputs:
%       ax (Axes): the axis handle that was plotted on
arguments
    SummaryData table
    Xvar {mustBeMember(Xvar, {'Cycle_Index','Charge_Throughput_Ah','Energy_Throughput_Wh','Absolute_Charge_Throughput_Ah','Absolute_Energy_Throughput_Wh',...
        'Equivalent_Full_Cycles','Equivalent_Charge_Slippage','Time_s','Time_d','Datenum_d','datenum_d','tsecs_start','tsecs_end','tsecs_cycle',...
        'Q_chg','q_chg','Q_dis','q_dis','CE','E_chg','e_chg','E_dis','e_dis','EE','DeltaV','deltaV','V_min','V_max','V_avg','I_min','I_max','I_avg',...
        'P_min','P_max','P_avg','T_min','T_max','T_avg'})}
    Yvar {mustBeMember(Yvar, {'Cycle_Index','Charge_Throughput_Ah','Energy_Throughput_Wh','Absolute_Charge_Throughput_Ah','Absolute_Energy_Throughput_Wh',...
        'Equivalent_Full_Cycles','Equivalent_Charge_Slippage','Time_s','Time_d','Datenum_d','datenum_d','tsecs_start','tsecs_end','tsecs_cycle',...
        'Q_chg','q_chg','Q_dis','q_dis','CE','E_chg','e_chg','E_dis','e_dis','EE','DeltaV','deltaV','V_min','V_max','V_avg','I_min','I_max','I_avg',...
        'P_min','P_max','P_avg','T_min','T_max','T_avg'})}
    LineSpec char = '.k'
    options.AxisHandle = [] % An axis handle to plot on
    options.AxisLimits = [] % Axis limits, [x_min x_max y_min y_max]
    options.Cycle = 'All' % 'All', a Cycle_Label from SummaryData, or a vector of cycle numbers.
    options.?matlab.graphics.chart.primitive.Line % Any property of Line
end

% Parse axis options
if isempty(options.AxisHandle)
    figure; hold on; box on; grid on; ax = gca;
else
    ax = options.AxisHandle; hold on; box on; grid on;
end
if isempty(options.AxisLimits)
    lims = [-Inf Inf -Inf Inf];
else
    lims = options.AxisLimits;
end
% Parse optional cycle input
if isstring(options.Cycle) || ischar(options.Cycle)
    if strcmp(options.Cycle, 'All')
        Cycles = unique(SummaryData.Cycle_Index);
        display_name = [];
    else
        validCycleLabels = unique(SummaryData.Cycle_Label);
        if ~any(strcmp(options.Cycle, validCycleLabels))
            error("Must specify one Cycle_Label out of: " + join(validCycleLabels, ", ") + ".")
        else
            Cycles = SummaryData.Cycle_Index(strcmp(SummaryData.Cycle_Label, options.Cycle), :);
            display_name = options.Cycle;
        end
    end
else
    assert(isvector(options.Cycle) && isnumeric(options.Cycle), 'The input options.Cycle must be a numeric vector.');
    Cycles = options.Cycle;
    display_name = [];
end
Cycles = any(SummaryData.Cycle_Index == Cycles', 2);
if ~isfield(options, 'DisplayName')
    % Use Cycle_Label as line label, or yvar if there's no Cycle_Label
    % filter
    if ~isempty(display_name)
        options.DisplayName = display_name;
    else
        options.DisplayName = get_tick_label(Yvar);
    end
end
% Remove from options struct and format options struct for input to plot
options = rmfield(options, {'AxisHandle', 'AxisLimits', 'Cycle'});

% Plot
options = namedargs2cell(options);
plot(ax, SummaryData{Cycles, Xvar}, SummaryData{Cycles, Yvar}, LineSpec, options{:})

% Decorations
axis(ax, lims)
xlabel(get_tick_label(Xvar))
ylabel(get_tick_label(Yvar))

% Output ax only if it's requested
if nargout > 0
    varargout = {ax};
else
    varargout = {};
end
end