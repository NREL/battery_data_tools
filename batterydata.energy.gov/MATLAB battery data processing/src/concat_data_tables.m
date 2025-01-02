function RawData = concat_data_tables(varargin)
% CONCAT_DATA_TABLES Concatenates battery data tables.
%   Concatenates tables of raw battery data. The concatenation process
%   combines the input tables in the order they were provided, updating
%   cycle count and time variables for each file in the sequence. Data
%   tables from various measurements may have new table variables (for
%   instance, impedance data variables); these variables are automatically
%   filled in with NaNs if they are missing from other tables.

% Sort tables in Datenum_d order. Check for datenum overlap.
dt_start = zeros(length(varargin),1);
dt_end = zeros(length(varargin),1);
for i = 1:length(varargin)
    data = varargin{i};
    dt_start(i) = data.Datenum_d(1);
    dt_end(i) = data.Datenum_d(end);
end
% useful for figuring out what's going on with your file starts/stops
dt_ = [dt_start - dt_start(1), dt_end - dt_start(1)]; 
dt_ = [dt_, dt_(:,2)-dt_(:,1)];
[dt_start, idx_sort] = sort(dt_start);
dt_end = dt_end(idx_sort);
n = 2*length(varargin);
dt = zeros(n,1);
dt(1:2:n-1) = dt_start;
dt(2:2:n) = dt_end;
if ~issorted(dt, 'ascend')
    error('Datenum across data tables overlaps, ensure raw data table do not overlap in time, which should be impossible')
end
% Calculate order to read tables
[~, idxSorted] = sort(dt_start);
for i = 1:length(varargin)
    data = varargin{idxSorted(i)};
    % disp(i)
    try
        check_IVt_data_formatting(data)
    catch
        warning("Data table number " + i + " failed formatting check using check_IVt_data_formatting.")
        err = lasterror();
        error(err.message)
    end
    if min(data.Cycle_Index) <= 0
        data.Cycle_Index = data.Cycle_Index + abs(min(data.Cycle_Index)) + 1;
    end
    if i == 1
        RawData = data;
    else
        % Check for variables that aren't either both tables, fill with NaN
        vars_xor = setxor(RawData.Properties.VariableNames, data.Properties.VariableNames);
        if ~isempty(vars_xor)
            vars_new = vars_xor(~contains(vars_xor, RawData.Properties.VariableNames));
            RawData{:, vars_new} = NaN;
            vars_new = vars_xor(~contains(vars_xor, data.Properties.VariableNames));
            data{:, vars_new} = NaN;
        end
        % Add cycle count, time, other offsets
        data.Cycle_Index = data.Cycle_Index + RawData.Cycle_Index(end);
        data.Time_s = data.Time_s + RawData.Time_s(end) + (data.Datenum_d(1) - RawData.Datenum_d(end))*(3600*24);
        % Append
        RawData = [RawData; data];
    end
end