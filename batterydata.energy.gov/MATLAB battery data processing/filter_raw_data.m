function data_filtered = filter_raw_data(ProcessedData, SummaryData, Cycle, Segment)
% FILTER_RAW_DATA Filters processed raw data by 'Cycle_Label' and
% 'Segment_Label'.
%   Filters processed raw I-V-t data by an input 'Cycle_Label' and an
%   optional 'Segment_Label'. Returns a cell array with each distinct
%   instance of the measurement being filtered. With no optional input
%   Required inputs:
%       ProcessedData (table): table of raw current-voltage-time timeseries
%           data. Must have been processed using the 'process_IVt_data'
%           function before using plot_measurements.
%       SummaryData (table): table of summarized battery performance
%           metrics, calculated using process_IVt_data.
%       Cycle: A Cycle_Label from SummaryData, 'All', or the numeric
%           indices of the cycles of interest
%   Optional inputs:
%       Segment (text): A Segment_Label to further filter each 'Cycle'
%   Output:
%       data_filtered (cell): cell array of data tables of each 'cycle'
arguments
    ProcessedData table
    SummaryData table
    Cycle % 'All', a Cycle_Label from SummaryData, or a vector of cycle numbers.
    Segment = '' % A Segment_Label from the 'Cycle' being filtered out
end
if isstring(Cycle) || ischar(Cycle)
    if strcmp(Cycle, 'All')
        if isempty(Segment)
            data_filtered = ProcessedData;
            return
        else
            Cycle = unique(ProcessedData.Cycle_Index);
        end
    else
        validCycleLabels = unique(SummaryData.Cycle_Label);
        if ~any(strcmp(Cycle, validCycleLabels))
            error("Must specify one Cycle_Label out of: " + join(validCycleLabels, ", ") + ".")
        else
            Cycle = SummaryData.Cycle_Index(strcmp(SummaryData.Cycle_Label, Cycle), :);
        end
    end
else
    assert(isvector(Cycle) && isnumeric(Cycle), 'The input Cycle must be a numeric vector.');
end
n = length(Cycle);
data_filtered = [];
for i_cycle = 1:n
    this_cycle = Cycle(i_cycle);
    mask_this_cycle = ProcessedData.Cycle_Index == this_cycle;
    Data_this_cycle = ProcessedData(mask_this_cycle, :);
    if ~isempty(Segment)
        mask_this_segment = strcmp(Data_this_cycle.Segment_Label, Segment);
        idx_segments = contiguous(double(mask_this_segment));
        if size(idx_segments, 1) > 1
            idx_segments = idx_segments{2,2};
            segment_starts = idx_segments(:,1);
            segment_stops = idx_segments(:,2);
        else
            segment_starts = 1;
            segment_stops = height(Data_this_cycle);
        end
    else
        segment_starts = 1;
        segment_stops = height(Data_this_cycle);
    end
    for i_segment = 1:length(segment_starts)
        Data_this_segment = Data_this_cycle(segment_starts(i_segment):segment_stops(i_segment), :);
        if isempty(data_filtered)
            data_filtered = {Data_this_segment};
        else
            data_filtered = [data_filtered; {Data_this_segment}];
        end
    end
end
end