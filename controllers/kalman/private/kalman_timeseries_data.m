function data = kalman_timeseries_data(value)
%KALMAN_TIMESERIES_DATA Convert numeric or timeseries data to a column.

if isa(value, 'timeseries')
    data = value.Data;
else
    data = value;
end

data = double(squeeze(data));
data = data(:);
end

