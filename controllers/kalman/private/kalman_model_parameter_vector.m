function values = kalman_model_parameter_vector(model)
%KALMAN_MODEL_PARAMETER_VECTOR Extract idnlgrey parameter values.

values = zeros(numel(model.Parameters), 1);
for k = 1:numel(model.Parameters)
    value = model.Parameters(k).Value;
    if iscell(value)
        value = value{1};
    end
    values(k) = double(value);
end
end

