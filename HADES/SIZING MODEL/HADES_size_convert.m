function converted_value = HADES_size_convert(val, from_unit, to_unit)
% HADES_size_convert: Converts pressure, temperature, and speed units.
% Inputs: 
% val (double): The numerical value to convert
% from_unit (string): Current unit (e.g., 'mmHg', 'degC', 'mph')
% to_unit (string): Target unit


    from_unit = lower(from_unit);
    to_unit = lower(to_unit);
    
    %% Pressure Conversions
    % Using bar as the internal pivot point
    pressure_map = struct('bar', 1, 'atm', 1.01325, 'psi', 14.5038, 'mmhg', 750.062);
    
    if isfield(pressure_map, from_unit) && isfield(pressure_map, to_unit)
        val_in_bar = val / pressure_map.(from_unit);
        converted_value = val_in_bar * pressure_map.(to_unit);
        return;
    end

    %% Speed Conversions
    speed_map = struct('ms', 1, 'fts', 3.28084);
    
    % Mapping common aliases
    if strcmp(from_unit, 'm/s'); from_unit = 'ms'; end
    if strcmp(to_unit, 'm/s'); to_unit = 'ms'; end
    if strcmp(from_unit, 'ft/s'); from_unit = 'fts'; end
    if strcmp(to_unit, 'ft/s'); to_unit = 'fts'; end

    if isfield(speed_map, from_unit) && isfield(speed_map, to_unit)
        val_in_ms = val / speed_map.(from_unit);
        converted_value = val_in_ms * speed_map.(to_unit);
        return;
    end

    %% Temperature Conversions
    switch from_unit
        case {'k', 'kelvin'}
            temp_k = val;
        case {'c', 'degc', 'celsius'}
            temp_k = val + 273.15;
        case {'f', 'degf', 'fahrenheit'}
            temp_k = (val - 32) * 5/9 + 273.15;
        case {'r', 'rankine'}
            temp_k = val * 5/9;
        otherwise
            temp_k = NaN;
    end
    
   
    if ~isnan(temp_k)
        switch to_unit
            case {'k', 'kelvin'}
                converted_value = temp_k;
            case {'c', 'degc', 'celsius'}
                converted_value = temp_k - 273.15;
            case {'f', 'degf', 'fahrenheit'}
                converted_value = (temp_k - 273.15) * 9/5 + 32;
            case {'r', 'rankine'}
                converted_value = temp_k * 9/5;
            otherwise
                error('Target unit not recognized.');
        end
        return;
    end

    error('Unit conversion from %s to %s is not supported.', from_unit, to_unit);
end