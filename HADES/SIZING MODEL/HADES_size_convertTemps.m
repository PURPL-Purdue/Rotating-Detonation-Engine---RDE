function T_out = HADES_size_convertTemps(T_in,unit_in,unit_out)

unit_in = upper(unit_in); unit_out = upper(unit_out);

switch unit_in
    case 'K', T_K = T_in;
    case 'C', T_K = T_in + 273.15;
    case 'F', T_K = (T_in-32)*5/9 + 273.15;
    case 'R', T_K = T_in*5/9;
    otherwise, error('Unsupported temperature unit');
end

switch unit_out
    case 'K', T_out = T_K;
    case 'C', T_out = T_K - 273.15;
    case 'F', T_out = (T_K-273.15)*9/5 + 32;
    case 'R', T_out = T_K*9/5;
end
end
