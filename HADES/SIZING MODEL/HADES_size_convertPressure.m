function p_out = HADES_size_convertPressure(p_in,unit_in,unit_out)

unit_in = lower(unit_in); unit_out = lower(unit_out);

switch unit_in
    case 'psia', p_pa = p_in*6894.76;
    case 'atm',  p_pa = p_in*101325;
    case 'bar',  p_pa = p_in*1e5;
    case 'mmhg', p_pa = p_in*133.322;
    otherwise, error('Unsupported pressure unit');
end

switch unit_out
    case 'psia', p_out = p_pa/6894.76;
    case 'atm',  p_out = p_pa/101325;
    case 'bar',  p_out = p_pa/1e5;
    case 'mmhg', p_out = p_pa/133.322;
    otherwise, error('Unsupported pressure unit');
end
end