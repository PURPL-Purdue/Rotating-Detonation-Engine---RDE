function a = get_speed_sound(params,V)
% Speed of sound a² = a0²       - 0.5*(gamma-1)*V²
%                   = gamma*R*T - 0.5*(gamma-1)*V²
  a = sqrt ( params.gamma * params.R * params.T - ...
             0.5 * ( params.gamma - 1.) * V.^2 ) ;
end