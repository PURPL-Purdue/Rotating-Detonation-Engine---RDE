function [phi, y_inj, x_inj] = MoC_Infill_Triangle(v_inj, v_cj, y_inj)
    % Base of the triangle is the unwrapped circumference
    
    % Height of the fresh mixture layer
    x_inj = y_inj * (v_cj / v_inj);
    
    % Angle of the bounding gas interface
    phi = asin(v_inj / v_cj);

    % Convert to degrees
    phi = phi * 180 / pi;
end