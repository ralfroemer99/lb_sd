function [delta] = invert_raw_to_delta(raw_thrust)
    pwm_inv = raw_thrust / 65535.0;
    current_thrust_inv = 1.00058 - 0.0000992573 * sqrt( 1.03503 * 10^8 - 9.0112 * 10^7 * pwm_inv);
    delta = current_thrust_inv - (0.033 * 9.8);
end

