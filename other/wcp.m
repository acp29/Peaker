function [R_s, C_m, R_m] = wcp (t, y, V_step, step_start, step_duration, sampInt)

% [R_s, C_m, R_m] = wcp (t, y, V_step, step_start, step_duration, sampInt)
%
% Measures whole cell properties. Specifically, this function returns the
% voltage clamp step estimates of series resistance, input resistance, cell
% membrane resistance, cell membrane capacitance, cell surface area and
% specific membrane resistance.
%   
% The series (or access) resistance is obtained my dividing the voltage step
% by the peak amplitude of the current transient (Ogden, 1994): R_s = V / I_p
%    
% The input resistance is obtained by dividing the voltage step by the average
% amplitude of the steady-state current (Barbour, 2014): R_in = V / I_ss
%    
% The cell membrane resistance is calculated by subtracting the series
% resistance from the input resistance (Barbour, 1994): R_m = R_in - R_s
%    
% The cell membrane capacitance is estimated by dividing the transient charge
% by the size of the voltage-clamp step (Taylor et al. 2012): C_m = Q / V
%    
% The cell surface area is estimated by dividing the cell capacitance by the
% specific cell capacitance, c (1.0 uF/cm^2; Gentet et al. 2000; Niebur, 2008):
% Area = C_m / c
%    
% The specific membrane resistance is calculated by multiplying the cell
% membrane resistance with the cell surface area: rho = R_m * Area
% Users should be aware of the approximate nature of determining cell
% capacitance and derived parameters from the voltage-clamp step method
% (Golowasch, J. et al., 2009)
%
% References:
%   Barbour, B. (2014) Electronics for electrophysiologists. Microelectrode
%     Techniques workshop tutorial.
%     www.biologie.ens.fr/~barbour/electronics_for_electrophysiologists.pdf
%    Gentet, L.J., Stuart, G.J., and Clements, J.D. (2000) Direct measurement
%     of specific membrane capacitance in neurons. Biophys J. 79(1):314-320
%    Golowasch, J. et al. (2009) Membrane Capacitance Measurements Revisited:
%     Dependence of Capacitance Value on Measurement Method in Nonisopotential
%     Neurons. J Neurophysiol. 2009 Oct; 102(4): 2161-2175.
%    Niebur, E. (2008), Scholarpedia, 3(6):7166. doi:10.4249/scholarpedia.7166
%     www.scholarpedia.org/article/Electrical_properties_of_cell_membranes
%     (revision #13938, last accessed 30 April 2018)
%    Ogden, D. Chapter 16: Microelectrode electronics, in Ogden, D. (ed.)
%     Microelectrode Techniques. 1994. 2nd Edition. Cambridge: The Company
%     of Biologists Limited.
%    Taylor, A.L. (2012) What we talk about when we talk about capacitance
%     measured with the voltage-clamp step method J Comput Neurosci.
%     32(1):167-175
    
    % Prepare variables from input arguments
    t0 = step_start / sampInt;
    l = step_duration / sampInt;
    
    % Set cursors and update measurements
    b1 = 1 + fix((step_start - 0.001) / sampInt);
    b2 = 1 + fix(t0 - 1);
    p1 = 1 + fix(t0);
    p2 = 1 + fix((step_start + 0.001) / sampInt);
    f1 = 1 + fix(t0);
    f2 = 1 + fix(t0 + l - 1);

    % Subtract baseline
    b = mean(y(b1:b2));
    y = y - b;
    
    % Calculate series resistance (R_s) from initial negative transient
    peak = min(y(p1:p2));
    R_s = V_step / peak;  % in ohm

    % Calculate charge delivered during the voltage clamp step
    Q = trapz(t(f1:f2), y(f1:f2));

    % Set new fit cursor positions
    f1 = 1 + fix(t0 + l - 1 - ( step_duration / 4 ) / sampInt);
    f2 = 1 + fix(t0 + l - 1);

    % Measure steady state current and calculate input resistance
    I = mean(y(f1:f2));
    R_in  = V_step / I;                       % in ohm
    
    % Calculate cell membrane resistance
    R_m = R_in - R_s;                           % in ohm

    % Calculate voltage-clamp step estimate of the cell capacitance
    C_m = (Q - I * step_duration) / V_step;   % in F

    
end
