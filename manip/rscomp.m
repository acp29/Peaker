%     Script File: rscomp
%
%     This script is my Matlab/Octave adaptation of Erwin Neher's Igor functions:
%       RunSeriesResComp
%       SeriesresistanceComp
%     These were freely available in in Proc02_Apr4Web.ipf from Erwin Neher's webpage:
%      http://www3.mpibpc.mpg.de/groups/neher/index.php?page=software
%      (last accessed: 01 July 2014)
%
%     The function replaces current traces by their series-compensated version;
%     the value at i is replaced by the average at i and i+1
%     R_s is in ohms, C_m in Farads, fraction is the fraction to be compensated
%     if R_s was 5 MOhm in the experiment and if it was 50% hardware compensated
%     then R_s = 2.5e6 has to be entered and f=1 for a complete overall compensation
%     The routine, similarly to that published by Traynelis J. Neurosc. Meth. 86:25,
%     compensates the frequency response, assuming a single R_s*C_m - time constant
%     (at constant V-hold) and a three component equivalent circuit for the pipette cell
%     assembly with  R_s, C_m, R_m
%
%     Theory: V_h = R_s*I+U_m;
%     I_r is membrane resistive current,
%     I is total current
%     I_r = I-I_c = I-C_m*dU_m/dt = I+C_m*R_s*dI/dt (because R_s*I+U_m = const.)
%     G_m=I_r/ U_m = (I+C_m*R_s*dI/dt)/ (V_h-R_s*I)
%     For complete correction (fraction = 1) : I_corr = V_h*(I+C_m*R_s*dI/dt)/ (V_h-R_s*I)
%
%     rscomp v1.0 (last updated: 01/06/2014)
%     Author: Andrew Charles Penn
%     https://www.researchgate.net/profile/Andrew_Penn/


% Clear current variables, load data and define matrix columns
clear;
close all;
format short g
diary('on');
cwd = pwd;
file = input(sprintf('\nData filename (including extension): '),'s');
if ~isempty(regexpi(file(end-2:end),'.gz'))
  [pathstr,filename,ext] = fileparts(file(end-2:end));
elseif ~isempty(regexpi(file(end-3:end),'.zip'))
  [pathstr,filename,ext] = fileparts(file(end-3:end));
else
  [pathstr,filename,ext] = fileparts(file);
end
if ~isempty(pathstr)
  chdir(pathstr);
end
[data,xdiff,xunit,yunit,names,notes] = ephysIO (strcat(filename,ext));
t = data(:,1);
Y = data;
Y(:,1)=[];
numpoints = size(Y,1);
numtraces = size(Y, 2);
if ~all(diff(t))
 error('The data points are not evenly spaced')
end
sampInt = diff(t(1:2));

% Request and calculate the values for unknown variables
disp(sprintf('\nNumber of traces: \n')); disp(numtraces);
R_s = input(sprintf('\nEnter series resistance (in Mohm): '));
R_s = R_s * 1e6;
tau = input(sprintf('\nEnter cell time constant (in ms): '));
if isempty(tau)
 C_m = input(sprintf('\nEnter cell capacitance (in pF): '));
 C_m = C_m * 1e-12;
 tau = R_s * C_m;
 disp(sprintf('\nCalculated cell time constant (in ms): \n')); disp(tau*1e3);
else
 tau = tau * 1e-3;
 C_m = tau / R_s;
 disp(sprintf('\nCalculated cell capacitance (in pF): \n')); disp(C_m*1e12);
end
V_hold = input(sprintf('\nEnter holding potential (in mV): '));
if isempty(V_hold)
 error('The holding potential must be specified')
end
V_hold = V_hold * 1e-3;
V_reversal = input(sprintf('\nEnter reversal potential (in mV, default is 0): '));
if isempty(V_reversal)
 V_reversal = 0;
end
voltage = V_hold - V_reversal;
fraction = input(sprintf('\nEnter the fraction to be compensated (default is 1): '));
if isempty(fraction)
 fraction = 1;
end
disp(sprintf('\n'));
diary('off');
tau_corr = R_s * C_m * (1 - fraction);

for j=1:numtraces
 % Loop through the traces in the data file
 y = Y(:,j);

 % First point: (we have to calculate this separately, because we need the value at i-1 below)
 denominator = voltage - R_s * fraction * y(1);
 if denominator ~= 0
  y(1) = y(1) * (voltage / denominator);
 end

 for i = 2:numpoints-1
  % this is the loop doing the correction for all other points
  % first calculate R_m for zero series resistance under the assumptions
  % that  U_m + U_Rs = const = voltage
  current = (y(i+1) + y(i)) / 2;  % The in between(mean) value
  derivative =  (y(i+1) - y(i)) / sampInt;
  denominator = current + tau * derivative;
  if denominator ~= 0
   R_m = (voltage - R_s * current) / denominator;  % calculate the true R_m
  else
    % Do nothing. Leave R_m as is
  end
  % Now calculate current for new series resitance
  denominator = (R_m + (1 - fraction) * R_s) * (1 + tau_corr / sampInt);
  if denominator ~= 0
   y(i) = tau_corr / (tau_corr + sampInt) * y(i-1) + voltage/denominator;
  else
   y(i) = y(i-1);  % old value
  end
 end

 Y(:,j)=y;

end

% Save series resistance compensated data
newfilename = strcat(filename,'_cp.atf');
output_data=[t,Y];
ephysIO (newfilename,output_data,xunit,yunit);

% Save rscomp summary
movefile('./diary',strcat(filename,'_rscomp.txt'));
chdir(cwd)
