%  Function File: fitter
%
%  peak = fitter(file,TC,s,R)
%  peak = fitter(file,TC,s,R,...,'wave',wave)
%  peak = fitter(file,TC,s,R,...,'peak',peak)
%  peak = fitter(file,TC,s,R,...,'times',times)
%  peak = fitter(file,TC,s,R,...,'include',include)
%  peak = fitter(file,TC,s,R,...,'exclude',exclude)
%  peak = fitter(file,TC,s,R,...,'latency',latency)
%  peak = fitter(file,TC,s,R,...,'hpf',hpf)
%  peak = fitter(file,TC,s,R,...,'lpf',lpf)
%  peak = fitter(file,TC,s,R,...,'channel',channel)
%  peak = fitter(file,TC,s,R,...,'config',config)
%  [peak,area] = fitter(...)
%  [peak,area,tau1] = fitter(...)
%  [peak,area,tau1,tau2] = fitter(...)
%
%  peak = fitter(file,TC,s,R) returns the amplitudes of peaks measured by
%    least-squares fitting of the events, where each event is modelled as
%    the sum of two exponentials:
%
%      f(t) = exp ( - t / tau_decay ) - exp ( - t / tau_rise )
%
%    The time constants must be provided (in units seconds) as a 2-element
%    row vector (TC). The sign (s) of the event peak deflections in the wave
%    file must be specified as '-' or '+'. The fits can be restrained close
%    to the initial time constant values defined in TC by providing a non-
%    negative scalar value R, which is the scale factor of a penalty term
%    added to the least-squares objective function:
%
%      sum (( y - f(peak,TC))^2) + R * sum((TC - TC0)^2)
%
%    where TC0 is the initial values of TC
%
%    The default value of R is 1e-6. Increase R to tighten restraints on
%    the time constants. For unrestrained fits, set R to 0. Fitting is
%    performed by Nelder-Mead Simplex optimization of log-transformed free
%    parameter values. Thus, only positive parameter space is searched
%    (even when fitting is unrestrained). The optimization problem is solved
%    using MATLABs native fminsearch function. The file extension must be
%    included in the filename of the file input argument. See the associated
%    input-output function ephysIO for details of supported input file formats.
%
%  peak = fitter(file,TC,s,R,...,'wave',wave) sets fitter to select
%    the wave number from the file. If none is specified, eventer analyses
%    the first wave. Note that wave numbers have zero indexing
%
%  peak = fitter(file,TC,s,R,...,'peak',peak) sets the initial peak amplitude
%    for the fitting.
%
%  peak = fitter(file,TC,s,R,...,'times',times) sets the times of the
%    events for fitting.
%
%  peak = fitter(file,TC,s,R,...,'include',include) sets the region of
%    the wave to perform fitting. The region must be specified as a
%    a 2-element row vector.
%
%  peak = fitter(file,TC,s,R,...,'exclude',exclude) sets exclusion zones,
%    which must be specified as a 2-column matrix, where the first and
%    second columns define the start and end times of the excluded zones
%    (in seconds). Linear interpolation occurs within each exclusion zone.
%
%  peak = fitter(file,TC,s,R,...,'latency',latency) sets fitter to add
%    a delay of magnitude equal to latency to each event time.
%
%  peak = fitter(file,TC,s,R,...,'hpf',hpf) sets the -3 dB cut-off (in
%    Hz) of the low-pass median filter, where the filtered wave is then
%    subtracted from the original wave. Thus, this is essentially a
%    high-pass filter. The algorithm implements a bounce correction
%    to avoid end effects. The default cut-off is 0.1 Hz.
%
%  peak = fitter(file,TC,s,R,...,'lpf',lpf) sets the -3 dB cut-off (in
%    Hz) of the low-pass binomial filter applied to the original wave.
%    The default cut-off is inf Hz.
%
%  peak = fitter(file,TC,s,R,...,'channel',channel) sets fitter to
%    select the recording channel from the file. If none is specified,
%    eventer analyses channel number 1. This argument is ignored for
%    loaded filetypes that do not support multiple recording channels.
%
%  peak = fitter(file,TC,s,R,...,'config',config) sets the configuration
%    of the recording wave to either 'VC' (for voltage clamp) or 'CC'
%    (for current-clamp). The default is neither (i.e. '').
%
%  [peak,area] = fitter(...) returns the area of each event.
%
%  [peak,area,tau1] = fitter(...) returns the rise time constant of
%    the events.
%
%  [peak,area,tau1,tau2] = fitter(...) returns the decay time constant
%    of the events.
%
%  [peak,area,tau1,tau2,modelWave] = fitter(...) returns the model wave
%
%  Dependencies: fminsearch, medianf, bounce, ephysIO.
%
%
%  fitter v1.0 (last updated: 17/04/2018)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/

function [MSE,peak,area,tau1,tau2,modelWave] = fitter(file,TC,s,R,varargin)

  % Initialize
  close all
  format short g

  % Error checking
  if nargin<2 || sum(size(TC))>3
    error('A vector defining the two initial time constant values must be specified');
  else
    if sum(sign(TC))~=2
      error('The time constants must be non-zero and non-negative');
    end
    if TC(1)>=TC(2)
      error('The first time constant (rise) must be smaller than the second time constant (decay)');
    end
  end

  if nargin<3
    error('The sign of the peaks must be specified in the + or - direction');
  else
    if s~='+' && s~='-'
      error('The sign of the peaks must be specified in the + or - direction');
    end
  end

  if nargin<4 || isempty(R)
    R = 1e-6;  % Set to zero for unconstrained
  else
    if R ~= abs(R)
      error('The time constant restraint factor must be non-negative');
    end
  end

  % Set additional options
  options = varargin;
  wave = 1+find(strcmp('wave',options));
  peak0 = 1+find(strcmp('peak',options));
  ET = 1+find(strcmp('times',options));
  incl = 1+find(strcmp('include',options));
  excl = 1+find(strcmp('exclude',options));
  latency = 1+find(strcmp('latency',options));
  hpf = 1+find(strcmp('hpf',options));
  lpf = 1+find(strcmp('lpf',options));
  channel = 1+find(strcmp('channel',options));
  config = 1+find(strcmp('config',options));
  if ~isempty(peak0)
    try
      peak0 = options{peak0};
    catch
      peak0 = [];
    end
  else
    peak0 = [];
  end
  if ~isempty(ET)
    try
      ET = options{ET};
    catch
      ET = [];
    end
  else
    ET = [];
  end
  if ~isempty(wave)
    try
      wave = options{wave};
    catch
      wave = 0;
    end
  else
    wave = 0;
  end
  wave = wave+1; % Add 1 since the input uses zero indexing
  wave = wave+1; % Add another 1 since first column in data matrix is time
  if ~isempty(incl)
    try
      incl = options{incl};
    catch
      incl = [];
    end
  else
    incl = [];
  end
  if ~isempty(excl)
    try
      excl = options{excl};
    catch
      excl = [];
    end
  else
    excl = [];
  end
  if ~isempty(latency)
    try
      latency = options{latency};
    catch
      latency = 0;
    end
  else
    latency = 0;
  end
  if ~isempty(hpf)
    try
      hpf = options{hpf};
    catch
      hpf = 0.1;
    end
  else
    hpf = 0.1;
  end
  if ~isempty(lpf)
    try
      lpf = options{lpf};
    catch
      lpf = inf;
    end
  else
    lpf = inf;
  end
  if ~isempty(channel)
    try
      channel = options{channel};
    catch
      channel = 'none';
    end
  else
    channel = 1;
  end
  if ~isempty(config)
    try
      config = options{config};
    catch
      config = '';
    end
  else
    config = '';
  end

  %%%%%%%%%%%%%%%%%%%%%%%% CUSTOM SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%
  if isempty(peak0)
    peak0 = ones(1,4)*1e-3;   % 1 mV
  end
  delay = 0.2;
  if isempty(ET)
    if wave == 2
      ET = [delay delay+1.0*1 delay+1.0*2 delay+1.0*3];
    elseif wave == 3
      ET = [delay delay+0.8*1 delay+0.8*2 delay+0.8*3];
    elseif wave == 4
      ET = [delay delay+0.6*1 delay+0.6*2 delay+0.6*3];
    elseif wave == 5
      ET = [delay delay+0.4*1 delay+0.4*2 delay+0.4*3];
    elseif wave == 6
      ET = [delay delay+0.2*1 delay+0.2*2 delay+0.2*3];
    elseif wave == 7
      ET = [delay delay+0.1*1 delay+0.1*2 delay+0.1*3];
    elseif wave == 8
      ET = [delay delay+0.08*1 delay+0.08*2 delay+0.08*3];
    elseif wave == 9
      ET = [delay delay+0.06*1 delay+0.06*2 delay+0.06*3];
    elseif wave == 10
      ET = [delay delay+0.04*1 delay+0.04*2 delay+0.04*3];
    elseif wave == 11
      ET = [delay delay+0.02*1 delay+0.02*2 delay+0.02*3];
    end
  end
  if isempty(incl)
    if wave == 2
      incl = [0.1 4.20];
    elseif wave == 3
      incl = [0.1 3.40];
    elseif wave == 4
      incl = [0.1 2.60];
    elseif wave == 5
      incl = [0.1 2.00];
    elseif wave == 6
      incl = [0.1 0.95];
    elseif wave == 7
      incl = [0.1 0.66];
    elseif wave == 8
      incl = [0.1 0.55];
    elseif wave == 9
      incl = [0.1 0.50];
    elseif wave == 10
      incl = [0.1 0.42];
    elseif wave == 11
      incl = [0.1 0.35];
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Load data
  cwd = pwd;
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
  [data,xdiff,xunit,yunit,names,notes] = ephysIO ({strcat(filename,ext),channel});
  if wave > size(data,2)
    error('wave number exceeds data dimensions')
  end
  filewave = sprintf('%s_ch%d_%s',filename,channel,names{wave});
  trace = data(:,wave);
  if s=='-'
    trace = trace*-1;
  end
  t = data(:,1);
  trace = filter1 (trace, t, hpf, lpf, 'median');

  % Interpolate across excluded regions
  if ~isempty(excl)
    excl_idx = uint32(excl/xdiff);
    for i=size(excl,1)
      trace(excl_idx(i,1):excl_idx(i,2)) = ...
        interp1q([t(excl_idx(i,1)),t(excl_idx(i,2))]',...
        [trace(excl_idx(i,1)),trace(excl_idx(i,2))]',...
        t(excl_idx(i,1):excl_idx(i,2)));
    end
  end

  % Blank stimulus artifact and setup baseline
  % Pre-event baseline is considered the median during the 1 ms period before the event
  i = uint32(ET/xdiff);
  n = numel(ET);
  base = zeros(numel(t),1);
  k = ET(2)-ET(1);
  B = ones(1,n);
  for j=1:n
    trace(i(j):i(j)+latency/xdiff-1) = NaN;
    B(j) = median(trace(i(j)-50:i(j)-1));
  end
  lim = TC(2)*3 +...
        (TC(2)*TC(1))/(TC(2)-TC(1))*log(TC(2)/TC(1)); % 3 times the decay tau after peak
  if k>lim
    base = interp1q([0 ET t(end)]',[trace(1) B trace(end)]',t);
  else
    trace = trace-B(1);
  end
  trace = trace-base;
  figure(1);
  plot(t,trace,'b');

  % Create vector of time constants
  tau1 = ones(1,n)*TC(1);
  tau2 = ones(1,n)*TC(2);

  % Correct event times for latency
  times = ET + latency;

  % Initialize variables
  times = times - incl(1);
  incl_idx = uint32(1+incl/xdiff);
  trace = trace(incl_idx(1):incl_idx(2));
  t = t(incl_idx(1):incl_idx(2));
  M = numel(t);
  base = zeros(M,1);

  % Convert the initial peak amplitude to the initial P scale
  % parameter (P0) in the sum-of-two-exponentials expression
  tpeak = (tau2.*tau1)./(tau2-tau1).*log(tau2./tau1);
  P0 = peak0 ./ (-exp(-tpeak./tau1)+exp(-tpeak./tau2));

  % Fit events
  [P,tau1,tau2,MSE] = optim_fit(P0,tau1,tau2,times,R,xdiff,base,trace);

  % Calculate and plot the model trace
  hold on;plot(t,base,'g');hold off;
  [modelWave,Y] = sum_events([log(P);log(tau1);log(tau2)],times,xdiff,base);
  figure(1);
  hold on; plot(t,modelWave,'r','linewidth',3);hold off;
  hold on; plot(t,Y,'r-');hold off;

  % Measure the amplitude of the peaks
  tpeak = (tau2.*tau1)./(tau2-tau1).*log(tau2./tau1);
  peak = sum2exp(P,tau1,tau2,tpeak);

  % Measure the area under each event (i.e. the integral)
  area = abs(P(1,:).*(tau1-tau2));

  % Save data measurements and figures
  if exist('fitter.output','dir')==0
    mkdir('fitter.output');
  end
  cd fitter.output
  if exist(filewave,'dir')==0
    mkdir(filewave);
  end
  chdir(filewave);
  if exist('fig','dir')==0
    mkdir('fig');
  end
  chdir('fig');
  saveas(gcf,'output.fig','fig');
  cd ..
  if exist('txt','dir')==0
    mkdir('txt');
  end
  chdir('txt');
  dlmwrite('peak.txt',peak','\t');
  dlmwrite('area.txt',area','\t');
  dlmwrite('rise.txt',tpeak','\t');
  dlmwrite('taus.txt',[tau1' tau2'],'\t');
  cd ..
  fid = fopen('inp.m','w');
  fprintf(fid,'MSE = %d;\n',MSE);
  fprintf(fid,'[peak,area,tau1,tau2] = ');
  fprintf(fid,'fitter(''%s'',[%d %d],''%s'',%d,',file,TC(1),TC(2),s,R);
  fprintf(fid,'''wave'',%i,',wave-2);
  fprintf(fid,'''peak'',[');
  for i = 1:numel(peak0)
    fprintf(fid,'%d',peak0(i));
    if i == numel(peak0)
      break
    end
    fprintf(fid,',');
  end
  fprintf(fid,'],');
  fprintf(fid,'''times'',[');
  for i = 1:numel(ET)
    fprintf(fid,'%d',ET(i));
    if i == numel(ET)
      break
    end
    fprintf(fid,',');
  end
  fprintf(fid,'],');
  if isempty(incl)
    fprintf(fid,'''include'',[],');
  else
    fprintf(fid,'''include'',[%d %d],',incl);
  end
  fprintf(fid,'''exclude'',[');
  for i = 1:size(excl,1)
    fprintf(fid,'%d %d',excl(i,:));
    if i == size(excl,1)
      break
    end
    fprintf(fid,';',excl(1,:));
  end
  fprintf(fid,'],');
  fprintf(fid,'''latency'',%d,',latency);
  fprintf(fid,'''hpf'',%d,',hpf);
  fprintf(fid,'''lpf'',%d,',lpf);
  fprintf(fid,'''channel'',%i,',channel);
  fprintf(fid,'''config'',''%s''',config);
  fprintf(fid,')');
  fclose(fid);
  cd ../..

end

function [P,tau1,tau2,MSE] = optim_fit(P,tau1,tau2,times,R,xdiff,base,trace)

  % Adjust event amplitude and timecourse to minimize objective function

  % Initialize parameters
  n = numel(times);

  % Natural log fransformation of tau
  % Required for optimizer to only search parameter space > 0
  p0 = log([P;tau1;tau2]);

  % Generate objective function
  objfunc = @(p,R) nansum((trace-sum_events(p,times,xdiff,base)).^2) +...
                 R * sum(sum((p(2:end,:)-p0(2:end,:)).^2,1)); % Penalty to restrain tau(s)
  minfunc = @(p) objfunc(p,R);

  % Optimize using Nelder-Mead simplex algorithm
  options=optimset('PlotFcns',@optimplotfval,...
                   'TolFun',1e-3,'TolX',1e-3,...
                   'MaxFunEval',inf,'MaxIter',inf);
  [p,fval,exitflag] = fminsearch(minfunc,p0,options);
  R = 0;
  SSE = objfunc(p,R);
  MSE = SSE/sum(~isnan(trace));

  % Prepare return values
  P = exp(p(1,:));
  tau1 = exp(p(2,:));
  tau2 = exp(p(3,:));

end

function [modelWave,Y] = sum_events(p,times,xdiff,base)

  % Virtual summation of model events

  % Model individual events
  P = exp(p(1,:));
  tau1 = exp(p(2,:));
  tau2 = exp(p(3,:));
  [modelEvents] = model_events(P,tau1,tau2,xdiff);

  % Initialize variables
  [m,n] = size(modelEvents);
  m = uint32(m);
  n = uint32(n);
  M = numel(base);
  modelWave = base;
  i = uint32(times/xdiff);
  Y = zeros(M,n);

  % Event summation
  for j=1:n
    Y(i(j):i(j)+min(m-1,M-i(j)),j) = ...
    Y(i(j):i(j)+min(m-1,M-i(j)),j) + ...
    modelEvents(1:min(m,M-i(j)+1),j);
  end
  Y(~any(Y,2),:) = NaN;
  modelWave = sum(Y,2)+base;

end


function [modelEvents] = model_events(P,tau1,tau2,xdiff)

  % Vectorized construction of model events

  % Initialize variables
  t = (0:xdiff:1)'; % Calculate model event up to 1 second
  m = numel(t);
  n = numel(P);
  P = ones(m,1)*P;
  t = (t*ones(1,n));
  tau1 = (ones(m,1)*tau1);
  tau2 = (ones(m,1)*tau2);

  % Model events using the sum of two exponentials
  modelEvents = sum2exp(P,tau1,tau2,t);

end

function [y] = sum2exp(P,tau1,tau2,t)

  f = @(P,tau1,tau2,t) P .* (-exp(-t./tau1)+exp(-t./tau2));
  y = f(P,tau1,tau2,t);

end
