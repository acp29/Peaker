%  Function File: eventer
%
%  peak = eventer(file,TC,s,SF)
%  peak = eventer(file,TC,s,SF,...,'exclude',exclude)
%  peak = eventer(file,TC,s,SF,...,'rmin',rmin)
%  peak = eventer(file,TC,s,SF,...,'hpf',hpf)
%  peak = eventer(file,TC,s,SF,...,'lpf',lpf)
%  peak = eventer(file,TC,s,SF,...,'taus',taus)
%  peak = eventer(file,TC,s,SF,...,'baseline',baseline)
%  peak = eventer(file,TC,s,SF,...,'win',win)
%  peak = eventer(file,TC,s,SF,...,'lambda',lambda)
%  peak = eventer(file,TC,s,SF,...,'config',config)
%  peak = eventer(file,TC,s,SF,...,'average',average)
%  peak = eventer(file,TC,s,SF,...,'export',format)
%  peak = eventer(file,TC,s,SF,...,'channel',channel)
%  peak = eventer(file,TC,s,SF,...,'wave',wave)
%  peak = eventer(file,TC,s,SF,...,'exmode',wave)
%  [peak,IEI] = eventer(...)
%
%  peak = eventer(file,TC,s,SF) returns the amplitudes of peaks detected by
%    an algorithm implementing FFT-based deconvolution of a wave composed
%    of spontaneously occurring EPSC- or EPSP-like waveforms [1]. The
%    algorithm uses a template modeled by the sum of two exponentials,
%    whose time constants must be provided in units seconds as a vector
%    (TC). The sign (s) of the event peak deflections in the wave file
%    must be specified as '-' or '+'. Event times are then defined where
%    there are maxima of delta-like waves exceeding a threshold, which
%    is set to the standard deviation of the noise of the deconvoluted wave
%    multiplied by a scale factor (SF). Events are then extracted as
%    episodic data and the factor parameter of least-squares fits with the
%    template define the peak amplitudes. See the associated input-output
%    function ephysIO for details of supported input file formats. The file
%    extension must be included in the filename of the file input argument.
%    The template kinetics is modeled by the difference of 2 exponentials:
%      f(t) = exp ( - t / tau_decay ) - exp ( - t / tau_rise )
%
%  peak = eventer(file,TC,s,SF,...,'exclude',exclude) sets exclusion zones,
%    which must be specified as a 2-column matrix, where the first and
%    second columns define the start and end times of the excluded zones
%    (in seconds). Events detected in these regions are discarded. By
%    default, there are no exclusion zones.
%
%  peak = eventer(file,TC,s,SF,...,'rmin',rmin) sets the correlation
%    coefficient for fitted model template and the event. Events with a
%    correlation coefficient (r) less than the value defined in rmin are
%    discarded [2]. The default rmin is 0.4, where r < 0.4 is considered
%    a very weak correlation.
%
%  peak = eventer(file,TC,s,SF,...,'hpf',hpf) sets the -3 dB cut-off (in
%    Hz) of the low-pass median filter applied to the deconvoluted wave,
%    where the filtered wave is then subtracted from the deconvoluted
%    wave. Thus, this is essentially a high-pass filter. The algorithm
%    implements a bounce correction to avoid end effects. The default
%    cut-off is 1 Hz.
%
%  peak = eventer(file,TC,s,SF,...,'lpf',lpf) sets the -3 dB cut-off (in
%    Hz) of the low-pass binomial filter applied to the deconvoluted wave.
%    The default cut-off is 200 Hz.
%
%  peak = eventer(file,TC,s,SF,...,'taus',taus) sets the number of time
%    constants after the peak of the template to use when fitting the
%    template to the detected events. The default is 0.4, which corresponds
%    to the time for the template to decay by 25 % of the peak. This
%    limited fit ensures that peak measurements are not compromised at high
%    event frequencies where there is event overlap.
%
%  peak = eventer(file,TC,s,SF,...,'baseline',baseline) sets the length
%    of time to use as the pre-event baseline for the template fit in
%    seconds.
%
%  peak = eventer(file,TC,s,SF,...,'win',win) sets the event window limits
%    in seconds for the conversion of the continuous wave to episodic
%    data. By default the limits are set at [-0.01 0.04].
%
%  peak = eventer(file,TC,s,SF,...,'lambda',lambda) sets the damping
%    factor used in the Levenberg-Marquardt ordinary non-linear least-
%    squares fitting procedures. The default is 1. The higher the value
%    of lambda, the more robust the fitting procedure is to the initial
%    values but the greater the number of iterations it will use.
%
%  peak = eventer(file,TC,s,SF,...,'config',config) sets the configuration
%    of the recording wave to either 'VC' (for voltage clamp) or 'CC'
%    (for current-clamp).The default is 'VC'.
%
%  peak = eventer(file,TC,s,SF,...,'average',average) sets eventer to
%    compute either the ensemble 'mean' or 'median' of the merged event
%    data. The default is 'mean'.
%
%  peak = eventer(file,TC,s,SF,...,'export',format) sets eventer to
%    export the episodic wave data of all detected events in either
%    Axon ('atf') or Igor ('itx') text file formats.
%
%  peak = eventer(file,TC,s,SF,...,'channel',channel) sets eventer to
%    select the recording channel from the file. If none is specified,
%    eventer analyses channel number 1. This argument is ignored for
%    loaded filetypes that do not support multiple recording channels.
%
%  peak = eventer(file,TC,s,SF,...,'wave',wave) sets eventer to select
%    the wave number from the file. If none is specified, eventer analyses
%    the first wave.
%
%  peak = eventer(file,TC,s,SF,...,'exmode',exmode) tells eventer what to
%    do with the first event proceeding each exclusion zone. For mode = 1
%    (default) IEIs are calculated for these events from the last event
%    preceeding the exclusion zone under the assumption that no events
%    occurred during the exclusion zone. For mode = 2, these events are
%    assigned an IEI of NaN. Note that events with NaN values are excluded
%    during the merge (i.e. that are in the ALL_events output directory).
%
%  [peak,IEI] = eventer(...) returns the interevent intervals preceding
%    each peak event. The first event from the start of the wave is assigned
%    an IEI of NaN.
%
%  See the example distributed with this function.
%
%  Dependencies: sinv, binomialf, bounce, hpfilter, lpfilter, lsqfit and
%                ephysIO.
%
%  References:
%   [1] Pernia-Andrade AJ, Goswami SP, Stickler Y, Frobe U, Schlogl A,
%    Jonas P. (2012) A deconvolution-based method with high sensitivity
%    and temporal resolution for detection of spontaneous synaptic currents
%    in vitro and in vivo. Biophys J. 103(7):1429-39.
%   [2] Jonas P, Major G and Sakmann (1993) Quantal components of unitary
%    EPSCs at the mossy fibre synapse on CA3 pyramidal cells of the rat
%    hippocampus. J Physiol. 472:615-663.
%
%  eventer v1.4 (last updated: 27/07/2016)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/


function eventer(file,TC,s,SF,varargin)

  if nargin<2 || sum(size(TC))>3
    error('A vector defining the time constants of the template has to be specified');
  else
    if sum(sign(TC))~=2
      error('The time constants must be non-zero and non-negative');
    end
    if TC(1)>=TC(2)
      error('The first time constant (rise) must be smaller than the second time constant (decay)');
    end
  end

  if nargin<3
    error('The sign of the peaks has to be specified');
  end

  if s~='+' && s~='-'
    error('The sign of the peaks has to be specified in the + or - direction');
  end

  if nargin<4
    error('The factor of standard deviations has to be specified to set the detection threshold');
  else
    if isempty(SF)
      SF = 4;
    end
  end

  % Set additional options
  options = varargin;
  excl = 1+find(strcmp('exclude',options));
  rmin = 1+find(strcmp('rmin',options));
  hpf = 1+find(strcmp('hpf',options));
  lpf = 1+find(strcmp('lpf',options));
  taus = 1+find(strcmp('taus',options));
  base = 1+find(strcmp('baseline',options));
  win = 1+find(strcmp('win',options));
  config = 1+find(strcmp('config',options));
  lambda = 1+find(strcmp('lambda',options));
  average = 1+find(strcmp('average',options));
  export = 1+find(strcmp('export',options));
  channel = 1+find(strcmp('channel',options));
  wave = 1+find(strcmp('wave',options));
  exmode = 1+find(strcmp('exmode',options));
  if ~isempty(excl)
    try
      excl = options{excl};
    catch
      excl = [];
    end
  else
    excl = [];
  end
  if ~isempty(rmin)
    try
      rmin = options{rmin};
    catch
      rmin = 0.4;
    end
  else
    rmin = 0.4;
  end
  % Default value of 0.4 corresponds to the time for the event to decay by 25% of the peak
  if ~isempty(taus)
    try
      taus = options{taus};
    catch
      taus = 0.4;
    end
  else
    taus = 0.4;
  end
  if ~isempty(hpf)
    try
      hpf = options{hpf};
    catch
      hpf = 1;
    end
  else
    hpf = 1;
  end
  if ~isempty(lpf)
    try
      lpf = options{lpf};
    catch
      lpf = 200;
    end
  else
    lpf = 200;
  end
  if ~isempty(base)
    try
      base = abs(options{base});
    catch
      base = 1e-03;
    end
  else
    base = 1e-03;
  end
  if ~isempty(win)
    try
      win = options{win};
    catch
      win = [-0.01 0.04];
    end
  else
    win = [-0.01 0.04];
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
  if ~isempty(lambda)
    try
      lambda = options{lambda};
    catch
      lambda = 1;
    end
  else
    lambda = 1;
  end
  if ~isempty(average)
    try
      average = options{average};
    catch
      average = 'mean';
    end
  else
    average = 'mean';
  end
  if ~isempty(export)
    try
      export = options{export};
    catch
      export = 'none';
    end
  else
    export = 'none';
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
  if ~isempty(wave)
    try
      wave = options{wave};
    catch
      wave = 1;
    end
  else
    wave = 1;
  end
  wave = wave+1; % Add 1 since first column in data matrix is time
  if ~isempty(exmode)
    try
      exmode = options{exmode};
    catch
      exmode = 1;
    end
  else
    exmode = 1;
  end


  % Error checking
  if size(win,1)~=1 && size(win,2)~=2
      error('win must be a vector defining the limits of the event window');
  end

  % Load data
  close all
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
  if numel(wave) == 1
    filewave = sprintf('%s_ch%d_%s',filename,channel,names{wave});
  else
    filewave = sprintf('%s_ch%d_%s',filename,channel,'_list');
  end
  if isempty(xunit)
    warning('xunit is undefined. X dimenion units assumed to be seconds.')
    xunit = 's';
  elseif ~strcmp(xunit,'s')
    error('expected xunit to be seconds (s)')
  end
  if isempty(yunit)
    if strcmp(config,'VC')
      warning('yunit assumed to be amps (A)')
      yunit = 'A';
    elseif strcmp(config,'CC')
      yunit = 'V';
      warning('yunit assumed to be volts (V)')
    else
      error('The recording mode or units must be defined')
    end
  else
    if ~strcmp(yunit,'A') && ~strcmp(yunit,'V')
      error('expected yunit to be amps (A) or volts (V)')
    end
    if ~isempty(config)
      warning('config will overide any units descriptions from the file')
    end
  end
  try
    if xdiff == 0
      warning('data array must have a constant sampling interval');
    end
  catch
    if strcmp(xdiff,'variable')
      warning('data array must have a constant sampling interval');
    end
  end

  % Assign variables from data
  N = size(data,1);
  RecordTime = range(data(:,1));
  sample_rate = round((N-1)/RecordTime);
  tau_rise = TC(1);
  tau_decay = TC(2);
  % Concatenate traces
  numwave = numel(wave);
  if numel(wave) > 1
    Trace = reshape(data(:,wave),numwave*N,1);
    T = (0:numwave*N-1)'*1/sample_rate;
    t = T;
  else
    Trace = data(:,wave);
    T = (0:N-1)'*1/sample_rate;
    t = data(:,1);
  end
  clear data

  % Create model template event
  Template = -exp(-T/tau_rise)+exp(-T/tau_decay);
  TemplatePeak = max(Template);
  Template = Template/TemplatePeak;
  if s=='-'
    Template = Template*-1;
  end

  % Perform fourier transform-based deconvolution (pad ends to mask end effects)
  Template_FFT = fft(padarray(Template,2*sample_rate,0,'post'));
  Trace_FFT = fft(bounce(Trace,sample_rate));
  DEC = real(ifft(Trace_FFT./Template_FFT));
  clear Trace_FFT Template_FFT
  DEC(1:sample_rate) = [];
  DEC(end-sample_rate+1:end) = [];

  %% Band-pass filter the deconvoluted trace (default is 1-200 Hz)
  %DEC = filter1 (DEC, t, hpf, inf, 'median'); % high pass median filter
  %DEC = filter1 (DEC, t, 0, lpf, 'binomial'); % low pass binomial filter
  
  % Band-pass filter the deconvoluted trace (default is 1-200 Hz)
  DEC = hpfilter(DEC,t,hpf);
  DEC = lpfilter(DEC,t,lpf);

  % Assign NaN to deconvoluted waves for values inside user-defined exclusion zones
  % Calculate actual recording time analysed (not including exclusion zones)
  if ~isempty(excl)
    excl_idx(:,1) = dsearchn(t,excl(:,1));
    excl_idx(:,2) = dsearchn(t,excl(:,2));
  end
  for i=1:size(excl,1)
    DEC(excl_idx(i,1):excl_idx(i,2)) = NaN;
  end
  AnalysedTime = (sum(~isnan(DEC))-1)/sample_rate;

  % Create all-point histogram with optimal bin width determined by modified Silverman's
  % rule of thumb on the putative noise component based on robust statistics.
  % References:
  % -Hardle, W. SMOOTHING TECHNIQUES, with implementations in S. Springer, New York. 1991
  % -Silverman, B.W. DENSITY ESTIMATION FOR STATISTICS AND DATA ANALYSIS.
  %  Monographs on Statistics and Applied Probability, London: Chapman and Hall, 1986.
  noise = DEC(~isnan(DEC));
  lim = abs(min(noise));
  noise(noise>lim) = [];
  Q2 = median(noise);
  Q1 = median(noise(noise<Q2));
  Q3 = median(noise(noise>Q2));
  IQR = abs(Q1-Q3);
  h = 0.79*IQR*numel(noise)^(-1/5);
  binwidth = 2*h;
  bins = ceil(range(DEC)/binwidth);
  [counts,x] = hist(DEC,bins);

  % Find the center of the noise peak
  peak_count = max(counts);
  mu = x(counts==max(counts));

  % Use linear interpolation between histogram bins to find half-maximum and
  % extrapolate the FWHM from the bottom half of the noise peak. Use this to
  % approximate the standard deviation of the noise peak.
  top = x(counts>max(counts)/2);
  idx = ones(1,2)*find(x==top(1));
  idx(1) = idx(1)-1;
  LL = interp1(counts(idx),x(idx),peak_count/2,'linear','extrap');
  FWHM = abs((LL-mu))*2;
  sigma = FWHM/(2*sqrt(2*log(2)));

  % Center and scale the noise component of the deconvoluted data using initial estimates
  DEC = (DEC-mu)/sigma;
  x = (x-mu)/sigma;

  % Least-squares Gaussian fitting of the noise peak using Levenberg-Marquardt
  % algorithm with the initial estimates for the fit parameters (data is scaled
  % for best performance)
  xdata = x((x>-4)&(x<1))';
  ydata = counts((x>-4)&(x<1))'/max(counts);
  fun1 = @(p,xdata)p(1)*p(3)*sqrt(2*pi)*normpdf(xdata,p(2),p(3));
  p0 = [1;0;1];
  optimoptions = struct;
  optimoptions.Algorithm = {'levenberg-marquardt',lambda};
  optimoptions.TolX = eps;
  optimoptions.TolFun = eps;
  [p,resnorm,residual,exitflag] = lsqfit(fun1,p0,xdata,ydata,[],[],optimoptions);
  clear xdata ydata
  if exitflag<1
    error('Noise peak fitting failed to reach tolerance. Try a higher lambda value.')
  end

  % Recenter and scale the noise component of the deconvoluted data using the
  % optimization result
  DEC = (DEC-p(2))/p(3);
  x = (x-p(2))/p(3);

  % Scan superthreshold deconvoluted wave for local maxima (vectorized)
  N = numel(Trace); % Number of sample points in cropped wave
  Event = zeros(N,1);
  y = DEC.*(DEC>SF);
  Event(2:N-1) = diff(sign(diff(y)))==-2;

  % Exclude events too close to the ends of the event wave
  samples_pre = round(abs(win(1))*sample_rate);
  samples_post = round(win(2)*sample_rate);
  Event_idx = find(Event);
  Event(Event_idx(Event_idx<=samples_pre)) = 0;
  Event(Event_idx(Event_idx>=N-samples_post)) = 0;

  % Create episodic data
  Event_idx = find(Event);
  idx_pre = Event_idx-samples_pre;
  idx_post = Event_idx+samples_post;
  n = sum(Event);
  t_events(:,1) = win(1):1/sample_rate:win(2);
  y_events = NaN(length(t_events),n);
  if n > 0
    for i=1:n
      clear idx
      idx = idx_pre(i):idx_post(i); idx=idx(:);
      y_events(:,i) = Trace(idx);
    end
  end
  Event_time = t(Event_idx);

  if n > 0
    % Fit template onto events using linear least squares
    SCALE = ones(1,n);
    OFFSET = zeros(1,n);
    % Time-to-peak = -log(tau_decay/tau_rise)/(1/tau_decay-1/tau_rise)
    if s=='+'
      idx = find(Template==max(Template));
    elseif s=='-'
      idx = find(Template==min(Template));
    end
    templateTime = idx+round(sample_rate*tau_decay*taus);
    numelBase = base*sample_rate;
    A = [zeros(numelBase,1); Template(1:templateTime)];
    l = length(A);
    A = [ones(l,1),A];
    start = dsearchn(t_events,base*-1);
    t_fit = t_events(start:start+numelBase+templateTime-1);
    y_fit = y_events(start:start+numelBase+templateTime-1,:);
    tzero = dsearchn(t_fit,0);
    for i=1:n
      % Perform linear least-squares fit using singular value decomposition on events
      kmin = sinv(A'*A)*A'*y_fit(:,i);
      OFFSET(i) = kmin(1);
      SCALE(i) = abs(kmin(2)); % Absolute value constraint: The template cannot be inverted
      y_template(:,i) = SCALE(i)*A(:,2)+OFFSET(i);
      temp = diag(rot90(corrcoef(y_template(:,i),y_fit(:,i)))); % Matlab/Octave compatibility
      r(i) = temp(1);
    end
    y_events = y_events-ones(length(t_events),1)*OFFSET;

    % Discard events that are poorly correlated with the fitted template (r < rmin)
    ridx = find(r<rmin);
    Event(ridx) = 0; %#ok<NASGU> Keep in case required for diagnostic tests
    Event_idx(ridx') = [];
    Event_time(ridx') = [];
    SCALE(ridx) = [];
    OFFSET(ridx) = []; %#ok<NASGU> Keep in case required for diagnostic tests
    y_events(:,ridx) = [];
    y_fit(:,ridx) = [];
    y_template(:,ridx) = [];
    n = numel(Event_idx);
    y_avg = mean(y_fit,2);

    % Create wave of fitted templates to overlay onto the event wave
    % Note that the baseline period is not plotted
    FITtrace = nan(N,1);
    l = length(t_fit(tzero:end));
    for i=1:n
      FITtrace(Event_idx(i):Event_idx(i)+(l-1)) = y_template(tzero:end,i);
      FITtrace(Event_idx(i)-1) = NaN;
    end
  end

  % Figure 1: All-points histogram of the deconvoluted data points
  clear h
  h1 = figure(1);
  bar(x,counts,1,'EdgeColor','b','Facecolor','b');
  hold on; plot(x,max(counts)*fun1([p(1);0;1],x),'r','linewidth',3);
  ylimits=ylim;
  plot(SF*ones(1,2),ylimits,'g','linewidth',2)
  hold off;
  ylim(ylimits); xlim([min(DEC),max(DEC)]);
  box('off'); grid('off');
  title('All-points histogram of the deconvoluted wave (truncated at -5 and +20 SD)');
  xlabel('Deconvoluted wave (SD)');
  ylabel('Number of points');
  xlim([-5,20]);  % Set x-axis limits to -5 and +20 SD of the baseline noise for clarity

  % Figure 2: Filtered deconvoluted wave (blue) and detection threshold (green)
  h2 = figure(2);
  plot(t,DEC,'b');
  hold on;
  plot(t,SF*ones(N,1),'g','linewidth',2);
  hold off;
  xlim([min(t),max(t)]);
  grid('off'); box('off');
  title('Deconvoluted wave');
  ylabel('Deconvoluted wave (SD)');
  xlabel('Time (s)');

  % Print result summary and escape from the function if no events were detected
  if n == 0
    % Print basic information
    fprintf(['--------------------------EVENTER---------------------------\n',...
             ' Automatic PSC/PSP detection using FFT-based deconvolution\n',...
             ' and event peak analysis by least-squares template fitting\n',...
             ' Version v1.0 Copyright 2014 Andrew Charles Penn                \n\n'])
    diary('on');
    fprintf('Filename: %s\n',filename);
    fprintf('Channel: %d\n',channel);
    fprintf('Wave number: %d\n',wave-1);
    fprintf('Wave name: %s\n',names{wave});
    fprintf('Total number of events detected: %d\n',n);
    fprintf('Duration of recording analysed (in s): %.1f\n',AnalysedTime);
    fprintf('High-pass filter cut-off on deconvoluted wave (at -3 dB, in Hz): %.3g\n',hpf);
    fprintf('Low-pass filter cut-off on deconvoluted wave (at -3 dB, in Hz): %.3g\n',lpf);
    fprintf('Vector of model template time constants (in s): [%.3g,%.3g]\n',TC);
    fprintf('Standard deviation of the noise of the deconvoluted wave (a.u.): %.3g\n',sigma*p(3));
    fprintf('Scale factor of noise standard deviations for threshold setting: %.3g\n',SF);
    fprintf('False positive event rate (in Hz): %.3g\n',(1-normcdf(SF))*sample_rate);
    fprintf('Sign of the event peaks: %s\n',s);
    fprintf('Sample rate of the recording (in kHz): %d\n',sample_rate*1e-03);
    fprintf('lsqfit exitflag for fitting the noise peak: %d\n',exitflag);
    format short g
    fprintf('Exclusion zones:\n');disp(excl);
    diary('off');

    % Save episodic data waves to file
    if exist('eventer.output','dir')==0
      mkdir('eventer.output');
    end
    cd eventer.output
    if exist(filewave,'dir')==0
      mkdir(filewave);
    end
    chdir(filewave);

    % Save detection parameters and event amplitudes and times to text file
    dlmwrite('_parameters',[TC';sigma*p(3);AnalysedTime;sample_rate])
    dlmwrite('_offset',win(1));
    if exist('summary.txt','file')~=0
      delete('summary.txt');
    end
    movefile('../../diary','summary.txt');
    cd ../..

    merge_data(average,s,win,export,optimoptions,cwd);

    return

  end

  % Calculate axes autoscaling for figure 3
  y_autoscale = 0.1*max(range(y_events));
  y_maxlim = max(max(y_events))+y_autoscale;
  y_minlim = min(min(y_events))-y_autoscale;
  % Figure 3: Identified events (black) aligned and overlayed with the mean event (red)
  h3 = figure(3);
  ylimits = [y_minlim,y_maxlim]; % Encoded y-axis autoscaling
  plot(t_events,y_events,'color',[0.85,0.85,0.85]);
  hold on;plot(t_events,mean(y_events,2),'-k');
  xlim(win);
  ylim([ylimits(1),ylimits(2)]);
  if strcmp(yunit,'A')
    title('Mean EPSC');
  elseif strcmp(yunit,'V')
    title('Mean EPSP');
  end
  xlabel('Time (s)');
  if strcmp(yunit,'A')
    ylabel('Current (A)');
  elseif strcmp(yunit,'V')
    ylabel('Voltage (V)');
  end
  box('off'); grid('off');
  hold off;

  % Figure 4: Overlay of model template (green) and mean event (red)
  h4 = figure(4);
  plot(t_fit,y_avg,'-g','linewidth',3);
  hold on;plot(t_fit,mean(y_template,2),'-r','linewidth',2);
  hold off;
  xlim([t_fit(1) t_fit(end)]);
  ylim([min(y_avg)-range(y_avg)*0.1 max(y_avg)+range(y_avg)*0.1]);
  xlabel('Time (s)');
  if strcmp(yunit,'A')
    ylabel('Current (A)');
  elseif strcmp(yunit,'V')
    ylabel('Voltage (V)');
  end
  grid('off'); box('off');
  title('Template and Mean Event Overlay')

  % Calculate axes autoscaling for figure 5
  y_autoscale = 0.1*range(Trace);
  y_maxlim = max(Trace)+y_autoscale;
  y_minlim = min(Trace)-y_autoscale;
  % Figure 5: Event wave overlaid with template fits
  h5 = figure(5);
  ylimits = [y_minlim y_maxlim]; % Encoded y-axis autoscaling
  plot(t,Trace,'-','color',[0.85,0.85,0.85]);
  hold on; plot(t,FITtrace,'-r'); hold off;
  xlim([min(t),max(t)]);
  ylim([ylimits(1), ylimits(2)]);
  grid('off'); box('off');
  if strcmp(yunit,'A')
    title('EPSC wave');
  elseif strcmp(yunit,'V')
    title('EPSP wave');
  end
  xlabel('Time (s)');
  if strcmp(yunit,'A')
    ylabel('Current (A)');
  elseif strcmp(yunit,'V')
    ylabel('Voltage (V)');
  end

  % Calculate statistics
  ET = Event_time;
  IEI = [NaN; diff(ET)];
  if exmode==1
    % Do nothing
  elseif exmode==2
    % Convert interevent intervals corresponding to events spanning exclusion zones to NaN
    if ~isempty(excl)
      for i=1:size(excl,1)
        temp = find(ET>=excl(i,2));
        if ~isempty(temp)
          IEI(temp(1)) = NaN;
        end
      end
    end
  else
    error('exclusion method not recognised')
  end
  Amplitude = mean(SCALE);

  % Print basic information
  fprintf(['--------------------------EVENTER---------------------------\n',...
           ' Automatic PSC/PSP detection using FFT-based deconvolution\n',...
           ' and event peak analysis by least-squares template fitting\n',...
           ' Version v1.0 Copyright 2014 Andrew Charles Penn                \n\n'])
  diary('on');
  fprintf('Filename: %s\n',filename);
  fprintf('Channel: %d\n',channel);
  fprintf('Wave number: %d\n',wave-1);
  fprintf('Wave name: %s\n',names{wave});
  fprintf('Total number of events detected: %d\n',n);
  fprintf('Duration of recording analysed (in s): %.1f\n',AnalysedTime);
  if strcmp(yunit,'A')
    fprintf('Mean event amplitude (pA): %.3g\n',Amplitude*1e+12);
  elseif strcmp(yunit,'A')
    fprintf('Mean event amplitude (mV): %.3g\n',Amplitude*1e+3);
  end
  fprintf('Event frequency (in Hz): %.3g\n',n/AnalysedTime);
  fprintf('High-pass filter cut-off on deconvoluted wave (at -3 dB, in Hz): %.3g\n',hpf);
  fprintf('Low-pass filter cut-off on deconvoluted wave (at -3 dB, in Hz): %.3g\n',lpf);
  fprintf('Vector of model template time constants (in ms): [%.3g,%.3g]\n',TC);
  fprintf('Standard deviation of the noise of the deconvoluted wave (a.u.): %.3g\n',sigma*p(3));
  fprintf('Scale factor of noise standard deviations for threshold setting: %.3g\n',SF);
  fprintf('False positive event rate (in Hz): %.3g\n',(1-normcdf(SF))*sample_rate);
  fprintf('Sign of the event peaks: %s\n',s);
  fprintf('Number of decay time constants used in the template fit: %.3g\n',taus);
  fprintf('Minimum acceptable correlation coefficient for the template fit: %.3g\n',rmin);
  fprintf('Episodic data window limits centred around each event: [%.3g,%.3g]\n',win);
  fprintf('Sample rate of the recording (in kHz): %d\n',sample_rate*1e-03);
  fprintf('lsqfit exitflag for fitting the noise peak: %d\n',exitflag);
  format short g
  fprintf('Exclusion zones:\n');disp(excl);
  diary('off');

  % Save episodic data waves to file
  event_data = cat(2,sample_rate^-1*[0:size(y_events,1)-1]',y_events);
  ensemble_mean = [sample_rate^-1*[0:size(y_events,1)-1]',mean(y_events,2)];
  if exist('eventer.output','dir')==0
    mkdir('eventer.output');
  end
  cd eventer.output
  if exist(filewave,'dir')==0
    mkdir(filewave);
  end
  chdir(filewave);
  state = warning('query');
  warning off
  ephysIO ('event_data.mat',event_data,xunit,yunit);
  ephysIO ('ensemble_mean.mat',ensemble_mean,xunit,yunit);
  warning(state)
  if ~strcmpi(export,'none')
    ephysIO (sprintf('event_data.%s',export),event_data,xunit,yunit);
    ephysIO (sprintf('ensemble_mean.%s',export),ensemble_mean,xunit,yunit);
  end

  % Save detection parameters and event amplitudes and times to text file
  dlmwrite('_parameters',[TC';sigma*p(3);AnalysedTime;sample_rate])
  dlmwrite('_offset',win(1));
  if exist('summary.txt','file')~=0
    delete('summary.txt');
  end
  movefile('../../diary','summary.txt');
  if exist('txt','dir')==0
    mkdir('txt');
  end
  chdir('txt');
  dlmwrite('times.txt',Event_time,'\t');
  dlmwrite('peak.txt',SCALE','\t');
  dlmwrite('IEI.txt',IEI,'\t');
  cd ..

  % Save figures and images
  if exist('fig','dir')==0
    mkdir('fig');
  end
  chdir('fig');
  saveas(h1,'output_histogram.fig','fig');
  saveas(h2,'output_decon.fig','fig');
  saveas(h3,'output_avg.fig','fig');
  saveas(h4,'output_template.fig','fig');
  saveas(h5,'output_wave.fig','fig');
  cd ..
  %if exist('img','dir')==0
  %  mkdir('img');
  %end
  %chdir('img');
  %if exist('eps','dir')==0
  %  mkdir('eps');
  %end
  %chdir('eps');
  %print(1,'output_histogram.eps','-depsc');
  %print(2,'output_decon.eps','-depsc');
  %print(3,'output_avg.eps','-depsc');
  %print(4,'output_template.eps','-depsc');
  %print(5,'output_wave.eps','-depsc');
  %cd ..
  %if exist('png','dir') == 0
  %  mkdir('png');
  %end
  %chdir('png');
  %print(1,'output_histogram.png','-dpng');
  %print(2,'output_decon.png','-dpng');
  %print(3,'output_avg.png','-dpng');
  %print(4,'output_template.png','-dpng');
  %print(5,'output_wave.png','-dpng');
  %cd ..
  cd ../..

  merge_data(average,s,win,export,optimoptions,cwd);

%----------------------------------------------------------------------

function merge_data(average,s,win,export,optimoptions,cwd)

  % Merge the event data for all the waves analysed in this experiment

  % Merge all data from eventer.output folder
  count = 0;
  cd eventer.output
  numtraces=0;
  dirlist = dir('.');
  dirlist = char(dirlist.name);
  numdir = size(dirlist,1);
  dirarray = mat2cell(dirlist,ones(numdir,1),size(dirlist,2));
  for i=1:numdir
    dirname{i} = strtrim(dirarray{i});
    if isdir(dirname{i}) && ~strcmp(dirname{i},'ALL_events') &&...
    ~strcmp(dirname{i},'.') &&...
    ~strcmp(dirname{i},'..')
      count = count+1;
      chdir(dirname{i});
      t = [];
      if exist('event_data.mat','file')
        [temp,xdiff,xunit,yunit] = ephysIO('event_data.mat');
        temp(:,1) = temp(:,1)+win(1);
        temp(abs(temp(:,1))<eps,1)=0;
        if isempty(t)
          t = temp(:,1);
        elseif count>1
          if numel(t)==numel(temp(:,1))
             if all(t==temp(:,1))~=1
               error('Inconsistent window dimensions. Cannot merge data.');
             end
          else
            error('Inconsistent window dimensions. Cannot merge data.');
          end
        end
        temp(:,1) = [];
        data{count} = temp;
        numtraces(count,1) = size(data{count},2);
      else
        numtraces(count,1) = 0;
      end
      parameters(:,count) = load('-ascii','_parameters');
      TC(count,:) = parameters(1:2,count)';
      sigma(count,1) = parameters(3,count);
      AnalysedTime(count,1) = parameters(4,count);
      sample_rate(count,1) = parameters(5,count);
      if sample_rate(count)~=sample_rate(1)
        error('Inconsistent sample rate. Cannot merge data.')
      end
      if exist('txt','dir')
        cd txt
        IEI{i,1} = load('-ascii','IEI.txt');
        peak{i,1} = load('-ascii','peak.txt');
        cd ..
      end
      cd ..
    else
      % The name in the directory list is not a folder suitable for the
      % merge process so do nothing
    end
  end
  numEvents = sum(numtraces);
  if exist('ALL_events','dir')==0
    mkdir('ALL_events');
  end
  cd('ALL_events')
  if numEvents > 1
    tau_rise = sum(TC(:,1).*numtraces/sum(numtraces));
    tau_decay = sum(TC(:,2).*numtraces/sum(numtraces));
    y = cell2mat(data);
    IEI = cell2mat(IEI);
    peak = cell2mat(peak);
    nanidx = isnan(IEI);
    peak(nanidx) = [];
    IEI(nanidx) = [];
    y(:,find(nanidx)) = [];
    numEvents = numel(peak);
    freq = 1/median(IEI);
    events = [sample_rate(1)^-1*[0:size(y,1)-1]',y];
    if strcmp(average,'mean')
      y_avg = mean(y,2);
    elseif strcmp(average,'median')
      y_avg = median(y,2);
    end
    ensemble_average = [sample_rate(1)^-1*[0:size(y_avg,1)-1]',y_avg];

    % Fit a sum of exponentials function to the ensemble average event
    p0 = [1,tau_rise,tau_decay];
    tpeak0 = p0(3)*p0(2)/(p0(3)-p0(2))*log(p0(3)/p0(2));
    idx = t>=0 & t<=tpeak0+p0(3);
    tdata = t(idx);
    fun2 = @(p,tdata)-exp(-tdata/p(1))+exp(-tdata/p(2));
    ypeak0 = fun2([p0(2),p0(3)],tpeak0);
    if s=='-'
      NF = ypeak0/min(y_avg(idx));  % Normalization factor
      ydata = y_avg(idx)*NF;
    elseif s=='+'
      NF = ypeak0/max(y_avg(idx));  % Normalization factor
      ydata = y_avg(idx)*NF;
    end
    fun3 = @(p,tdata)p(1)*(-exp(-tdata/p(2))+exp(-tdata/p(3)));
    [p,resnorm,residual,exitflag] = lsqfit(fun3,p0,tdata,ydata,[],[],optimoptions);
    tpeak = p(3)*p(2)/(p(3)-p(2))*log(p(3)/p(2));
    fitAmplitude = abs((fun3(p,tpeak))/NF);
    fitIntegral = abs(p(1)*(p(3)-p(2))/NF);
    fit = fun3(p,tdata)/NF;
    residuals = y_avg(idx)-fit;

    % Calculate axes autoscaling for figure 6
    y_autoscale = 0.1*max(range(y));
    y_maxlim = max(max(y))+y_autoscale;
    y_minlim = min(min(y))-y_autoscale;
    % Figure 6: Identified events aligned and overlayed with the ensemble average event
    h6 = figure(6);
    ylimits = [y_minlim y_maxlim]; % Encoded y-axis autoscaling
    plot(t,y,'color',[0.85,0.85,0.85]);
    hold on;
    plot(t,y_avg,'-b');
    if s=='-'
      plot(tdata,fit,'r-','linewidth',2);
    elseif s=='+'
      plot(tdata,fit,'r-','linewidth',2);
    end
    xlim(win);
    ylim([ylimits(1), ylimits(2)]);
    if strcmp(yunit,'A')
      title('Merged EPSC data: ensemble average (blue) and fit (red)');
    elseif strcmp(yunit,'V')
      title('Merged EPSP data: ensemble average (blue) and fit (red)');
    end
    xlabel('Time (s)');
    if strcmp(yunit,'A')
      ylabel('Current (A)');
    elseif strcmp(yunit,'V')
      ylabel('Voltage (V)');
    end
    box('off'); grid('off');
    hold off;

    % Save data
    save('ensemble_average.txt','ensemble_average','-ascii','-tabs');
    save('fit.txt','fit','-ascii','-tabs');
    save('residuals.txt','residuals','-ascii','-tabs');
    ephysIO ('event_data.mat',events,xunit,yunit);
    ephysIO ('ensemble_average.mat',ensemble_average,xunit,yunit);
    if ~strcmpi(export,'none')
      ephysIO (sprintf('event_data.%s',export),events,xunit,yunit);
      ephysIO (sprintf('ensemble_average.%s',export),ensemble_average,xunit,yunit);
    end
    dlmwrite('_parameters',[tau_rise;tau_decay;sigma;AnalysedTime]);
    dlmwrite('_offset',win(1));
    if exist('img','dir')==0
      mkdir('img');
    end
    cd('img')
    if exist('fig','dir')==0
      mkdir('fig');
    end
    cd('fig')
    saveas(h6,'output_avg.fig','fig');
    cd ..
    %if exist('png','dir')==0
    %  mkdir('png');
    %end
    %cd('png')
    %print(6,'output_avg.png','-dpng');
    %cd ..
    %if exist('eps','dir')==0
    %  mkdir('eps');
    %end
    %cd('eps')
    %print(6,'output_avg.eps','-depsc');
    %cd ..
    cd ..
    if exist('txt','dir')==0
      mkdir('txt');
    end
    cd('txt')
    dlmwrite('IEI.txt',IEI,'\t');
    dlmwrite('peak.txt',peak,'\t');
    cd ..
  end
  sigma = sum(sigma.*AnalysedTime/sum(AnalysedTime));
  AnalysedTime = sum(AnalysedTime);

  % Print basic information
  fprintf('------------------------------------------------------------\n');
  diary('on');
  fprintf('Number of analyses merged: %d\n',count);
  fprintf('Total recording time analysed (in s): %.1f\n',AnalysedTime);
  fprintf('Total number of events: %d\n',numEvents);
  fprintf('Event frequency (in Hz): %.3g\n',numEvents/AnalysedTime)
  if numEvents > 1
    fprintf('Ensemble average: %s\n',average);
    if strcmp(yunit,'A')
      fprintf('Amplitude of the model EPSC fit (pA): %.3g\n',fitAmplitude*1e+12);
      fprintf('Integral (charge) of the model EPSC fit (fC): %.4g\n',fitIntegral*1e+15);
      fprintf('Rise time constant of the model EPSC fit (ms): %.3g\n',p(2)*1e+03);
      fprintf('Decay time constant of the model EPSC fit (ms): %.3g\n',p(3)*1e+03);
    elseif strcmp(yunit,'V')
      fprintf('Amplitude of the model EPSP fit (mV): %.3g\n',fitAmplitude*1e+03);
      fprintf('Rise time constant of the model EPSP fit (ms): %.3g\n',p(2)*1e+03);
      fprintf('Decay time constant of the model EPSP fit (ms): %.3g\n',p(3)*1e+03);
    end
      fprintf('lsqfit exitflag for fitting the ensemble average event: %d\n',exitflag);
  else
    fprintf('Note: Not enough events for analysis\n');
  end
  fprintf('Standard deviation of the noise of the deconvoluted waves (a.u.): %.3g\n',sigma);
  diary('off');
  movefile('./diary','summary.txt');
  cd ../..
  chdir(cwd)
