%     Script file: peaker
%
%     Searches for peaks in smoothed episodic trace data that exceed a user-
%     defined amplitude threshold. The local preceding baseline for each
%     peak is determined and used to calculate relative peak amplitides and
%     rising edge statistics.
%
%     The amplitude threshold is set by left-button mouse clicking the cursor
%     on the preview graph at the desired y-value(s) and then pressing RETURN.
%     From a single coordinate, a constant threshold is set across the whole
%     trace. If multiple coordinates are entered then a dynamic threshold is
%     constructed by linear interpolation. If no coordinates are entered,
%     then a constant threshold is automatically set to 20 percent of the
%     maximum peak amplitude. This latter feature requires the trace data
%     to have little or no constant offset. Finally, a trace can be rejected
%     by right-button mouse clicking the cursor on the preview graph.
%
%     Peaker can tolerate moderate levels of noise, particularly with
%     appropriate use of data smoothing. As a general rule of thumb, for a
%     new data set iteratively increase the amount of smoothing until the
%     script reliably detects the baseline of each peak. This will occur
%     when the noise is removed from the rising edge of the peaks. During
%     the execution of the valley algorithm, the constant baseline offset
%     is assumed from the baseline point of the first peak.
%
%     The valley algorithm is crucial for finding peaks in the presense
%     of noise. By default, the valley threshold is set to 0.5. This means
%     that a baseline value exceeding 50 percent of the smallest adjacent
%     candidate peak is discarded. The largest peak between geniune baseline
%     points is then selected. Usually, the default threshold is suitable.
%     However, for instances where geniune baseline points between peaks
%     exceeds this threshold, such as with extensive temporal summation,
%     then the valley threshold should be raised accordingly.
%
%     Peak and baseline coordinates are plotted as crosshairs.
%     Lower and upper bounds for risetimes are plotted as plus signs.
%     Half-amplitude decay time point is plotted as an empty circle.
%     The peak in the first derivative is plotted as an asterix.
%
%     When there is more than one trace, peaker can operate in batch mode.
%     In batch mode peaker exports data tables of the peak times, baseline
%     times, absolute amplitudes, relative amplitudes, risetimes, slopes
%     half-amplitude decay times, half-width and 20 percent risetime and
%     50 percent decay time points. In the data tables, each row corresponds
%     to a trace and each column to a peak. The table values have units
%     without prefixes. The value NaN is returned if the number of peaks
%     detected in a particular trace is not the same as for the first
%     trace or if the trace was rejected. The indices of traces that are
%     rejected are listed in the file named 'rejected'. The indices of
%     traces where peaker fails to find the same number of peaks as the
%     specified peak number reference are listed in the file named 'failed'.
%     The script can be rerun in normal mode to repeat peak detection on
%     individual traces and automatically update the data tables. When
%     analysis of a particular file is complete, it is strongly advised to
%     clear all global variables by entering at the command prompt (>>):
%
%     >> clear
%
%     To abandon viewing the output in the Octave terminal press 'q'. This
%     will not affect the saved output files. To terminate the script
%     prematurely, press 'CTRL' + 'C'. Be sure to clear the global variables
%     and check your working directory if you do this.
%
%     This script requires the following functions and their dependencies:
%     'ndiff', 'binomialf', 'sma', 'fpeaks' and 'edge'.
%
%     See the example distributed with this script.
%
%     Bibliography:
%     IGOR Pro Version 4, Vol II Users Guide, 2000. WaveMetrics, Inc.
%     Axograph 4.8 User Manual, 2002. Axon Instruments, Inc.
%
%     peaker v1.6 (last updated: 05/06/2016)
%     Author: Andrew Charles Penn
%     https://www.researchgate.net/profile/Andrew_Penn/
%


% Clear current variables, load data and define matrix columns
close all
diary('off');
if exist('diary','file')
 delete('diary');
end
diary('on');
format short g
 if exist('filename')==1 && exist(strcat('./peaker.output/',filename,'/tables/_parameters'))~=0
   clear BasePointsX Slopes gy slope BasePointsY i stl D Tr idx stl_ref G V intervals subtracted_baseline L batch syl N button k syl_ref Nref clamp l t P dsign m tl PeakPointsX dydt_ref n tl_ref PeakPointsY dydtlimits p smooth txtfilename PeakRelative edgeX imgfilename y Rise20 edgeY edge_dydx edge_d2ydx2 SlopePlots reject yl RisePlots gX ridx yl_ref RiseTimes gY s ylimits gx PT BT A R RT S R20 SR flist parameters scan j rlist n sample_rate HTF Fc scale_factor i_ref rise50 Decay50 decayX decayY HalfMax HalfWidth HW D50 DT DecayTime;
   disp(sprintf(strcat('\nData matrix filename (including extension): \n'))); disp(filename);
   if exist('data')~=1
   % if exist(strcat(filename,'.txt.gz')) == 2
   %  gunzip(strcat(filename,'.txt.gz'));
   % end
   % data=load('-ascii', strcat(filename,'.txt'));
   % gzip(strcat(filename,'.txt'));
   % delete(strcat(filename,'.txt'));
     data=ephysIO(filename);
   end
   cd(strcat('./peaker.output/',filename,'/tables/'));
   parameters=load('-ascii','_parameters');
   tl_ref=load('-ascii','_tl_ref');
   cd ../../..
 elseif exist('filename')==0 && exist('_parameters') && exist('_tl_ref')~=0
   parameters=load('-ascii','_parameters');
   tl_ref=load('-ascii','_tl_ref');
   cd ..
   [pathstr,name,ext]=fileparts(pwd);
   cd ../..
   filename=strcat(name,ext);
   [data,xdiff]=ephysIO(filename);
 else
   clear;
   filename=input(sprintf('\nData matrix filename (including extension): '),'s');
   % if exist(strcat(filename,'.txt.gz'),'file') == 2
   %  gunzip(strcat(filename,'.txt.gz'));
   % end
   % data=load('-ascii', strcat(filename,'.txt'));
   % gzip(strcat(filename,'.txt'));
   % delete(strcat(filename,'.txt'));
   [data,xdiff]=ephysIO(filename);
   if exist(strcat('./peaker.output/',filename,'/tables/'),'dir')
    cd(strcat('./peaker.output/',filename,'/tables/'))
    parameters=load('-ascii','_parameters');
    tl_ref=load('-ascii','_tl_ref');
    cd ../../..
   end
 end
% Timebase conversion if offset provided in file '_offset'
if exist('_offset','file')
  offset=load('-ascii','_offset');
  t=data(:,1)+offset;
else
  t=data(:,1);
end
y=data; y(:,1)=[];
n=size(y,2);
 if n > 1
  disp(sprintf('\nNumber of traces in this file: ')); disp(n);
  Tr=input(sprintf('\nSelect the trace number for this analysis (default is all): '));
   if isempty(Tr)
    batch=1;
    y=y(:,1);
   elseif length(Tr) == 1 && (Tr > 0) && (Tr <= n) && (Tr == round(Tr))
    batch=0;
    y=y(:,Tr);
    i=1; n=1;
   elseif length(Tr) > 1
    batch=1;
    n=length(Tr);
   else
    error('The trace number must be a finite positive integer <= the above value');
   end
 elseif n == 1
  batch=0;
  Tr=1; i=1; n=1;
 end

 for i=1:n
  if i > 1
   clear tl yl s l L m gx gy gX gY button reject ylimits dydt G N D P BasePointsX BasePointsY PeakPointsX PeakPointsY PeakRelative idx Rise20 Slopes RiseTimes RisePlots SlopePlots k m yl_ref stl_ref syl_ref dydtlimits i_ref rise50 Decay50 decayX decayY HalfMax HalfWidth DecayTime;
   disp('---------------------------------------------------------------------------');
   diary('on');
   disp(sprintf('\nData matrix filename (excluding extension): ')); disp(filename);
  end
  if length(Tr) > 1
   i_ref=i;
   i=Tr(i); %#ok<FXSET, the i value is saved as i_ref and is restored at the end of the loop >
   clamp=parameters(1);
   p=0;
   dsign=parameters(3);
   V=parameters(4);
   Nref=0;
   subtracted_baseline=parameters(6);
   reject=1;
  end
  if batch == 1
   disp(sprintf('\nTrace number: \n')); disp(i);
   y=data; y(:,1)=[];
   y=y(:,i);
  end
 if i > 1
  disp(sprintf('\nThe settings used are the same as those documented in the parameters file'));
 end

% Perform linear int1erpolation on non-evenly spaced datasets
l=numel(t);
if exist('sample_rate') == 0
 sample_rate=0.001*(l-1)/(max(t)-min(t));
end
 if xdiff == 0
  disp(sprintf('Input must consist of data sampled at evenly spaced time points.\nLinear interpolation will precede smoothing and differentiation.'));
  if i == 1
   if exist(strcat('./peaker.output/',filename,'/tables/_sample_rate')) == 0
    disp(sprintf('Current average sampling rate: \n')); disp(sample_rate);
    old_sample_rate=sample_rate;
    sample_rate=input(sprintf('\nEnter the desired sampling frequency of the recording (in kHz): '));
     if isempty(sample_rate)
      sample_rate=0.001*(l-1)/(max(t)-min(t));disp(sample_rate);   % Unit kHz
     end
   else
    cd(strcat('./peaker.output/',filename,'/tables/'));
    sample_rate=load('_sample_rate')/1000;
    old_sample_rate=sample_rate;
    cd ../../..
   end
  elseif i > 1
   disp(sprintf('\nThe sampling frequency of the recording used for linear interpolation (in kHz): \n')); disp(sample_rate);
  end
  if exist('tl_ref') == 0
   tl_ref=linspace(min(t),max(t),l);
   tl_ref=tl_ref(:);
   tl=tl_ref;
  elseif exist('tl_ref') == 1
   tl=tl_ref;
  end
  state = warning('query');
  warning off
  yl=interp1q(t,y,tl);
  warning(state)
  yl=yl(:);
  yl_ref=yl;
  scale_factor=sample_rate/old_sample_rate;
 else
  if exist('parameters') == 1
   if length(parameters) == 7
    scale_factor=parameters(7);
   end
  end
  if exist('scale_factor') == 0
   disp(sprintf('\nThe sampling frequency of the recording (in kHz): \n')); disp(sample_rate);
%  scale_factor=input(sprintf('\nEnter scaling factor for the desired sampling frequency (default is 1): '));
   scale_factor=1;
  end
   if isempty(scale_factor) || (scale_factor == 1)
    scale_factor=1;
    tl=t;
    yl=y;
    tl_ref=tl;
    yl_ref=yl;
%   disp(sprintf('\nNo upsampling was implemented.\n\nThe sampling frequency of the recording (in kHz): \n'));
   elseif ~isempty(scale_factor)
     if scale_factor < 1
      error('The scaling factor is for upsampling only');
     end
    l=1+scale_factor*(l-1);
    tl_ref=linspace(min(t),max(t),l);
    tl_ref=tl_ref(:);
    tl=tl_ref;
    state = warning('query');
    warning off
    yl=interp1(t,y,tl,'spline');
    yl_ref=yl;
    warning(state)
    disp(sprintf('\nUpsampling achieved with spline interpolation.\n\nThe new sampling frequency of the recording (in kHz): \n'));
   end
 sample_rate=0.001*(l-1)/(max(tl)-min(tl));
%disp(sample_rate);
 end

 if i == 1 && (length(Tr) <= 1)
  % Prompts for input parameters
  % Type of recording
  if exist('parameters') == 1
   clamp=parameters(1);
  else
   clamp=input(sprintf('\nAre the recordings voltage (0) or current (1) clamp? (default is 0): '));
    if isempty(clamp)
     clamp=0;
    end
    if (clamp ~= 0) && (clamp ~=1)
     error('The type of recording must be specified with a logical value');
    end
  end
  % Smoothing method
  % Currently, The option to use a boxcar filter is disabled.
  %if batch == 1
  % smooth=input(sprintf('\n(1) Binomial (binomialf)\n(2) Boxcar (sma)\nSelect which smoothing algorithm to use (default is binomial): '));
  %  if isempty(smooth)
  %   smooth=1;
  %  end
  %elseif batch == 0
  % smooth=input(sprintf('\n(1) Binomial (binomialf)\n(2) Boxcar (sma)\nSelect which smoothing algorithm to use (default is binomial): '));
  %   if isempty(smooth)
  %    smooth=1;
  %   end
  %end
  % Number of points for smoothing
  %if smooth == 1
   Fc=input(sprintf('\nInput the cut-off frequency of the binomial smoothing filter (-3 dB, in kHz): '));
   if isempty(Fc)
    p=[];
   else
    HTF=Fc/0.71;
    p=1/log2((cos(pi*(2*HTF/sample_rate)/2)^-2));
    p=round(p);
   end
  %elseif smooth == 2
  % Fc=input(sprintf('\nInput the cut-off frequency of the boxcar smoothing filter (-3 dB, in kHz): '));
  % p=0.5*(0.443/Fc*sample_rate-1);
  % p=round(p);
  %end
   if isempty(p)
    p=0;
   end
   if p > 0
     if isinf(p) || ~all(size(p) == 1) || p<0 || p~=abs(p) || p ~= round(p)
      if smooth == 1
       error('If nonzero, the filter order must be a positive integer');
      elseif smooth == 2
       error('If nonzero, the box size must be an odd number > 2 ');
      end
     end
   end
  % Set sign of peak deflections
  if exist('parameters') == 1 && parameters(3) ~= 0
   dsign=parameters(3);
  else
   dsign=input(sprintf('\nAre the peaks positive (1) or negative (-1) deflections? (default is auto): '));
  end
  if isempty(dsign)
   dsign=0;
  end
  if dsign == 1
   s='+';
  elseif dsign == -1
   s='-';
  end
  % Set the number of peaks expected per trace when in batch mode
  if (batch == 1) && isempty(Tr)
   Nref=input(sprintf('\nInput the number of peaks expected per trace: '));
    if isempty(Nref)==1 || isinf(Nref) || ~all(size(Nref) == 1) || Nref<0 || Nref~=round(Nref)
     error('Nref must be a non-negative integer');
    end
  end
  % Specify whether the input data is baseline-subtracted
  if (batch == 1) && isempty(Tr)
   subtracted_baseline=input(sprintf('\nIs the input data baseline subtracted? (1=yes, 0=no, default is yes): '));
    if isempty(subtracted_baseline)
     subtracted_baseline=1;
    end
    if (subtracted_baseline ~= 0) && (subtracted_baseline ~=1)
     error('The answer to the question must be specified with a logical value');
    end
  else
   subtracted_baseline=0;
  end
  % Set valley threshold
  V=input(sprintf('\nInput the fractional valley threshold (default is 0.5): '));
   if isempty(V)
	V=0.5;
   end
  % Set event feature scanning time window
  if (batch == 1) && (Nref == 1) && (subtracted_baseline == 1) && isempty(Tr) && ~isempty(regexp(filename,'event','once'))
   scan=input(sprintf('\nSet event feature scanning time window (in ms): '));
   scan=scan*1e-3;
    if isempty(scan)==1 || isinf(scan) || ~all(size(scan) == 1) || scan<0
     error('Scan time must be a non-negative number');
    end
  end
 end
 if exist('parameters')==1 && batch==0
  disp(sprintf('\nOther settings used are the same as those documented in the parameters file'));
 end

% Plot preview graph for setting threshold level(s)
if i == 1
  figure(2);
  clf(2);
end
figure(1);
clf(1);
if subtracted_baseline == 1
 close(1)
end
 if p>0
   %if smooth == 1
    [syl, stl, HTF]=binomialf(yl,tl,p,'on');   % Binomial smoothing with bounce end correction
    Fc=0.001*0.71*HTF;
    disp(sprintf('\nActual cut-off of the binomial filter (-3 dB, in kHz): \n')); disp(Fc);
   %elseif smooth == 2
   % [syl]=sma(yl,p,'mean','on');       % Box smoothing with bounce end correction
   % stl=tl;
   % Fc=0.443/(2*p+1)*sample_rate;
   % disp(sprintf('\nActual cut-off of the boxcar filter (-3 dB, in kHz): \n')); disp(Fc);
   %end
  stl_ref=stl;
  syl_ref=syl;
 elseif p==0
  syl=yl;
  stl=tl;
  stl_ref=stl;
  syl_ref=syl;
 end
dydt_ref=ndiff(syl,stl);
gX=[]; gY=[]; button=1;
if subtracted_baseline == 0 & length(Tr) <= 1
 y_autoscale=0.1*(max(max(yl))-min(min(yl)));
 y_maxlim=max(max(yl))+y_autoscale;
 y_minlim=min(min(yl))-y_autoscale;
 ylimits=[y_minlim y_maxlim]; % Encoded y-axis autoscaling
 plot(tl,yl,'-','color',[0.75,0.75,0.75]);
 hold on;
 plot(stl_ref,syl_ref,'k');
 plot([min(tl),max(tl)],[0, 0],'-','color',[0.9,0.9,0.9]);
 plot([0, 0],[ylimits(1), ylimits(2)],'-','color',[0.9,0.9,0.9]);
 hold off; grid('off');xlim([min(tl),max(tl)]); ylim([ylimits(1), ylimits(2)]); box('off');
 diary('off');
 disp(sprintf('\nSet threshold coordinates by clicking the mouse cursor on the preview graph'));
 diary('on');
 while (button ~= 13) & (button ~= 3)
 clear gx gy button
 [gx, gy, button]=ginput(1);
  if exist('m','var') == 0
   m = 1;
  elseif exist('m','var') == 1
   m = m + 1;
  end
  if button == 1
   gX(m,:) = gx;
   gY(m,:) = gy;
   hold on; plot(gx, gy,'ko','markersize',10); hold off;
  end
 end
end
if ~isempty(strfind(filename,'event')) & (batch == 1) & (subtracted_baseline == 1) & isempty(Tr)
  clear idx
  idx=find((tl_ref >= -scan) .* (tl_ref <= scan));
  tl=tl_ref(idx); yl=yl_ref(idx);
  l=numel(tl);
  clear idx
  idx=find((stl_ref >= -scan) .* (stl_ref <= scan));
  stl=stl_ref(idx); syl=syl_ref(idx);
end

if (length(Tr) <= 1)
 if button == 3
  reject = 1;
 else
  reject = 0;
 end
end

if reject == 0
 if dsign == 0
  % Automatic sign determination of the peak deflections
   if length(gY) > 1
    clear idx
    idx=find((tl_ref >= min(gX)) .* (tl_ref <= max(gX)));
    tl=tl_ref(idx); yl=yl_ref(idx);
    l=numel(tl);
    clear idx
    idx=find((stl_ref >= min(gX)) .* (stl_ref <= max(gX)));
    stl=stl_ref(idx); syl=syl_ref(idx);
    dydt=dydt_ref(idx(2):idx(end-1));
   else
    dydt=dydt_ref;
   end
   if isempty(find(dydt/max(abs(dydt)) == -1))
    s='+';
   elseif isempty(find(dydt/max(abs(dydt)) == 1))
    s='-';
   end
 end
 if dsign == 1
  s='+';
 elseif dsign == -1
  s='-';
 end


 % Create vector of threshold levels
 if length(gY) > 1
  G=[gX, gY];
  G=sortrows(G,1);
  L=interp1q(G(:,1),G(:,2),stl);
 elseif length(gY) == 1
  L=ones(l,1)*gY;
 elseif isempty(gY) == 1
  if s == '+'
   L=0.20*max(syl);
  elseif s == '-'
   L=0.20*min(syl);
  end
 L=ones(l,1)*L;
 end
 if s == '+'
  L=L/max(syl);
 elseif s == '-'
  L=L/min(syl);
 end

 % Plot threshold levels on preview graph
 if subtracted_baseline == 0
  hold on;
   if s == '+'
    plot(stl,L*max(syl),'k-');
   elseif s == '-'
    plot(stl,L*min(syl),'k-');
   end
  hold off;
 end
end

% Run fpeaks peak detection algorithm and report results
if reject == 0
  try
   [N, P, D]=fpeaks(stl,syl,0,s,L,V);
  catch
   N=0;
   BasePointsX=[];
   BasePointsY=[];
   PeakPointsX=[];
   PeakPointsY=[];
   PeakRelative=[];
  end
 disp(sprintf('\nNumber of peaks: \n')); disp(N);
 if N > 0
  BasePointsX=P.tBase;
  BasePointsY=P.yBase;
  PeakPointsX=P.tPeaks;
  PeakPointsY=P.yPeaks;
  PeakRelative=P.yAmpl;
  if (subtracted_baseline == 1) && ~isempty(strfind(filename,'event')) && N > 0
   clear ridx
   ridx=find((BasePointsX > 0) + (PeakPointsX < 0));
   BasePointsX(ridx)=[];
   BasePointsY(ridx)=[];
   PeakPointsX(ridx)=[];
   PeakPointsY(ridx)=[];
   PeakRelative(ridx)=[];
   N=length(PeakPointsX);
   if (N > 0) && (Nref == 1)
    clear idx
    idx=dsearchn(PeakPointsX,0);
    BasePointsX=BasePointsX(idx);
    BasePointsY=BasePointsY(idx);
    PeakPointsX=PeakPointsX(idx);
    PeakPointsY=PeakPointsY(idx);
    PeakRelative=PeakRelative(idx);
   end
   N=length(PeakPointsX);
  end
 disp(sprintf('\nPeak Time (in ms): \n')); disp(PeakPointsX*1e3);
  if clamp == 0
   disp(sprintf('\nPeak Current (in pA): \n')); disp(PeakPointsY*1e12);
  elseif clamp == 1
   disp(sprintf('\nPeak Voltage (in mV): \n')); disp(PeakPointsY*1e3);
  end
 disp(sprintf('\nBaseline Time (in ms): \n')); disp(BasePointsX*1e3);
  if clamp == 0
   disp(sprintf('\nBaseline Current (in pA): \n')); disp(BasePointsY*1e12);
  elseif clamp == 1
   disp(sprintf('\nBaseline Voltage (in mV): \n')); disp(BasePointsY*1e3);
  end
  if clamp == 0
   disp(sprintf('\nRelative Peak Current (in pA): \n')); disp(PeakRelative*1e12);
  elseif clamp == 1
   disp(sprintf('\nRelative Peak Voltage (in mV): \n')); disp(PeakRelative*1e3);
  end
 end
elseif reject == 1
 try
   [N, P, D]=fpeaks(stl,syl,0,s,L,V);
 catch
 end
 N = 0;
end
if batch == 0
 Nref=N;
end

% Run edge algorithm to calculate rising edge statistics and calculate half-width
if N == Nref && N > 0
 for k=1:Nref
  clear idx edgeX edgeY edge_dydx edge_d2ydx2 slope intervals P
  try
  idx=find((D(:,1) >= BasePointsX(k)) .* (D(:,1) <= PeakPointsX(k)));
     edgeX=D(idx,1);
     edgeY=D(idx,2);
    if s == '-'
     edgeY=edgeY*-1;   % make edge rising for negative deflections
    end
  [slope, intervals, P]=edge(edgeX,edgeY,0.2,0.8);	% By default this script calculates 20-80 percent rise time
  Rise20(k,:)=P(1,1);
  RiseTimes(k,:)=intervals(1);
   if s == '+'
    RisePlots(k*2-1:k*2,2)=P(:,2);
   elseif s == '-'
    RisePlots(k*2-1:k*2,2)=P(:,2)*-1;
   end
  RisePlots(k*2-1:k*2,1)=P(:,1);
  catch
   Rise20(k,:)=NaN;
   Slopes(k,:)=NaN;
   RiseTimes(k,:)=NaN;
   RisePlots(k*2-1:k*2,1)=NaN;
   RisePlots(k*2-1:k*2,2)=NaN;
  end
  [slope, intervals, P]=edge(edgeX,edgeY,0.5,1);
  Rise50(k,:)=P(1,1);
  if k < N
   idx=find((D(:,1) >= PeakPointsX(k)) .* (D(:,1) <= BasePointsX(k+1)));
  else
   idx=find((D(:,1) >= PeakPointsX(k)));
  end
  %decayX=D(idx,1);
  %decayY=D(idx,2);
  decayX=tl_ref(tl_ref > PeakPointsX(k));
  decayY=syl_ref(tl_ref > PeakPointsX(k));
  HalfMax(k,:)=(PeakRelative(k) * 0.5) + BasePointsY(k);
  if s == '+'
   if HalfMax(k) < min(decayY)
    Decay50(k,:)=NaN;
    HalfWidth(k,:)=NaN;
    DecayTime(k,:)=NaN;
   else
    Decay50(k,:)=decayX(find(decayY <= HalfMax(k),1));
    HalfWidth(k,:)=Decay50(k)-Rise50(k);
    DecayTime(k,:)=Decay50(k)-PeakPointsX(k);
   end
  elseif s == '-'
   if HalfMax(k) > max(decayY)
    Decay50(k,:)=NaN;
    HalfWidth(k,:)=NaN;
    DecayTime(k,:)=NaN;
   else
    Decay50(k,:)=decayX(find(decayY >= HalfMax(k),1));
    HalfWidth(k,:)=Decay50(k)-Rise50(k);
    DecayTime(k,:)=Decay50(k)-PeakPointsX(k);
   end
  end

 % Find the peak of the first derivative
 try
  idx=find((D(:,1) >= BasePointsX(k)) .* (D(:,1) <= PeakPointsX(k)));
  edgeX=D(idx,1);
  edgeY=D(idx,2);
  if s == '-'
   edgeY=edgeY*-1;   % make edge rising for negative deflections
  end
  [edge_dydx,edgeY,edgeX]=ndiff(edgeY,edgeX);
  %%% OPTION %%%
  %%% Code for the shortest latency peak of the first derivative
  %%[edge_d2ydx2,edge_dydx,edgeX]=ndiff(edge_dydx,edgeX);
  %%Slopes(k,:)=edge_dydx(find(edge_d2ydx2 <= 0,1));
  %%SlopePlots(k,1)=edgeX(find(edge_d2ydx2 <= 0,1));
  % Code for the maximum peak of the first derivative
  Slopes(k,:)=max(edge_dydx);
  SlopePlots(k,1)=edgeX(find(edge_dydx==max(edge_dydx)));
 catch
  Slopes(k,:)=NaN;
  SlopePlots(k,1)=NaN;
 end
 if s == '+'
  SlopePlots(k,2)=Slopes(k,:);
 elseif s == '-'
  SlopePlots(k,2)=Slopes(k,:)*-1;
 end

 end
 disp(sprintf('\nPeak 20-80 percent rise time (in ms): \n')); disp(RiseTimes*1e3);
 if clamp == 0
  disp(sprintf('\nPeak rising slope (in nA/ms): \n')); disp(Slopes*1e6);
 elseif clamp == 1
  disp(sprintf('\nPeak rising slope (in mV/ms): \n')); disp(Slopes);
 end
 disp(sprintf('\nHalf-amplitude decay time (in ms): \n')); disp(DecayTime*1e3);
 disp(sprintf('\nHalf-width (in ms): \n')); disp(HalfWidth*1e3);
 disp(sprintf('\n'));

 % Create data tables when running in batch mode
 if (batch == 1) && (N == Nref)
    for j=1:Nref
     PT(i,j)=PeakPointsX(j);
     BT(i,j)=BasePointsX(j);
     A(i,j)=PeakPointsY(j);
     R(i,j)=PeakRelative(j);
     RT(i,j)=RiseTimes(j);
     S(i,j)=Slopes(j);
     SR(i,j)=SlopePlots(j,1);
     R20(i,j)=Rise20(j);
     D50(i,j)=Decay50(j);
     DT(i,j)=DecayTime(j);
     HW(i,j)=HalfWidth(j);
    end
   flist(i,1)=0;
   rlist(i,1)=0;
 elseif (batch == 1) && (N ~= Nref)
	for j=1:Nref
     PT(i,j)=NaN;
     BT(i,j)=NaN;
     A(i,j)=NaN;
	 R(i,j)=NaN;
     RT(i,j)=NaN;
     S(i,j)=NaN;
     SR(i,j)=NaN;
     R20(i,j)=NaN;
     D50(i,j)=NaN;
     DT(i,j)=NaN;
     HW(i,j)=NaN;
    end
  if reject==0
   flist(i,1)=1;
   rlist(i,1)=0;
  end
 end

 % Plot graphs
 figure(2);
 subplot(2,1,1);
 y_autoscale=0.1*(max(dydt_ref)-min(dydt_ref));
 y_maxlim=max(dydt_ref)+y_autoscale;
 y_minlim=min(dydt_ref)-y_autoscale;
 dydtlimits=[y_minlim y_maxlim];
 plot([min(tl_ref), max(tl_ref)],[0, 0],'-','color',[0.9,0.9,0.9]);
 hold on;
 plot(stl_ref(2:end-1),dydt_ref,'-','color',[0.75,0.75,0.75]);
 if ~isempty(regexp(filename,'event','once'))
  plot([0, 0],[dydtlimits(1), dydtlimits(2)],'-','color',[0.9,0.9,0.9]);
 end
 plot(SlopePlots(:,1),SlopePlots(:,2),'k*','markersize',10);
 ylim([dydtlimits(1), dydtlimits(2)]);hold off;
  if clamp == 0
   ylabel('Smoothed first derivative (A/s)');
  elseif clamp == 1
   ylabel('Smoothed first derivative (V/s)');
  end
 xlabel('Time (s)'); grid('off'); xlim([min(tl_ref),max(tl_ref)]); box('off');
 subplot(2,1,2);
 y_autoscale=0.1*(max(yl_ref)-min(yl_ref));
 y_maxlim=max(yl_ref)+y_autoscale;
 y_minlim=min(yl_ref)-y_autoscale;
 ylimits=[y_minlim y_maxlim];
 plot(tl_ref,yl_ref,'-','color',[0.75,0.75,0.75]);
 hold on;
 plot(BasePointsX,BasePointsY,'ko','markersize',10);
 plot(PeakPointsX,PeakPointsY,'ko','markersize',10);
 plot(BasePointsX,BasePointsY,'k+','markersize',10);
 plot(PeakPointsX,PeakPointsY,'k+','markersize',10);
 plot(RisePlots(:,1),RisePlots(:,2),'k+','markersize',10);
 plot(Decay50,HalfMax,'ko','markersize',10);
 ylim([ylimits(1), ylimits(2)]); hold off;
 if clamp == 0
   ylabel('Current (A)');
  elseif clamp == 1
   ylabel('Voltage (V)');
  end
 xlabel('Time (s)');grid('off'); xlim([min(tl_ref),max(tl_ref)]); box('off');


 % Save data
 diary('off');
  if exist('peaker.output','dir') == 0
    mkdir('peaker.output');
  end
 cd peaker.output
  if exist(filename,'dir') == 0
    mkdir(filename);
  end
 cd(filename);
  if exist('img','dir') == 0
   mkdir('img');
  end
  if exist('txt','dir') == 0
   mkdir('txt');
  end
  if batch == 0
   imgfilename=strcat('Tr', num2str(Tr), '_', filename);
   txtfilename=strcat('Tr', num2str(Tr), '_', filename,'.txt');
  elseif  batch == 1
   imgfilename=strcat('Tr', num2str(i), '_', filename);
   txtfilename=strcat('Tr', num2str(i), '_', filename,'.txt');
  end
 cd txt
  if exist(txtfilename) ~= 0
   delete(txtfilename);
  end
 movefile('../../../diary',txtfilename);
 cd ../img
 if exist('eps','dir') == 0
   mkdir('eps');
 end
 cd eps
  if exist(strcat(imgfilename,'.eps')) ~= 0
   delete(strcat(imgfilename,'.eps'));
  end
 print(2,strcat(imgfilename,'.eps'),'-depsc');
 cd ..
 if exist('png','dir') == 0
   mkdir('png');
 end
 cd png
  if exist(strcat(imgfilename,'.png')) ~= 0
   delete(strcat(imgfilename,'.png'));
  end
 print(2,strcat(imgfilename,'.png'),'-dpng');
 cd ..
 cd ../../..

elseif N ~= Nref || N == 0
 % Plot graphs
 figure(2);
 clf(2);
 subplot(2,1,1);
 y_autoscale=0.1*(max(dydt_ref)-min(dydt_ref));
 y_maxlim=max(dydt_ref)+y_autoscale;
 y_minlim=min(dydt_ref)-y_autoscale;
 dydtlimits=[y_minlim y_maxlim]; % Encoded y-axis autoscaling
 plot([min(tl_ref), max(tl_ref)],[0, 0],'-','color',[0.9,0.9,0.9]);
 hold on;
 plot(stl_ref(2:end-1),dydt_ref,'-','color',[0.75,0.75,0.75]);
 if ~isempty(regexp(filename,'event','once'))
  plot([0, 0],[dydtlimits(1), dydtlimits(2)],'-','color',[0.9,0.9,0.9]);
 end
 ylim([dydtlimits(1),dydtlimits(2)]); hold off;
  if clamp == 0
   ylabel('Smoothed first derivative (A/s)');
  elseif clamp == 1
   ylabel('Smoothed first derivative (V/s)');
  end
 xlabel('Time (s)'); grid('off'); xlim([min(tl_ref),max(tl_ref)]); box('off');
 subplot(2,1,2);
 y_autoscale=0.1*(max(yl_ref)-min(yl_ref));
 y_maxlim=max(yl_ref)+y_autoscale;
 y_minlim=min(yl_ref)-y_autoscale;
 ylimits=[y_minlim y_maxlim]; % Encoded y-axis autoscaling
 plot(tl_ref,yl_ref,'-','color',[0.75,0.75,0.75]); ylim([ylimits(1),ylimits(2)]);
 if clamp == 0
  ylabel('Current (A)');
 elseif clamp == 1
  ylabel('Voltage (V)');
 end
 xlabel('Time (s)');grid('off'); xlim([min(tl_ref),max(tl_ref)]); box('off');
 if reject == 1
  disp(sprintf('\nTrace rejected\n'));
  flist(i,1)=0;
  rlist(i,1)=1;
 elseif reject == 0
  disp(sprintf('\nThe number of peaks detected did not match the specified peak number reference\n'));
  flist(i,1)=1;
  rlist(i,1)=0;
 end
 if batch == 1
  for j=1:Nref
     PT(i,j)=NaN;
	 BT(i,j)=NaN;
     A(i,j)=NaN;
     R(i,j)=NaN;
     RT(i,j)=NaN;
     S(i,j)=NaN;
     SR(i,j)=NaN;
     R20(i,j)=NaN;
     D50(i,j)=NaN;
     DT(i,j)=NaN;
     HW(i,j)=NaN;
  end
 end
   diary('off');
 if exist('peaker.output','dir') == 0
  mkdir('peaker.output');
 end
cd peaker.output
 if exist(filename,'dir') == 0
  mkdir(filename);
 end
cd(filename);
 if batch == 0
  imgfilename=strcat('Tr', num2str(Tr), '_', filename);
  txtfilename=strcat('Tr', num2str(Tr), '_', filename,'.txt');
 elseif  batch == 1
  imgfilename=strcat('Tr', num2str(i), '_', filename);
  txtfilename=strcat('Tr', num2str(i), '_', filename,'.txt');
 end
 if exist('txt','dir') == 0
  mkdir('txt');
 end
cd txt
 if exist(txtfilename) ~= 0
  delete(txtfilename);
 end
movefile('../../../diary',txtfilename);
cd ..
 if exist('img','dir') == 0
  mkdir('img');
 end
cd img
if exist('eps','dir') == 0
  mkdir('eps');
end
cd eps
if exist(strcat(imgfilename,'.eps')) ~= 0
 delete(strcat(imgfilename,'.eps'));
end
print(2,strcat(imgfilename,'.eps'),'-depsc');
cd ..
if exist('png','dir') == 0
 mkdir('png');
end
cd png
if exist(strcat(imgfilename,'.png')) ~= 0
 delete(strcat(imgfilename,'.png'));
end
print(2,strcat(imgfilename,'.png'),'-dpng');
cd ..
cd ../../..
end
if length(Tr) > 1
 i=i_ref; %#ok<FXSET, this resets i to the reference value saved earlier in the loop>
end
end
cd peaker.output
cd(filename);
 if (batch == 1) && isempty(Tr)
  mkdir('tables');
  cd tables
  save -ascii _ID filename
  dlmwrite('_sample_rate',sample_rate*1000)
  parameters=[clamp p dsign V Nref subtracted_baseline scale_factor]';
  dlmwrite('_parameters',parameters)
  dlmwrite('_offset',offset);
  dlmwrite('_tl_ref',tl_ref)
  cd ..
 end
 if (batch == 1) && isempty(Tr)
  cd tables
  save('peak_time.txt','PT','-ascii','-tabs');
  save('baseline_time.txt','BT','-ascii','-tabs');
  save('absolute.txt','A','-ascii','-tabs');
  save('relative.txt','R','-ascii','-tabs');
  save('risetime.txt','RT','-ascii','-tabs');
  save('slope.txt','S','-ascii','-tabs');
  save('steepest_rise.txt','SR','-ascii','-tabs');
  save('rise20.txt','R20','-ascii','-tabs');
  save('decay50.txt','D50','-ascii','-tabs');
  save('decay_time.txt','DT','-ascii','-tabs');
  save('halfwidth.txt','HW','-ascii','-tabs');
  dlmwrite('_results',[flist rlist],'\t');
  dlmwrite('failed.txt',find(flist),'\t');
  dlmwrite('rejected.txt',find(rlist),'\t');
  close all
  cd ..
 elseif exist('tables','dir') ~= 0
  cd tables
  if exist('peak_time.txt') ~= 0
   PT=load('-ascii','peak_time.txt');
   if size(PT,2) == Nref
    for j=1:size(PT,2)
     PT(Tr,j)=PeakPointsX(j);
    end
   elseif size(PT,2) ~= Nref
    for j=1:size(PT,2)
     PT(Tr,j)=NaN;
    end
   end
   delete('peak_time.txt');
   save('peak_time.txt','PT','-ascii','-tabs');
  end
  if exist('baseline_time.txt') ~= 0
   BT=load('-ascii','baseline_time.txt');
   if size(BT,2) == Nref
    for j=1:size(BT,2)
     BT(Tr,j)=BasePointsX(j);
    end
   elseif size(BT,2) ~= Nref
    for j=1:size(BT,2)
     BT(Tr,j)=NaN;
    end
   end
  delete('baseline_time.txt');
  save('baseline_time.txt','BT','-ascii','-tabs');
  end
  if exist('absolute.txt') ~= 0
   A=load('-ascii','absolute.txt');
   if size(A,2) == Nref
    for j=1:size(A,2)
     A(Tr,j)=PeakPointsY(j);
    end
   elseif size(A,2) ~= Nref
    for j=1:size(A,2)
     A(Tr,j)=NaN;
    end
   end
   delete('absolute.txt');
   save('absolute.txt','A','-ascii','-tabs');
  end
  if exist('relative.txt') ~= 0
   R=load('-ascii','relative.txt');
   if size(R,2) == Nref
    for j=1:size(R,2)
     R(Tr,j)=PeakRelative(j);
    end
   elseif size(R,2) ~= Nref
    for j=1:size(R,2)
     R(Tr,j)=NaN;
    end
   end
   delete('relative.txt');
   save('relative.txt','R','-ascii','-tabs');
  end
  if exist('risetime.txt') ~= 0
   RT=load('-ascii','risetime.txt');
   if size(RT,2) == Nref
    for j=1:size(RT,2)
     RT(Tr,j)=RiseTimes(j);
    end
   elseif size(RT,2) ~= Nref
    for j=1:size(RT,2)
     RT(Tr,j)=NaN;
    end
   end
   delete('risetime.txt');
   save('risetime.txt','RT','-ascii','-tabs');
  end
  if exist('slope.txt') ~= 0
   S=load('-ascii','slope.txt');
   if size(S,2) == Nref
    for j=1:size(S,2)
     S(Tr,j)=Slopes(j);
    end
   elseif size(S,2) ~= Nref
    for j=1:size(S,2)
     S(Tr,j)=NaN;
    end
   end
   delete('slope.txt');
   save('slope.txt','S','-ascii','-tabs');
  end
  if exist('steepest_rise.txt') ~= 0
   SR=load('-ascii','steepest_rise.txt');
   if size(SR,2) == Nref
    for j=1:size(SR,2)
     SR(Tr,j)=SlopePlots(j,1);
    end
    elseif size(SR,2) ~= Nref
     for j=1:size(SR,2)
      SR(Tr,j)=NaN;
     end
    end
   delete('steepest_rise.txt');
   save('steepest_rise.txt','SR','-ascii','-tabs');
  end
  if exist('rise20.txt') ~= 0
   R20=load('-ascii','rise20.txt');
   if size(R20,2) == Nref
    for j=1:size(R20,2)
     R20(Tr,j)=Rise20(j);
    end
   elseif size(R20,2) ~= Nref
    for j=1:size(R20,2)
     R20(Tr,j)=NaN;
    end
   end
   delete('rise20.txt');
   save('rise20.txt','R20','-ascii','-tabs');
  end
  if exist('decay50.txt') ~= 0
   D50=load('-ascii','decay50.txt');
   if size(D50,2) == Nref
    for j=1:size(D50,2)
     D50(Tr,j)=Decay50(j);
    end
   elseif size(D50,2) ~= Nref
    for j=1:size(D50,2)
     D50(Tr,j)=NaN;
    end
   end
   delete('decay50.txt');
   save('decay50.txt','D50','-ascii','-tabs');
  end
  if exist('decay_time.txt') ~= 0
   DT=load('-ascii','decay_time.txt');
   if size(DT,2) == Nref
    for j=1:size(DT,2)
     DT(Tr,j)=DecayTime(j);
    end
   elseif size(DT,2) ~= Nref
    for j=1:size(DT,2)
     DT(Tr,j)=NaN;
    end
   end
   delete('decay_time.txt');
   save('decay_time.txt','DT','-ascii','-tabs');
  end
  if exist('halfwidth.txt') ~= 0
   HW=load('-ascii','halfwidth.txt');
   if size(HW,2) == Nref
    for j=1:size(HW,2)
     HW(Tr,j)=HalfWidth(j);
    end
   elseif size(HW,2) ~= Nref
    for j=1:size(HW,2)
     HW(Tr,j)=NaN;
    end
   end
   delete('halfwidth.txt');
   save('halfwidth.txt','HW','-ascii','-tabs');
  end
  if exist('_results') ~= 0
   results=load('-ascii','_results');
    if (reject == 0) && (N == Nref) && (Nref ~= 0)
     results(Tr,1)=0;
    else
     results(Tr,1)=1;
    end
    if reject == 0
     results(Tr,2)=0;
    elseif reject == 1
     results(Tr,1)=0;
     results(Tr,2)=1;
    end
   delete('_results');
   dlmwrite('_results',results,'\t');
   delete('failed.txt');
   dlmwrite('failed.txt',find(results(:,1)),'\t');
   delete('rejected.txt');
   dlmwrite('rejected.txt',find(results(:,2)),'\t');
  end
  cd ..
 end
cd ../..
disp(sprintf('\nRemember to clear variables from the workspace before analysing other files\n'));
