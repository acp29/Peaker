%     Script File: avgtraces
%
%     When running stand-alone, this script calculates the mean trace.
%     If executed after running peaker, the baseline offset is subtracted
%     and the trace data can be optionally aligned to features of the
%     first peak. All the peaker-accepted trace data is plot overlaid.
%
%     If data contains 1 peak/trace, the event filter becomes available.
%     When the event filter is selected, peak-aligned events whose
%     decay cannot be fit with a single exponential using the Chebyshev
%     algorithm are excluded. The fit parameters of the real events are
%     then optimized using Nelder-Mead simplex minimization of the SSE.
%     and simplex restarts are performed to ensure an optimal solution.
%     The effectiveness of these methods is evident from the mean decay
%     trend line, which represents the mean of the fits to the individual
%     event decays.
%
%     Event filtering can also be achieved by providing maximum cut-offs
%     for event rise and decay time. For example, this can be used to
%     exclude distal synaptic conductances that are strongly filtered
%     by the dendritic arborization and are thus under poor voltage
%     control. By default, the cut-offs are set automatically using
%     Tukey's box plot rules (k=1.5). To disable cut-off filtering,
%     enter 'Inf' as the cut-off.
%
%     Compared to peak amplitudes, the integrals of synaptic events are
%     less attenuated by dendritic filtering. The integral of each event
%     is calculated from the preceding baseline to the point where the
%     event has decayed by 3 time constants (i.e. by more than 95 %).
%     Note that when the latter point exceeds the user-defined maximum
%     time point for decay fitting, the event is excluded. This is to
%     prevent closely following events contributing to the integral
%     calculation. The time integrals of current and voltage have units
%     coulomb (C) and weber (Wb) respectively. The median of the peak
%     integrals is reported.
%
%     This script is distributed with the script 'peaker'. Note that
%     peaker MUST be run prior to the execution of this script.
%
%     Bibliography:
%     Brant (1990) J Am Stat Assoc 85(412): 1083-90
%     Magee and Cook (2000) Nat Neurosci 3: 895-903
%     Malachowski, Clegg and Redford (2007) J Microsc 228(3): 282-95
%     Press et al. (1992) Numerical Recipes in C. Cambridge University Press
%     Williams and Mitchell (2008) Nat Neurosci 11: 790-798
%
%     avgtraces v2.00 (last updated: 01/06/2013)
%     Author: Andrew Charles Penn
%     https://www.researchgate.net/profile/Andrew_Penn/


if exist('diary','file') == 2
 diary('off');
 delete('diary');
end

if exist('filename','var')
 if exist(strcat('peaker.output/',filename,'/tables'))
  cd(strcat('peaker.output/',filename,'/tables'))
  if exist('diary','file') == 2
   diary('off');
   delete('diary');
  end
 end
end

clear
close all
format short g

if exist('_ID') == 0
% Ensure that you are in the "tables" directory if you wish to reanalyse peaker output when the global scope is clear
% The following code is to perform arithmetic mean of all traces
% This function of avgtraces is independent of peaker
 filename=input(sprintf('\nData matrix filename (excluding extension): '),'s');
 %if exist(strcat(filename,'.txt.gz')) == 2
 % gunzip(strcat(filename,'.txt.gz'));
 %end
 %data=load('-ascii', strcat(filename,'.txt'));
 data = ephysIO(filename)
 x=data(:,1);
 y=data; y(:,1)=[];
 Y=mean(y,2);
 Y_var=std(y,0,2).^2;
 % Plot average trace
 try
  figure(1,'visible','on');
 catch
  figure(1);
 end
 clf;
 y_autoscale=0.05*(max(Y)-min(Y)); y_maxlim=max(Y)+y_autoscale; y_minlim=min(Y)-y_autoscale; % Encoded y-axis autoscaling
 plot(x,Y,'k');
 box('off');
 grid('off'); xlim([min(x),max(x)]); ylim([y_minlim y_maxlim]);
 % Save data
 output_data=cat(2,x,Y,Y_var);
% output_data=cat(2,x,Y);
 newfilename=strcat(filename,'_avg.txt');
 dlmwrite(newfilename,output_data,'\t');
 gzip(newfilename);
 delete(newfilename);
 gzip(strcat(filename,'.txt'));
 delete(strcat(filename,'.txt'));
else

if exist('_sample_rate') == 0
 error('The sample file does not exist in the current directory. Ensure that you are in the "tables" directory.')
end

% Load data
ID=load('_ID');
sample_rate=load('_sample_rate');
parameters=load('-ascii','_parameters');
clamp=parameters(1);
s=parameters(3);
if length(parameters) == 7
 scale_factor=parameters(7);
end
unitTime=1/sample_rate;
filename=char(ID);
R20=load('-ascii','rise20.txt');
D50=load('-ascii','decay50.txt');
DT=load('-ascii','decay_time.txt');
HW=load('-ascii','halfwidth.txt');
A=load('-ascii','absolute.txt');
R=load('-ascii','relative.txt');
BT=load('-ascii','baseline_time.txt');
PT=load('-ascii','peak_time.txt');
RT=load('-ascii','risetime.txt');
SP=load('-ascii','slope.txt');
SR=load('-ascii','steepest_rise.txt');
cd ../../..
if exist('diary') == 2
 diary('off');
 delete(diary);
end
diary('on');
%if exist(strcat(filename,'.txt.gz')) == 2
% gunzip(strcat(filename,'.txt.gz'));
%end
%data=load('-ascii', strcat(filename,'.txt'));
%gzip(strcat(filename,'.txt'));
%delete(strcat(filename,'.txt'));
data=ephysIO(filename);
if ~isempty(regexp(filename,'_events','once'))
 if exist('events.output','dir') ~= 0
  cd events.output
 end
 if exist(filename,'dir') ~= 0
  chdir(filename);
 end
 ET=load('-ascii',strrep(filename,'_events','_times.txt'));
 cd ../..
end
if exist('./txt/IEI.txt')
  % Quickfix for analysing events detected with eventer
  [tok,rem]=strtok(filename,'.');
  filename=strcat(tok,'_events',rem);
  cd txt
  IEI=load('-ascii','IEI.txt');
  if isnan(IEI(1))
    IEI(1)=[];
    R20(1)=[];
    D50(1)=[];
    DT(1)=[];
    HW(1)=[];
    A(1)=[];
    R(1)=[];
    BT(1)=[];
    PT(1)=[];
    RT(1)=[];
    SP(1)=[];
    SR(1)=[];
    data(:,2)=[];
  end
  cd ..
  ET=cumsum(IEI);
end


% Upsample data if applicable
% This feature is currently inactive
if exist('scale_factor')
 if scale_factor ~= 1
%  l=size(data,1);
%  l=1+scale_factor*(l-1);
%  t=linspace(min(data(:,1)),max(data(:,1)),l);
%   for i=1:size(data,2)-1
%    state = warning('query');
%    warning off
%    y(:,i)=interp1q(data(:,1),data(:,i+1),t);
%    warning(state)
%   end
  t=data(:,1);
  y=data; y(:,1)=[];
 elseif scale_factor == 1
  t=data(:,1);
  y=data; y(:,1)=[];
 end
else
 t=data(:,1);
 y=data; y(:,1)=[];
end
if exist('_offset','file')
  offset=load('-ascii','_offset');
  t=t+offset;
end

Tr=input(sprintf('\nSelect upto which trace number you wish to analyse (default is all): '));
 if isempty(Tr)
  N=size(y,2);
  Tr=N;
 elseif Tr==round(Tr) && Tr > 1
   if exist('ET')==1
    if Tr<size(y,2)
      ET=ET(1:Tr);
    end
   end
  y=y(:,1:Tr);
  N=size(y,2);
  R20=R20(1:Tr,:);
  D50=D50(1:Tr,:);
  DT=DT(1:Tr,:);
  HW=HW(1:Tr,:);
  A=A(1:Tr,:);
  R=R(1:Tr,:);
  PT=PT(1:Tr,:);
  BT=BT(1:Tr,:);
  RT=RT(1:Tr,:);
  SP=SP(1:Tr,:);
  SR=SR(1:Tr,:);
 else
  error('The trace number must be a finite positive integer > 1');
 end
edx=1:Tr; edx=edx(:);

avg=input(sprintf('\nWould you like to plot the average trace (1) or not (0, default is yes)?: '));
if isempty(avg)
 avg=1;
elseif (avg ~= 0) && (avg ~= 1)
 error('The answer should be a logical value');
end

if size(R,2) > 1
 aln=input(sprintf('\nWould you like to perform alignment (1) or not (0, default is no)?: '));
  if isempty(aln)
   aln=0;
  elseif (aln ~= 0) && (aln ~= 1)
   error('The answer should be a logical value');
  end
elseif size(R,2) == 1
 aln=input(sprintf('\nWould you like to perform alignment (1) or not (0, default is yes)?: '));
  if isempty(aln)
   aln=1;
  elseif (aln ~= 0) && (aln ~= 1)
   error('The answer should be a logical value');
  end
end

% Baseline subtraction
y_zeroed=y;
 for i=1:N
  if isnan(A(i,1)) == 0
   B(i,:)=A(i,:)-R(i,:);
   y_zeroed(:,i)=y(:,i)-B(i,1);
   A(i,:)=A(i,:)-B(i,1);
   B(i,:)=B(i,:)-B(i,1);
  elseif isnan(A(i,1)) == 1
   y_zeroed(:,i)=NaN;
  end
 end

% Make data evenly spaced if this is not the case
isDiscrete=~any(round(diff(t)*10e9)-mean(round(diff(t)*10e9)));
 for i=1:N
  clear yl
  if isDiscrete == 0
    if ~exist('tl')
     l=1+sample_rate*(max(t)-min(t));
     tl=linspace(min(t),max(t),l+1);
     tl=tl(:);
    end
   state = warning('query');
    warning off
     yl=interp1q(t,y_zeroed(:,i),tl);
    warning(state)
   y_zeroed(:,i)=yl(:);
  elseif isDiscrete == 1
   tl=t;
   tl=tl(:);
  end
 end

% Align events to the peak for decay analysis
ref=PT(:,1);
t_offset=max(ref)-ref;
y_shift_samples=round(t_offset/unitTime);
t_zero=max(ref);
t_shift_samples=round(t_zero/unitTime);
y_decays=y_zeroed;
tl_decays=tl;
for i=1:N
 if isnan(y_shift_samples(i)) == 0
  y_decays(:,i)=circshift(y_decays(:,i),y_shift_samples(i));
 elseif isnan(y_shift_samples(i)) == 1
  y_decays(:,i)=NaN;
 end
end
if ~isempty(regexp(filename,'_events','once'))
  tl_decays=circshift(tl_decays,t_shift_samples);
else
  tl_decays=tl;
end

% Remove traces rejected during peaker execution
for i=1:N
 if all(isnan(y_decays(:,i))) == 1
  reject(i,:)=1;
 elseif all(isnan(y_decays(:,i))) == 0
  reject(i,:)=0;
 end
end
y_decays(:,find(reject))=[];
y_zeroed(:,find(reject))=[];
R20(find(reject),:)=[];
D50(find(reject),:)=[];
DT(find(reject),:)=[];
HW(find(reject),:)=[];
A(find(reject),:)=[];
R(find(reject),:)=[];
PT(find(reject),:)=[];
BT(find(reject),:)=[];
RT(find(reject),:)=[];
SP(find(reject),:)=[];
SR(find(reject),:)=[];
edx(find(reject),:)=[];
 if exist('ET')==1
  ET(find(reject))=[];
 end
n=length(edx);

% Decay fitting feature
decay_fit=input(sprintf('\nWould you like to fit an exponential to the decays (1=yes, 0=no, default is no)?: '));
 if isempty(decay_fit)
  decay_fit=0;
 end
if decay_fit == 1
 % Fit decay-time constants
 template=mean(y_decays,2);
 if ~isempty(regexp(filename,'_events','once'))
  idxPk=dsearchn(tl_decays,0);
 else
  idxPk=dsearchn(tl_decays,t_zero);
 end
 L=input(sprintf('\nEnter the fraction of the peak to fit exponential from (default is 1): '));
 if isempty(L)
  L=1;
 end
 L=1-L;
 hi_time=input(sprintf('\nEnter the time period after peak to fit the exponential to (in ms): '));
 if isempty(hi_time)==0
  if ~isempty(regexp(filename,'_events','once'))
   hi=dsearchn(tl_decays,hi_time*1e-3);
  else
   hi=dsearchn(tl_decays,t_zero+hi_time*1e-3);
  end
 elseif isempty(hi_time)==1
  hi=length(tl_decays);
 end
 template_edge=template(idxPk:hi);
 [slope,intervals,P]=edge(tl_decays(idxPk:hi),template_edge,L,1);
 lo=dsearchn(tl_decays,P(1,1));
 ti=tl_decays(lo:hi);
 yi=y_decays(lo:hi,:);
 EXPR='p(1)+p(2)*exp(-ti/p(3))';
 F=inline(EXPR,'ti','p');
 EXPR=regexprep('sum((y_inp-(E)).^2)','E',EXPR);
 Fmin=inline(EXPR,'p','ti','y_inp');
 %simSize=input(sprintf('\nEnter Simplex size for decay fitting (default is 1e-15): '));
 %if isempty(simSize)
  simSize=1e-12;
 %end
 y_fit=yi*NaN;
 if clamp == 0
   SF=1e12;
 elseif clamp == 1
   SF=1e3;
 end
  for i=1:n
   clear y_inp S wf p pmin
   try
   [fit,S]=chebexp(ti,SF*yi(:,i),1);
    if (sign(S.tau) == 1) && (sign(S.a) == s)
      % Nelder-Mead Simplex optimization
      y_inp=SF*yi(:,i);
	  p=[S.offset S.a S.tau];p=p(:);
      try
       [pmin]=fminsearch(Fmin,p,[0,simSize,inf,inf,1,0],[],ti,y_inp);     % For Octave
	   [pmin]=fminsearch(Fmin,pmin,[0,simSize,inf,inf,1,0],[],ti,y_inp);  % Optimization restart in case of anomalous step
      catch
       [pmin]=fminsearch(Fmin,p,[],ti,y_inp);                             % For Matlab
       [pmin]=fminsearch(Fmin,pmin,[],ti,y_inp);                          % Optimization restart in case of anomalous step
      end
      wf=feval(F,ti,pmin);
      fit_SSE(i,1)=sum((y_inp-wf).^2);
       if fit_SSE(i) <= S.normr^2
        TC(i,1)=pmin(3);
	    y_fit(:,i)=1/SF*wf;
       elseif S.normr^2 < fit_SSE(i)
        TC(i,1)=S.tau;
        y_fit(:,i)=1/SF*feval(F,ti,p);
        fit_SSE(i,1)=S.normr^2;
       end
    else
     TC(i,1)=NaN;
     y_fit(:,i)=NaN;
     fit_SSE(i,1)=NaN;
    end
   catch
    TC(i,1)=NaN;
    y_fit(:,i)=NaN;
    fit_SSE(i,1)=NaN;
   end
  end

 % Calculate integral (if blanked stimulus artifact does not overlap with the event))
 if exist('ET')==1
  IEI=cat(1,diff(ET),inf);
 end
 I=ones(n,1);
 for i=1:n
  if isnan(TC(i))
   I(i,:)=0;
  else
   if exist('IEI')==1
    if (P(1,1)+TC(i)*3 > t_zero+hi_time*1e-3) || (IEI(i) < t_zero+hi_time*1e-3)
	 I(i,:)=0;
	else
	 lo_idx=dsearchn(tl,BT(i));
     hi_idx=dsearchn(tl,P(1,1)+TC(i)*3);
     I(i,:)=trapz(tl(lo_idx:hi_idx),y_zeroed(lo_idx:hi_idx,i));
     %Plot cumulative peak integral for diagnostic purposes
     %figure(4); hold on; plot(tl(lo_idx:hi_idx),cumtrapz(tl(lo_idx:hi_idx),y_zeroed(lo_idx:hi_idx,i)));hold off;
    end
   else
    lo_idx=dsearchn(tl,BT(i));
    hi_idx=dsearchn(tl,P(1,1)+TC(i)*3);
    I(i,:)=trapz(tl(lo_idx:hi_idx),y_zeroed(lo_idx:hi_idx,i));
    %Plot cumulative peak integral for diagnostic purposes
    %figure(4); hold on; plot(tl(lo_idx:hi_idx),cumtrapz(tl(lo_idx:hi_idx),y_zeroed(lo_idx:hi_idx,i)));hold off;
   end
  end
 end
 % Exclusion of decay and integral values for which decay fitting or integral calculations failed
 exclude=find(I==0);
 y_decays(:,exclude)=[];
 y_zeroed(:,exclude)=[];
 R20(exclude,:)=[];
 D50(exclude,:)=[];
 DT(exclude,:)=[];
 HW(exclude,:)=[];
 A(exclude,:)=[];
 R(exclude,:)=[];
 PT(exclude,:)=[];
 BT(exclude,:)=[];
 RT(exclude,:)=[];
 SP(exclude,:)=[];
 SR(exclude,:)=[];
 TC(exclude,:)=[];
 I(exclude,:)=[];
 y_fit(:,exclude)=[];
 fit_SSE(exclude)=[];
 excluded=edx(exclude);
 edx(exclude)=[];
  if exist('ET')==1
   ET(exclude)=[];
  end
 tl_decays=tl_decays(idxPk:hi);
 y_decays=y_decays(idxPk:hi,:);

 % Set cut-offs for event filtering by Tukey's box plot rules
 % Applies only to the first event of a multi-event set of traces
 % Rise-time (20-80 percent)
 RTmax=input(sprintf('\nEnter optional cut-off for 20-80 percent event rise-time (in ms, default is auto): '));
 logRT=log(RT(:,1));  % Make logaritmic transformation to reduce distribution skew. Impossible to have rise time < 0
  if isempty(RTmax)
   Q2=median(logRT);
   if round(n/2) == n/2  % if n is even
    Q1=median(logRT(find(logRT < Q2)));
    Q3=median(logRT(find(logRT > Q2)));
   elseif round(n/2) ~= n/2  % if n is odd
    Q1=median(logRT(find(logRT <= Q2)));
    Q3=median(logRT(find(logRT >= Q2)));
   end
   IQR=Q3-Q1;
   k=1.5;  % For outliers and extreme values
   %k=3;   % For extreme values only
   TukeysL=Q1-k*IQR;
   TukeysU=Q3+k*IQR;
   RTmin=exp(TukeysL)*1e+3;
   RTmax=exp(TukeysU)*1e+3;
   %RTmin
   RTmax
  elseif sign(RTmax)==-1 || RTmax==0
   error('The cut-off must be a non-zero positive value');
  else
   RTmin=0;
  end
  RTmin=RTmin*1e-3;
  RTmax=RTmax*1e-3;
 % Decay time constant
 TCmax=input(sprintf('\nEnter optional cut-off for event decay tau (in ms, default is auto): '));
  logTC=log(TC); % Make logaritmic transformation to reduce distribution skew. Impossible to have decay tau < 0
  if isempty(TCmax)
   n=numel(logTC);
   Q2=median(logTC);
   if round(n/2) == n/2  % if n is even
    Q1=median(logTC(find(logTC < Q2)));
    Q3=median(logTC(find(logTC > Q2)));
   elseif round(n/2) ~= n/2  % if n is odd
    Q1=median(logTC(find(logTC <= Q2)));
    Q3=median(logTC(find(logTC >= Q2)));
   end
   IQR=Q3-Q1;
   k=1.5;  % For outliers and extreme values
   %k=3;   % For extreme values only
   TukeysL=Q1-k*IQR;
   TukeysU=Q3+k*IQR;
   TCmin=exp(TukeysL)*1e+3;
   TCmax=exp(TukeysU)*1e+3;
   %TCmin
   TCmax
  elseif sign(TCmax)==-1 || TCmax==0
   error('The cut-off must be a non-zero positive value');
  else
   TCmin=0;
  end
 TCmin=TCmin*1e-3;
 TCmax=TCmax*1e-3;
 clear exclude
 exclude=find((RT(:,1) > RTmax) + (TC > TCmax));
 %exclude=find((RT < RTmin) + (RT > RTmax) + (TC < TCmin) + (TC > TCmax));
 y_decays(:,exclude)=[];
 y_zeroed(:,exclude)=[];
 R20(exclude,:)=[];
 D50(exclude,:)=[];
 DT(exclude,:)=[];
 HW(exclude,:)=[];
 A(exclude,:)=[];
 R(exclude,:)=[];
 PT(exclude,:)=[];
 BT(exclude,:)=[];
 RT(exclude,:)=[];
 SP(exclude,:)=[];
 SR(exclude,:)=[];
 TC(exclude,:)=[];
 I(exclude,:)=[];
 y_fit(:,exclude)=[];
 fit_SSE(exclude)=[];
 excluded=sort(cat(1,excluded,edx(exclude)));
 edx(exclude)=[];
 nreal=size(edx,1);
  if exist('ET')==1
   ET(exclude)=[];
  end
 if (RTmax ~= inf) || (TCmax ~= inf)
  event_filter=1;
 else
  event_filter=0;
 end
 for j=2:size(R,2)
  TC(:,j)=NaN;
  I(:,j)=NaN;
 end
else
 nreal=n;
 event_filter=0;
end


% Summary after decay fitting and event filtering
if decay_fit == 1
 if event_filter == 1
  disp(sprintf('\nTraces that could not be fit and/or were excluded by event filter: \n'));
 elseif event_filter == 0
  disp(sprintf('\nTraces that could not be fit: \n'));
 end
 if isempty(excluded) == 1
  disp(0);
 else
  disp(excluded);
 end
 disp(sprintf('\nNumber of successful traces: \n')); disp(nreal);
 if (clamp == 0) && (size(R,2) > 1)
  disp(sprintf('\nMedian absolute peak amplitude (in pA): \n')); disp(median(A,1)'*1e12);
 elseif (clamp == 1) && (size(R,2) > 1)
  disp(sprintf('\nMedian absolute peak amplitude (in mV): \n')); disp(median(A,1)'*1e3);
 end
 if clamp == 0
  disp(sprintf('\nMedian relative peak amplitude (in pA): \n')); disp(median(R,1)'*1e12);
 elseif clamp == 1
  disp(sprintf('\nMedian relative peak amplitude (in mV): \n')); disp(median(R,1)'*1e3);
 end
 disp(sprintf('\nCoefficient of variation of the relative event amplitude: \n')); disp(std(R,1)'./abs(mean(R,1)'));
 if clamp == 0
  disp(sprintf('\nMedian peak integral (in fC): \n')); disp(median(I,1)'*1e15);
 elseif clamp == 1
  disp(sprintf('\nMedian peak integral (in mWb): \n')); disp(median(I,1)'*1e3);
 end
 disp(sprintf('\nMedian 20-80 percent risetime (in ms): \n')); disp(median(RT,1)'*1000);
 if clamp == 0
  disp(sprintf('\nMedian initial rising slope (in nA/ms): \n')); disp(median(SP,1)'*1e6);
 elseif clamp == 1
  disp(sprintf('\nMedian initial rising slope (in mV/ms): \n')); disp(median(SP,1)');
 end
 disp(sprintf('\nMedian half-amplitude decay time (in ms): \n')); disp(median(DT,1)'*1000);
 disp(sprintf('\nMedian half-width (in ms): \n')); disp(median(HW,1)'*1000);
 disp(sprintf('\nMedian decay time constant (in ms): \n')); disp(median(TC,1)'*1000);
 if exist('ET')==1
  ET=ET+R20(:,1); % redefines event time reference to the event 20 percent rise point
  IEI=diff(ET);
  freq=1/median(IEI);
  %IEI=padarray(IEI,1,NaN,'post');
  IEI=padarray(IEI,1,NaN,'pre');
  disp(sprintf('\nMedian instantaneous event frequency (per second): \n')); disp(freq); disp(sprintf('\n'));
 else
  disp(sprintf('\n'));
 end
elseif decay_fit == 0
 disp(sprintf('\nNumber of traces: \n')); disp(nreal);
  if avg == 1
   if (clamp == 0) && (size(R,2) > 1)
    disp(sprintf('\nMedian absolute peak amplitude (in pA): \n')); disp(median(A,1)'*1e12);
   elseif (clamp == 1) && (size(R,2) > 1)
    disp(sprintf('\nMedian absolute peak amplitude (in mV): \n')); disp(median(A,1)'*1e3);
   end
   if clamp == 0
    disp(sprintf('\nMedian relative peak amplitude (in pA): \n')); disp(median(R,1)'*1e12);
   elseif clamp == 1
    disp(sprintf('\nMedian relative peak amplitude (in mV): \n')); disp(median(R,1)'*1e3);
   end
   disp(sprintf('\nCoeffcient of variation of the event amplitude: \n')); disp(std(R,1)'./abs(mean(R,1)'));
   disp(sprintf('\nMedian 20-80 percent risetime (in ms): \n')); disp(median(RT,1)'*1000);
   if clamp == 0
    disp(sprintf('\nMedian initial rising slope (in nA/ms): \n')); disp(median(SP,1)'*1e6);
   elseif clamp == 1
    disp(sprintf('\nMedian initial rising slope (in mV/ms): \n')); disp(median(SP,1)');
   end
 disp(sprintf('\nMedian half-amplitude decay time (in ms): \n')); disp(median(DT,1)'*1000);
 disp(sprintf('\nMedian half-width (in ms): \n')); disp(median(HW,1)'*1000);
   if exist('ET')==1
    ET=ET+R20(:,1); % redefines event time reference to the event 20 percent rise point
    IEI=diff(ET);
    IEI=padarray(IEI,1,NaN,'post');
    freq=1/median(IEI);
    disp(sprintf('\nMedian instantaneous event frequency (Hz): \n')); disp(freq); disp(sprintf('\n'));
   else
    disp(sprintf('\n'));
   end
  end
end

% Average decays
y_avgdecay=mean(y_decays,2);

% Plot decays and fits
 try
  figure(3,'visible','on');
 catch
  figure(3);
 end
 clf;
if decay_fit == 1
 y_autoscale=0.05*(max(max(y_decays))-min(min(y_decays)));
 y_maxlim=max(max(y_decays))+y_autoscale;
 y_minlim=min(min(y_decays))-y_autoscale;
 plot(tl_decays,y_decays,'o','markersize',1,'color',[0.9,0.9,0.9],'markerfacecolor',[0.9,0.9,0.9]);hold on;
 if decay_fit == 1
  plot(ti,y_fit,'-','color',[0.75,0.75,0.75]);
  end
 plot(tl_decays,y_avgdecay,'k')
 if decay_fit == 1
  plot(ti,mean(y_fit,2),'k-');
 end
 hold off;box('off');
 if clamp == 0
  ylabel('Current (A)');
 elseif clamp == 1
  ylabel('Voltage (V)');
 end
 xlabel('Time (s)'); grid('off');
 if ~isempty(regexp(filename,'_events','once'))
  xlim([0,hi_time*1e-3]);
 else
  xlim([t_zero,t_zero+hi_time*1e-3]);
 end
 ylim([y_minlim y_maxlim]);
else
 close(3)
end

% Realign traces to rise
if aln == 1
 aln_method=input(sprintf('1) Risetime \n2) Peak \n3) Steepest rise\nSelect which reference point you want to align the events to (default is 1): '));
 if isempty(aln_method)
  aln_method=1;
 end
 if aln_method == 1
  ref=R20(:,1);
 elseif aln_method == 2
  ref=PT(:,1);
 elseif aln_method == 3
  ref=SR(:,1);
 end
 t_offset=max(ref)-ref;
 y_shift_samples=round(t_offset/unitTime);
 t_zero=max(ref);
 t_shift_samples=round(t_zero/unitTime);
 y_shifted=y_zeroed;
 tl_shifted=tl;
 for j=1:size(R,2)
  PT(:,j)=PT(:,j)+t_offset;
  BT(:,j)=BT(:,j)+t_offset;
  R20(:,j)=R20(:,j)+t_offset;
  D50(:,j)=D50(:,j)+t_offset;
 end
elseif aln == 0
 t_offset=zeros(N,1);
 y_shift_samples=zeros(N,1);
 t_zero=0;
 t_shift_samples=0;
 y_shifted=y_zeroed;
 tl_shifted=tl;
end
for i=1:nreal
 y_shifted(:,i)=circshift(y_shifted(:,i),y_shift_samples(i));
end
if ~isempty(regexp(filename,'_events','once'))
 tl_shifted=circshift(tl_shifted,t_shift_samples);
 PT=PT-t_zero;
 BT=BT-t_zero;
 R20=R20-t_zero;
 D50=D50-t_zero;
else
 tl_shifted=tl;
end
y_shifted(1:max(y_shift_samples),:)=[];
tl_shifted(1:max(y_shift_samples))=[];
if sign(t_shift_samples) == -1
y_shifted(find(tl_shifted == max(tl_shifted)):end,:)=[];
tl_shifted(find(tl_shifted == max(tl_shifted)):end)=[];
elseif sign(t_shift_samples) == 1
 y_shifted(1:find(tl_shifted == min(tl_shifted)),:)=[];
 tl_shifted(1:find(tl_shifted == min(tl_shifted)))=[];
end

% Average realigned traces and calculate variance
y_mean=mean(y_shifted,2);
y_var=std(y_shifted,0,2).^2;
y_median=median(y_shifted,2);
mean_trace=cat(2,tl_shifted,y_mean);
mean_var=cat(2,tl_shifted,y_mean,y_var);
median_trace=cat(2,tl_shifted,y_median);
if decay_fit == 1
 mean_var(hi+1:end,:)=[];
 mean_var(1:idxPk-1,:)=[];
end
 % Calculate peak-scaled variance
 SF=R(:,1)'/mean(R(:,1));
 y_sqdiff=(y_shifted-y_mean*SF).^2;
 for i=1:nreal
  y_scaled(:,i)=ones(size(y_shifted,1),1)*1/SF(i);
 end
 y_scaled=y_scaled.*y_shifted;
 y_scaled_var=sum(y_sqdiff,2)/(nreal-1);
 mean_var_scaled=cat(2,tl_shifted,y_mean,y_scaled_var);
 if decay_fit == 1
  mean_var_scaled(hi+1:end,:)=[];
  mean_var_scaled(1:idxPk-1,:)=[];
 end

% Plot realigned traces
try
 figure(1,'visible','on');
catch
 figure(1);
end
clf;
y_autoscale=0.05*(max(max(y_shifted))-min(min(y_shifted)));
y_maxlim=max(max(y_shifted))+y_autoscale;
y_minlim=min(min(y_shifted))-y_autoscale;
if avg == 1
 plot(tl_shifted,y_shifted,'-','markersize',1,'color',[0.85,0.85,0.85],'markerfacecolor',[0.85,0.85,0.85]);hold on;
 plot(tl_shifted,y_mean,'k-');
% For diagnistic purposes
% plot(PT,A,'m+','markersize',10);
% plot(BT,A-R,'m+','markersize',10);
% plot(R20,A-R+(0.2*R),'g+','markersize',10);
% plot(D50,A-R+(0.5*R),'c+','markersize',10);
elseif avg == 0
 plot(tl_shifted,y_shifted,'k-');hold on;
% For diagnistic purposes
% plot(PT,A,'m+','markersize',10);
% plot(BT,A-R,'m+','markersize',10);
% plot(R20,A-R+(0.2*R),'g+','markersize',10);
% plot(D50,A-R+(0.5*R),'c+','markersize',10);
end
hold off;box('off');
if clamp == 0
 ylabel('Current (A)');
elseif clamp == 1
 ylabel('Voltage (V)');
end
xlabel('Time (s)'); grid('off'); xlim([min(tl),max(tl)]); ylim([y_minlim y_maxlim]);
 try
 figure(2,'visible','on');
 catch
  figure(2);
 end
 clf;
 y_autoscale=0.05*(max(max(y_scaled))-min(min(y_scaled)));
 y_maxlim=max(max(y_scaled))+y_autoscale;
 y_minlim=min(min(y_scaled))-y_autoscale;
 if avg == 1
  plot(tl_shifted,y_scaled,'-','markersize',1,'color',[0.85,0.85,0.85],'markerfacecolor',[0.85,0.85,0.85]);
  hold on; plot(tl_shifted,y_mean,'k-');
 elseif avg == 0
  plot(tl_shifted,y_scaled,'k-');
 end
 hold off; box('off');
 if clamp == 0
  ylabel('Current (A)');
 elseif clamp == 1
  ylabel('Voltage (V)');
 end
 xlabel('Time (s)'); grid('off'); xlim([min(tl),max(tl)]); ylim([y_minlim y_maxlim]);


% Save data
diary('off');
if exist('./txt/IEI.txt')
  % Quickfix for analysing events detected with eventer
  filename=regexprep(filename,'_events','');
end
%if aln == 1
% aligned_data=cat(2,tl_shifted,y_shifted);
% aligned_filename=strcat(filename,'_shifted.txt');
% if exist(aligned_filename) ~= 0
%  delete(aligned_filename);
% end
% save(aligned_filename,'aligned_data','-ascii','-tabs');
% gzip(aligned_filename);
% delete(aligned_filename);
%end
% scaled_data=cat(2,tl_shifted,y_scaled);
% scaled_filename=strcat(filename,'_scaled.txt');
% if exist(scaled_filename) ~= 0
%  delete(scaled_filename);
% end
% save(scaled_filename,'scaled_data','-ascii','-tabs');
% gzip(scaled_filename);
% delete(scaled_filename);
%if decay_fit == 1
% decay_data=cat(2,tl_decays,y_decays);
% decay_filename=strcat(filename,'_decays.txt');
% if exist(decay_filename) ~= 0
%  delete(decay_filename);
% end
% save(decay_filename,'decay_data','-ascii','-tabs');
% gzip(decay_filename);
% delete(decay_filename);
% fits_data=cat(2,ti,y_fit);
% fits_filename=strcat(filename,'_fits.txt');
%  if exist(fits_filename) ~= 0
%   delete(fits_filename);
%  end
% save(fits_filename,'fits_data','-ascii','-tabs');
% gzip(fits_filename);
% delete(fits_filename);
%end
if exist('avgtraces.output','dir') == 0
 mkdir('avgtraces.output');
end
cd avgtraces.output
if exist(filename,'dir') == 0
 mkdir(filename);
end
cd(filename);
save('mean_trace.txt','mean_trace','-ascii','-tabs');
%gzip('mean_trace.txt');
%delete('mean_trace.txt');
save('median_trace.txt','median_trace','-ascii','-tabs');
%save('mean_var.txt','mean_var','-ascii','-tabs');
%gzip('mean_var.txt');
%delete('mean_var.txt');
%save('mean_var_scaled.txt','mean_var_scaled','-ascii','-tabs');
%gzip('mean_var_scaled.txt');
%delete('mean_var_scaled.txt');
if event_filter == 0
 movefile('../../diary','summary.txt');
elseif event_filter == 1
 movefile('../../diary','event_filter_summary.txt');
end
if event_filter == 0
 if exist('tables','dir') == 0
  mkdir('tables');
 end
 cd tables
elseif event_filter == 1
 if exist('event_filter_tables','dir') == 0
  mkdir('event_filter_tables');
 end
 cd event_filter_tables
end
save('absolute.txt','A','-ascii','-tabs');
save('relative.txt','R','-ascii','-tabs');
save('risetime.txt','RT','-ascii','-tabs');
save('slope.txt','SP','-ascii','-tabs');
save('decay_time.txt','DT','-ascii','-tabs');
save('halfwidth.txt','HW','-ascii','-tabs');
if decay_fit == 1
 dlmwrite('excluded.txt',excluded,'\t');
 save('tau_decay.txt','TC','-ascii','-tabs');
 save('integral.txt','I','-ascii','-tabs');
end
dlmwrite('successful.txt',edx,'\t');
if (exist('ET')==1)
% save('event_times.txt','ET','-ascii','-tabs');
 save('interevent_intervals.txt','IEI','-ascii','-tabs');
end
dlmwrite('_offset',tl_shifted(1));
cd ..
if decay_fit== 1
 figure(3);
 print(3,'decays.png','-dpng');
 print(3,'decays.eps','-depsc');
end
figure(2);
print(2,'scaled_output.png','-dpng');
print(2,'scaled_output.eps','-depsc');
figure(1);
print(1,'output.png','-dpng');
print(1,'output.eps','-depsc');
%if aln == 1
 % Save aligned data traces
 aligned_data=cat(2,tl_shifted,y_shifted);
 %if exist('output.txt.gz','file') ~= 0
 % gunzip('output.txt.gz');
 %end
 %if exist('output.txt','file') ~= 0
 % delete('output.txt');
 %end
 %save('output.txt','aligned_data','-ascii','-tabs');
 %gzip('output.txt');
 %delete('output.txt');
 state = warning('query');
 warning off
 if clamp == 0
  ephysIO('output.mat',aligned_data,'s','A');
 elseif clamp == 1
  ephysIO('output.mat',aligned_data,'s','V');
 end
 warning(state)

 % Save aligned and scaled data traces
 scaled_data=cat(2,tl_shifted,y_scaled);
 %if exist('scaled_output.txt.gz','file') ~= 0
 % gunzip('scaled_output.txt.gz');
 %end
 %if exist('scaled_output.txt') ~= 0
 % delete('scaled_output.txt');
 %end
 %save('scaled_output.txt','scaled_data','-ascii','-tabs');
 %gzip('scaled_output.txt');
 %delete('scaled_output.txt');
 state = warning('query');
 warning off
 if clamp == 0
  ephysIO('scaled_output.mat',scaled_data,'s','A');
 elseif clamp == 1
  ephysIO('scaled_output.mat',scaled_data,'s','V');
 end
 warning(state)
%end
cd ../..
cd(strcat('peaker.output/',filename))
if exist('img')
 tar('img.tar','img');
 gzip('img.tar');
 delete('img.tar');
 rmdir('img','s');
end
cd ../..
if ~isempty(regexp(filename,'_events','once'))
 meet(filename,clamp);
end
end

disp(sprintf('\nRemember to clear variables from the workspace before analysing other files\n'));
