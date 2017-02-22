%     Script File: nsfa
%
%     Performs non-stationary fluctuation analysis (NSFA) of events
%     detected by peaker and aligned by avgtraces. Execute this script
%     when the current working directory is the tables subfolder for
%     the avgtraces output of the relevant analysis.
%
%     The ensemble mean and variance are calculated from aligned events,
%     where peak scaling of the ensemble mean for calculation of the
%     ensemble variance is optional. The traces can then be low-pass
%     filtered and manually inspected to discard traces contaminated
%     with noise. Plots and statistics for the time stability of the
%     peak amplitude, rise time and decay time are made [1]. This also
%     includes a frequency histogram to confirm the assumed unimodal
%     distribution of the decay times. Tukey's box plot rules are then
%     used to exclude traces with extreme and outlying noise.
%
%     To achieve accurate peak scaling, a timecourse-matched template
%     is derived by resampling the ensemble mean to have equivalent
%     halfwidth for each event. This template is amplitude scaled and
%     fit by ordinary least squares methods (OLS) including an offset
%     parameter. The amplitude and offset parameters only are then used
%     to scale the ensemble mean current to each event to calculate the
%     ensemble variance for NSFA.
%
%     For conventional NSFA, an option is available to calculate the
%     ensemble variance from successive, overlapping traces [3,4]. This
%     method is more robust to drift but the uncertainties are larger.
%     Ensemble current mean and variance samples are allocated to
%     uniformly spaced bins of the mean and averaged. Bins corresponding
%     to mean current values below the standard deviation of the pre-
%     event background noise are excluded before then fitting a parabola
%     to the binned data by weighted least squares using errors determined
%     by autocovariance analysis as decribed [1-4]. Standard deviations
%     for the binned variances and the single-channel paramaters are
%     calculated from the diagonals of the error-covariance and the
%     parameter-covariance matrices respectively [4].
%
%     For solving linear least squares problems, instability from using
%     normal equations is avoided by calculating the inverse of square
%     matrices indirectly by singular value decomposition. The option
%     is available to limit the fit to a user-defined fraction of the
%     mean-variance relationship.
%
%     The output includes a data matrix of columns corresponding to the,
%     binned ensemble mean, the binned ensemble variance, the standard
%     deviation of the binned ensemble variance, the bin sizes and the
%     fit values.
%
%     This function requires the 'binxy', 'sinv' and 'pkscale' functions.
%     Note that peaker and avgtraces MUST be run prior to executing this
%     script.
%
%     See the example distributed with this script.
%
%     References:
%     [1] Hartveit and Veruki (2007) Nat Protoc 2(2): 434-48
%     [2] Hartveit and Veruki (2006) J Physiol 574(3): 751-85
%     [3] Heinemann & Conti (1992) Methods Enzymol 207: 131-48
%     [4] Steffan and Heinemann (1997) J Neurosci Methods 78: 51-63
%
%     nsfa v1.0 (last updated: 03/08/2013)
%     Author: Andrew Charles Penn
%     https://www.researchgate.net/profile/Andrew_Penn/


if exist('diary','file') == 2
 diary('off');
 delete('diary');
end

if exist('filename','var')
 if regexp(filename,'_Tr._')
  folderALL=regexprep(filename,'Tr.','ALL');
 elseif regexp(filename,'_Tr.._')
  folderALL=regexprep(filename,'Tr..','ALL');
 else
   folderALL=filename;
 end
 if exist(strcat('avgtraces.output/',folderALL,'/tables'))
  cd(strcat('avgtraces.output/',folderALL,'/tables'))
  save -ascii _ID filename
  if exist('diary','file') == 2
   diary('off');
   delete('diary');
  end
 elseif exist('tables','dir')
  cd tables
  save -ascii _ID filename
 else
  save -ascii _ID filename
 end
end

clear
close all
format short g
warning('off','MATLAB:nearlySingularMatrix')

if exist('_ID','file')
 ID=load('_ID');
 filename=char(ID);
 delete('_ID');
end

% Load data
% Ensure that you are in the 'tables' directory if you wish to run nsfa on
% avgtraces output
try
 if exist('amplitude.txt','file');
  R=load('-ascii','amplitude.txt');
 end
 if exist('relative.txt','file');
  R=load('-ascii','relative.txt');
 end
 RT=load('-ascii','risetime.txt');
 DT=load('-ascii','decay_time.txt');
 HW=load('-ascii','halfwidth.txt');
 if exist('interevent_intervals.txt')
  IEI=load('-ascii','interevent_intervals.txt');
 else
  IEI=[];
 end
 if exist('_offset','file')
  offset=load('-ascii','_offset');
 end
catch
 error(['Ensure that you are in the "tables" directory if you wish to '...
        'run nsfa on avgtraces output']);
end
cd ..
if exist('diary','file') == 2
 diary('off');
 delete('diary');
end
diary('on');
scale=input(sprintf(['\n(1) Conventional \n(2) Peak-scaled\nSelect '...
                       'which analysis to perform (default is 1): ']));
%if exist('output.txt.gz','file')
% gunzip('output.txt.gz');
%end
%data=load('-ascii','output.txt');
%gzip('output.txt');
%delete('output.txt');
data=ephysIO('output.mat');
if isempty(scale)
 scale=1;
end
if exist('offset','var')
  t=data(:,1)+offset;
else
  t=data(:,1);
end
y=data;
y(:,1)=[];
if mean(R) > 0
elseif mean(R) < 0
 y=y*-1;
end
data=cat(2,t,y);
Fc=input(sprintf(['\nSet the -3 dB cut-off of the binomial filter '...
                  '(in kHz, default is none): ']));
if isempty(Fc)
else
 y=filter1(y,t,0,Fc*1000,'binomial');
end

% Discard traces with missing data
ridx=[];
ridx=vertcat(find(isnan(IEI)),find(isnan(RT)),find(isnan(HW)),find(isnan(DT)));
ridx=unique(ridx);
R(ridx)=[];
RT(ridx)=[];
DT(ridx)=[];
HW(ridx)=[];
if ~isempty(IEI)
 IEI(ridx)=[];
end
y(:,ridx)=[];
data(:,ridx+1)=[];
Nr=size(R,1);
disp(sprintf('\nNumber of events discarded for missing data: '));
disp(size(ridx,1));

% Subtract residual y-offset and define window size for peak scaling
rise_time=mean(RT);
winMin=rise_time*-3-0.001;
% Subtract the baseline amplitude from a one millisecond pre-event
% time segment. The median is used in this step to ensure robust
% baseline determination in the presence of outlying noisy samples.
y_baseline=y(1:dsearchn(t,winMin+0.001),:);
y=y-repmat(median(y_baseline,1),size(y,1),1);
data=cat(2,t,y);
numDT=input(sprintf(['\nEnter the number of decay half-times to study '...
                     'peak fluctuations (default is 7): ']));
if isempty(numDT)
 numDT=7;
end
decay_time=mean(DT);
y_avg=mean(y,2);
peak_time=t(find(y_avg == max(y_avg)));
winMax=peak_time+decay_time*numDT;


% Discard traces with overlapping events
% Find indices of events that have another event detected in the decay
clear ridx
ridx=find(IEI < (peak_time+numDT*decay_time));
% Include indices of those events whose decay overlaps with the next event
ridx=vertcat(ridx,ridx-1);
ridx(ridx==0)=[];
R(ridx)=[];
RT(ridx)=[];
DT(ridx)=[];
HW(ridx)=[];
if ~isempty(IEI)
 IEI(ridx)=[];
end
y(:,ridx)=[];
data(:,ridx+1)=[];
Nr=size(R,1);
disp(sprintf('\nNumber of overlapping events discarded: '));
disp(size(ridx,1));

% Peak scaling of traces that were baseline subtracted and aligned by
% avgtraces. To obtain more accurate peak scaling, the algorithm
% here involves fitting a scaled average template whose timecourse
% has been scaled using halfwidth data to match each event. The
% fitting includes a baseline offset coefficient.
y_avg=mean(y,2);
if scale == 1
 % Set scale factor to one for conventional NSFA
 yOffset(:,1)=zeros(Nr,1);
 ySF(:,1)=ones(Nr,1);
 tSF(:,1)=ones(Nr,1);
 for i=1:Nr
  B(:,i)=pkscale(y_avg,t,1,0,1);
 end
elseif scale == 2
 % Calculate scale factor for template matching
 tSF=HW/mean(HW);
 % Calculation of TIMECOURSE scaled average matrix for peak scaling
 for i=1:Nr
  B(:,i)=pkscale(y_avg,t,1,0,tSF(i));
 end
end
winIDX=find(((t>winMin) + (t<winMax)) -1);
t=t(winIDX);
y=y(winIDX,:);
B=B(winIDX,:);
y_avg=y_avg(winIDX);
for i=1:Nr
 if scale == 2
  % Fitting of timecourse matched template by ordinary least squares.
  % Since the order of the polynomial regression is low, it is not
  % necessary to scale the data. The inverted matrix is calculated
  % indirectly by singular value decomposition.
  A(:,1)=ones(size(B(:,i)));
  A(:,2)=B(:,i);
  kmin=sinv(A'*A)*A'*y(:,i);
  yOffset(i,1)=kmin(1);
  ySF(i,1)=kmin(2);
  clear A kmin
 end
 % Calculation of AMPLITUDE and TIMECOURSE scaled average matrix B
 B(:,i)=pkscale(B(:,i),t,ySF(i),yOffset(i),1);
 % Calculation of residuals
 F(:,i)=y(:,i)-B(:,i);
end


% Manual inspection of peak scaling
manual=input(sprintf(['\n(1) Perform manual inspection of traces'...
                      '\n(2) Load previous manual inspection results'...
                      '\n(3) Skip this step'...
                      '\nSelect an option (default is 1): ']));
if isempty(manual)
 manual=1;
 disp(sprintf(['\nUse left and right cursor keys to navigate the traces'...
               '\nToggle a trace for deletion by pressing the t key.'...
               '\nPress the RETURN key when the inspection is complete.'...
               '\nPress escape to quite and skip manual inspection.']));
end
if manual == 1
 figure(1);
 clf
 clear ridx
 ridx=[];
 i=0;
 try
  while i <= Nr
   i=i+1;
   clear gx gy button
   y_autoscale=0.1*(max(max(y(:,i)))-min(min(y(:,i))));
   y_maxlim=max(max(y(:,i)))+y_autoscale;
   y_minlim=min(min(y(:,i)))-y_autoscale;
   plot(t,y(:,i),'color',[0.75 0.75 0.75]);hold on;plot(t,B(:,i),'k');
   hold off;
   xlim([min(t),max(t)]);
   ylim([min(y_minlim),max(y_maxlim)]);
   box('off');
   if any(ridx==i)
    title(sprintf(strcat(['Trace ',num2str(i),' will be discarded'])));
   else
    title(sprintf(strcat(['Trace ',num2str(i)])));
   end
   if mean(R) > 0
    ylabel('Current (A)');
   elseif mean(R) < 0
    ylabel('Current (-A)');
   end
   xlabel('Time (s)');
   button=0;
   while (button ~= 32) && (button ~= 27) && (button ~= 13)...
         && (button ~= 120) && (button ~= 122)
    [gx, gy, button]=ginput(1);
   end
   if button == 27
    % Escape manual inspection by pressing the escape button
    error('This error message is to escape the manual inspection')
   end
   if button == 13
    % Escape manual inspection by pressing the escape button
    error('This error message is to escape the manual inspection')
   end
   if button == 32
    % Toggle traces for deletion by pressing the space bar
    if any(ridx==i)
     ridx(ridx==i)=[];
     i=i-1;
    else
     ridx=vertcat(ridx,i);
     i=i-1;
    end
   elseif button == 120
    % Continue to next trace by pressing the x key
    if i < Nr
    elseif i == Nr
     i=i-1;
    end
   elseif button == 122
    % Return to the previous trace by pressing z key
    if i > 1
     i=i-2;
    elseif i == 1
     i=i-1;
    end
   end
  end
 catch
  close(1)
 end
end
rem=pwd;
if ismac
 while ~isempty(rem)
  [tok,rem]=strtok(rem,'/');  %#ok<STTOK> No TEXTCAN in Octave
 end
elseif ispc
 while ~isempty(rem)
  [tok,rem]=strtok(rem,'\');  %#ok<STTOK> No TEXTCAN in Octave
 end
end
%if ~exist('filename','var')
% filename=regexprep(tok,'ALL','TrX');
%end
if exist('../../nsfa.output','dir') == 0
 mkdir('../../nsfa.output');
end
if exist(strcat('../../nsfa.output/',tok),'dir') == 0
 mkdir(strcat('../../nsfa.output/',tok));
end
if manual == 1
 dlmwrite('ridx.txt',ridx)
 movefile('ridx.txt',strcat('../../nsfa.output/',tok,'/ridx.txt'));
elseif manual == 2
 try
  ridx=[];
  ridx=load('-ascii',strcat('../../nsfa.output/',tok,'/ridx.txt'));
 catch
  disp(sprintf(['\nInspection results file ridx.txt is empty or '...
                'could not be found']));
 end
end
try
 if manual < 3
  R(ridx)=[];
  RT(ridx)=[];
  DT(ridx)=[];
  HW(ridx)=[];
  if ~isempty(IEI)
   IEI(ridx)=[];
  end
  y(:,ridx)=[];
  B(:,ridx)=[];
  F(:,ridx)=[];
  tSF(ridx)=[];
  ySF(ridx)=[];
  yOffset(ridx)=[];
  data(:,ridx+1)=[];
  Nr=size(R,1);
  disp(sprintf('\nNumber of manually discarded traces: '));
  disp(size(ridx,1));
 end
catch
 error('The time window is not the same as the loaded inspection')
end

% Discard 'Bad' traces defined by their outlying and extreme root
% mean square deviation of the fluctuations around the timecourse-
% matched template
RMSD=sqrt(mean(F.^2,1))';
Q2=mean(RMSD);
if round(Nr/2) == Nr/2  % if Nr is even
 Q1=median(RMSD(find(RMSD < Q2)));
 Q3=median(RMSD(find(RMSD > Q2)));
elseif round(Nr/2) ~= Nr/2  % if Nr is odd
 Q1=median(RMSD(find(RMSD <= Q2)));
 Q3=median(RMSD(find(RMSD >= Q2)));
end
IQR=Q3-Q1;
k=1.5;  % For outliers and extreme values
%k=3;   % For extreme values only
TukeysU=Q3+k*IQR;
clear ridx
ridx=find((RMSD > TukeysU));
R(ridx)=[];
RT(ridx)=[];
DT(ridx)=[];
if ~isempty(IEI)
 IEI(ridx)=[];
end
y(:,ridx)=[];
data(:,ridx+1)=[];
B(:,ridx)=[];
F(:,ridx)=[];
tSF(ridx)=[];
ySF(ridx)=[];
yOffset(ridx)=[];
RMSD(ridx)=[];
Nr=size(R,1);
disp(sprintf('\nNumber of noisy traces discarded: '));
disp(size(ridx,1));

% Discard traces whose events have superthreshold rise times
clear ridx
RTmax=input(sprintf(['\nInput a threshold for the maximum event '...
                     'rise time (in ms, default is inf): ']));
if isempty(RTmax)
 RTmax=inf;
end
clear ridx
ridx=find(RT > 1e-3*RTmax);
R(ridx)=[];
RT(ridx)=[];
DT(ridx)=[];
HW(ridx)=[];
if ~isempty(IEI)
 IEI(ridx)=[];
end
y(:,ridx)=[];
data(:,ridx+1)=[];
if exist('B','var')
 B(:,ridx)=[];
end
F(:,ridx)=[];
tSF(ridx)=[];
ySF(ridx,:)=[];
yOffset(ridx)=[];
Nr=size(R,1);
disp(sprintf('\nNumber of slow rising events discarded: '));
disp(size(ridx,1));
disp(sprintf('\nTotal number of successful traces: '));
disp(Nr);


% Linear fits and statistics for stability plots
idx=(1:Nr)';
P=polyfit(idx,R,1);
R_fit=polyval(P,idx);
P=polyfit(idx,RT,1);
RT_fit=polyval(P,idx);
P=polyfit(idx,DT,1);
DT_fit=polyval(P,idx);
if mean(R) > 0
 param=cat(2,R,RT,DT);
elseif mean(R) < 0
 param=cat(2,-1*R,RT,DT);
end
for j=1:3
 r=corr(idx,param(:,j),'type','Spearman');
 if j == 1
  disp(sprintf('\nCorrelation statistics for peak amplitudes: '));
 elseif j == 2
  disp(sprintf('\nCorrelation statistics for peak rise times: '));
 elseif j == 3
  disp(sprintf('\nCorrelation statistics for peak half-decay times: '));
 end
 if abs(r) < 0.2
  disp(sprintf('Very Weak correlation'));
 elseif (abs(r) > 0.2) && (abs(r) <= 0.4)
  disp(sprintf('Weak correlation'));
 elseif (abs(r) > 0.4) && (abs(r) <= 0.7)
  disp(sprintf('Modest correlation'));
 elseif (abs(r) > 0.7) && (abs(r) <= 0.9)
  disp(sprintf('Strong correlation'));
 elseif (abs(r) > 0.9) && (abs(r) <= 1.0)
  disp(sprintf('Very Strong correlation'));
 end
 disp(sprintf('Spearman Rank Correlation Coefficient (r): '));
 disp(r);
 t_stat=abs(r)*sqrt((Nr-2)/(1-r^2));
 p_value=2*(1-tcdf(t_stat,Nr-2));
 if p_value > 0.05
  disp(sprintf('Statistically not significant'));
 elseif (p_value <= 0.05) && (p_value > 0.01)
  disp(sprintf('Statistically significant'));
 elseif (p_value <= 0.01) && (p_value > 0.001)
  disp(sprintf('Statistically very significant'));
 elseif (p_value <= 0.001)
  disp(sprintf('Statistically extremely significant'));
 end
 disp(sprintf('Two-tailed P-value from Studentized r coefficient: '));
 disp(p_value);
end


% Half-decay time, which is related to channel open time. Since
% this parameter can effect estimates of the single channel
% parameters obtained by NSFA, it is important that it is the
% same across the experimental conditions.
% See Benke et al (2001) J Neurosci
disp(sprintf('\nHalf-decay time (+/- SD, in ms): '));
disp(sprintf(strcat([char(32),...
                     num2str(1e+3*mean(DT),3),...
                     ' +/- ',...
                     num2str(1e+3*sqrt(var(DT)),2)])));


% Repeat peak scaling since traces were discarded and the average trace
% to be scaled has changed. As before, to obtain more accurate peak
% scaling, the algorithm here involves fitting a scaled average template
% whose timecourse has been scaled using halfwidth data to match each
% event. The fitting includes a baseline offset coefficient.The only
% difference from earlier is that upon completion, the determined scale
% factors are used to scale the amplitude ONLY for subsequent non-
% stationary fluctuation analysis (NSFA).
if scale == 1
 y_avg=mean(y,2);
elseif scale == 2
 clear B F k ySF
 y=data;
 y(:,1)=[];
 t=data(:,1);
 y_avg=mean(y,2);
 % Calculation of TIMECOURSE scaled average matrix for peak scaling
 for i=1:Nr
  B(:,i)=pkscale(y_avg,t,1,0,tSF(i));
 end
end
winIDX=find(((t>winMin) + (t<winMax)) -1);
t=t(winIDX);
y=y(winIDX,:);
B=B(winIDX,:);
y_avg=y_avg(winIDX);
for i=1:Nr
 if scale == 2
  % Fitting of timecourse matched template by ordinary least squares.
  % Since the order of the polynomial regression is low, it is not
  % necessary to scale the data. The inverted matrix is calculated
  % indirectly by singular value decomposition.
  A(:,1)=ones(size(B(:,i)));
  A(:,2)=B(:,i);
  kmin=sinv(A'*A)*A'*y(:,i);
  yOffset(i,1)=kmin(1);
  ySF(i,1)=kmin(2);
  clear A kmin
 end
 % Calculation of AMPLITUDE and OFFSET ONLY scaled average matrix B
 B(:,i)=pkscale(y_avg,t,ySF(i),yOffset(i),1);
 % Calculation of fluctuations for NSFA analysis
 F(:,i)=y(:,i)-B(:,i);
end


if scale == 1
 opt=input(sprintf(['\n(1) Independent traces\n(2) Overlapping '...
                    'difference traces\nSelect which method to '...
                    'calculate the ensemble variance (default is 1): ']));
 if isempty(opt)
  opt=1;
 end
end
if ~exist('opt','var')
 opt=1;
end
if opt == 1
 % Calculates the ensemble variance from independent records. This method
 % is more powerful since the errors are only proportional to 2/(Nr-1).
 % However, it is sensitive to drift.
 y_var=1/(Nr-1)*sum(F.^2,2);
elseif opt == 2
 % Calculates the ensemble variance from successive, overlapping traces.
 % This method is not available for peak-scaled NSFA. It is more robust to
 % drift. However, the errors are larger being proportional to 3/Nr.
 y_var=1/(2*(Nr-1))*sum(((y(:,1:end-1)-y(:,2:end))).^2,2);
end

% Crop data to decay
pkidx=find(y_avg == max(y_avg));
t_decay=t;
y_avg_decay=y_avg;
y_var_decay=y_var;
F_decay=F;
t_decay(1:pkidx,:)=[];
F_decay(1:pkidx,:)=[];
y_avg_decay(1:pkidx,:)=[];
y_var_decay(1:pkidx,:)=[];

% Variance-mean binning
topEdge=max(y_avg_decay);
botEdge=min(y_avg_decay);
numBins=input(sprintf('\nEnter the number of bins: '));
[binAvg, binVar, binSize]=binxy(y_avg_decay,y_var_decay,numBins);
% Remove empty bins
binAvg(binSize==0)=[];
binVar(binSize==0)=[];
if ~all(binSize)
 disp(sprintf(['\nNote that some of the bins are empty.'...
               '\nThe fit will proceed on unequally spaced bins.']));
end
binSize(binSize==0)=[];
numBins=size(binAvg,1);

% Autocovariance analysis and error determination of the binned data.
% See Heinemann and Conti (1992) Methods in Enzymol and Steffan and
% Heinemann (1997) J Neurosci Methods
% The data is scaled to avoid the limitations in accuracy caused when
% calculations result in numbers approaching the machine precision
binAvg=1e+12*binAvg;
binVar=1e+24*binVar;
F_decay=1e+12*F_decay;
l=size(F_decay,1);
C=zeros(l,l);
% m lines/rows
for m=1:l
 % n columns
 for n=1:l
  if opt == 1
   C(m,n)=(1/(Nr-1))*sum(F_decay(m,:).*F_decay(n,:));
  elseif opt == 2
   C(m,n)=1/(2*(Nr-1))*sum(F_decay(m,:).*F_decay(n,:));
  end
 end
end
cumBinSize=cumsum(binSize);
if opt == 1
 % For variance calculated from independent records
 K=2/(Nr-1);
elseif opt == 2
 % For variance calculated from overlapping difference records
 K=3/Nr;
end
% m lines/rows
for m=1:numBins
 % Make a new temporary copy of the unbinned, error-covariance matrix
 temp=C;
 if m > 1
  % When the line/row (m) to be parsed is greater than 1, remove all
  % previous lines/rows from the temporary matrix
  temp(1:cumBinSize(m-1),:)=[];
 end
 % n columns
 for n=1:numBins
  % After parsing the column (n) for the current line/row, remove the
  % column
  %omega(m,n) = K * (binSize(m) * binSize(n))^-1 ...
  %             * sum(sum(resize(temp,binSize(m),binSize(n)).^2));
  omega(m,n) = K * (binSize(m) * binSize(n))^-1 ...
               * sum(sum(temp(1:binSize(m),1:binSize(n)).^2));
  temp(:,1:binSize(n))=[];
 end
end
binVarSD=sqrt(diag(omega));

% Use the variance calculated in a one millisecond pre-event time
% segment to determine the level of baseline noise.
clear y_baseline
y_baseline=y(1:dsearchn(t,winMin+0.001),:);
noiseVar=mean(var(y_baseline,0,1));
noiseSD=sqrt(noiseVar);
disp(sprintf('\nVariance of the pre-event baseline noise (in pA^2): '));
disp(sprintf(strcat([char(32), num2str(1e+24*noiseVar,3)])));

% The standard deviation of the baseline noise represents the minimum
% mean current value at which non-stationary fluctuations will be
% resolvable. Remove bins corresponding to amplitudes lower than this
% value.
%clear ridx
%ridx=find(binAvg<1e+12*noiseSD);
%disp(sprintf(['\nNumber of excluded bins below the baseline '...
%              'noise level: ']))
%disp(sprintf(strcat([char(32),num2str(numel(ridx))])));
%binAvg(ridx)=[];
%binVar(ridx)=[];
%binVarSD(ridx)=[];
%binSize(ridx)=[];

% Enter fit constraints. The feature to constrain the baseline variance
% is currently disabled. It is not recommended to constrain the fit in
% this manner since the baseline variance cannot be known with absolute
% certainty. As a result, constraining to a given baseline variance can
% give inaccurate estimates for single channel parameters, and under-
% estimate their uncertainties.
%
%constrain=input(sprintf(['\nTo constrain the fit, enter the '...
%                         'baseline variance (in pA^2): ']));
constrain=0;
if isempty(constrain)
 constrain=0;
end
fraction=input(sprintf(['\nEnter up to what fraction of the '...
                        'relationship to fit (default is 1): ']));
if isempty(fraction)
 fraction=1;
end
if fraction < 0 || fraction > 1
 error(['The entered value of the must be a fraction expressed as a '...
       'decimal'])
end
fidx=find(binAvg<fraction*y_avg_decay(1)*1e+12);
binAvg_ref=1e-12*binAvg;
binVar_ref=1e-24*binVar;
binVarSD_ref=1e-24*binVarSD;
binAvg=binAvg(fidx);
binVar=binVar(fidx);
binVarSD=binVarSD(fidx);
%omega=resize(omega,size(binVar,1));
omega=omega(1:size(binVar,1),1:size(binVar,1));
if fraction < 1
 disp(sprintf('\nNumber of terminal points excluded from the fit: '));
 disp(size(binAvg_ref,1)-size(binAvg,1));
end


% Perform weighted linear least squares fit taking into account the errors
% determined above and the user-defined constraints. Solution using normal
% equations where X is the data matrix.
%
% Create the data matrix (X), either with (unconstrained) or without
% (constrained), the first column of ones.
%  X(:,1)=ones(L,1);
%  X(:,2)=binAvg;
%  X(:,3)=binAvg.^2;
%  if constrain > 0
%   X(:,1)=[];
%  end
%
% Create diaganol matrix of weights. This is equivalent to inverting the
% error-covariance matrix (omega) with all off-diaganol elements set to 0.
%  w=diag(omega); W=diag(1./w);
%
% Obtain the parameter-covariance matrix for estimates of the
% parameter uncertainties. Note that these could be slightly
% biased when errors are correlated.
%  ksi = (X' * W * X)^-1
%
% Obtain unbiased estimates of the parameters by Weighted Least Squares
% from determined errors (WLSD). Subtract the constrained baseline
% variance from binVar values (for constrained fit only)
%  p = ksi * X' * W * (binVar - constrain)
%
% Calculate the variances of the parameters
%  pVar = diag(ksi)
%
% Convert parameter variances to standard deviations
%  pSD = sqrt(Var)
%
% The inverse (or inv) function used normally for the matrix inversion
% (^-1) can be unstable. Instead use the function sinv included with PEAK,
% which calculates the inverted matrix indirectly by singular value
% decomposition. Note that in the previous calculations, the data was
% scaled to avoid the limitations in accuracy caused as numbers approach
% the machine precision during matrix arithmetic. The following operations
% are equivalent to executing calculations for OLS after dividing the
% elements of the data matrix (X) as well as the individual data points
% (binVar) by the square root of the binned variance (i/e. the standard
% deviations)
%
%
L=length(binAvg);
X(:,1)=ones(L,1);
X(:,2)=binAvg;
X(:,3)=binAvg.^2;
if constrain > 0
 X(:,1)=[];
end
w=diag(omega);
W=diag(1./w);
ksi=sinv(X'*W*X);
p=ksi*X'*W*(binVar-constrain);
pVar=diag(ksi);
if constrain > 0
 p=vertcat(constrain,p);
 pVar=vertcat(0,pVar);
end
pSD=sqrt(pVar);
% Calculate residuals of the fit
residuals=binVar-polyval(flipud(p),binAvg);
SSE=sum(residuals.^2);
SST=sum((binVar-mean(binVar)).^2);
Rsq=1-(SSE/SST);
warning('off','MATLAB:nearlySingularMatrix')

% Display fit equation
if p(2) > 0 && p(3) > 0
 disp(sprintf(strcat(['\nf(x)= ' num2str(p(1)) ' + ' num2str(abs(p(2)))...
                      '*x + ' num2str(abs(p(3))) '*x.^2'])));
elseif p(2) < 0 && p(3) > 0
 disp(sprintf(strcat(['\nf(x)= ' num2str(p(1)) ' - ' num2str(abs(p(2)))...
                      '*x + ' num2str(abs(p(3))) '*x.^2'])));
elseif p(2) > 0 && p(3) < 0
 disp(sprintf(strcat(['\nf(x)= ' num2str(p(1)) ' + ' num2str(abs(p(2)))...
                      '*x - ' num2str(abs(p(3))) '*x.^2'])));
elseif p(2) < 0 && p(3) < 0
 disp(sprintf(strcat(['\nf(x)= ' num2str(p(1)) ' - ' num2str(abs(p(2)))...
                      '*x - ' num2str(abs(p(3))) '*x.^2'])));
end

disp(sprintf('\nCoefficient of determination (R-sq) of the fit: '));
disp(sprintf(strcat([char(32),num2str(Rsq,3)])));

% Fitted baseline variance
if constrain == 0
 disp(sprintf('\nFitted baseline variance (+/- SD, in pA^2): '));
 disp(sprintf(strcat([char(32),num2str(p(1),3),' +/- ',num2str(pSD(1),...
                     2)])));
elseif constrain > 0
 disp(sprintf('\nConstrained baseline variance (in pA^2): '));
 disp(sprintf(strcat([char(32),num2str(p(1),3)])));
end
if mean(R) > 0
 % Mean weighted single-channel current
 i=p(2);
elseif mean(R) < 0
 % Mean weighted single-channel current
 i=-1*p(2);
end
disp(sprintf('\nWeighted-mean single channel current (+/- SD, in pA): '));
disp(sprintf(strcat([char(32),num2str(i,3),' +/- ',num2str(pSD(2),2)])));

% The total number of channels, or when peak scaling this corresponds
% average number of channels open at the peak
N=abs(-1/p(3));

% The calculation of the variance of N from the parameter-covariance matrix
NVar=N^4*pVar(3);

% The conversion of the variance of N to standard deviation
NSD=sqrt(NVar);

% The maximum binned current
Imax=max(binAvg);

% The maximum peak open probability, or when peak scaling this corresponds
% to the observed over expected mean peak current
if constrain == 0
 Po=max(binAvg)/(p(2)*N);
elseif constrain > 0
 Po=NaN;
end

% The calculation of the variance of Po from error and parameter
% covariance matrices (with) error propagation)
if constrain == 0
 PoVar=(Imax^2/((i*N)^2)) * ((pVar(2)/(p(2)^2)) + (NVar/(N^2))...
       - ((2*N*ksi(2,3))/p(2)));
elseif constrain > 0
 PoVar=NaN;
end

% The conversion of the variance of Po to standard deviation
PoSD=sqrt(PoVar);

if scale == 1
 disp(sprintf('\nTotal number of channels (N):'));
 disp(sprintf(strcat([char(32),num2str(N,3),' +/- ',num2str(sqrt(NSD),...
                     2)])));
 disp(sprintf('\nMaximum peak open probability (Po): '));
 disp(sprintf(strcat([char(32),num2str(Po,3),' +/- ',num2str(PoSD,2)])));
elseif (scale == 2) && (fraction == 1)
 disp(sprintf('\nAverage number of channels open at peak (+/- SD):'));
 disp(sprintf(strcat([char(32),num2str(N,3),' +/- ',num2str(sqrt(NSD),...
                      2)])));
end
disp(sprintf('\n'));
binAvgFit=linspace(binAvg(1),binAvg(end),101)';
binVarFit=polyval(flipud(p),binAvgFit);

% Figures
figure(1); clf;
plot(t,y,'color',[0.8,0.8,0.8]);hold on; plot(t,y_avg,'k-'); hold off;
% Encoded y-axis autoscaling
y_autoscale=0.05*(max(max(y))-min(min(y)));
y_maxlim=max(max(y))+y_autoscale;
y_minlim=min(min(y))-y_autoscale;
xlim([min(t) max(t)]); ylim([y_minlim y_maxlim]); box('off'); grid('off');
xlabel('Time (s)');
if mean(R) > 0
 ylabel('Current (A)');
elseif mean(R) < 0
 ylabel('Current (-A)');
end

figure(2);
% Inverted colour map of squared fluctuations
imagesc(F.^2); colormap('gray'); cmap=colormap;...
     cmap=flipud(cmap); colormap(cmap);
set(gca,'tickdir','out','box','off');
% The following code is deprecated
%pcolor(y); colormap('gray');shading('flat');
%set(gca,'tickdir','out','linewidth',0.25,'ticklength',[0.005 0.025],...
%    'ydir','reverse','ytick',[],'box','off','xaxislocation','bottom',...
%    'xminortick','on');
xlabel('Event Number');ylabel('<--- Time (samples)');
h=colorbar;set(h,'location','eastoutside','tickdir','out','box','on',...
               'linewidth',0.25,'dataaspectratio',[0.2 1 1]);
set(get(h,'ylabel'),'string','Squared Difference Current (A^2)');


figure(3);
imagesc(1e-24*C); colormap('gray'); cmap=colormap;...
     cmap=flipud(cmap); colormap(cmap);
set(gca,'tickdir','out','box','off');
% The following code is deprecated
%pcolor(1e-24*C); colormap('gray');shading('flat');
%set(gca,'tickdir','out','linewidth',0.25,'ticklength',[0.005 0.025],...
%    'ydir','reverse','ytick',[],'xtick',[],'box','off','xaxislocation',...
%    'bottom');
xlabel('Time (samples) --->');ylabel('<--- Time (samples)');
h=colorbar;
set(h,'location','eastoutside','tickdir','out','box','on','linewidth',...
    0.25,'dataaspectratio',[0.2 1 1]);
set(get(h,'ylabel'),'string','Variance (A^2)');

figure(4);
subplot(2,1,1);
plot(t_decay,y_avg_decay,'k-');
% Encoded y-axis autoscaling
y_autoscale=0.05*(max(y_avg_decay)-min(min(y_avg_decay)));
y_maxlim=max(y_avg_decay)+y_autoscale;
y_minlim=min(y_avg_decay)-y_autoscale;
xlim([min(t_decay) max(t_decay)]);
ylim([y_minlim y_maxlim]); box('off'); grid('off');
ylabel('Current (A)');
subplot(2,1,2);
plot(t,y_var,'k-');
% Encoded y-axis autoscaling
y_autoscale=0.05*(max(y_var_decay)-min(y_var_decay));
y_maxlim=max(y_var_decay)+y_autoscale;
y_minlim=min(y_var_decay)-y_autoscale;
xlim([min(t_decay) max(t_decay)]);
ylim([y_minlim y_maxlim]); box('off'); grid('off');
ylabel('Variance (A^2)'); xlabel('Time (s)');

% Ensembl mean-variance plot with the offsets of the pre-event baseline
% noise subtracted
figure(5);
if max(binVarFit) >= max(1e+24*binVar_ref)
 % Encoded y-axis autoscaling
 y_autoscale=0.2*(max(1e-24*binVarFit));
 y_maxlim=max(1e-24*binVarFit)+y_autoscale;
elseif max(binVarFit) < max(1e+24*binVar_ref)
 % Encoded y-axis autoscaling
 y_autoscale=0.2*(max(binVar_ref));
 y_maxlim=max(binVar_ref)+y_autoscale;
end
if max(binAvgFit) >= 1e+12*y_avg_decay(1)
 % Encoded y-axis autoscaling
 x_autoscale=0.2*(max(1e-12*binAvgFit));
 x_maxlim=max(1e-12*binAvgFit)+x_autoscale;
elseif max(binAvgFit) < 1e+12*y_avg_decay(1)
 % Encoded y-axis autoscaling
 x_autoscale=0.2*(max(y_avg_decay)-0);
 x_maxlim=max(y_avg_decay)+x_autoscale;
end
hold on;
%h=errorbar(binAvg_ref,binVar_ref,binVarSD_ref,'~.k');
h=errorbar(binAvg_ref,binVar_ref,binVarSD_ref);
set(h,'linestyle','none'); set(h,'marker','o');
plot(1e-12*binAvgFit,1e-24*binVarFit,'-k')
hold off;
xlim([0 x_maxlim]); ylim([0 y_maxlim]); box('off'); grid('off');
ylabel('Ensemble variance (A^2)');
if mean(R) > 0
 xlabel('Ensemble mean (A)');
elseif mean(R) < 0
 xlabel('Ensemble mean (-A)');
end

figure(6);
subplot(2,2,1);
if mean(R) > 0
 % Encoded y-axis autoscaling
 y_autoscale=0.05*(max(R)-0);
 y_maxlim=max(R)+y_autoscale;
 plot(idx,R,'ok');
 hold on; plot(idx,R_fit,'k-'); hold off; box('off');
 ylabel('Amplitude (A)'); xlabel('Event Number');
 xlim([0 idx(end)]); ylim([0 y_maxlim]);
elseif mean(R) < 0
 % Encoded y-axis autoscaling
 y_autoscale=0.05*(max(-1*R)-0); y_maxlim=max(-1*R)+y_autoscale;
 plot(idx,-1*R,'ok');
 hold on; plot(idx,-1*R_fit,'k-'); hold off; box('off');
 ylabel('Amplitude (-A)'); xlabel('Event Number');
 xlim([0 idx(end)]); ylim([0 y_maxlim]);
end
subplot(2,2,2);
% Encoded y-axis autoscaling
y_autoscale=0.05*(max(RT)-0);
y_maxlim=max(RT)+y_autoscale;
plot(idx,RT,'ok');
hold on; plot(idx,RT_fit,'k-'); hold off; box('off');
ylabel('Rise time (s)');
xlabel('Event Number');
xlim([0 idx(end)]);
ylim([0 y_maxlim]);
subplot(2,2,3);
% Encoded y-axis autoscaling
y_autoscale=0.05*(max(DT)-0);
y_maxlim=max(DT)+y_autoscale;
y_minlim=min(DT)+y_autoscale;
plot(idx,DT,'ok');
hold on; plot(idx,DT_fit,'k-'); hold off; box('off');
ylabel('Decay time (s)');
xlabel('Event Number');
xlim([0 idx(end)]);
ylim([0 y_maxlim]);
subplot(2,2,4);
hist(DT,15,1,'facecolor',[0.75 0.75 0.75]);
box('off'); xlabel('Decay time (s)'); ylabel('Frequency');

% Save data
diary('off');
movefile('diary','../../nsfa.output/diary');
cd('../../nsfa.output');
if exist(tok,'dir') == 0
 mkdir(tok);
end
cd(tok);
movefile('../diary','./summary.txt');
binned_data=cat(2,binAvg_ref,binVar_ref,binVarSD_ref,binSize,...
                 1e-24*polyval(flipud(p),1e+12*binAvg_ref));
save('binned_data.txt','binned_data','-ascii','-tabs');
if exist('./eps','dir') == 0
 mkdir('./eps');
end
cd('eps');
print(1,'overlay.eps','-depsc');
print(2,'matrix.eps','-depsc');
print(3,'autocovariance.eps','-depsc');
print(4,'ensemble_traces.eps','-depsc');
print(5,'mean_var.eps','-depsc');
print(6,'stability.eps','-depsc');
cd ..
if exist('./png','dir') == 0
 mkdir('./png');
end
cd('png');
print(1,'overlay.png','-dpng');
print(2,'matrix.png','-dpng');
print(3,'autocovariance.png','-dpng');
print(4,'ensemble_traces.png','-dpng');
print(5,'mean_var.png','-dpng');
print(6,'stability.png','-dpng');
cd('../../..');
disp(sprintf('\nRemember to clear variables from the workspace before analysing other files\n'));
