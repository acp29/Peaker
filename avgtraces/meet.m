%     Function File: [N] = meet (filename, clamp)
%
%     Average the Mean Event of Every Trace
%
%     This function requires prior execution of 'peaker' and 'avgtraces'.
%
%     meet v1.0 (last updated: 30/06/2013)
%     Author: Andrew Charles Penn
%     https://www.researchgate.net/profile/Andrew_Penn/


function [N] = meet (filename, clamp)

if nargin < 2
 error('Invalid number of input arguments');
end

if nargin > 2
 error('Invalid number of input arguments');
end

if (clamp ~= 0) && (clamp ~=1)
 error('The type of recording must be specified with a logical value');
end

if regexp(filename,'_Tr._')
 Tr=eval(filename(regexp(filename,'_Tr')+3));
elseif regexp(filename,'_Tr.._')
 Tr=eval(filename((regexp(filename,'_Tr')+3):(regexp(filename,'_Tr')+4)));
end
cd('avgtraces.output')
for i=1:Tr
 if regexp(filename,'_Tr._')
  filename=regexprep(filename,'Tr.',strcat('Tr',num2str(i)));
 elseif regexp(filename,'_Tr.._')
  filename=regexprep(filename,'Tr..',strcat('Tr',num2str(i)));
 end
 if exist(filename,'dir')
  cd(filename)
  if exist('event_filter_summary.txt','file') ~= 0
   fid=fopen('event_filter_summary.txt','r');
  elseif exist('summary.txt','file') ~= 0
   fid=fopen('summary.txt','r');
  end
  temp=fgetl(fid);
  temp_str=sscanf(temp,'%c');
  while ~strcmp(temp_str,'Number of traces: ') && ~strcmp(temp_str,'Number of successful traces: ')
   temp=fgetl(fid);
   temp_str=sscanf(temp,'%c');
  end
  temp=fgetl(fid);
  temp=fgetl(fid);
  n(i,1)=sscanf(temp,'%u');
  gunzip('mean_trace.txt.gz');
  trace_array{i}=load('-ascii','mean_trace.txt'); % Add mean trace data to cell array
  startValues(i,:)=trace_array{i}(1,1);
  endValues(i,:)=trace_array{i}(end,1);
  gzip('mean_trace.txt');
  delete('mean_trace.txt');
  if exist('output.txt.gz','file') ~= 0
   gunzip('output.txt.gz');
   output_array{i}=load('-ascii','output.txt'); % Add output data to cell array
   output_array{i}(:,1)=[];
   gzip('output.txt');
   delete('output.txt')
  end
%  if exist('scaled_output.txt.gz','file') ~= 0
%   gunzip('scaled_output.txt.gz');
%   scaled_output_array{i}=load('-ascii','scaled_output.txt'); % Add scaled output data to cell array
%   scaled_output_array{i}(:,1)=[];
%   gzip('scaled_output.txt');
%   delete('scaled_output.txt')
%  end
  if exist('event_filter_tables','dir') ~= 0
   cd event_filter_tables
  elseif exist('tables','dir') ~= 0
   cd tables
  end
   DT{i}=load('-ascii','decay_time.txt');
   HW{i}=load('-ascii','halfwidth.txt');
   R{i}=load('-ascii','relative.txt');
   RT{i}=load('-ascii','risetime.txt');
   SP{i}=load('-ascii','slope.txt');
   IEI{i}=load('-ascii','interevent_intervals.txt');
   if exist('tau_decay.txt','file') ~= 0
    TC{i}=load('-ascii','tau_decay.txt');
   end
   if exist('integral.txt','file') ~= 0
    I{i}=load('-ascii','integral.txt');
   end
  cd ../..
 else
  startValues(i,:)=NaN;
  endValues(i,:)=NaN;
  n(i,1)=0;
 end
end
lo=max(startValues);
hi=min(endValues);
N=sum(n);

for i=1:Tr
 if ~isempty(trace_array{i})
  idx(:,i)=find(trace_array{i}(:,1) >= lo & trace_array{i}(:,1) <= hi);
  y(:,i)=trace_array{i}(idx(:,i),2);
  scale_factor(:,i)=n(i)/N * ones(size(y(:,1)));
  if ~exist('x','var')
   x(:,1)=trace_array{i}(idx(:,i),1);
  end
  if exist('output_array') ~= 0
   output_array{i}=output_array{i}(idx(:,i),:);
  end
%  if exist('scaled_output_array') ~= 0
%   scaled_output_array{i}=scaled_output_array{i}(idx(:,i),:);
%  end
 end
end
y_avg=sum(scale_factor.*y,2);
mean_trace=[x y_avg];
Tr=size(find(n),1);
if regexp(filename,'_Tr._')
 newfilename=regexprep(filename,'Tr.','ALL');
elseif regexp(filename,'_Tr.._')
 newfilename=regexprep(filename,'Tr..','ALL');
end
if ~exist(newfilename,'dir')
 mkdir(newfilename);
end
cd(newfilename);
if ~exist('tables','dir')
mkdir('tables');
end
cd('tables');
R=vertcat(R{1,:});
RT=vertcat(RT{1,:});
SP=vertcat(SP{1,:});
DT=vertcat(DT{1,:});
HW=vertcat(HW{1,:});
IEI=vertcat(IEI{1,:});
FRQ=1./IEI;
save('amplitude.txt','R','-ascii','-tabs');
save('risetime.txt','RT','-ascii','-tabs');
save('slope.txt','SP','-ascii','-tabs');
save('decay_time.txt','DT','-ascii','-tabs');
save('halfwidth.txt','HW','-ascii','-tabs');
save('interevent_interval.txt','IEI','-ascii','-tabs');
save('frequency.txt','FRQ','-ascii','-tabs');
if exist('TC','var') ~= 0
 TC=vertcat(TC{1,:});
 save('tau_decay.txt','TC','-ascii','-tabs');
end
if exist('I','var') ~= 0
 I=vertcat(I{1,:});
 save('integral.txt','I','-ascii','-tabs');
end
cd ..
save('mean_trace.txt','mean_trace','-ascii','-tabs');
fid=fopen('summary.txt','w');
fputs(fid,sprintf(cstrcat(num2str(N),' events averaged from ',num2str(Tr),' traces')));
if clamp == 0
 fputs(fid,sprintf(cstrcat('\n\nMedian peak amplitude (in pA): \n\n',num2str(1e12*median(R)))));
 if exist('I','var') ~= 0
  fputs(fid,sprintf(cstrcat('\n\nMedian peak integral (in fC): \n\n',num2str(1e15*median(I)))));
 end
 fputs(fid,sprintf(cstrcat('\n\nMedian 20-80 percent risetime (in ms): \n\n',num2str(1e3*median(RT)))));
 if exist('TC','var') ~= 0
  fputs(fid,sprintf(cstrcat('\n\nMedian decay time constant (in ms): \n\n',num2str(median(1e3*TC)))));
 end
 fputs(fid,sprintf(cstrcat('\n\nMedian instantaneous event frequency (Hz): \n\n',num2str(median(FRQ)))));
 fputs(fid,sprintf(cstrcat('\n\nMedian interevent interval (s): \n\n',num2str(median(IEI)))));
 fputs(fid,sprintf(cstrcat('\n\nMedian initial rising slope (in nA/ms): \n\n',num2str(median(1e6*SP)))));
 fputs(fid,sprintf(cstrcat('\n\nMedian half-amplitude decay time (in ms): \n\n',num2str(median(1e3*DT)))));
 fputs(fid,sprintf(cstrcat('\n\nMedian half-width (in ms): \n\n',num2str(1e3*median(HW)))));
elseif clamp == 1
 fputs(fid,sprintf(cstrcat('\n\nMedian peak amplitude (in mV): \n\n',num2str(1e3*median(R)))));
 if exist('I','var') ~= 0
  fputs(fid,sprintf(cstrcat('\n\nMedian peak integral (in mWb): \n\n',num2str(1e3*median(I)))));
 end
 fputs(fid,sprintf(cstrcat('\n\nMedian 20-80 percent risetime (in ms): \n\n',num2str(1e3*median(RT)))));
 if exist('TC','var') ~= 0
  fputs(fid,sprintf(cstrcat('\n\nMedian decay time constant (in ms): \n\n',num2str(median(1e3*TC)))));
 end
 fputs(fid,sprintf(cstrcat('\n\nMedian instantaneous event frequency (Hz): \n\n',num2str(median(FRQ)))));
 fputs(fid,sprintf(cstrcat('\n\nMedian interevent interval (s): \n\n',num2str(median(IEI)))));
 fputs(fid,sprintf(cstrcat('\n\nMedian initial rising slope (in mV/ms): \n\n',num2str(median(SP)))));
 fputs(fid,sprintf(cstrcat('\n\nMedian half-amplitude decay time (in ms): \n\n',num2str(1e3*median(DT)))));
 fputs(fid,sprintf(cstrcat('\n\nMedian half-width (in ms): \n\n',num2str(1e3*median(HW)))));
end
fclose(fid);
figure(5);
y_autoscale=0.05*(max(max(y_avg))-min(min(y_avg)));
y_maxlim=max(max(y_avg))+y_autoscale;
y_minlim=min(min(y_avg))-y_autoscale;
plot(x,y_avg,'k-');
xlim([min(x),max(x)]); ylim([y_minlim y_maxlim]);
xlabel('Time (s)'); grid('off'); box('off');
if clamp == 0
 ylabel('Current (A)');
elseif clamp == 1
 ylabel('Voltage (V)');
end
print(5,'output.png','-dpng');
print(5,'output.eps','-depsc');
close(5)
if exist('output_array') ~= 0
 % Save aligned data
 y_shifted=horzcat(output_array{:});
 aligned_data=cat(2,x,y_shifted);
 if exist('output.txt.gz','file') ~= 0
  gunzip('output.txt.gz');
 end
 if exist('output.txt','file') ~= 0
  delete('output.txt');
 end
 save('output.txt','aligned_data','-ascii','-tabs');
 gzip('output.txt');
 delete('output.txt');
end
%if exist('scaled_output_array') ~= 0
%%Save scaled data
% y_scaled=horzcat(scaled_output_array{:});
% scaled_data=cat(2,x,y_scaled);
% if exist('scaled_output.txt.gz','file') ~= 0
%  gunzip('scaled_output.txt.gz');
% end
% if exist('scaled_output.txt') ~= 0
%  delete('scaled_output.txt');
% end
% save('scaled_output.txt','scaled_data','-ascii','-tabs');
% gzip('scaled_output.txt');
% delete('scaled_output.txt');
%end
cd ../..

