%     Script file: rmoffset
%
%     Removes the constant y-offset calculated as the mean of values between
%     lower and upper limits defined graphically by left-clicking the mouse
%     cursor then right-clicking or pressing the 'RETURN' key. The resulting
%     data matrix is exported as an ascii text tile (.txt) with the filename
%     ending '_cp' if it does not exist already.
%
%     rmoffset v1.1 (last updated: 22/02/2017)
%     Author: Andrew Charles Penn
%     https://www.researchgate.net/profile/Andrew_Penn/


% Clear current variables, load data and define matrix columns
clear;
close all;
format short g
filename=input(sprintf('\nData matrix filename (including extension): '),'s');
cwd = pwd;
if ~isempty(regexpi(filename(end-2:end),'.gz'))
  [pathstr,filename,ext] = fileparts(filename(end-2:end));
elseif ~isempty(regexpi(filename(end-3:end),'.zip'))
  [pathstr,filename,ext] = fileparts(filename(end-3:end));
else
  [pathstr,filename,ext] = fileparts(filename);
end
if ~isempty(pathstr)
  chdir(pathstr);
end
[data,xdiff,xunit,yunit,names,notes] = ephysIO (strcat(filename,ext));
t=data(:,1);
y=data;
y(:,1)=[];
n=size(y,2);
l=numel(t);

% Plot preview graph
figure(1);
clf;
y_autoscale=0.05*(max(max(y))-min(min(y))); y_maxlim=max(max(y))+y_autoscale; y_minlim=min(min(y))-y_autoscale; % Encoded y-axis autoscaling
plot(t,y,'k-');grid('off');xlim([min(t),max(t)]); ylim([y_minlim y_maxlim]); box('off');
disp(sprintf('\nSet coordinates by clicking the mouse cursor on the preview graph\n'));
button=1;
while button ~= 3 & button ~=13
 clear gx gy
  if exist('gX') == 0
   m = 1;
  elseif exist('gX') == 1
   if length(gX) == 1
    m = 2;
   elseif length(gX) == 2
    clear gX
    m = 1;
   end
  end
 [gx, gy, button]=ginput(1);
  if button == 1
   gX(m,1)=gx;
   hold on; plot(gx,gy,'ko'); hold off;
  elseif button == 3
  end
  if exist('gX') == 1
   if length(gX) == 2 && button ~=3 && button ~=13
	gX=sort(gX);
    limits=dsearchn(t,gX);
    clf
    y_autoscale=0.05*(max(max(y(limits(1):limits(2),:)))-min(min(y(limits(1):limits(2),:)))); y_maxlim=max(max(y(limits(1):limits(2),:)))+y_autoscale; y_minlim=min(min(y(limits(1):limits(2),:)))-y_autoscale; % Encoded y-axis autoscaling
    plot(t,y,'k-');grid('off');xlim([gX(1),gX(2)]); ylim([y_minlim y_maxlim]); box('off');
   end
  end
end

% Calculate and subtract y-offset
output_data=data;
if exist('limits')
  for i=1:n
   clear y_offset
   y_offset=mean(y(limits(1):limits(2),i));
   output_data(:,i+1)=y(:,i)-y_offset;
  end
else
 output_data=data;
end
clf(1);
y_autoscale=0.05*(max(max(output_data(:,2:end)))-min(min(output_data(:,2:end)))); y_maxlim=max(max(output_data(:,2:end)))+y_autoscale; y_minlim=min(min(output_data(:,2:end)))-y_autoscale; % Encoded y-axis autoscaling
plot(output_data(:,1),output_data(:,2:n+1),'k-');grid('off');xlim([min(t),max(t)]); ylim([y_minlim y_maxlim]); box('off');

% Save data
cploc=regexp(filename,'_cp');
newfilename=strcat(filename,'_cp',ext);
ephysIO(newfilename,output_data,xunit,yunit,names,notes);
