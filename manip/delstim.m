%     Script file: delstim
%
%     Interpolates between lower and upper limits defined graphically by
%     left-clicking the mouse cursor then right-clicking or pressing the
%     'RETURN' key. The resulting data matrix is exported as an ascii text
%     tile (.txt) with the filename ending '_cp' if it does not exist
%     already.
%
%     delstim v1.1 (last updated: 22/02/2017)
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

% Perform linear interpolation on non-evenly spaced datasets
if xdiff == 0
 disp(sprintf('\nInput must consist of data sampled at evenly spaced time points.\nLinear interpolation will precede stimulus artifact deletion.'));
 tl=linspace(min(t),max(t),l);
 tl=tl(:);
 for i=1:n
  yl(:,i)=interp1q(t,y(:,i),tl,'linear');
 end
elseif xdiff > 0
  tl=t;
  yl=y;
end

% Plot preview graph
figure(1);
clf;
y_autoscale=0.05*(max(max(yl))-min(min(yl))); y_maxlim=max(max(yl))+y_autoscale; y_minlim=min(min(yl))-y_autoscale; % Encoded y-axis autoscaling
plot(tl,yl,'k-');grid('off');xlim([min(tl),max(tl)]); ylim([y_minlim y_maxlim]); box('off');
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
    limits=dsearchn(tl,gX);
    y_autoscale=0.05*(max(max(yl(limits(1):limits(2),:)))-min(min(yl(limits(1):limits(2),:)))); y_maxlim=max(max(yl(limits(1):limits(2),:)))+y_autoscale; y_minlim=min(min(yl(limits(1):limits(2),:)))-y_autoscale; % Encoded y-axis autoscaling
    clf
    plot(tl(limits(1):limits(2)),yl(limits(1):limits(2),:),'k-');grid('off'); xlim([gX(1),gX(2)]);ylim([y_minlim y_maxlim]); box('off');
   end
  end
end

% Interpolate between coordinates and replot data
output_data=[tl yl];
if exist('limits')
  for i=1:n
   output_data(limits(1):limits(2),i+1)=interp1(tl([limits(1) limits(2)]),yl([limits(1) limits(2)],i),tl(limits(1):limits(2)),'linear','extrap');
  end
end
clf(1);
y_autoscale=0.05*(max(max(output_data(:,2:end)))-min(min(output_data(:,2:end)))); y_maxlim=max(max(output_data(:,2:end)))+y_autoscale; y_minlim=min(min(output_data(:,2:end)))-y_autoscale; % Encoded y-axis autoscaling
plot(output_data(:,1),output_data(:,2:n+1),'k-');grid('off');xlim([min(tl),max(tl)]); ylim([y_minlim y_maxlim]); box('off');

% Save data
cploc=regexp(filename,'_cp');
newfilename=strcat(filename,'_cp',ext);
ephysIO(newfilename,output_data,xunit,yunit,names,notes);
