%     Script file: viewer
%
%     View trace data. Press z and x to scroll through the traces
%     and left mouse click to break from scrolling and examine the
%     trace (e.g. zoom). Press escape at any time to return to the
%     command window.
%
%     viewer v1.0 (last updated: 22/02/2017)
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

% View traces
figure(1);
clf
i=0;
try
 disp(sprintf('Press ESCAPE key to return to the command window'))
 while i <= n
  i=i+1;
  clear gx gy button
  y_autoscale=0.1*(max(max(y(:,i)))-min(min(y(:,i))));
  y_maxlim=max(max(y(:,i)))+y_autoscale;
  y_minlim=min(min(y(:,i)))-y_autoscale;
  plot(t,y(:,i),'color','k');
  xlim([min(t),max(t)]);
  ylim([min(y_minlim),max(y_maxlim)]);
  box('off');
  title(sprintf(strcat(['Trace ',num2str(i)])));
  xlabel(strcat(names{1},' (',xunit,')'));
  ylabel(strcat(names{i+1},' (',yunit,')'));
  button=0;
  while (button ~= 32) && (button ~= 27) && (button ~= 13)...
        && (button ~= 120) && (button ~= 122) && (button ~= 1)
   [gx, gy, button]=ginput(1);
  end
  if button == 1
    break
  end
  if button == 27
   % Escape manual inspection by pressing the escape button
   error('This error message is to escape the manual inspection')
  end
  if button == 13
   % Escape manual inspection by pressing the escape button
   error('This error message is to escape the manual inspection')
  end
  if button == 120
   % Continue to next trace by pressing the x key
   if i < n
   elseif i == n
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


