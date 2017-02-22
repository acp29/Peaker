%     Function File: [y] = pkscale (y, t, ySF, yOffset, tSF)
%
%     Scales the values of the y vector by the factor ySF with an offset
%     equal to yOffset. Also, the option is available to 'stretch' or
%     'squeeze' the y vector by the factor tSF by resampling the y vector
%     using linear interpolation and antialiasing filtering and reassigning
%     the values to the time vector by cropping (upsampling) or NaN padding
%     (downsampling) the ends around the zero time point in the t vector.
%
%     For example, the following command will output the y vector values
%     scaled up by a factor of 5 with an offset of 3 and scaled in time
%     by a factor of 0.5 (and therefore squeezed):
%
%     [y] = pkscale (y, t, 5, 3, 0.5)
%
%     pkscale v1.0 (last updated: 19/08/2013)
%     Author: Andrew Charles Penn
%     https://www.researchgate.net/profile/Andrew_Penn/


function [y] = pkscale (y, t, ySF, yOffset, tSF)

if nargin < 5
 error('Invalid number of input arguments');
end

if all(size(t) == 1) || ~any(size(t) == 1) || length(t) ~= length(y)
 error('t and y must be vectors of the same size');
end

if prod(size(ySF),2) ~= 1
 error('The y scale factor must be scalar');
end

if (prod(size(yOffset),2) ~= 1) && any((size(yOffset) ~= size(y)))
 error('The y offset must be scalar or a vector the same size as y');
end

if prod(size(tSF),2) ~= 1
 error('The t scale factor must be scalar');
end

% Assess sampling characteristics of input with precision of 10e-9
isDiscrete=~any(round(diff(t)*10e9)-mean(round(diff(t)*10e9)));
 if isDiscrete == 0
  warning('non-discrete',...
          'Input must consist of data sampled at evenly spaced points');
 end

% Set all input vectors as column vectors where applicable
t=t(:); y=y(:);

% Scale the amplitude of the y vector
y=y*ySF;

% Add offset to the y vector
y=y+yOffset;

% Scale the time course of the y vector
n=size(t,1);
n_scaled=round(n*tSF);
if tSF > 1
 % Upsampling by linear interpolation and antialias filtering
 t_scaled=linspace(min(t),max(t),n_scaled)';
 y_scaled=interp1q(t,y,t_scaled);
 [y_scaled]=binomialf(y_scaled,t_scaled,1,'on');
elseif tSF < 1
 % Downsampling by antialias filtering and linear interpolation
 order=1+round(1/tSF);
 [y]=binomialf(y,t,order,'on');
 t_scaled=linspace(min(t),max(t),n_scaled)';
 y_scaled=interp1q(t,y,t_scaled);
elseif tSF == 1
 t_scaled=t;
 y_scaled=y;
end
y=y_scaled;

% Resize scaled trace around the zero time point
post_size=numel(find(t_scaled>0));
if tSF > 1
 % Crop ends of the y vector
 post_crop=round(post_size/n_scaled*(n_scaled-n));
 y(n_scaled-post_crop+1:end)=[];
 y(1:n_scaled-n-post_crop)=[];
elseif tSF < 1
 % Pad ends of the y vector with zero values
 post_pad=round(post_size/n_scaled*(n-n_scaled));
 y=padarray(y,post_pad,0,'post');
 y=padarray(y,(n-n_scaled-post_pad),0,'pre');
end

