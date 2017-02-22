%     Function File: [Y, X] = medianf (y, x, r)
%
%     Calculates the median over a sliding window of y-values of 2r+1
%     number of points, where r is the filter rank. In order to prevent
%     exhausting the memory, this function segments the data before
%     implementing an efficient median algorithm block-by-block and
%     recompiling the data. The algorithm performs only one sort for
%     every data block and stores the original indices to recall sorted
%     data values for the sliding window. Bounce end-effect correction
%     is used by default.
%
%     This function requires the 'sma' and 'bounce' functions.
%
%
%     Bibliography:
%     Moore and Jorgenson (1993) Anal Chem 6: 188-191
%     Friedrichs (1995) J Biomol NMR 5: 147-153
%
%     medianf v1.0 (last updated: 05/04/2012)
%     Author: Andrew Charles Penn
%     https://www.researchgate.net/profile/Andrew_Penn/


function [Y, X] = medianf (y, x, r)

if nargin ~= 3
 error('Invalid number of input arguments');
end

if all(size(x) == 1) || ~any(size(x) == 1) || any(size(x)~=size(y))
 error('x and y must be vectors of the same size');
end

if isinf(r) || ~all(size(r) == 1) || r<=0 || r~=round(r)
 error('r must be a nonnegative integer');
end

% Assess sampling characteristics of input with precision of 10e-9
isDiscrete=~any(round(diff(x)*10e9)-mean(round(diff(x)*10e9)));
if isDiscrete == 0
 warning('non-discrete',...
         'Input must consist of data sampled at evenly spaced points');
end

% Set all input vectors as column vectors and calculate vector size
x=x(:); y=y(:);

% Calculate total number of y-points to average within sliding box
p=2*r+1;

% Bounce end-effect correction
[y]=bounce(y,r);
X=x;
N=length(y);

% Block-by-block median filter algorithm
Y=[];
Y_block_matrix=zeros(100,p);
blockSize=100+2*r;
while N > blockSize
 y_block=y(1:blockSize);
 [y_block_sorted, idx] = sort(y_block);
 for i=1:100
  idx_window=find(idx >= i & idx <= 2*r+i);
  y_block_matrix(i,:)=y_block_sorted(idx_window);
 end
 Y_block=y_block_matrix(:,r+1);
 Y=cat(1,Y,Y_block);
 y(1:100)=[];
 N=length(y);
end

% Median filter algorithm on final block
clear y_block_sorted idx idx_window
y_block_matrix=zeros(N-2*r,p);
y_block=y(1:N);
[y_block_sorted, idx] = sort(y_block);
for i=1:N-2*r
 idx_window=find(idx >= i & idx <= 2*r+i);
 y_block_matrix(i,:)=y_block_sorted(idx_window);
end
Y_block=y_block_matrix(:,r+1);
Y=cat(1,Y,Y_block);
