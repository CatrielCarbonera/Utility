function [ wt ] = LWT( x, n , WT_steps, norm ,rounding) 
% LWT is the linear wavelet transform also known as CDF(2,2)
%  This code is based on the lifting scheme
% s - smooth samples at various scales, also storage for the accumulated detail
% coefficients
% d - detail at various scale

% limit the maximum number of transform stages 
steps=WT_steps;
if log2(length(x))<WT_steps
    steps=log2(length(x));
end
s=x;    
wt=[];
for j=1:steps
    [d, s] = LWT_stage(s, rounding);   % compute LWT_stage
    wt = [wt d];             % pack coefficients into the output stream
end

% also pack the last stage smooth coeffs
wt=[wt s];


function [d, s]= LWT_stage(y, rounding)
% s - smooth coefficients
% d - detail coefficients
leny=length(y);

for i=1:(leny/2)
% dual lifting
    if i==1
       % front boundary treatment equivalent to inserting a zero
       d(1) = y(1)-rounding( ( 0 + y(2))); %/2);
    else
       d(i) = y(2*i-1) - rounding((y(2*i) + y(2*i-2))/2);
    end
end

for i=1:(leny/2)
% primal lifting
    if i<leny/2
        s(i) = y(2*i) + rounding( (d(i) + d(i+1))/4);
    else
        % end boundary treatment equivalent to inserting a zero
        s(i) = y(2*i) + rounding((d(i) + 0 )/2); %4);
    end
end
