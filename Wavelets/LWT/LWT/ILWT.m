function wt = ILWT(x, n, WT_steps, norm, rounding)
% s - smooth coefficients
% d - detail coefficients
% n is the data length, but it is not used in this implementation yet

lenx=length(x);

% limit the maximum number of transform stages 
steps=WT_steps;
if lenx<2^WT_steps
    steps=log2(lenx);
end

offset=lenx/2^(steps-1); % offset from the END of the data array
wt=x;
for j=1:steps
    wt=wt((lenx-offset+1):lenx);
    [s] = ILWT_stage(wt, rounding);  % compute LWT_stage
    wt = [x(1:(lenx-offset)) s];             % pack coefficients into the output stream
    offset=2*offset;
end

% also pack the last stage smooth coeffs

function [yr] = ILWT_stage(wt, rounding)
% s - smooth coefficients
% d - detail coefficients

len=length(wt);
d=wt(1:(len/2));
s=wt((len/2+1):len);

%% inverse transform
yr=zeros(1,length(wt));

for i=1:(len/2-1)
    yr(2*i) = s(i) - rounding((d(i)+d(i+1))/4);
end

% front boundary treatment: no division corresponds to symmetrical
% reflection over the first element, division by 2 correspond to prepending
% with zero
yr(1)= d(1) + rounding(yr(2) ); %/2);
% back boundary treatment: more symmetrical treatment is achieved by the
% factor of 2
yr(len)= s(len/2)-rounding(d(len/2)/ 2); %4);

for i=1:(len/2-1)
    yr(2*i+1) = d(i+1) + rounding((yr(2*i)+yr(2*i+2))/2);
end

return;
