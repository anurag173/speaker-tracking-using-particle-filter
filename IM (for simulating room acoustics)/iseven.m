function [res] = iseven(val)
%function [res] = iseven(val)
%  What do you think this one does???
%  Results in 1 if val is even, 0 otherwise!...

valby2 = ceil(val/2);
if val==2*valby2
   res=1;
else
   res=0;
end

