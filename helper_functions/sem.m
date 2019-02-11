function [y,x] = nanstderr(varargin)

x = varargin{1};

[d1 d2] = size(x);

%determine if x is 1d vector
if d1==1 | d2==1
    %answer is a scalar
    x = x(~isnan(x));
    y = nanstd(x) / sqrt(length(x));
else
    if nargin == 1
        dim = 1;
    else
        dim = varargin{2};
    end
    y = nanstd(x,[],dim) ./ sqrt(sum(~isnan(x),dim));
end

end


