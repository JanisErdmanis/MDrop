function [p,t] = meshgen() %,boxsize)
    function f = fd(p)
    f = p(:,1).^2/1+p(:,2).^2/1+p(:,3).^2/1-1;
end


    
    [p,t]=distmeshsurface(@fd,@huniform,0.2,[-1.1 -1.1 -1.1
 1.1 1.1 1.1]);
end

function h=huniform(p,varargin)
    h=ones(size(p,1),1);
end

