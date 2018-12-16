function [dpdf,meanRT,phit] = ddmpdf_Mat(T,x,A,T0,x0,z,sigma)

import EVC.DDM.*;

if any(x<0) || any(x>1)
   error; 
end

s = sigma;
nT = length(T);
dpdf = zeros(length(T0),nT);
errtol = 1e-4;
for k = 1:nT
    t = T(k)-T0;
    dpdf(t < 0, k) = 0;
    if(any(t >= 0))
        dpdf(t >= 0, k) = wfpt_Mat(t(t >= 0), x(t >= 0), A(t >= 0)./s(t >= 0), 2.*z(t >= 0)./s(t >= 0), (x0(t >= 0)+z(t >= 0))./s(t >= 0), errtol);
    end
end



flipped = zeros(size(A));
flipped(A<0) = 1;
x0(A<0) = -x0(A < 0);
A(A < 0) = -A(A < 0);

zz = z./A;
aa = (A./s).^2;
xx = x0./A;

ER = 1 ./ (1 + exp(2.*zz.*aa)) - ...
     (1 - exp(-2.*xx.*aa))./(exp(2.*zz.*aa)-exp(-2.*zz.*aa));

phit = nan(size(x));
phit(x == 0) = ER(x == 0);
phit(x ~= 0) = 1 - ER(x ~= 0);

phit(logical(flipped)) = 1 - phit(logical(flipped));

% Compute expected value, careful to normalize
dt = T(2)-T(1);
I = sum(trapz(dpdf,2).*dt',2);
I = repmat(I, 1, size(dpdf,2));
cpdf = (1./I) .* dpdf;

meanRT = sum(trapz(repmat(T,size(cpdf,1),1).*cpdf,2).*dt,2);