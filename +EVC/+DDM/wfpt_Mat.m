% Modified by Sam Feng
%
% original code from Navarro Fuss 2009

function p=wfpt_Mat(t,x,v,a,z,err)

tt=t./(a.^2); % use normalized time
w=z./a; % convert to relative start point

v(x == 1) = -v(x == 1);
w(x == 1) = 1-w(x == 1);


% calculate number of terms needed for large t
idx = find(pi.*tt.*err<1); % if error threshold is set low enough
kl(idx) = sqrt(-2*log(pi*tt(idx)*err)./(pi^2*tt(idx))); % bound
kl(idx) = max([kl(idx);1./(pi*sqrt(tt(idx)))]); % ensure boundary conditions met

idx2 = find(pi.*tt.*err>=1);    % if error threshold set too high
kl(idx2) = 1./(pi*sqrt(tt(idx2)));  % set to boundary condition

% calculate number of terms needed for small t
idx = find(2*sqrt(2*pi.*tt)*err<1); % if error threshold is set low enough
ks(idx) = 2+sqrt(-2*tt(idx).*log(2*sqrt(2*pi*tt(idx))*err)); % bound

% if error threshold was set too high
ks(2*sqrt(2*pi.*tt)*err>=1)=2; % minimal kappa for that case


% compute f(tt|0,1,w)
p = zeros(1,length(t)); %initialize density
K = nan(1, length(t));

tmp.ks = ks;
tmp.kl = kl;
tmp.w = w;
tmp.tt = tt;
tmp.p = p;
tmp.K = K;

for i = 1:length(tmp.p)
    
    ks = tmp.ks(i);
    kl = tmp.kl(i);
    w = tmp.w(i);
    tt = tmp.tt(i);
    p = tmp.p(i);

    if ks<kl % if small t is better (i.e., lambda<0)
        K=ceil(ks); % round to smallest integer meeting error
        for k=-floor((K-1)/2):ceil((K-1)/2) % loop over k
            p=p+(w+2*k)*exp(-((w+2*k)^2)/2/tt); % increment sum
        end
        p=p/sqrt(2*pi*tt^3); % add constant term

    else % if large t is better...
        K=ceil(kl); % round to smallest integer meeting error
        for k=1:K
            p=p+k*exp(-(k^2)*(pi^2)*tt/2)*sin(k*pi*w); % increment sum
        end
        p=p*pi; % add constant term
    end
    
    tmp.p(i) = p;

end

p = tmp.p;
w = tmp.w;

% convert to f(t|v,a,w)
p=p.*exp(-v.*a.*w -(v.^2).*t/2)./(a.^2); 

