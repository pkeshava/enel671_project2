function   [e,W]= LMS_algorithm(u,d,mu,delta,M);
%% REQUIRES
    
    % u: tap-input vector
    % d: desired signal vector
    % mu: step size parameter
    % delta: delay of channel
    % M: filter order
% Returns    
    % W: Tap weight vector
    % e: estimation error (Nx1)
    
W = zeros(M,1);
N = length(u);
u = u(:);
d = d(:);
for n=M:N
    % Define tap input vector with length of n-M+1
    u_vec = u(n:-1:n-M+1);
    e(n) = d(n-delta)-W'*u_vec;
    W = W + mu*u_vec*e(n);
end
e = e(:);
end