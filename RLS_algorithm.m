function [e,W] = RLS_algorithm(M,u,d,pre,delta)
%% REQUIRES
    
    % u: tap-input vector
    % d: desired signal vector
    % mu: step size parameter
    % delta: delay of channel
    % M: filter order
% Returns   
    % W: Tap weight vector
    % e: estimation error 
    
W = zeros(M,1);
N = length(u);
u = u(:);
d = d(:);
lamda = 1;
P = pre^-1*eye(M);
for n=M:N
    % Define tap input vector with length of n-M+1 = 11
    u_vec = u(n:-1:n-M+1);
    % Compute Kalman gain vector
    Kal = P*u_vec/(lamda + u_vec'*P*u_vec);
    % Update the inverse correlation matrix
    P = lamda^(-1)*(P - Kal*u_vec'*P);
    % Compute a priori error and update weight vector
    e(n) = d(n-delta)-W'*u_vec;
    W = W + Kal*e(n);

end
e = e(:);
end