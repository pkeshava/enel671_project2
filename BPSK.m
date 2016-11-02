function   a = BPSK(N)
%% REQUIRES 
    % N: data vector length
%RETURNS
    % a: randomly generated BPSK data
%%
a = round(rand(1,N));
for i=1:N
    if a(i) == 0
        a(i) = -1;
    end
end
end