 function u= filterinput(a,h)
 
% Takes values of autocorrelation matrix and a noise matrix and calculates
% the overall autocorellation along with the eigenvalues.
u = zeros(size(a,2),4);
for i = 1:4
    %define output of channel as a filtered output of the rational transfer
    %function defined by the channel transfer function
    x= filter(h(i,:),1,a);
    % add noise
    % output of channel/input of filter defined as channel output +
    % additive white gaussian noise
    u(:,i) = awgn(x,40); 
end