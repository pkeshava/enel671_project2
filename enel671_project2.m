%% ENEL 671 Project2
% Pouyan Keshavarzian
% FALL 2016
% https://github.com/pkeshava/enel671_Project2
% LMS algorithm
    % Filter order: M
    % Machine Precision Inialization: 0.01
    % vector data length: N
    % channel delay + adaptive filter delay: delta
    % number of runs: K
    
clear all
close all
clc
%% Effect of Eignevalue Spread
M = 11;
delta = 0.01;
N = 600;
h = [0.2194 1.0 0.2194;0.2798 1.0 0.2798;0.3365 1.0 0.3365;0.3887 1.0 0.3887];
delta_d = (M-1)/2 + (length(h(1,:))-1)/2;
K = 350;
e = zeros(1,N);
%P = delta*eye(M);
%W = zeros(M,1);
for k=1:K
      a = BPSK(N);
      u = filterinput(a,h);
      [e,W]= RLS_algorithm(M,u,a,delta,delta_d); 
      ed(:,k)= e.^2;
end
    % averaging error over the number of runs K
    MSEe=sum(ed,2)/K;
    figure(1),semilogy(1:N,MSEe,'LineWidth',2),legend('Ch1','Ch2','Ch3','Ch4'),grid on
    xlabel('Time(n)'),ylabel('Mean Squared Error(MSE)'),hold on 

    % Calculate u(n)

    % Run data vector through LMS algorithm for each channel and get calculated
    % error and weight vector
    