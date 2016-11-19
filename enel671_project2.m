%% ENEL 671 Project2
% Pouyan Keshavarzian
% FALL 2016
% https://github.com/pkeshava/enel671_Project2

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
for k=1:K
      a = BPSK(N);
      u = filterinput(a,h); 
      [e1,W1] = RLS_algorithm(M,u(:,1),a,delta,delta_d);
      ed1(:,k) = e1.^2;
      MSEE1 = sum(ed1,2)/K;
      [e2,W2] = RLS_algorithm(M,u(:,2),a,delta,delta_d);
      ed2(:,k) = e2.^2;
      MSEE2 = sum(ed2,2)/K;
      [e3,W3] = RLS_algorithm(M,u(:,3),a,delta,delta_d);
      ed3(:,k) = e3.^2;
      MSEE3 = sum(ed3,2)/K;
      [e4,W4] = RLS_algorithm(M,u(:,4),a,delta,delta_d);
      ed4(:,k) = e4.^2;
      MSEE4 = sum(ed4,2)/K;
end 
% Plot learning curves

for MSEE = [MSEE1 MSEE2 MSEE3 MSEE4]
    figure(1)
    semilogy(1:N,MSEE,'LineWidth',2)
    legend('Channel 1','Channel 2','Channel 3','Channel 4')
    grid on
    xlabel('Time (s)');
    ylabel('Mean Squared Error'); 
    title('Effect of Eigenvalue Spread');
    hold on
end
hold off
%% Effect of Filter Order

N = 500;
K = 200;
delta = 0.01;
for M = [9 11 21]
    delta_d =(M-1)/2 + (length(h(1,:))-1)/2;
    for k=1:K
        a = BPSK(N);
        u = filterinput(a,h);
        [e2_2,W2_2] = RLS_algorithm(M,u(:,2),a,delta,delta_d);
        ed2_2(:,k) = e2_2.^2;
        MSEE2_2 = sum(ed2_2,2)/K;
    end
    figure(2)
    semilogy(1:N,MSEE2_2,'LineWidth',2)
    legend('M = 9','M = 11','M = 21')
    grid on
    xlabel('Time (s)');
    ylabel('Mean Squared Error');
    title('Effect of Filter Order');
    hold on
end
hold off
%% Study of the Amplitude Spectra of the Adaptive Filter

clear all 
clc
N = 500;
K = 400;
M = 11;
h = [0.2194 1.0 0.2194;0.2798 1.0 0.2798;0.3365 1.0 0.3365;0.3887 1.0 0.3887];
delta_d = 5;
e = zeros(1,N);
Wtap_ave = zeros(M,N-M+1);
Wtap = zeros(M,N-M+1);

delta = 0.01;


for k=1:K
    
      a = BPSK(N);
      d = a; 
      % function below creats channel outputs and also adds white noise
      u = filterinput(a,h);
      u = u(:,1);
      d = d(:);
      e = zeros(1,N);     
      lamda = 1;
      P = delta^-1*eye(M);
      W = zeros(M,1);
      
for n=M:N
    
    
    % Define tap input vector with length of n-M+1 = 11
    u_vec = u(n:-1:n-M+1);
    % Compute Kalman gain vector
    Kal = P*u_vec/(lamda + u_vec'*P*u_vec);
    % Compute a priori error and update weight vector
    e(n) = d(n-delta_d)-W'*u_vec;
    W = W + Kal*e(n);
    % Update the inverse correlation matrix
    P = lamda^(-1)*(P - Kal*u_vec'*P);
    Wtap(:,n-M+1)= W;
    
end
   Wtap_ave = Wtap_ave + Wtap;
end

Wtap_ave = Wtap_ave/K;
WW = zeros(M,N);
WW(:,M:end) = Wtap_ave;
figure(3)
plot(WW(5,:),'r','LineWidth',2.5);
legend('Ch1')
grid on
xlabel('Time(n)')
ylabel('Tap-weight coefficient #5')

figure(4); 
stem(W,'color','r','LineWidth',2.5)
grid on
xlabel('Filter order(M)')
ylabel('Steady state values of tap-weight coefficients')

Ch = h(1,:);
W_f = 0;
H_f = 0;
freq = linspace(0,1,2000);
z = exp(j * 2 * pi * freq); 
%=================Magnitude spectrum of Adaptive filter W(f)===============
for i = 1:M
    X = W(i)*(z).^(-i); 
    W_f = X + W_f; 
end
   Magnitude_W_f = abs(W_f);
%==================Magnitude spectrum for Channel H(f)=====================
for j = 1:length(h(1,:))
    Y = Ch(j)*(z).^(-j);
    H_f = Y + H_f;
end
   Magnitude_H_f = abs(H_f);
   Magnitude_Cascade = ( Magnitude_W_f.* Magnitude_H_f);
%========================Plotings========================================== 
figure(5), plot( Magnitude_W_f,'LineWidth',2.5) 
grid on
hold on
plot( Magnitude_H_f,'LineWidth',2.5)
plot( Magnitude_Cascade,'LineWidth',2.5)
legend('Adaptive filter:|W(f)|','Second Channel:|H(f)|','Cascade:|W(f)H(f)|')
xlabel('Frequency(f)')
ylabel('Magnitude spectrum')

