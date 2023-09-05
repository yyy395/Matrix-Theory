clear; close all;
N = 8;               % 阵元个数        
M = 3;               % 信源数目
J = 5;
J_N=N-J+1;
theta = [-pi/6,0,pi/3];  % 待估计角度
snr = 20;            % 信噪比
K = 1024;             % 快拍数
fs = 1000;
Ts = 0.001;
T1 = Ts*(K-1);
T = 0:Ts:T1;
dd = 0.5;            % 阵元间距 
d = 0:dd:(N-1)*dd;
S = sin(100*pi*T);
A = exp(-1j*2*pi*d'*sin(theta));
X = zeros(N,K);
for i = 1:N
    for j = 1:length(theta)
        X(i,:)=X(i,:)+S*A(i,j);
    end
end
X1 = awgn(X,snr);
% 计算协方差矩阵
R = X1*X1'/N;                                   % 原始协方差矩阵
Rf = zeros(J_N, J_N);                           %运用空间平滑算法重新构造协方差矩阵
for i = 1:J
    Rf = Rf+R(i:i+J_N-1,i:i+J_N-1);
end
Rf = Rf/J;
[V,D] = eig(Rf);     %特征值分解
Uw=V(:,1:J_N-M);      %噪声子空间
P = zeros(1,181);
w = -pi/2:pi/180:pi/2;
theta1 = -90:1:90;
for i = 1:length(w)       % 遍历每个角度，计算空间谱
    a = exp(-1j*2*pi*d(1:J_N)'*sin(w(i)));
    P(i) = 1/(a'*(Uw*Uw')*a);
end
P = abs(P);
Pmax=max(P);
P = 10*log10(P/Pmax);   %功率转化为dB
plot(theta1,P);
title('空间平滑算法MUSIC');
xlabel('入射角/(degree)');
ylabel('功率/(dB)');