clear;
close all;
N = 8;               % 阵元个数        
M = 1;               % 信源数目
theta = -pi/6;     % 单个信号待估计角度
%M = 2;               % 信源数目
%theta = [-pi/6,pi/3];  % 两个信号待估计角度
snr = 20;            % 信噪比
K = 1024;             % 快拍数
fs = 1000;
Ts = 0.001;
T1 = Ts*(K-1);
T = 0:Ts:T1;
dd = 0.5;            % 阵元间距 
d = 0:dd:(N-1)*dd;
S = sin(100*pi*T);      %信源信号
%S2 = cos(100*pi*T);    %不相干信号
A = exp(-1j*2*pi*d'*sin(theta)*50/1500);   %方向向量
X = zeros(N,K);
for i = 1:N
    X(i,:)=S*A(i,1);                  %单个信号
    %X(i,:)=S*A(i,1)+S*A(i,2);         %两个相干信号
    %X(i,:)=S*A(i,1)+S2*A(i,2);       %两个不相干信号
end
X1 = awgn(X,snr);   %添加噪声，使得信噪比为20dB
R = X1*X1'/K;       %计算协方差矩阵
[V,D] = eig(R);     %特征值分解
Uw=V(:,1:N-M);      %噪声子空间
P = zeros(1,181);
w = -pi/2:pi/180:pi/2;
theta1 = -90:1:90;
for i = 1:length(w)       % 遍历每个角度，计算空间谱
    a = exp(-1j*2*pi*d'*sin(w(i))*50/1500);
    P(i) = 1/(a'*(Uw*Uw')*a);
end
P = abs(P);
[Pmax,index]=max(P);
P = 10*log10(P/Pmax);   %归一化功率转化为dB
plot(theta1,P);
title('MUSIC');
xlabel('入射角/(degree)');
ylabel('功率/(dB)');
grid on;

