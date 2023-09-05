clear; close all;
N = 8;               % 阵元个数        
M = 1;               % 信源数目
theta = -pi/6;  % 待估计角度
snr = -5;            % 信噪比
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
    X(i,:)=X(i,:)+S*A(i);
end
X1 = awgn(X,snr);
I=eye(N);
I1=rot90(I);
X2=conj(X1);
Y=I1*X2;
% 计算协方差矩阵
R2 = X1*X1'/K;
R3 = Y*Y'/K;
R = R2+R3;
[V,D] = eig(R);     %特征值分解
Uw=V(:,1:N-M);
D1 = diag(D);
P = zeros(1,181);
w = -pi/2:pi/180:pi/2;
theta1 = -90:1:90;
for i = 1:N-M
    Uw(:,i)=Uw(:,i)*(D1(i)^(1));
end
for i = 1:length(w)       % 遍历每个角度，计算空间谱
    a = exp(-1j*2*pi*d'*sin(w(i)));
    P(i) = 1/(a'*Uw*Uw'*a);
end
P = abs(P);
Pmax=max(P);
P = 10*log10(P/Pmax);   %功率转化为dB
plot(theta1,P);
xlabel('入射角/(degree)');
ylabel('功率/(dB)');
grid on;