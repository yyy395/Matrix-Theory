clear; close all;
N = 8;               % 阵元个数        
M = 1;               % 信源数目
theta = -pi/6;  % 待估计角度
snr = 20;            % 信噪比
K = 3;             % 快拍数
fs = 1000;
Ts = 0.001;
T1 = Ts*(K-1);
T = 0:Ts:T1;
dd = 0.5;            % 阵元间距 
d = 0:dd:(N-1)*dd;
S = sin(100*pi*T);
A = exp(-1j*2*pi*d'*sin(theta));   %方向向量
X = zeros(N,K);
for i = 1:N
    X(i,:)=S*A(i,1);
end
X1 = awgn(X,snr);   %添加噪声，使得信噪比为20dB
R = X1*X1'/K;
[V,D] = eig(R);     %特征值分解
Uw=V(:,1:N-M);      %噪声子空间
syms z;
pz = z.^([0:N-1]');
pz1 = (z^(-1)).^([0:N-1]);
fz = z.^(N-1)*pz1*Uw*Uw'*pz;
a=sym2poly(fz);
zx=roots(a);
zx1=zx';
[as,ad]=(sort(abs((abs(zx1)-1))));
angleest=asin((angle(zx1(ad(1)))/pi))*180/pi;
disp(angleest);