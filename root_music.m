clear; close all;
N = 8;               % ��Ԫ����        
M = 1;               % ��Դ��Ŀ
theta = -pi/6;  % �����ƽǶ�
snr = 20;            % �����
K = 3;             % ������
fs = 1000;
Ts = 0.001;
T1 = Ts*(K-1);
T = 0:Ts:T1;
dd = 0.5;            % ��Ԫ��� 
d = 0:dd:(N-1)*dd;
S = sin(100*pi*T);
A = exp(-1j*2*pi*d'*sin(theta));   %��������
X = zeros(N,K);
for i = 1:N
    X(i,:)=S*A(i,1);
end
X1 = awgn(X,snr);   %���������ʹ�������Ϊ20dB
R = X1*X1'/K;
[V,D] = eig(R);     %����ֵ�ֽ�
Uw=V(:,1:N-M);      %�����ӿռ�
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