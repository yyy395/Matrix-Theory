clear; close all;
N = 8;               % ��Ԫ����        
M = 1;               % ��Դ��Ŀ
theta = -pi/6;  % �����ƽǶ�
snr = -5;            % �����
K = 1024;             % ������
fs = 1000;
Ts = 0.001;
T1 = Ts*(K-1);
T = 0:Ts:T1;
dd = 0.5;            % ��Ԫ��� 
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
% ����Э�������
R2 = X1*X1'/K;
R3 = Y*Y'/K;
R = R2+R3;
[V,D] = eig(R);     %����ֵ�ֽ�
Uw=V(:,1:N-M);
D1 = diag(D);
P = zeros(1,181);
w = -pi/2:pi/180:pi/2;
theta1 = -90:1:90;
for i = 1:N-M
    Uw(:,i)=Uw(:,i)*(D1(i)^(1));
end
for i = 1:length(w)       % ����ÿ���Ƕȣ�����ռ���
    a = exp(-1j*2*pi*d'*sin(w(i)));
    P(i) = 1/(a'*Uw*Uw'*a);
end
P = abs(P);
Pmax=max(P);
P = 10*log10(P/Pmax);   %����ת��ΪdB
plot(theta1,P);
xlabel('�����/(degree)');
ylabel('����/(dB)');
grid on;