clear;
close all;
N = 8;               % ��Ԫ����        
M = 1;               % ��Դ��Ŀ
theta = -pi/6;     % �����źŴ����ƽǶ�
%M = 2;               % ��Դ��Ŀ
%theta = [-pi/6,pi/3];  % �����źŴ����ƽǶ�
snr = 20;            % �����
K = 1024;             % ������
fs = 1000;
Ts = 0.001;
T1 = Ts*(K-1);
T = 0:Ts:T1;
dd = 0.5;            % ��Ԫ��� 
d = 0:dd:(N-1)*dd;
S = sin(100*pi*T);      %��Դ�ź�
%S2 = cos(100*pi*T);    %������ź�
A = exp(-1j*2*pi*d'*sin(theta)*50/1500);   %��������
X = zeros(N,K);
for i = 1:N
    X(i,:)=S*A(i,1);                  %�����ź�
    %X(i,:)=S*A(i,1)+S*A(i,2);         %��������ź�
    %X(i,:)=S*A(i,1)+S2*A(i,2);       %����������ź�
end
X1 = awgn(X,snr);   %���������ʹ�������Ϊ20dB
R = X1*X1'/K;       %����Э�������
[V,D] = eig(R);     %����ֵ�ֽ�
Uw=V(:,1:N-M);      %�����ӿռ�
P = zeros(1,181);
w = -pi/2:pi/180:pi/2;
theta1 = -90:1:90;
for i = 1:length(w)       % ����ÿ���Ƕȣ�����ռ���
    a = exp(-1j*2*pi*d'*sin(w(i))*50/1500);
    P(i) = 1/(a'*(Uw*Uw')*a);
end
P = abs(P);
[Pmax,index]=max(P);
P = 10*log10(P/Pmax);   %��һ������ת��ΪdB
plot(theta1,P);
title('MUSIC');
xlabel('�����/(degree)');
ylabel('����/(dB)');
grid on;

