clear; close all;
N = 8;               % ��Ԫ����        
M = 3;               % ��Դ��Ŀ
J = 5;
J_N=N-J+1;
theta = [-pi/6,0,pi/3];  % �����ƽǶ�
snr = 20;            % �����
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
    for j = 1:length(theta)
        X(i,:)=X(i,:)+S*A(i,j);
    end
end
X1 = awgn(X,snr);
% ����Э�������
R = X1*X1'/N;                                   % ԭʼЭ�������
Rf = zeros(J_N, J_N);                           %���ÿռ�ƽ���㷨���¹���Э�������
for i = 1:J
    Rf = Rf+R(i:i+J_N-1,i:i+J_N-1);
end
Rf = Rf/J;
[V,D] = eig(Rf);     %����ֵ�ֽ�
Uw=V(:,1:J_N-M);      %�����ӿռ�
P = zeros(1,181);
w = -pi/2:pi/180:pi/2;
theta1 = -90:1:90;
for i = 1:length(w)       % ����ÿ���Ƕȣ�����ռ���
    a = exp(-1j*2*pi*d(1:J_N)'*sin(w(i)));
    P(i) = 1/(a'*(Uw*Uw')*a);
end
P = abs(P);
Pmax=max(P);
P = 10*log10(P/Pmax);   %����ת��ΪdB
plot(theta1,P);
title('�ռ�ƽ���㷨MUSIC');
xlabel('�����/(degree)');
ylabel('����/(dB)');