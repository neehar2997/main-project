tic;
clc
clear all

%square wave
ta = 20;
t = linspace(0,ta,ta*1e2);
ti = 4; %frequency of square wave
n = ta/ti;
u = (100/150)*0.5*(square(2*pi*(1/ti)*t,1)+1);
U = u;
subplot(2,1,1)
plot(t,12*u)
title('incident');
xlabel('time(ns)');
ylabel('amplitude(v)');

%transmission line parameters
zo = 100;
zs = 50;
zl = 200;
T = 0.2;
T1 = T;
rl = (zl-zo)/(zl+zo);
rs = (zs-zo)/(zs+zo);
p = length(t)/(n);

%vector
N = 20;
vec = zeros(N,1);
positions = [1:4];
vec(positions) = 1;
A1 = repmat(vec,2000);
v = A1(:,1);
B1 = v(1:ta*1e2)';
for k = 1:n
    A = (k-1)*p + 1;
    B = k*p;
    for j = 1:(ti/T1)/2
        u1 = (rl^j)*(rs^(j-1))*B1;
        for i = A : B
            if t(i)<T
                a(i) = u(i);
            else if t(i) >T
                    a(i) = U(i) + u1(i);
                end
            end
        end
        T = T+0.2;
        u2 = (rl^j)*(rs^j)*B1;
        for i = A:B
            if t(i) < T
                b(i) = a(i);
            else if t(i) > T
                    b(i) = U(i) + u2(i);
                end
            end
        end
        T = T+0.2;
        u = b;
    end
    u = U;
end
subplot(2,1,2)
plot(t,12*b)
title('after multiple reflections')
xlabel('time(ns)')
ylabel('ampitude(v)');
x = 4.2;
v = b-U;
for i = 1:length(t)
    if t(i) - x <=2e-3
        s(i) = b(i)-U(i);
    end
end

%square wave generation
Fs = 1/x;
hold on
p = linspace(0,ta*pi,length(t));
s1 = square(Fs*p); %sampling frequency
for i = i:length(s1)
    if s1(i) == 1
        s2(i) = -1;
    else if s1(i) == -1
            s2(i) = 1;
        end
    end
end
hold on
plot(p/(pi),s1)
 
%sampling action
for i = 1:length(s1)
    if i == length(t)
        v1(i) = v(i);
    else if s1(i+1) + s1(i) == 0
        v1(i+1) = v(i+1);
        end
    end
end
v1 = 8*v1;
figure
plot(t,v1)
title('sampled sequence');
xlabel('n');
ylabel('x[n]');
 k=nonzeros(v1);
 k1=repelem(k,x*1e2);
 k2=find(v1~=0);
 k4=[rl;k1];
 k5 = zeros(1,ta*1e2);
 k5(1:1,1:size(k4)) = k4;
 k6 = (k5+1)*1e15;
 figure
 plot(t,k6)
 title('reconstructed signal');
 xlabel('time(ns)');
 ylabel('x(t)');
 toc;
 

%a = sin((2*pi*(10/ti)*t));
%figure
%plot(t,a)
%a = a.*k5;
%figure
%plot(t,a)