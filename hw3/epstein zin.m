clc
clear all
close all
%% form of the value function, guess and verify
syms w psi alpha mu j(w) c v1 v2 rra

% what is the form of j
dsolve(diff(j,2)/diff(j,1)*w == 1/(1-alpha))
% verify
j=w^((alpha - 2)/(alpha - 1))
v1=diff(j,w,1)
v2=diff(j,w,2)
rra=-w*v2/v1
simplify(rra)

%% plot value function and policy
close all; clear all; clc
cd('C:\Users\arun\SkyDrive\workspaces\hw3');
rawData1 = importdata('policy1a');
g1a=rawData1;
rawData1 = importdata('value1a');
v1a=rawData1;
rawData1 = importdata('policy1b');
g1b=rawData1;
rawData1 = importdata('value1b');
v1b=rawData1;

rawData1 = importdata('policy2a');
g2a=rawData1;
rawData1 = importdata('value2a');
v2a=rawData1;
rawData1 = importdata('policy2b');
g2b=rawData1;
rawData1 = importdata('value2b');
v2b=rawData1;

rawData1 = importdata('policy3a');
g3a=rawData1;
rawData1 = importdata('value3a');
v3a=rawData1;
rawData1 = importdata('policy3b');
g3b=rawData1;
rawData1 = importdata('value3b');
v3b=rawData1;

rawData1 = importdata('policy4a');
g4a=rawData1;
rawData1 = importdata('value4a');
v4a=rawData1;
rawData1 = importdata('policy4b');
g4b=rawData1;
rawData1 = importdata('value4b');
v4b=rawData1;


rawData1 = importdata('asset');
a=rawData1;
rawData1 = importdata('gamma');
gamma=rawData1;
rawData1 = importdata('shock');
s=rawData1;
%%
figure
for i=0:2
plot(a,[v1a(end-i,:)])
hold on
end
hold off

figure
for i=0:2
plot(a,[g1a(end-i,:)])
hold on
end
plot(a,a,'-g')
hold off

figure
for i=0:2
plot(a,[v1b(end-i,:)])
hold on
end
hold off

figure
for i=0:2
plot(a,[g1b(end-i,:)])
hold on
end
plot(a,a,'-g')
hold off

%%
figure
for i=0:2
plot(a,[v2a(end-i,:)])
hold on
end
hold off

figure
for i=0:2
plot(a,[g2a(end-i,:)])
hold on
end
plot(a,a,'-g')
hold off

figure
for i=0:2
plot(a,[v2b(end-i,:)])
hold on
end
hold off

figure
for i=0:2
plot(a,[g2b(end-i,:)])
hold on
end
plot(a,a,'-g')
hold off

%%
figure
for i=0:2
plot(a,[v3a(end-i,:)])
hold on
end
hold off

figure
for i=0:2
plot(a,[g3a(end-i,:)])
hold on
end
plot(a,a,'-g')
hold off

figure
for i=0:2
plot(a,[v3b(end-i,:)])
hold on
end
hold off

figure
for i=0:2
plot(a,[g3b(end-i,:)])
hold on
end
plot(a,a,'-g')
hold off

%%
figure
for i=0:2
plot(a,[v4a(end-i,:)])
hold on
end
hold off

figure
for i=0:2
plot(a,[g4a(end-i,:)])
hold on
end
plot(a,a,'-g')
hold off

figure
for i=0:2
plot(a,[v4b(end-i,:)])
hold on
end
hold off

figure
for i=0:2
plot(a,[g4b(end-i,:)])
hold on
end
plot(a,a,'-g')
hold off

%%
figure
plot(a,[g(end,:)])
hold on
plot(a,a)

figure
plot(a,[v(end,:)])

%%
beta=0.9
r=0.1
rho=0.7
u=beta*(1+r)
alpha=1/(rho-1)
s=u^alpha/(1+u^alpha)
%%
syms x a y 