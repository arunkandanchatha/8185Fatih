%% plot value function and policy
close all; clear all; clc
cd('C:\Users\arun\SkyDrive\git\8185Fatih\hw3');
rawData1 = importdata('EZ_policyr2e2');
g1a=rawData1;
rawData1 = importdata('EZ_valuer2e2');
v1a=rawData1;

rawData1 = importdata('EZ_policyr2e0p1');
g2a=rawData1;
rawData1 = importdata('EZ_valuer2e0p1');
v2a=rawData1;

rawData1 = importdata('EZ_policyr10e2');
g3a=rawData1;
rawData1 = importdata('EZ_valuer10e2');
v3a=rawData1;

rawData1 = importdata('EZ_policyr10e0p1');
g4a=rawData1;
rawData1 = importdata('EZ_valuer10e0p1');
v4a=rawData1;

rawData1 = importdata('EZ_policyCRRA2');
g5a=rawData1;
rawData1 = importdata('EZ_valueCRRA2');
v5a=rawData1;

rawData1 = importdata('EZ_policyCRRA10');
g6a=rawData1;
rawData1 = importdata('EZ_valueCRRA10');
v6a=rawData1;

rawData1 = importdata('CRRA_policy2');
g5b=rawData1;
rawData1 = importdata('CRRA_value2');
v5b=rawData1;

rawData1 = importdata('CRRA_policy10');
g6b=rawData1;
rawData1 = importdata('CRRA_value10');
v6b=rawData1;

rawData1 = importdata('asset');
a=rawData1;
rawData1 = importdata('gamma');
gamma=rawData1;
rawData1 = importdata('shock');
s=rawData1;

figure
subplot(2,1,1)
for i=1
plot(a,[v1a(end-i,:)])
hold on
end
title('EZ RRA=2 EIS=2\newline Value')

subplot(2,1,2)
for i=1
plot(a,[g1a(end-i,:)])
hold on
end
plot(a,a,'-g')
title('Policy')
hold off

figure
subplot(2,1,1)
for i=1
plot(a,[v2a(end-i,:)])
hold on
end
title('EZ RRA=2 EIS=0.1\newline Value')

subplot(2,1,2)
for i=1
plot(a,[g2a(end-i,:)])
hold on
end
plot(a,a,'-g')
title('Policy')
hold off

figure
subplot(2,1,1)
for i=1
plot(a,[v3a(end-i,:)])
hold on
end
title('EZ RRA=10 EIS=2\newline Value')

subplot(2,1,2)
for i=1
plot(a,[g3a(end-i,:)])
hold on
end
plot(a,a,'-g')
title('Policy')
hold off

figure
for i=1
plot(a,[g5b(end-i,:)-g5a(end-i,:)])
end
title('Policy function difference, RRA=10')

figure
for i=1
plot(a,[g6b(end-i,:)-g6a(end-i,:)])
end
title('Policy function difference, RRA=10')
