%% plot value function and policy
close all; clear all; clc
cd('C:\Users\arun\SkyDrive\workspaces\hw4q1');
rawData1 = importdata('policyR2E2');
g1a=rawData1;

rawData1 = importdata('asset');
a=rawData1;

figure
plot(a,[g1a(end-1,:)])
hold on
plot(a,a,'-g')
xlim([-4 20])
