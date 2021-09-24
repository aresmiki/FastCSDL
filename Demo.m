%% demo 1
%% Code by Liu He
%    aresmiki@163.com
%    SWJTU
%%
clc
clear all
close all
K=30;
Fs=10000;
N=10000;
t=(0:1:N-1)/Fs;
%% load mix-signal, pure-signal
load('Simdata.mat');%odata
load('Simpure.mat');%x
%% 学习
Lq=20;
[dic1,xcofe1,rec1]=CSC_MOD(odata,Lq,K,100); %LoCOMP
figure
plot(t,x,'LineWidth',1)
hold on
plot(t,rec1,'-.','LineWidth',1)
xlabel('Time [s]','fontsize',12)
ylabel('Amplitude','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gca,'Box','on');
set(gcf,'position',[200,300,300,120]);
ylim([-1,2])
set(gca,'Ytick',[-1,0,1,2]);
l1=legend('Simulation fault signal','Reconstructed signal')
set(l1,'Fontname', 'Times New Roman','FontSize',12)
set(l1,'Box','off');

figure
plot(t,x,'LineWidth',1)
hold on
plot(t,rec1,'-.','LineWidth',1)
xlabel('Time [s]','fontsize',12)
ylabel('Amplitude','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gca,'Box','on');
set(gcf,'position',[200,300,400,150]);
xlim([0.2,0.206])
ylim([-1,2])
set(gca,'Ytick',[-1,0,1,2]);
l1=legend('Simulation fault signal','Reconstructed signal')
set(l1,'Fontname', 'Times New Roman','FontSize',12)
set(l1,'Box','off');
%%  UC-DLA 使用了SHIFT-INVARIANCE工具包 请用本地路径替换
addpath(genpath('/Users/heliu/Learning Algorithms /SHIFT-INVARIANCE-master'))
np=Lq;
close all
% UC-DLA
randn('seed',1)
Y=odata;
Ya=[];
npq=1000;
for i=1:floor(N/npq)
    Ya(:,i)=Y(((i-1)*npq+1):i*npq+0);
end
L=1;
m = size(Ya,1)+1-np;
[Dconvsu, Xuconvsu] = uconvdlasu(Ya,round(K/size(Ya,2)), L, np, m); %仿真信号10
Ruconvsu=Dconvsu*Xuconvsu;
Ry=Ruconvsu(:);

figure
plot(t(1:length(x)),x,'LineWidth',1)
hold on
for i=1:size(Ya,2)-1
    plot(ones(100,1)*i*npq*(1/Fs),-2+0.04:0.04:2,'r-.','LineWidth',1)
end
xlabel('Time [s]','fontsize',12)
ylabel('Amplitude','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gca,'Box','on');
set(gcf,'position',[200,300,400,150]);
ylim([-1,2])
set(gca,'Ytick',[-1,0,1,2]);


figure
plot(t(1:length(x)),odata,'LineWidth',1)
hold on
for i=1:size(Ya,2)-1
    plot(ones(100,1)*i*npq*(1/Fs),-2+0.04:0.04:2,'r-.','LineWidth',1)
end
xlabel('Time [s]','fontsize',12)
ylabel('Amplitude','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gca,'Box','on');
set(gcf,'position',[200,300,300,120]);
ylim([-1,2])
set(gca,'Ytick',[-1,0,1,2]);


figure
plot(t(1:length(x)),x,'LineWidth',1)
hold on
plot(t(1:length(Ry)),Ry,'-.','LineWidth',1)
xlabel('Time [s]','fontsize',12)
ylabel('Amplitude','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gca,'Box','on');
set(gcf,'position',[200,300,400,150]);
ylim([-1,2])
set(gca,'Ytick',[-1,0,1,2]);
l1=legend('Simulation fault signal','Reconstructed signal')
set(l1,'Fontname', 'Times New Roman','FontSize',12)
set(l1,'Box','off');
xlim([0.2,0.206])

figure
plot(t(1:length(x)),x,'LineWidth',1)
hold on
plot(t(1:length(Ry)),Ry,'-.','LineWidth',1)
xlabel('Time [s]','fontsize',12)
ylabel('Amplitude','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gca,'Box','on');
set(gcf,'position',[200,300,400,150]);
ylim([-1,2])
set(gca,'Ytick',[-1,0,1,2]);
l1=legend('Simulation fault signal','Reconstructed signal')
set(l1,'Fontname', 'Times New Roman','FontSize',12)
set(l1,'Box','off');

%% cbpdndl方法
% 使用了sporco工具包，请用本地路径替换
addpath(genpath('/Users/heliu/Learning Algorithms /sporco-m0.0.9'))

clc
npq=10000;
Y=odata;
Ya=[];
N=length(Y);
for i=1:floor(N/npq)
    Ya(:,i)=Y(((i-1)*npq+1):i*npq+0);
end

% Set up cbpdndl parameters
lambda = 0.5;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 250;
opt.rho = 10*lambda;
opt.sigma = size(Ya,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
n=20;
L=1;
% Construct initial dictionary
D0 = zeros(n,1,L);
D0(1:n,1,:) = randn(n,1,L);
% Do dictionary learning
[D, X, optinf] = cbpdndl(D0, Ya, lambda, opt);
% Compute reconstruction
Ry = ifft2(sum(bsxfun(@times, fft2(D, size(X,1), size(X,2)), fft2(X)),3), ...
           'symmetric');
       
figure
plot(t(1:length(x)),x,'LineWidth',1)
hold on
plot(t(1:length(Ry(:))),Ry(:),'-.','LineWidth',1)
xlabel('Time [s]','fontsize',12)
ylabel('Amplitude','fontsize',12)
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gca,'Box','on');
set(gcf,'position',[200,300,400,150]);
ylim([-1,2])
set(gca,'Ytick',[-1,0,1,2]);
l1=legend('Simulation fault signal','Reconstructed signal')
set(l1,'Fontname', 'Times New Roman','FontSize',12)
set(l1,'Box','off');
%%  时间统计  
ct1=tic
for i=1:10
    [dic1,xcofe1,rec1]=CSC_MOD(odata,Lq,K,100); %LoCOMP
end
ct2=toc(ct1)
(ct2)/10
%%
Y=odata;
Ya=[];
npq=400;
for i=1:floor(N/npq)
    Ya(:,i)=Y(((i-1)*npq+1):i*npq+0);
end
np=20;
L=1;
m = size(Ya,1)+1-np;
ct1=tic;
for i=1:10
[Dconvsu, Xuconvsu] = uconvdlasu(Ya,K/size(Ya,2), L, np, m); %仿真信号10
end
ct2=toc(ct1)
(ct2)/10
%%

Y=odata;
Ya=[];
npq=1000;
for i=1:floor(N/npq)
    Ya(:,i)=Y(((i-1)*npq+1):i*npq+0);
end
np=20;
L=1;
m = size(Ya,1)+1-np;
ct1=tic;
for i=1:10
[Dconvsu, Xuconvsu] = uconvdlasu(Ya,K/size(Ya,2), L, np, m); %仿真信号10
end
ct2=toc(ct1)
(ct2)/10
%%
Y=odata;
Ya=[];
npq=10000;
for i=1:floor(N/npq)
    Ya(:,i)=Y(((i-1)*npq+1):i*npq+0);
end
np=20;
L=1;
m = size(Ya,1)+1-np;
ct1=tic;
for i=1:10
[Dconvsu, Xuconvsu] = uconvdlasu(Ya,K/size(Ya,2), L, np, m); %仿真信号10
end
ct2=toc(ct1)/10
Ruconvsu=Dconvsu*Xuconvsu;
Ry=Ruconvsu(:);
%%
npq=10000;
Y=odata;
Ya=[];
N=length(Y);
for i=1:floor(N/npq)
    Ya(:,i)=Y(((i-1)*npq+1):i*npq+0);
end

% Set up cbpdndl parameters

opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 250;
opt.rho = 10*lambda;
opt.sigma = size(Ya,3);
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
n=20;
L=1;
% Construct initial dictionary
D0 = zeros(n,1,L);
D0(1:n,1,:) = randn(n,1,L);
% Do dictionary learning
lambda = 0.5;
t1=tic;
for i=1:10
    [D, X, optinf] = cbpdndl(D0, Ya, lambda, opt);
end
toc(t1)/10

lambda = 0.2;
t1=tic;
for i=1:10
    [D, X, optinf] = cbpdndl(D0, Ya, lambda, opt);
end
toc(t1)/10
