%% clear all the previous data and record
clear;
clc;
syms z
% input file name
[a,Fs] = audioread('sample.wav'); 
% or you can generate a 
white=wgn(length(a),1,0);
filename = 'newwhitenoise.wav';
% audiowrite(filename,test2,Fs);
audiowrite(filename,white,Fs);
[x,Fs]=audioread('newwhitenoise.wav');

system=inputparam(300,900,1500,2500,0,0,0,20,0,Fs);
%      inputparam(FLC,FB1,FB2,FHC,GL,GB1,GB2,GB3,GH,FS)
%results=simplify(system); % designed filter based on pass band cutoff and target gain
[fzn, fzd] = numden(system);
% extract the coefficient of the transfer function
fzc = sym2poly(fzn);
fzb = sym2poly(fzd);
%fvtool(fzc,fzb,'Fs',Fs)
freqz(fzc,fzb)
%% time domain convolution with impulse response conv with original signal
% get the impulse response based on the design parameters
h=impz(fzc,fzb);
% calculate the convolution of the system wit
fil1=conv(x,h);
y=fil1(length(h):length(fil1)); %
NFFT=length(y);
FIL=fft(y,NFFT);
F = ((0:1/NFFT:1-1/NFFT)*Fs).';
magnitudeY = abs(FIL);

%% draw the filtered signal via filter function (Method 2)
% test2=filter(fzc,fzb,x);
% NFFT=length(test2);
% TEST=fft(test2,NFFT);
% F = ((0:1/NFFT:1-1/NFFT)*Fs).';
% TE = abs(TEST);   
% plot(F,TE)
% grid on
% hold on

%% plot original signal before filtering
NFFT=length(x)
X=fft(x,NFFT);
F = ((0:1/NFFT:1-1/NFFT)*Fs).';
magnitudeX = abs(X);        % Magnitude of the FFT
%phaseX = unwrap(angle(X));  % Phase of the FFT
%% plot the  in freqeuncy domain
figure
plot(F,magnitudeX)
xlim([0 10000])
grid on
figure
%xlim([0 20000])
plot(F,magnitudeY)
xlim([0 10000])
grid on
%% plot in log scale:
% figure
% set(gca, 'XScale', 'log')
% xlim([0 20000])
% hold on
% plot(F,20*log10(magnitudeY))

%% plot and check the ratio of the gain after the filter
% figure
% H=magnitudeY./magnitudeX;
% ax=1:1:length(H);
% plot(ax,H)
% xlim([0 2000])
%% write into a new file 
filename = 'filteredsignal1.wav';
% audiowrite(filename,test2,Fs);
audiowrite(filename,fil1,Fs);
%% support function
function mainfunction=inputparam(FLC,FB1,FB2,FHC,GL,GB1,GB2,GB3,GH,FS)
syms z
L1=solveLength(FLC,FS)
B1=solveLength(FB1,FS)
B11=solveLength(FB2,FS)
H1=solveLength(FHC,FS)
mainfunction=Equalizer(L1,B1,B11,H1,GL,GB1,GB2,GB3,GH);
end

function result1=Drrs(L1,L2)
syms z
  result1=(1-z^(-L1))/(1-z^(-1))*(1-z^(-L2))/(1-z^(-1));
end

function y=getodd(x)
y = 2*floor(x/2)+1;
end

function bandLength=solveLength(fc,fs)
bandLength=getodd(1/(fc*2/fs));
%bandLength=getodd(1/(fc/fs));
end

function result2=Equalizer(L1,B1,B11,H1,GL,GB1,GB2,GB3,GH)
syms z
L2=getodd(L1/sqrt(2))
H2=getodd(H1/sqrt(2))
B2=getodd(B1/sqrt(2))
B22=getodd(B11/sqrt(2))
kH=1/(H1*H2);
kL=1/(L1*L2);
k1=1/(B1*B2);
k2=1/(B11*B22);
% change the log scale gain input to linear scale 
GL=10^(GL/20);
GB1=10^(GB1/20);
GB2=10^(GB2/20);
GB3=10^(GB3/20);
GH=10^(GH/20);
result2=GH*z^(-[L1+L2+B1+B2+B11+B22]/2+1)*[z^(-[H1+H2]/2+1)-kH*Drrs(H1,H2)]...
        +GB3*[kH*z^(-[L1+L2+B1+B2+B11+B22]/2+1)*Drrs(H1,H2)-k2*z^(-[H1+H2+B1+B2+L1+L2]/2+1)*Drrs(B11,B22)]...
        +GB2*[k2*z^(-[H1+H2+B1+B2+L1+L2]/2+1)*Drrs(B11,B22)-k1*z^(-[H1+H2+L1+L2+B11+B22]/2+1)*Drrs(B1,B2)]...
        +GB1*[k1*z^(-[H1+H2+L1+L2+B11+B22]/2+1)*Drrs(B1,B2)-kL*z^(-[H1+H2+B1+B2+B11+B22]/2+1)*Drrs(L1,L2)]...
        +GL*kL*z^(-[H1+H2+B1+B2+B11+B22]/2+1)*Drrs(L1,L2);
end


