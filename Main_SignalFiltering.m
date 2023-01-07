% filter the response signals with a highpass filter to remove the rigid body motion of the rocket
% add uncorrelated Gaussian sequences to the signals to simulate the measurement noise
clc;
clear;
close all

MCn = 200;
SNR = 40; 

datafile = '.\Response\Dis_MC1';
load (datafile);
[Nt,No] = size(dis); % Nt signal length, No signal channel
dt = t(2)-t(1); % sampling interval
fs = 1/dt;      % sampling rate

% time duration of the signal
ts = 72;
if ts>t(end)
    ts = t(end);
end

% highpass filter design
half_fs = fs/2;
A_pass = 3; % dB of ripple in the passband
A_stop = 40; % dB of attenuation in the stopband
freq_pass = 2/half_fs; % the normalized passband edge frequencies
freq_stop = 1/half_fs; % the normalized stopband edge frequencies
[filter_order,freq_n]=buttord(freq_pass,freq_stop,A_pass,A_stop);
[filter_num,filter_den]=butter(filter_order,freq_n,'high');
freqz(filter_num,filter_den,1000,fs);

% fitering and adding measurement noise
for MCii = 1:MCn
    MCii
    filename = strcat('.\Response\Dis_MC',num2str(MCii),'.mat');
    load(filename);
    
    dis = dis(1:(ts*fs+1),:);
    
    for jj = 1:No
        temp = filter(filter_num,filter_den,dis(:,jj)); % 滤波
        dis(:,jj) = awgn(temp,SNR,'measured');
    end
    
    % short-time Fourier transform
    NFFT = 512;
    if MCii == 1
        for ii = 1:No
            figure
            [~,F,T,P] = spectrogram(dis(:,ii),NFFT,NFFT/2,NFFT,fs,'yaxis');
            imagesc(T, F, 10*log10(P)) % add eps like pspectrogram does
            axis xy
            ylabel('Frequency/Hz')
            xlabel('Time/s')
            h = colorbar;
            h.Label.String = 'Power per frequency/dB·Hz^-^1';
            set(gca,'Fontsize',14,'Linewidth',1)
        end
    end
    
    filedir = strcat('.\Response_SNR',num2str(SNR));
    if ~exist(filedir,'dir')
        mkdir(filedir);
    end
    filename = strcat(filedir,'\Dis_MC',num2str(MCii),'.mat');
    save(filename,'dis','t')
end