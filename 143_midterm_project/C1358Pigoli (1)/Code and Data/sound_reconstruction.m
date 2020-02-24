%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script reconstructs the path of sound change from the changes     %
% in the log-spectrograms.                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%% Reconstruction of the path from a word pronounced by a speaker and its projection in a different language   %%%%%%%%%%%%%

%% load the log-spectrogram path

load('Spect_proj.mat');  % Projection of the digit five from the first speaker of French language into Portuguese

%% load original data and extract phase information

load PrunedLanguageDataset
window_size=160;      %window length
h=0.5*window_size;    % hop length
fs=1/16000; % sample frequency
sampling_frequency_InHz = 16000;
window_length_InMS = 10;
window_size = sampling_frequency_InHz  / (1000 / window_length_InMS  );

% Set Gaussian window of size equal to the specified duration of the window
% length in ms.
w = window(@gausswin,window_size);

[S,f,T,PSD] = spectrogram(PO_A_7_M_c, w, [],window_size, 16000);
dim=size(S);
ph=angle(S);
out=0;

%% Sound reconstruction
for u=1:1:11
    eval(['F2=SurfaceCubicInterpolator(Spect_path' num2str(u) ',dim(2));'])

    F=sqrt(exp(F2/10));

    B=F.*(cos(ph)+1i*sin(ph));

    [output,t]=myistft(B,h,window_size,fs);
    output=2*((output-min(output))/(max(output)-min(output)))-1;
    out=[out output];
end



audiowrite('SI2P_path_digit7.wav',out(2:length(out)),16000);


%%
%%%%%%%%%%%%%%% Reconstruction of the path between two different speakers %%%%%%%%%%%%%

%% load the log-spectrogram path
load('Spect_path.mat'); % Path between a French  speaker and a Portuguese speaker for the word five.

%% load original data and extract phase information

load PrunedLanguageDataset


window_size=160;      %window length
h=0.5*window_size;    % hop length
fs=1/16000; % sample frequency
sampling_frequency_InHz = 16000;
window_length_InMS = 10;
window_size = sampling_frequency_InHz  / (1000 / window_length_InMS  );

% Set Gaussian window of size equal to the specified duration of the window
% length in ms.
w = window(@gausswin,window_size);

[S,f,T,PSD] = spectrogram(PO_B_5_F_c, w, [],window_size, 16000);
dim=size(S);
ph=angle(S);
out=0;

%% Sound reconstruction


for u=1:1:11
    eval(['F2=SurfaceCubicInterpolator(Spect_path' num2str(u) ',dim(2));'])
    F=sqrt(exp(F2/10));
    B=F.*(cos(ph)+1i*sin(ph));

    [output,t]=myistft(B,h,window_size,fs);
    output=2*((output-min(output))/(max(output)-min(output)))-1;
    out=[out output];
end



audiowrite('F2P_spk_digit5.wav',out(2:length(out)),16000);





