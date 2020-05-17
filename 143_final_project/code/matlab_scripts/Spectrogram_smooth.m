function [PSD] = Spectro_smooth(signal_sequence)

sampling_frequency_InHz = 16000;
window_length_InMS = 10;
p_order=18;
smoothing_par = .5;

[F,Tspc,PSD] = GetSpectrogram(signal_sequence, sampling_frequency_InHz, window_length_InMS );

SmoothedPSD  = smoothn(10*log(PSD), smoothing_par);

PSD_collection.F= F;
PSD_collection.Tspc= Tspc;
PSD_collection.SmoothedPSD= SurfaceCubicInterpolator(SmoothedPSD,100);

PSD= PSD_collection;

