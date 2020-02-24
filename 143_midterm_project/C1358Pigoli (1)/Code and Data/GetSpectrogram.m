function [F,T,PSD] = GetSpectrogram(signal_sequence, sampling_frequency_InHz, window_length_InMS )

window_size = sampling_frequency_InHz  / (1000 / window_length_InMS  );

    % Set Gaussian window of size equal to the specified duration of the window
    % length in ms.
    w = window(@gausswin,window_size);
    %Calculate spectrogram of the signal sequence using the prespecified
    %window, assume that the length of the fft transformation will equal
    %the window size and using the known sampling frequency of the signal
    
   [S,F,T,PSD] = spectrogram(signal_sequence , w, [],window_size, sampling_frequency_InHz );
   
   %F being the signal frequencies, (hz)
   %S the STFT of the signal
   %T the time which we sampled at (seconds)
   %PSD the power spectral density (PSD) of each segment. 
   %For real x, PSD contains the one-sided modified periodogram estimate of the PSD of each segment

   
   
