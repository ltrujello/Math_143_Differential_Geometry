The statistical analysis of acoustic phonetic data: exploring differences between spoken Romance languages
D. Pigoli, P. Z. Hadjipantelis, J. S. Coleman and J. A. D. Aston
Appl. Statist., 67 (2018), 1103 -- 1145


Files descriptions:

licenses.txt: copyright notes for the code attached to the paper. 

Preprocessing:

GetSpectrogram.m : Matlab function which computes the spectrograms.

Preprocess_script.m : Matlab script to run the preprocessing on the linguistic dataset.

PrunedLanguageDataset.mat : Matlab object which containd the raw audio recordings.

Spectrogram_smooth.m: Matlab function which computes a smooth surface from the raw log-spectrogram.

SurfaceCubicInterpolator.m: Matlab function for spline interpolation.

dctn.m: Matalb function for the discrete cosine transform.

idctn.m: Matalb function for the inverse discrete cosine transform.

smoothn.m: Matlab function for robust spline smoothing of surfaces 

PACE-WARP_Customized: Matlab and C++ code for the PACE warping of surfaces.


Analysis of the smooth and aligned log-spectrograms:

WarpedPSD.RData: R object containing the smooth and aligned log-spectrograms

SpeakerInfo.txt: language, word and gender of the speakers for the log-spectrograms in WarpedPSD.RData.

acoustic_functions.R: R functions to estimate and test the components of the sound process.

acoustic_data_script.R: R script to perform the analysis described in the paper.


Sound reconstruction:

sound_reconstruction.m: Matlab script to reconstruct the sound.
myistfit.m: Matlab function for the inverse short-time Fourier transform.
Spect_path.mat: Matlab object containing an example of the changing path of the log-spectrograms between two speakers of different languages.
Spect_proj.mat: Matlab object containing an example of the changing path of the log-spectrograms between a speaker and its projection in a different language.

F2P_spk_digit5.wav: reconstructed sound for the path from a French speaker and a Portuguese speaker for the digit 5.
F2P_spk_digit7.wav: reconstructed sound for the path from a French speaker and a Portuguese speaker for the digit 7.
F2P_path_digit5.wav: reconstructed sound for the path from a French speaker and its projection into Portuguese for the digit 5.
F2P_path_digit7.wav: reconstructed sound for the path from a French speaker and its projection into Portuguese for the digit 7.


Analysis in the Supplementary Material:

Analysis_SuppMatt.R : R script to run the analysis with the alternative preprocessing methods.
SVRF_WarpedPSD_SuppMat.RData: smoothed and aligned log-spectrogram.
F2P_spk_digit5_SuppMat.wav : reconstructed sound for the path from a French speaker and a Portuguese speaker for the digit 5.
F2P_proj_digit5_SuppMat.wav : reconstructed sound for the path from a French speaker and its projection into Portuguese for the digit 5.



John A. D. Aston
Department of Pure Mathematics and Mathematical Statistics
Statistical Laboratory
University of Cambridge
Wilberforce Road
Cambridge
CB3 0WB
UK
E-mail: j.aston@statslab.cam.ac.uk




