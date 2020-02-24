% This scripts performs the estimation  of the smoothed log-spectrograms from the audio recording
% in PrunedLanguageDataset.mat and save the result in Phonetic.mat 

%% Smoothing of the log-spectrograms

save('Phonetic.mat')
load PrunedLanguageDataset; %Read in the data;
myvars= whos;               %Get the variable names

for i=1:length(myvars) 	    %For all variables
    i
	%Set the output variable names
	v = genvarname(['Analysis_Object_' myvars(i).name]);
    %v = genvarname(['Spectrogram' myvars(i).name]);
	%Give the sequence the full LPC/PSD treatment
	%eval([v ' = Analysis_Obj_Creator(eval(myvars(i).name)); '])
   
    eval([v ' = Spectrogram_smooth(eval(myvars(i).name)); '])
   save('Phonetic.mat',v,'-append')
end


%% Time warping

% The smoothed log-spectrograms in Phonetic.mat are warped together and
% the warped log-spectrograms are saved in WarpedPSD.mat, while warping functions 
% and their inverse are saved in Warping_functions.mat

clear;
%Read File Names
load('Phonetic');
ListOfVars = whos;
n = length(ListOfVars);

%Create Warped PSD-gram structure based
%on Digit and Speaker Sex
WarpedPSDs_DigitGender =struct;
WarpingdPSDs_DigitGender =struct;
unWarpingdPSDs_DigitGender =struct;
%Populate with empty warped PSD
for  u=1:1:n
    u
	NameOfVariable =  ListOfVars(u).name(17:(end-2));
	eval(['WarpedPSDs_DigitGender.' NameOfVariable '= zeros(81,100);']);
    eval(['WarpingPSDs_DigitGender.' NameOfVariable '= zeros(1,100);']);
    eval(['unWarpingPSDs_DigitGender.' NameOfVariable '= zeros(1,100);']);
end

%Set up the warper 
%Uncomment to recompile the first time
%cd ../PACE-WARP_Customized/
%mex -Ieigen-eigen-43d9075b23ef/ rthik_E_2D.cpp 
 % 
 addpath('PACE-WARP_Customized/')

%For each specific Digit
for DigitBase =1:1:10 
        DigitBase
        
        Y_dbF= {};   DigitVariableNamesF = {};    
        Y_dbM= {};   DigitVariableNamesM = {};         
        
        %Construct the unwarped female utterences dataset 
        k=1;
        for u=1:1:n
            u
            if (ListOfVars(u).name(22:end-4) == num2str(DigitBase) )
                if (ListOfVars(u).name(end-2) == 'F')		
                    ddd = eval(ListOfVars(u).name); 
                    NameOfVariable =  ListOfVars(u).name(17:(end-2)); 
                    Y_dbF{k} = ddd.SmoothedPSD;
                    DigitVariableNamesF{k} =  NameOfVariable;
                    k= 1+k;            
                end
            end
        end
                      
        %Construct the unwarped female utterences dataset 
        k=1;
        for u=1:1:n
            u
            if (ListOfVars(u).name(22:end-4) == num2str(DigitBase) )
                if (ListOfVars(u).name(end-2) == 'M')		
                    ddd = eval(ListOfVars(u).name); 
                    NameOfVariable =  ListOfVars(u).name(17:(end-2)); 
                    Y_dbM{k} = ddd.SmoothedPSD;
                    DigitVariableNamesM{k} =  NameOfVariable;
                    k= 1+k;            
                end
            end
        end
  
        %Warp the Y's together.
        [X_dbF,aligned_dbF,h_dbF,hinv_dbF]=WFPCA_2D(Y_dbF,linspace(0,1,100)); 
        [X_dbM,aligned_dbM,h_dbM,hinv_dbM]=WFPCA_2D(Y_dbM,linspace(0,1,100)); 

        %Save the warped frequency slice on the appropriate PSD position 
        %for every name find the corresponding index
        for u=1:1:length(DigitVariableNamesF ) 
            u
            eval(['WarpedPSDs_DigitGender.' char( DigitVariableNamesF(u)) ' = aligned_dbF{' num2str(u) '};']);
            eval(['WarpingPSDs_DigitGender.' char( DigitVariableNamesF(u)) ' = h_dbF(' num2str(u) ',:);']);
            eval(['unWarpingPSDs_DigitGender.' char( DigitVariableNamesF(u)) ' = hinv_dbF(' num2str(u) ',:);']);
        end        
        for u=1:1:length(DigitVariableNamesM )
            u
            eval(['WarpedPSDs_DigitGender.' char( DigitVariableNamesM(u)) ' = aligned_dbM{' num2str(u) '};']);
            eval(['WarpingPSDs_DigitGender.' char( DigitVariableNamesM(u)) ' = h_dbM(' num2str(u) ',:);']);
            eval(['unWarpingPSDs_DigitGender.' char( DigitVariableNamesM(u)) ' = hinv_dbM(' num2str(u) ',:);']);
        end
        
end 

save( 'WarpedPSD.mat', 'WarpedPSDs_DigitGender');
save('Warping_functions.mat','WarpingPSDs_DigitGender','unWarpingPSDs_DigitGender');

