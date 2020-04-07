import os
from numpy import genfromtxt
myPath = '/Users/luketrujillo/Desktop/Senior_year/Spring/Seminar in Diff Geo/143_midterm_project/' #Path containing this notebook; change this

def word_to_surface(word, myPath, accent='us'): #given a word, download the .mp3 from oxford and return the signal and frequency.
    myPath = '/Users/luketrujillo/Desktop/Senior_year/Spring/Seminar in Diff Geo/143_midterm_project/'
    
    ## Setting paths for later; keeps things organized.
    wavPath = myPath + "wavs/"                      #full path for the the .wav files
    psdPath = myPath + "psd_data/"                  #full path for the psd data 
    matlabScriptsPath = myPath + "matlab_scripts/" #full path to get to the matlab scripts
    
    ## Now we get the words from the web and compute the psd
    url = 'https://ssl.gstatic.com/dictionary/static/sounds/oxford/{}--_{}_1.mp3'.format(word, accent) #warning: the oxford website could change their website formats for .mp3s in the future; may be a bug source in the future
    wavFileName = wavPath + '{}-{}'.format(word, accent)
    
    print(wordFileName)
    if not os.path.isfile(wordFileName + '.wav'):            #if we haven't made the .wav before, then make it
        wget.download(url, wordFileName + '.mp3')            #download the .mp3
        sound = AudioSegment.from_mp3(wordFileName + '.mp3') 
        sound.export(wordFileName + '.wav', format="wav")    #the .wav is now created. 
    
    psdCSV = psdPath + word + "_smooth_psd.csv" #our psd data
    if not os.path.isfile(psdCSV):              #if we haven't made the .csv before, then make it
        wordWavPath = wordFileName + ".wav"     
        psd_smoother(word, wordWavPath, matlabScriptsPath) #call helper function to compute psd
    
    word_psd = genfromtxt(psdCSV, delimiter=',')

    ## finally, we plot the surface
    xdata = [np.linspace(0, 1, 81), np.linspace(0, 1, 100)] #our x-y coordinates
    i,j = np.meshgrid(*xdata, indexing = "ij")              #creates the 81x100 grid
    zdata_s = csaps(xdata, word_psd, xdata, smooth=0.998)   #using csaps to smooth our surface
    
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('none')
    c = [s['color'] for s in plt.rcParams['axes.prop_cycle']]
    ax.scatter(j, i, zdata_s, linewidths=0.5, color=c[0], alpha=0.5)
    ax.plot_surface(j, i, zdata_s, color=c[1], linewidth=0, alpha=1)
    ax.view_init(elev=15., azim=290)