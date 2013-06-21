/**
*	@file domain.cpp
*	@class Domain 
*	
*	@author Chad R. Befus, Chris Sanden, Cody Rioux
*	@version 0.1.3:Bravo
*
*	License
*	-----------
*	Copyright (c) 2010 Chad R. Befus, Chris Sanden
*	
*	CAMEL (Content-based Audio and Music Extraction Library) is unrestricted, open 
*	source "software", with respect, but not limited, to modification, use, 
*	publishing, and distribution, subject to the following conditions:
*	
*	1. This copyright notice is maintained in its current form throughout CAMEL.
*	
*	2. Any publication or distribution credits the usage of CAMEL appropriately.
*	
*	This software is provided "as is" without warranty. The authors are in no way
*	liable for any misuse or damages arising from the use, modification or 
*	distribution of it or any part thereof.
*
*	Description
*	-----------
*	@brief Extract domain values from a specified location in a given PCM text file
*
*	Calculate the domain (ex: Time, Frequency, Peak, ect.) values for a given window
*	within a PCM text file (as per instructions for input file found in the manual.
*	This includes differing between the various fft spectrum types and applying any
*	windowing functions.
*
*	Required Dependencies
*	-----------
*	std::vector
*	std::string
*	fileVector.h
*	configFile.h
*
*	Example Usage
*	-----------
*	#include "domain.h"
*	#include <vector>
*	using namespace std;
*
*	int main(){	 
*		Domain domain;
*		domain.setFileName("myPcmFile.txt");
*		domain.setWindowSize(1024);
*		domain_.setStartPosition(0);
*		domain_.setEndPosition(10000);
*		domain.setup();	
*		vector<double> domainVec = domain.getFrequencyDomainWindow(0);
*		for(int i = 0; i < 1024; ++i){
*			cout << domainVec[i] << endl;
*		}
*	}
*
*
*	More examples to be found in example folder
*
*	@todo
*	consider why we preset start/end position but ask for window size and startingAts
*	Switch namespace std to using localised scope
*
**/

#include "domain.h"

#include <string>
#include <vector>
#include <math.h>
///for fftw make sure you have the fftw3-dev package installed
#include <fftw3.h>

#include "fileVector.h"
#include "configFile.h"

using namespace std;


/** 
*	Domain
*	@brief The Constructor: Sets private members fileName_ and windowSize_
*	@param fileName A string: the name of a PCM file
*	@param windowSize An integer: the number of PCM values to include in a window
*	@post fileName_ and windowSize_ are set
*	@pre the file associated with fileName exists and windowSize is non-negative
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
Domain::Domain(string fileName, int windowSize)
	:fileName_(fileName), 
	windowSize_(windowSize){}


/** 
*	setup
*	@brief sets up the dependensies of the domain class -> pcmVec
*	@param
*	@post private member pcmVec_ is initialised
*	@pre file associated with fileName_ exists and is proper format
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void Domain::setup(){
	pcmVec_.setFileName(fileName_);
	pcmVec_.setup();
}


/** 
*	setFileName
*	@brief accessor to set the fileName private member
*	@param fileName the new name of the file object
*	@post private member fileName_ is changed
*	@pre fileName parameter names a existing file
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void Domain::setFileName(string fileName){
	fileName_ = fileName;
}


/** 
*	setWindowSize
*	@brief accessor to set private member windowSize_
*	@param windowSize the new positive integer windowSize in pcm value count
*	@post windowSize_ has changed
*	@pre windowSize is a positive integer
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void Domain::setWindowSize(int windowSize){
	windowSize_ = windowSize;
}


/** 
*	setStartPosition
*	@brief accessor to set the private member begin_
*	@param startPosition the positive integer representing where to begin extraction
*	@post begin_ has changed
*	@pre startPosition is set with a positive integer
*	@return
*	@par Method If startPosition is set to a value greater then the end_ position then
*		startPosition is set to 0.
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void Domain::setStartPosition(int startPosition){
	if(startPosition > end_) startPosition = 0;
	begin_ = startPosition;
}


/** 
*	setEndPosition
*	@brief accessor to set end_ private member
*	@param endPosition a positive integer representing where to stop extraction
*	@post end_ has changed
*	@pre endPosition is set with a positive integer
*	@return
*	@par Method If endPosition is greater than the length of the file we set end_ to the
*		length of the file.  Otherwise set it to endPosition parameter value.
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void Domain::setEndPosition(int endPosition){
	if(pcmVec_.size() < endPosition) endPosition = pcmVec_.size();
	end_ = endPosition;
}


/** 
*	getStartPosition
*	@brief accessor for retrieving the begin_ private member
*	@param
*	@post
*	@pre begin_ is set
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
int Domain::getStartPosition(){
	return begin_;
}


/** 
*	getEndPosition
*	@brief accessor for retrieving end_ private member
*	@param
*	@post
*	@pre end_ is set
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
int Domain::getEndPosition(){
	return end_;
}


/** 
*	getTimeDomainWindow
*	@brief returns a window (of windowSize_) of time domain values starting at startPosition
*	@param startPosition An int: defines the starting position of the window to extract
*	@post
*	@pre windowSize_ is set, startPosition is a positive integer
*	@return a vector if int: windowSize_ time domain values
*	@par Method Starting at startPosition for the number of set windowSize_ values or until the
*		end of file, grav the value from pcmVec_ and push it onto our return Vector
*	@note The first domain offered in CAMEL is the time domain.  This is the natural domain of 
*		the PCM audio file and therefore no //transformations are applied.  Time domain 
*		signals are values of amplitude stored at some sampling rate (44100Hz, 22050Hz, 
*		ect.) which is constant over time.
*	@remarks
*	@warning
*	@todo
**/
vector<int> Domain::getTimeDomainWindow(int startPosition){
	vector<int> pcmWindow;
	for(int i = 0; i < windowSize_ && i < end_ - startPosition; i++){
		pcmWindow.push_back(pcmVec_[startPosition+i]);
	}

	return pcmWindow;
}

//\subsubsection{Frequency Domain}
//+


/** 
*	getFrequencyDomainWindow
*	@brief returns a vector of the frequency domain associated with a window of pcm starting at parameter value
*	@param startPosition An int: the begining position for our window of extraction
*	@post
*	@pre windowSize_ is set, startPosition is positive, end_ is set
*	@return A vector of double: vector of frequency values over a window of PCM
*	@par Method Set up a vector of 0's as pcm values so that if we do not have enough PCM values
		(AKA we are at the end of the file) we are pre padded with 0's.  Overwrite the 0's with
		actual pcm values when possible. Send the vector of PCM values to the FrequencyDomain
		calculating function.  Apply the windowing function.
*	@note	Spectral features are calculated over the frequency domain.  The frequency domain is 
*		a Fourier transformation from the time domain. In CAMEL the Fourier transform is 
*		provided for us through the use of the FFTW Library \cite{frigo2005}.
*	@remarks
*	@warning
*	@todo Should be called Spectrum Domain... Frequency Domain is before spectrum function
**/
vector<double> Domain::getFrequencyDomainWindow(int startPosition){
	vector<int> pcmWindow;
	for(int i = 0; i < windowSize_; i++){
		pcmWindow.push_back(0);
	}
	for(int i = 0; i < windowSize_ && i < end_ - startPosition; i++){
		pcmWindow[i] = pcmVec_[startPosition+i];
	}
	vector<double> windowVec = calcFrequencyDomainWindow(pcmWindow);
	windowVec = applyWindowFunc(windowVec);

	return windowVec;
}


/** 
*	getPeakDomainWindow
*	@brief gets a window of peak domain values representing a window of PCM values
*	@param startPosition an int: represents the start of the PCM window for extraction
*	@post
*	@pre windowSize_ is set, startPosition is positive.
*	@return A vector of double: vector of peak domain values for the specified window
*	@par Method Calculate the Frequency Domain over the pcm window, calculate peak domain over that.
*	@note The peak domain is a calculation of amplitude and spectral peaks over the spectral domain.
*	@remarks
*	@warning
*	@todo
**/

vector<double> Domain::getPeakDomainWindow(int startPosition){

	vector<double> freqVec = getFrequencyDomainWindow(startPosition);
	vector<double> windowVec = calcPeakDomainWindow(freqVec);

	return windowVec;
}

/** 
*	calcFrequencyDomainWindow
*	@brief Calculates the Frequency Domain over a vector of pcm values and applies spectrum
*	@param pcmWindow A vector of int: pcm values
*	@post
*	@pre windowSize_ is set, pcmWindow is of size windowSize_
*	@return A vector of double: representing the spectral domain of the pcm window
*	@par Method Copy the vector over to an array (since FFTW uses arrays).  Call calcFFT to get
*		The frequency values. Call applySpectrumFunc to get the spectrum over FFT values.
*	@note	
*	@remarks
*	@warning
*	@todo should be called calcSpectrumDomainWindow... not Freq (which is a sub part)
**/
vector<double> Domain::calcFrequencyDomainWindow(vector<int> pcmWindow){

	float pcmVals[windowSize_];
	for(int i = 0; i < windowSize_; i++){
		pcmVals[i] = static_cast<double>(pcmWindow[i]);
	}
	vector<double> fftVals = calcFFT(pcmVals);
	vector<double> specVec = applySpectrumFunc(fftVals);
	
	return specVec;
}


/** 
*	calcPeakDomainWindow
*	@brief Calculate the Peak Domain over a vector of spectrum values
*	@param windowVec a vector of double: spectrum values
*	@post
*	@pre windowSize_ is set, windowVec is of length windowSize_
*	@return A vector of double: the peak domain values for the given windowVec
*	@par Method Given a window of spectral values $w$, it is formulated via:
*
*		\begin{equation}
*			\forall_{i=1}^{(n/2)-2} w_{i} = 
*				\frac{1}{4} * (w_{i-1} - w_{i+1}) * 
*				(\frac{\frac{1}{2}*(w_{i-1}-w_{i+1})}{w_{i-1}-2(w_{i}+w_{i+1})})
*		\end{equation}	
*
*		and	
*
*		\begin{equation}
*			\forall_{i=n/2}^{n} w_{i} = 
*				\frac{SR}{n}(i+(\frac{\frac{1}{2}(w_{i-1}-w_{i+1})}{w_{i-1}-2(w_{i}+w_{i+1})}))
*		\end{equation}	
*
*		where $w_{i} > th$ and $w_{i} > w_{i-1}$ and $w_{i} > w_{i+1}$ and
*
*		\begin{equation}
*			\forall_{i=0}^{n} w_{i} = 0 
*		\end{equation}	
*
*		where $w_{i} < th$ or $w_{i} < w_{i-1}$ or $w_{i} < w_{i+1}$.  
*		Note that $th$ is a desired threshold for deciding peak and $SR$ is the sampling rate.
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> Domain::calcPeakDomainWindow(vector<double> windowVec){

	int N = static_cast<int>(windowSize_/2.0);
	vector<double> result(windowSize_);
    	double q = settingsFile_.read<double>("PCM_SAMPLE_RATE", 44100.0) / (windowSize_/2.0);
    	double threshold = settingsFile_.read<double>("PEAK_THRESHOLD_PERCENT", 10.0);
	double max = 0.0;

	for(int i = 0; i < N; i++){
		if(windowVec[i] > max){
			max = windowVec[i];
		}
	}

    	threshold *= .01 * max;
    	result[0] = 0;
    	result[N] = 0;
    	for(int i = 1; i < N; i++){
        	if(windowVec[i] >= threshold){
            		if(windowVec[i] > windowVec[i - 1] && i + 1 < N && windowVec[i] > windowVec[i + 1]){
				double y = windowVec[i-1];
				double y2 = windowVec[i];
				double y3 = windowVec[i+1];
				double p = .5 * (y - y3) / (windowVec[i - 1] - 2 * y2 + windowVec[i + 1]);

                		result[N + i] = q * (i + p);
                		result[i] = y2 - .25 * (y - y3) * p;
            		}else{
                		result[i] = 0;
                		result[N + i] = 0;
            		}
        	}else{
            		result[i] = 0;
            		result[N + i] = 0;
        	}
    	}	  
	return result;
}

/** 
*	applySpectrumFunc
*	@brief
*	@param name descr
*	@post
*	@pre
*	@return
*	@par Method 
*	@note	Once the frequency domain has been calculated via a Fourier 
*		Transform it can be converted into one of several spectums 
*		via a spectral function.  We implement four such functions.
*		Note that for our descriptions below in commenting 
*		we designate $w$ to be the window on which we are working 
*		and $n$ to be its length.
*	@remarks
*	@warning
*	@todo
**/

vector<double> Domain::applySpectrumFunc(vector<double> fftVec){
	
	vector<double> specVec;
	int specType = settingsFile_.read<int>("SPECTRUM_TYPE",4);
	switch(specType){
		case 1:{
			specVec = applyLogMagnitude(fftVec);
			break;
		}
		case 2:{
			specVec = applyPower(fftVec);
			break;
		}
		case 3:{
			specVec = applyLogPower(fftVec);
			break;
		}
		default: {
			specVec = applyMagnitude(fftVec);
			break;
		}
	}

	//we could normalize here!

	return specVec;
}



/** 
*	applyLogMagnitude
*	@brief Applies the log magnitude funtion to a vector of fft values
*	@param fftVec A vector of doubles representing windowSize_ fft values
*	@post
*	@pre windowSize_ is set
*	@return a vector of doubles applying the log magnitude funtion to the fft values
*	@par Method For the log magnitude spectrum the function is very similar to that 
*	of the magnitude spectrum however we take the log of the values in first half of 
*	the window.  Given a window of frequency domain values $w$ we calculate the log 
*	magnitude spectrum via:
*
*	\begin{equation}
*		\forall_{i=0}^{(n/2)-1} w_{i} = 
*			\frac{log(\frac{1}{n}\sqrt{w_{ip}^2 + w_{n-ip}^2}) + off}{off}
*	\end{equation}
*	and
*	\begin{equation}
*		\forall_{i=n/2}^{n} w_{i} = ip * (SR/n)
*	\end{equation}
*
*	where $ip = i$ if the DC component is requested and $p=i+1$ otherwise. Also, $SR$
*	is the sampling rate of the audio and $off$ is the offset of the decible scale, 
*	defaulted to 96.0.
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> Domain::applyLogMagnitude(vector<double> fftVec){

	bool blnWithDC = settingsFile_.read<bool>("INCLUDE_DC_COMPONENT",false);
	vector<double> lmagVec(windowSize_);
	int intHalfN = windowSize_ / 2;
	double weightPerSample = settingsFile_.read<double>("PCM_SAMPLE_RATE", 44100.0) / static_cast<double>(windowSize_);
	
	int iPlus = 0;
	for(int i = 0; i < intHalfN; ++i){ 
     		if(!blnWithDC && iPlus == 0){
                    ++iPlus;
     		}

		double temp = pow(fftVec[iPlus],2) + pow(fftVec[windowSize_ - iPlus],2);
		if (temp > settingsFile_.read<double>("SMALL_NUMBER_LIMIT", 2e-42)){
     			temp = logf(sqrtf(temp) / static_cast<double>(windowSize_));
		} else {
			temp = settingsFile_.read<double>("LOG_LIMIT_DB", -96.0);
		}
       		
   		lmagVec[i] = (temp + settingsFile_.read<double>("DB_SCALE_OFFSET", 96.0)) / settingsFile_.read<double>("DB_SCALE_OFFSET", 96.0);
 		lmagVec[intHalfN + i] = iPlus * weightPerSample; 
		++iPlus;
	}

	return lmagVec;
}



/** 
*	applyPower
*	@brief Applies the power function to a vector of FFT values
*	@param fftVec A vector of doubles representing windowSize_ fft values
*	@post
*	@pre windowSize_ is set
*	@return a vector of doubles applying the power funtion to the fft values
*	@par Method The power spectrum is quite similar to the magnitude spectrum 
*	however instead of taking the square root of the numerator we square the 
*	denominator.  The formula for the power spectrum, given a window of frequency 
*	domain values $w$ is:
*
*	\begin{equation}
*		\forall_{i=n/2}^{n} w_{i} = ip * (SR/n)
*	\end{equation}
*	and
*	\begin{equation}
*		\forall_{i=0}^{(n/2)-1} w_{i} = 
*			\frac{1}{n^2}(w_{ip}^2 + w_{n-ip}^2)
*	\end{equation}
*
*	Note that $ip = i$ if the DC component is wanted and ip=i+1 otherwise. where 
*	$SR$ is the sampling rate of the audio. 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> Domain::applyPower(vector<double> fftVec){

	bool blnWithDC = settingsFile_.read<bool>("INCLUDE_DC_COMPONENT",false);
	vector<double> powerVec(windowSize_);
	int intHalfN = windowSize_ / 2;
	double weightPerSample = settingsFile_.read<double>("PCM_SAMPLE_RATE", 44100.0) / static_cast<double>(windowSize_);
	
	int iPlus = 0;
	for(int i = 0; i < intHalfN; ++i){
		if(!blnWithDC && iPlus == 0){
			++iPlus;
		}
		powerVec[i] = (pow(fftVec[iPlus],2) + pow(fftVec[windowSize_ - iPlus], 2)) / pow(windowSize_,2);
		powerVec[intHalfN + i] = iPlus * weightPerSample;
		++iPlus;
	}

	return powerVec;
}

/** 
*	applyLogPower
*	@brief Applies the log power function to a window of fft values
*	@param fftVec A vector of doubles representing windowSize_ fft values
*	@post
*	@pre windowSize_ is set
*	@return a vector of doubles applying the log power funtion to the fft values
*	@par Method Exactly the same way that the Log Magnitude spectrum added 
*	in the log function to the Magnitude spectrum, the Log Power spectrum 
*	adds a log function into the Power Spectrum.  The function for the Log 
*	Power Spectrum, given a window of frequency domain values $w$ is: 
*
*	\begin{equation}
*		\forall_{i=0}^{(n/2)-1} w_{i} = 
*			\frac{log(\frac{1}{n^2}(w_{ip}^2 + w_{n-ip}^2))+off}{off}	
*	\end{equation}
*	and
*	\begin{equation}
*		\forall_{i=n/2}^{n} w_{i} = ip * (SR/n)	
*	\end{equation}
*
*	Note that $ip = i$ if the DC component is wanted and ip=i+1 otherwise.  
*	Also note that $SR$ is the sampling rate of the audio and $off$ is the offset 
*	of the decible scale, defaulted to 96.0.
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> Domain::applyLogPower(vector<double> fftVec){
	
	bool blnWithDC = settingsFile_.read<bool>("INCLUDE_DC_COMPONENT",false);
	vector<double> lpowerVec(windowSize_);
	int intHalfN = windowSize_ / 2;
	double weightPerSample = settingsFile_.read<double>("PCM_SAMPLE_RATE", 44100.0) / static_cast<double>(windowSize_);

	int iPlus = 0;
	for(int i = 0; i < intHalfN; ++i){
		if(!blnWithDC && iPlus == 0){
			++iPlus;
		}
		double temp = pow(fftVec[iPlus],2) + pow(fftVec[windowSize_ - iPlus],2);
		if (temp > settingsFile_.read<double>("SMALL_NUMBER_LIMIT", 2e-42)){
			temp = logf(temp / pow(windowSize_,2));
		}else {
			temp = settingsFile_.read<double>("LOG_LIMIT_DB", -96.0);
		}
		
		lpowerVec[i] = (temp + settingsFile_.read<double>("DB_SCALE_OFFSET", 96.0)) / settingsFile_.read<double>("DB_SCALE_OFFSET", 96.0); 
		lpowerVec[intHalfN + i] = iPlus * weightPerSample;
		++iPlus;
	}

	return lpowerVec;
        
}


/** 
*	applyMagnitude
*	@brief Applies the magnitude function to a window of FFT values
*	@param fftVec A vector of doubles representing windowSize_ fft values
*	@post
*	@pre windowSize_ is set
*	@return a vector of doubles applying the magnitude funtion to the fft values
*	@par Method The most popular and simplest spectrum is the magnitude 
*	spectrum.  It is calculated given a window of the frequency domain $w$ 
*	via:
*
*	\begin{equation}
*		\forall_{i=0}^{(n/2)-1} w_{i} = 
*			\frac{1}{n}\sqrt{w_{ip}^2 + w_{n-ip}^2}
*	\end{equation}
*	and
*	\begin{equation}
*		\forall_{i=n/2}^{n} w_{i} = ip * (SR/n)
*	\end{equation}
*
*	where $ip = i$ if the DC component is requested and $p=i+1$ otherwise. 
*	Also, $SR$ is the sampling rate of the audio.  
*	@note	
*	@remarks
*	@warning
*	@todo
**/

vector<double> Domain::applyMagnitude(vector<double> fftVec){

	bool blnWithDC = settingsFile_.read<bool>("INCLUDE_DC_COMPONENT",false);
	vector<double> magVec(windowSize_);
	int intHalfN = windowSize_ / 2;
	double weightPerSample = settingsFile_.read<double>("PCM_SAMPLE_RATE", 44100.0) / static_cast<double>(windowSize_);	

	int iplus = 0;
	for(int i = 0; i < intHalfN; ++i){
		if(!blnWithDC && iplus == 0){
			++iplus;
		}
		magVec[i] = sqrt(pow(fftVec[iplus],2) + pow(fftVec[windowSize_ - iplus],2)) / windowSize_;
		magVec[intHalfN + i] = iplus * weightPerSample;
		++iplus;
	}

	return magVec;
}


/** 
*	calcFFT
*	@brief calculates the FFT values of a pcm array
*	@param pcmVals an array of Floats of size windowSize
*	@post
*	@pre windowSize_ is set
*	@return a vector of double representing the FFT values over pcmVals parameter
*	@par Method Simply use the FFTW library to calculate FFT over an array of floats.
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> Domain::calcFFT(float pcmVals[]){
	fftwf_plan ffpPlan;

	size_t allocSize = windowSize_ * sizeof(float);

	float *fltReturnArr = (float *)fftwf_malloc(allocSize);
	
	//FFT Optimization Options: FFTW_ESTIMATE,FFTW_PATIENT, FFTW_MEASURE
	ffpPlan = fftwf_plan_r2r_1d(windowSize_, pcmVals, fltReturnArr, FFTW_R2HC, FFTW_ESTIMATE); 
    
	fftwf_execute_r2r(ffpPlan, pcmVals, fltReturnArr);

	vector<double> fftVec(windowSize_);
	for(int i = 0; i < windowSize_; i++){ //because fftw requires arrays not vectors... 
		fftVec[i] = static_cast<double>(fltReturnArr[i]);
	}

	fftwf_destroy_plan(ffpPlan);
    	fftwf_free(fltReturnArr);
	return fftVec;
}

/* +\subsection{Windowing Functions}
+
+ */

/** 
*	applyWindowFunc
*	@brief  Switch through to the correct window function
*	@param windowVec a vector of double holding windowSize_ spectral values
*	@post
*	@pre settingsFile_ is set
*	@return
*	@par Method Once a domain window has been calculated we apply a windowing 
*	function.  Windowing functions wieght the importance of values at different 
*	locations within the window by different amounts.  We implement 8 such functions 
*	and describe them in this section. Note that for our descriptions here we designate 
*	$w$ to be the window on which we are working and $n$ to be its length.
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> Domain::applyWindowFunc(vector<double> windowVec){

	int windowType = settingsFile_.read<int>("SPECTRUM_WINDOW", 8);
	switch (windowType) {
		case 1:{
			windowVec = applyHann(windowVec);
			break;
		}
		case 2:{
			windowVec = applyHamming(windowVec);
			break;
		}
		case 3:{
			windowVec = applyTriangular(windowVec);
			break;
		}
		case 4:{
			windowVec = applyBartlett(windowVec);
			break;
		}
		case 5:{
			windowVec = applyBartlett_hann(windowVec);
			break;
		}
		case 6:{
			windowVec = applyBlackman(windowVec);
			break;
		}
		case 7:{
			windowVec = applyBlackman_harris(windowVec);
			break;
		}
		default:{
			windowVec = applyRectangular(windowVec);
			break;
		}
	}
	return windowVec;
}



/** 
*	applyRectangular
*	@brief  Just returns the windowVec as is
*	@param windowVec a vector of double representing windowSize_ spectrum values
*	@post
*	@pre
*	@return the windowVec parameter
*	@par Method The first window function, called a Rectangular Windowing function, 
*	is just to apply equal wieghting to all the values.  As such no actual formula 
*	is applied and the result is the same as the input.
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> Domain::applyRectangular(vector<double> windowVec){

	return windowVec;
}


/** 
*	applyHamming
*	@brief
*	@param windowVec a vector of double representing windowSize_ spectrum values
*	@post
*	@pre
*	@return a vector of double with the window wieghtings applied to the spectrum values
*	@par Method Similar to the Hann Window, the Hamming window also applies a cosine 
*	like wieghting.  The function to calculate a Hamming window over values $w$ is:
*
*	\begin{equation}
*		\forall_{i=0}^{n} w_{i} = w_{i} * 0.53836 - (\frac{0.46164 * cos(2\pi*i)}{n-1})
*	\end{equation}
*
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> Domain::applyHamming(vector<double> windowVec){

	double PI = settingsFile_.read<double>("PI", 3.14159265358979323846);

    	for (int i = 0; i < windowVec.size(); i++){
        	windowVec[i] *= 0.53836 - (0.46164 * cosf(2.0 * PI * static_cast<double>(i) / (windowVec.size() - 1)));
	}

	return windowVec;
}


/** 
*	applyHann
*	@brief
*	@param windowVec a vector of double representing windowSize_ spectrum values
*	@post
*	@pre
*	@return a vector of double with the window wieghtings applied to the spectrum values
*	@par Method Also called the Hanning Window the Hann window applies wieghtings in a 
*	cosine like fashion where in elements which are closer to the center of the window 
*	are assigned higher importance.  The function to calculate the Hann window over a 
*	window of values $w$ is:
*
*	\begin{equation}
*		\forall_{i=0}^{n} w_{i} = w_{i} * \frac{1}{2}(\frac{1-cos(2\pi*i)}{n-1})
*	\end{equation}
*
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> Domain::applyHann(vector<double> windowVec){

	double PI = settingsFile_.read<double>("PI", 3.14159265358979323846);

    	for (int i = 0; i < windowVec.size(); i++){
        	windowVec[i] *= 0.5 * (1.0 - cosf(2.0 * PI * static_cast<double>(i) / (windowVec.size() - 1)));
	}
	
	return windowVec;
}


/** 
*	applyBartlett
*	@brief
*	@param windowVec a vector of double representing windowSize_ spectrum values
*	@post
*	@pre
*	@return a vector of double with the window wieghtings applied to the spectrum values
*	@par Method The Bartlett windowing function, like the Hann and Hamming windowing 
*	functions, applies a greater wieght to the values at the center of the window than 
*	the extremes.  However, unlike the cosine shape of the previous functions Bartlett 
*	windowing functions are triangular in shape, meaning that the progression in wieght 
*	to and from the center location linearly increases and decreases respectively.  The 
*	function to calculate a Bartlett window over values $w$ is:
*
*	\begin{equation}
*		\forall_{i=0}^{n} w_{i} =
*			 w_{i} * \frac{2}{n-1} * (\frac{n-1}{2} - \mid i - \frac{n-1}{2} \mid)
*	\end{equation}
*
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> Domain::applyBartlett(vector<double> windowVec){

    	for (int i = 0; i < windowVec.size(); i++){
        	windowVec[i] *= 2.0 / (windowVec.size() - 1) * ((windowVec.size() - 1) / 2.0 - fabsf(i - (windowVec.size() - 1) / 2.0));
	}

	return windowVec;
}


/** 
*	applyTriangular
*	@brief
*	@param windowVec a vector of double representing windowSize_ spectrum values
*	@post
*	@pre
*	@return a vector of double with the window wieghtings applied to the spectrum values
*	@par Method Note that the Bartlett window above calculates a function which will 
*	have end points which are values of 0.  The triangular window is a Bartlett window 
*	which evaluates to have non-zero end points.  To accomplish this the only difference 
*	between the two functions is in a Triangular window we do not subtract 1 from $n$ in 
*	the two fractions. The function to calculate a Triangular Window over values $w$ is:
*
*	\begin{equation}
*		\forall_{i=0}^{n} w_{i} = 
*			w_{i} * \frac{2}{n} * (\frac{n}{2} - \mid i - \frac{n-1}{2} \mid)
*	\end{equation}
*
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> Domain::applyTriangular(vector<double> windowVec){

    	for (int i = 0; i < windowVec.size(); i++){
        	windowVec[i] *= 2.0 / windowVec.size() * (windowVec.size() / 2.0 - fabsf(i - (windowVec.size() - 1) / 2.0));
	}

	return windowVec;
}


/** 
*	applyBartlett_hann
*	@brief
*	@param windowVec a vector of double representing windowSize_ spectrum values
*	@post
*	@pre
*	@return a vector of double with the window wieghtings applied to the spectrum values
*	@par Method The Bartlett Hann window is a cross between the triangular shape of 
*	Bartlett and cosine shape of Hann.  It is calculated over a window of values $w$ via:
*
*	\begin{equation}
*		\forall_{i=0}^{n} w_{i} = 
*			w_{i} * 0.62 - 0.48 \mid \frac{i}{n-1} - \frac{1}{2} \mid 
*				- 0.38 cos(\frac{2\pi i}{n-1})
*	\end{equation}
*
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> Domain::applyBartlett_hann(vector<double> windowVec){

	double PI = settingsFile_.read<double>("PI", 3.14159265358979323846);

	double a0 = 0.62;
	double a1 = 0.5; //0.48?
	double a2 = 0.38;
	double term1 = 0.0;
	double term2 = 0.f;

    	for (int i = 0; i < windowVec.size(); i++){

        	term1 = a1 * fabsf(i / (windowVec.size() - 1) - 0.5);
        	term2 = a2 * cosf(2.0 * PI * static_cast<double>(i) / (windowVec.size() - 1));

        	windowVec[i] *= a0 - term1 - term2;
    	}

	return windowVec;
}


/** 
*	applyBlackman
*	@brief
*	@param windowVec a vector of double representing windowSize_ spectrum values
*	@post
*	@pre
*	@return a vector of double with the window wieghtings applied to the spectrum values
*	@par Method The Blackman window function is calculated over a window of values $w$ via:
*
*	\begin{equation}
*		\forall_{i=0}^{n} w_{i} = 
*			w_{i} * 0.42 - \frac{1}{2}cos(\frac{2\pi i}{n-1}) - 0.08 cos(\frac{4 \pi i}{n-1})
*	\end{equation}
*
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> Domain::applyBlackman(vector<double> windowVec){

	double PI = settingsFile_.read<double>("PI", 3.14159265358979323846);

       	double a0 = 0.42;
    	double a1 = 0.5;
     	double a2 = 0.08;
    	double term1 = 0.0;
     	double term2 = 0.0;

    	for (int i = 0; i < windowVec.size(); i++) {
    
        	term1 = a1 * cosf(2.0 * PI * static_cast<double>(i) / (windowVec.size() - 1));
        	term2 = a2 * cosf(4.0 * PI * static_cast<double>(i) / (windowVec.size() - 1));

        	windowVec[i] *= a0 - term1 + term2;
    	}

	return windowVec;
}


/** 
*	applyBlackman_harris
*	@brief
*	@param windowVec a vector of double representing windowSize_ spectrum values
*	@post
*	@pre
*	@return a vector of double with the window wieghtings applied to the spectrum values
*	@par Method The Blackman Harris window function is calculated over a window of values $w$ via:
*
*	\begin{equation}
*		\forall_{i=0}^{n} w_{i} = 
			w_{i} * 0.35875 - 0.48829cos(\frac{2 \pi i}{n-1}) 
				+ 0.14128cos(\frac{4 \pi i}{n-1}) - 0.01168cos(\frac{6 \pi i}{n-1})
*	\end{equation}
*
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> Domain::applyBlackman_harris(vector<double> windowVec){

	double PI = settingsFile_.read<double>("PI", 3.14159265358979323846);

       	double a0 = 0.35875;
      	double a1 = 0.48829;
      	double a2 = 0.14128;
      	double a3 = 0.01168;
   	double term1 = 0.0;
      	double term2 = 0.0;
       	double term3 = 0.0;

    	for (int i = 0; i < windowVec.size(); i++) {

        	term1 = a1 * cosf(2.0 * PI * i / (windowVec.size() - 1));
        	term2 = a2 * cosf(4.0 * PI * i / (windowVec.size() - 1));
        	term3 = a3 * cosf(6.0 * PI * i / (windowVec.size() - 1));

        	windowVec[i] *= a0 - term1 + term2 - term3;
    	}

	return windowVec;
}


