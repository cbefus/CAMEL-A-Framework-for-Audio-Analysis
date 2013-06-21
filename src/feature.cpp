/**
*	@file feature.cpp
*	@class Feature 
*	
*	@author Chad R. Befus, Chris Sanden
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
*	@brief Given a vector of appropriate domain values extracts a feature from them
*
*	Given a window of values this class can be called to extract a representative 
*	statistic over those values.
*	
*
*	Required Dependencies
*	-----------
*	std::vector
*	std::string
*	std::math.h
*	domain.h
*	configFile.h
*
*	Example Usage
*	-----------
	//should really be done from the featureExtract.cpp level...
*	#include "feature.h"
*	#include <vector>
*	using namespace std;
*
*	int main(){	 
*		Feature f;
*		f.setWindowSize(1024);
*		double value = f.calcWindowMean(domainWindowVector); //note see domain.h for how to get a domainwindovector
*		cout << value << endl;
*	}
*
*
*	More examples to be found in example folder
*
*	@todo replace namespace std with localised std
*	
*
**/
#include "feature.h"

#include <math.h>
#include <vector>

#include "fileVector.h"
#include "configFile.h"

using namespace std;

/** 
*	Feature
*	@brief constructor - inits the private member windowSize_
*	@param windowSize an int representing the window size of a frame
*	@post windowSize_ is set
*	@pre windowSize is a positive int
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
Feature::Feature(int windowSize)
	:windowSize_(windowSize){}

/** 
*	setWindowSize
*	@brief accessor to setting the windowSize_ private member
*	@param windowSize an int representing window size of a frame for extraction
*	@post windowSize_ is set
*	@pre windowSize is a positive integer
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void Feature::setWindowSize(int windowSize){
	windowSize_ = windowSize;
}




/** 
*	calcWindowMean
*	@brief calculates the mean over the parameter vector
*	@param windowVec a vector of ints representing time domain values
*	@post
*	@pre
*	@return a double representing the mean of windowVec
*	@par Method The mean, or average, is a statistical measure of a window 
*	and the balancing point of the data within.
*	
*	\begin{equation}
*		M(w) = \frac{1}{n} \sum_{i=0}^{n-1}w_{i}
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowMean(vector<int> windowVec){

	double result = 0.0;

	for(int i = 0; i < windowVec.size(); i++){
		result += windowVec[i];
	}

	if ( result == 0.0  ) {
		return 0.0;
	}
	
	return result / windowVec.size();

}

/** 
*	calcWindowVariance
*	@brief calculates the variance over the parameter vector
*	@param windowVec a vector of ints representing time domain values
*	@post
*	@pre
*	@return a double representing the variance of windowVec
*	@par Method The variance of a window is a measure of the distribution 
*	of data within that window away from the mean.
*
*	\begin{equation}
*		V(w) = \frac{1}{n} \sum_{i=0}^{n-1} (w_{i} - M(w))^2
*	\end{equation}
*	@note	
*	@remarks This is the variation of a finite (and known) population, not a sample.  
*	For this reason we do not divide by n-1... other libraries may do so.
*	@warning
*	@todo
**/
double Feature::calcWindowVariance(vector<int> windowVec){
	double result = 0.0;
	double mean = calcWindowMean(windowVec);

	for(int i = 0; i < windowVec.size(); i++){
		result += pow(windowVec[i] - mean, 2);
	}

	if ( result == 0.0  ) { 
		return 0.0;
	}

	return result / (windowVec.size() - 1);

}

/** 
*	calcWindowStandardDeviation
*	@brief calculates the standard deviation over the parameter vector
*	@param windowVec a vector of ints representing time domain values
*	@post
*	@pre
*	@return a double representing the standard deviation of windowVec
*	@par Method Standard Deviation is another measure of data distribution 
*	within the window away from the mean.
*
*	\begin{equation}
*		StD(w) = \sqrt{V(w)}
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowStandardDeviation(vector<int> WindowVec){
	
	double variance = calcWindowVariance(WindowVec);
	return sqrt(variance);
}

/** 
*	calcWindowAverageDeviation
*	@brief calculates the average deviation over the parameter vector
*	@param windowVec a vector of ints representing time domain values
*	@post
*	@pre
*	@return a double representing the average deviation of windowVec
*	@par Method Average Deviation is a third method of measuring data 
*	distribution within a window and away from the mean.
*
*	\begin{equation}
*		AvD = \frac{1}{n} \sum_{i=0}^{n-1} \mid w_{i} - M(w) \mid
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowAverageDeviation(vector<int> windowVec){
	double result = 0.0;
	double mean = calcWindowMean(windowVec);

	for(int i = 0; i < windowVec.size(); i++){
		result += fabs(windowVec[i] - mean);
	}

	if ( result == 0.0 ) { 
		return 0.0;
	}

	return result / windowVec.size();

}

/** 
*	calcWindowSkewness
*	@brief calculates the skewness over the parameter vector
*	@param windowVec a vector of ints representing time domain values
*	@post
*	@pre
*	@return a double representing the skewness of windowVec
*	@par Method Skewness is a measure of the data distribution which 
*	tends furthest from the mean.
*
*	\begin{equation}
*		Sk(w) = \frac{1}{n} \sum_{i=0}^{n-1} (\frac{w_{i} - M(w)}{StD(w)})^3
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSkewness(vector<int> windowVec){
	double result = 0.0;
	double mean = calcWindowMean(windowVec);
	double standardDeviation = calcWindowStandardDeviation(windowVec);

	for(int i = 0; i < windowVec.size(); i++){
		result += pow( (windowVec[i] - mean) / standardDeviation, 3);
	}

	if ( result == 0.0 || result != result) { //result = 0 or nan
		return 0.0;
	}

	return result / windowVec.size();

}

/** 
*	calcWindowKurtosis
*	@brief calculates the kurtosis over the parameter vector
*	@param windowVec a vector of ints representing time domain values
*	@post
*	@pre
*	@return a double representing the Kurtosis of windowVec
*	@par Method Kurtosis is much like skewness but deals even more so 
*	with measuring the extremes of the data distribution.
*
*	\begin{equation}
*		Ku(w) = (\frac{1}{n} \sum_{i=0}^{n-1} (\frac{w_{i} - M(w)}{StD(w)})^4) - 3
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowKurtosis(vector<int> windowVec){
	double result = 0.0;
	double mean = calcWindowMean(windowVec);
	double standardDeviation = calcWindowStandardDeviation(windowVec);

	for(int i = 0; i < windowVec.size(); i++){
		result += pow( (windowVec[i] - mean) / standardDeviation, 4);
	}

	if ( result == 0.0  || result != result) { //result = 0 or nan
		return 0.0;
	}

	return (result / windowVec.size()) - 3.0;

}

/** 
*	calcWindowZCR
*	@brief calculates the zero crossing rate over the parameter vector
*	@param windowVec a vector of ints representing time domain values
*	@post
*	@pre
*	@return a double representing the ZCR of windowVec
*	@par Method Is calcuated over the time domain and represents the frequency 
*	at which the signal crosses over from positive to negative or vice versa.  
*	ZCR has been associated with several different aspects of music including 
*	``the dominant frequency'' \cite{cord2008}. 
*
*	\begin{equation}
*		ZCR(w) = \frac{1}{n} \mid\mid w_{i} * w_{i+1} < 0 \mid\mid
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowZCR(vector<int> windowVec){
	int result = 0;

	double prevVal = windowVec[0];
 
	for (int i = 1; i < windowVec.size(); i++){
		double currVal = windowVec[i];
		if(currVal * prevVal < 0){
			++result;
		}
		prevVal = currVal;
	}

	if ( result == 0 ) {
		return 0.0;
	}

	return static_cast<double>(result)/(windowVec.size());
}

/** 
*	calcWindowRMS
*	@brief calculates the root mean squared (energy) over the parameter vector
*	@param windowVec a vector of ints representing time domain values
*	@post
*	@pre
*	@return a double representing the RMS of windowVec
*	@par Method RMS energy has been attributed to a good indication of loudness as 
*	well as a good statistic on which to conduct high level MIR tasks such as 
*	segmentation \cite{cord2008}.
*
*	\begin{equation}
*		RMS(w) = \sqrt{\frac{1}{n}\sum_{i=1}^{n-1}w_{i}^2}
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowRMS(vector<int> windowVec){
	double result = 0.0;
	
	for (int i = 0; i < windowVec.size(); i++){
		result += pow(windowVec[i],2);
	}

	if ( result == 0.0 ) {
		return 0.0;
	}

	result = sqrt(result/(windowVec.size()));

	return result;
}

/** 
*	calcWindowNonZeroCount
*	@brief calculates the non zero count over the parameter vector
*	@param windowVec a vector of ints representing time domain values
*	@post
*	@pre
*	@return a double representing the non zero count of windowVec
*	@par Method The Non Zero Count is just a count of the number of 
*	elements in a window which have some non-zero value.  It is a good 
*	measure of silence vs sound.
*
*	\begin{equation}
*		NZ(w) = \mid\mid w_{i} \neq 0 \mid\mid
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowNonZeroCount(vector<int> windowVec){
	int nonZeroCntr = 0;
	for(int i = 0; i < windowVec.size(); i++){
		if(windowVec[i] != 0){
			++nonZeroCntr;
		}
	}
	return static_cast<double>(nonZeroCntr);
}


/** 
*	calcWindowSpectralCentroid
*	@brief calculates the centroid over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre windowSize_ is set
*	@return a double representing the spectral centroid of the parameter vector
*	@par Method Spectral Centroid measures the balancing point of the spectrum, 
*	that is, the frequency where the energy below matches the energy above.  It 
*	is considered a measure of brightness and spectral shape \cite{cord2008}. 
*
*	\begin{equation}
*		SCe(w) = \frac{\sum_{i=0}^{(n/2)-1}w_{i} * w_{(n/2)+i}}{\sum_{i=0}^{(n/2)-1}w_{i}}
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralCentroid(vector<double> windowVec){
	
	double freqsByAmpsSum = 0;
	double ampsSum = 0;
	for(int j = 0; j < windowSize_/2.0; j++){
		freqsByAmpsSum += windowVec[j] * windowVec[static_cast<int>(windowSize_/2.0)+j];
		ampsSum += windowVec[j];
	}

	if(freqsByAmpsSum == 0.0 || ampsSum == 0.0){
		return 0.0;
	}
	return freqsByAmpsSum / ampsSum;
}

/** 
*	calcWindowSpectralVariance
*	@brief calculates the variance over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre windowSize_ is set
*	@return a double representing the spectral variance of the parameter vector
*	@par Method Spectral centroid is the first statistic required when attempting 
*	to measure spectral shape however several others are also required to gain an 
*	accurate understanding of the shape.  Spectral Variance describes how uniform 
*	or variable the sprectral shape is around the Spectral Centroid.
*
*	\begin{equation}
*		SVa(w) = \frac{\sum_{i=0}^{(n/2)-1}((w_{i} * w_{(n/2)+i})-SCe(w))^2}{\sum_{i=0}^{(n/2)-1}w_{i}}
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralVariance(vector<double> windowVec){
	
	double specMean = calcWindowSpectralCentroid(windowVec);
	double ampsSum = 0;
	double varianceSum = 0;

	for(int j = 0; j < windowSize_/2.0; j++){
		ampsSum += windowVec[j];
		varianceSum+= powf((windowVec[static_cast<int>(windowSize_/2.0)+j]  - specMean) * windowVec[j], 2);
	}
	
	if( varianceSum == 0.0 || ampsSum == 0.0 ){
		return 0.0;
	}

	return varianceSum / ampsSum;
}


/** 
*	calcWindowSpectralStandardDeviation
*	@brief calculates the standard deviation over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre 
*	@return a double representing the spectral standard deviation of the parameter vector
*	@par Method An alternative measurement of Spectral Variance is Spectral Standard 
*	Deviation, which is really just the square root of the first..
*
*	\begin{equation}
*		SStD(w) = \sqrt{SVa(w)}
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralStandardDeviation(vector<double> windowVec){
	double specVariance = calcWindowSpectralVariance(windowVec);
	return sqrt(specVariance);
}


/** 
*	calcWindowSpectralAverageDeviation
*	@brief calculates the average deviation over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre windowSize_ is set
*	@return a double representing the spectral average deviation of the parameter vector
*	@par Method Similar to Spectral Variance, Spectral Average Deviation takes only 
*	the absolute value of the difference between each of the values and the Spectral Centroid.
*
*	\begin{equation}
*		SAvD(w) = \frac{\sum_{i=0}^{(n/2)-1} \mid (w_{i} * w_{(n/2)+i}) - SCe(w) \mid}{\sum_{i=0}^{(n/2)-1}w_{i}}
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralAverageDeviation(vector<double> windowVec){
	
	double specMean = calcWindowSpectralCentroid(windowVec);

	double ampsSum = 0;
	double deviationSum = 0;

	for(int j = 0; j < windowSize_/2.0; j++){
		ampsSum += windowVec[j];
		deviationSum += fabs((windowVec[j]*windowVec[static_cast<int>(windowSize_/2.0)+j]) - specMean);
	}

	if( deviationSum == 0.0 || ampsSum == 0.0 ){
		return 0.0;
	}

	return deviationSum / ampsSum;
}

/** 
*	calcWindowSpectralSkewness
*	@brief calculates the skewness over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre windowSize_ is set
*	@return a double representing the spectral skewness of the parameter vector
*	@par Method Spectral Skewness is another measure of spectral shape.  Skewness 
*	is a measurement of the spectral distribution which is further away from the 
*	Spectral Centroid.
*
*	\begin{equation}
*		SSk(w) = \frac{ \sum_{i=0}^{(n/2)-1}((w_{i} * w_{(n/2)+i}-SCe(w)) / SStD(w))^3 }{ \sum_{i=0}^{(n/2)-1}w_{i} }
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralSkewness(vector<double> windowVec){
	
	double specMean = calcWindowSpectralCentroid(windowVec);
	double specStdDev = calcWindowSpectralStandardDeviation(windowVec);
	
	double ampsSum = 0;
	double temp = 0;
	double skewnessSum = 0;

	for(int j = 0; j < windowSize_/2.0; j++){
		ampsSum += windowVec[j];
		skewnessSum += pow(((windowVec[j]*windowVec[static_cast<int>(windowSize_/2.0)+j]) - specMean)/specStdDev, 3);
	}
	
	if( skewnessSum == 0.0 || ampsSum == 0.0 ){
		return 0.0;
	}

	return skewnessSum / ampsSum;
}


/** 
*	calcWindowSpectralKurtosis
*	@brief calculates the kurtosis over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre windowSize_ is set
*	@return a double representing the spectral kurtosis
*	@par Method Spectral Kurtosis is similar to Spectral Skewness only it focuses 
*	on measuring the extremes of the spectral distribution even more so than Spectral Skewness.
*
*	\begin{equation}
*		SKu(w) = (\frac{ \sum_{i=0}^{(n/2)-1}((w_{i} * w_{(n/2)+i}-SCe(w)) / SStD(w))^4 }{ \sum_{i=0}^{(n/2)-1}w_{i} }) - 3
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralKurtosis(vector<double> windowVec){
	
	double specMean = calcWindowSpectralCentroid(windowVec);
	double specStdDev = calcWindowSpectralStandardDeviation(windowVec);

	double ampsSum = 0;
	double temp = 0;
	double kurtSum = 0;

	for(int j = 0; j < windowSize_/2.0; j++){
		ampsSum += windowVec[j];
		kurtSum += pow(((windowVec[j]*windowVec[static_cast<int>(windowSize_/2.0)+j]) - specMean)/specStdDev, 4);
	}

	if( kurtSum == 0.0 || ampsSum == 0.0 ){
		return 0.0;
	}

	return (kurtSum / ampsSum) - 3.0;
}


/** 
*	calcWindowSpectralIrregularityK
*	@brief calculates the irregularity K over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre windowSize_ is set
*	@return a double representing the irregularity K
*	@par Method Spectral Irregularity is a measurement of spectral distribution 
*	and has been used as a measure of the noisyness of a spectral shape.  Two 
*	different Irregularity functions are provided in CAMEL, named Irregularity K 
*	and Irregularity J.
*
*	\begin{equation}
*		SIk(w) = \sum_{i=1}^{(n/2)-2} \mid w_{i} - \frac{w_{i-1} + w_{i} + w_{i+1}}{3} \mid
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralIrregularityK(vector<double> windowVec){
	
	double irregSum = 0;

	for(int j = 0; j < (windowSize_/2.0)-1; j++){
		irregSum += fabs(windowVec[j] - (windowVec[j-1] + windowVec[j] + windowVec[j+1]) / 3);
	}

	return irregSum;
}


/** 
*	calcWindowSpectralIrregularityJ
*	@brief calculates the irregularity J over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre
*	@return a double representing the irrgularity J
*	@par Method 
*	\begin{equation}
*		SIj(w) = \frac{ \sum_{i=0}^{(n/2)-2}(w_{i} - w_{i+1})^2 }{ \sum_{i=0}^{(n/2)-2} w_{i}^2}
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralIrregularityJ(vector<double> windowVec){
	
	double irregNumSum = 0;
	double irregDenSum = 0;

	for(int j = 0; j < (windowSize_/2.0); j++){
		irregNumSum += pow(windowVec[j] - windowVec[j+1], 2);
		irregDenSum += pow(windowVec[j], 2);
	}

	if( irregNumSum == 0.0 || irregDenSum == 0.0 ){
		return 0.0;
	}
	
	return irregNumSum/irregDenSum;
}

/** 
*	calcWindowSpectralFlatness
*	@brief calculates the flatness over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre windowSize_ is set, settingsFile is set
*	@return a double representing the spectral flatness of the parameter vector
*	@par Method Spectral flatness measures the deviation of the spectrum of each 
*	window from a flat line \cite{cord2008}.
*
*	\begin{equation}
*		SFl(w) = 10 * Log_{10}(\frac{\prod_{i=0}^{(n/2)-1}w_{i}^{1/\frac{n}{2}}}{ \sum_{i=0}^{(n/2)-1}w_{i} })
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralFlatness(vector<double> windowVec){
	
	double flatNumSum;
	
	if (windowVec[0] == 0.0){
		flatNumSum = 1.0;
	}else{
		flatNumSum = windowVec[0];
	}

	double flatDenSum = windowVec[0];
	double temp = 0;

	for(int j = 0; j < (windowSize_/2.0); j++){
		if(windowVec[j] != 0.0){
			flatNumSum *= windowVec[j];
			flatDenSum += windowVec[j];
		}
	}
	
	flatNumSum *= powf(flatNumSum, 1.0/(windowSize_/2.0));
	flatDenSum /= (windowSize_/2.0);
	
	double smallNumLimit = settingsFile_.read<double>("SMALL_NUMBER_LIMIT", 2e-42);
	if(flatNumSum < smallNumLimit){
		flatNumSum = smallNumLimit;
	}

	if(flatDenSum < smallNumLimit){
		flatDenSum = smallNumLimit;
	}	

	return	10 * log10(flatNumSum / flatDenSum);


}


/** 
*	calcWindowSpectralTonality
*	@brief calculates the tonality over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre 
*	@return a double representing the spectral tonality of the parameter vector
*	@par Method Spectral Tonality is an altered measure of spectral flatness.
*
*	\begin{equation}
*		STo(w) = SFl(w) / -60 ~where~STo(w) \geq 1 ~otherwise~return~1
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralTonality(vector<double> windowVec){
	double specFlatness = calcWindowSpectralFlatness(windowVec);

	if( specFlatness / -60.0 < 1.0){
		return specFlatness / -60.0;
	} else {
		return 1.0;
	}
}


/** 
*	calcWindowSpectralMin
*	@brief calculates the min over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre windowSize_ is set
*	@return a double representing the spectral min of the parameter vector
*	@par Method The Spectral Min is the smallest spectral value in a given window.
*
*	\begin{equation}
*		SMin(w) = Min_{i=0}^{(n/2)-1} (w{i})
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralMin(vector<double> windowVec){
	
	double minVal = windowVec[0];
		
	for(int j = 1; j < (windowSize_/2.0); j++){
		if(windowVec[j] < minVal){
			minVal = windowVec[j];
		}
	}

	return minVal;
}


/** 
*	calcWindowSpectralMax
*	@brief calculates the max over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre windowSize_ is set
*	@return a double representing the spectral max of the parameter vector
*	@par Method The Spectral Max is the greates spectral value in a given window.
*
*	\begin{equation}
*		SMax(w) = Max_{i=0}^{(n/2)-1} (w{i})
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralMax(vector<double> windowVec){
	
	double maxVal = 0.0;
	maxVal = windowVec[0];
	
	for(int j = 1; j < (windowSize_/2.0); j++){
		if(windowVec[j] > maxVal){
			maxVal = windowVec[j];
		}
	}
	return maxVal;
}


/** 
*	calcWindowSpectralCrest
*	@brief calculates the crest over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre
*	@return a double representing the spectral crest of the parameter vector
*	@par Method Spectral Crest is a measurement of spectral distribution from the Spectral Centroid.
*
*	\begin{equation}
*		SCr(w) = \frac{SMax(w)}{SCe(w)}
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralCrest(vector<double> windowVec){

	double specMean = calcWindowSpectralCentroid(windowVec);
	double specMax = calcWindowSpectralMax(windowVec);

	if( specMax == 0.0 || specMean == 0.0 ){
		return 0.0;
	}

	return specMax / specMean;
}

/** 
*	calcWindowSpectralSlope
*	@brief calculates the slope over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre windowSize_ is set
*	@return a double representing the spectral slope of the parameter vector
*	@par Method Also known as Spectral Gradient, Spectral Slope is a measure of 
*	how quickly the waveform tends towards the higher frequencies or a measure 
*	of how little energy exists at the higher frequencies.
*
*	\begin{equation}
*		SSlo(w) = \frac{\frac{1}{\sum_{i=0}^{(n/2)-1}w_{(n/2)+i}} * (\frac{n}{2}\sum_{i=0}^{(n/2)-1} w_{i} * w_{(n/2)+i}- \sum_{i=0}^{(n/2)-1} w_{(n/2)+i} * \sum_{i=0}^{(n/2)-1} w_{i})}{\frac{n}{2}\sum_{i=0}^{(n/2)-1}w_{i}^{2} - (\sum_{i=0}^{(n/2)-1}w_{i})^2}
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralSlope(vector<double> windowVec){

	double freqsSum = 0;
	double ampsSum = 0;
	double freqsByAmpsSum = 0;
	double freqsSquaredSum = 0;

	for(int j = 0; j < (windowSize_/2.0); j++){
		freqsSum += windowVec[j];
		ampsSum += windowVec[static_cast<int>(windowSize_/2.0)+j];
		freqsByAmpsSum += windowVec[j] * windowVec[static_cast<int>(windowSize_/2.0)+j];
		freqsSquaredSum += pow(windowVec[j],2);
	}

	double slopeNum = (1.0/ampsSum)*((windowSize_/2.0)*freqsByAmpsSum-freqsSum*ampsSum);
	double slopeDen = (windowSize_/2.0)*freqsSquaredSum - freqsSum * freqsSum;

	if( slopeNum == 0.0 || slopeDen == 0.0){
		return 0.0;
	}

	return slopeNum / slopeDen;
}

/** 
*	calcWindowSpectralSpread
*	@brief calculates the spread over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre windowSize_ is set
*	@return a double representing the spectral spread of the parameter vector
*	@par Method Spectral spread indicates whether or not the spectrum is centered 
*	near the spectral centroid or spread out over the entire spectrum.  As described 
*	by Cord \cite{cord2008} this potentially enables discrimination between pure 
*	and noisy sounds.
*
*	\begin{equation}
*		SSp(w) = \sqrt{ \frac{ \sum_{i=0}^{(n/2)-1}( i - SCe(w) )^{2} * w_{i} }{ \sum_{i=0}^{(n/2)-1}w_{i} } }
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralSpread(vector<double> windowVec){
	
	double specMean = calcWindowSpectralCentroid(windowVec);
	double numSum = 0.0;
	double denSum = 0.0;

	for(int j = 0; j < windowSize_/2.0; j++){
		numSum += pow(j-specMean, 2) * windowVec[j];
		denSum += windowVec[j];
	}

	if( numSum == 0.0 || denSum == 0.0){
		return 0.0;
	}

	return sqrt(numSum / denSum);
}

/** 
*	calcWindowSpectralRolloff
*	@brief calculates the rolloff over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre windowSize_ is set, settingsFile is set
*	@return a double representing the spectral rolloff of the parameter vector
*	@par Method Spectral Rolloff is a measure of spectral shape with a focus on 
*	the skewness of the spectrum \cite{cord2008}.
*
*	\begin{equation}
*		SRo(w) = \mid\mid \sum_{i=0}^{(n/2)-1}w_{i} < (\sum_{i=0}^{(n/2)-1}w_{i}) * P\mid\mid * \frac{SR}{n/2}
*	\end{equation}
*
*	where $P$ is a percentile (note Ong \cite{ong2007} sets it to 0.95, while Cord
*	 \cite{cord2008} sets it to 0.9 and jAudio sets it to 0.85) and $SR$ is the sample rate.
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralRolloff(vector<double> windowVec){

	double percentile = settingsFile_.read<double>("SPECTRAL_ROLLOFF_THRESH", 0.95); 
	double sampleRate = settingsFile_.read<double>("PCM_SAMPLE_RATE", 44100.0);	
	double pivot = 0.0;
		
	for(int j = 0; j < windowSize_/2.0; j++){
		pivot += windowVec[j];
	}
	pivot *= percentile;

	double temp = 0;
	int cntr = 0;		
	while(temp < pivot){
		temp += windowVec[cntr++];
	}
	return cntr * (sampleRate/(windowSize_/2.0));
}

/** 
*	calcWindowSpectralHPS
*	@brief calculates the Harmonic Product over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre windowSize_ is set, settingsFile is set
*	@return a double representing the harmonic product spectrum of the parameter vector
*	@par Method 
*	@note actually works correctly as opposed to other extraction libraries
*	@remarks
*	@warning
*	@todo
**/

double Feature::calcWindowSpectralHPS(vector<double> windowVec){
	
	double sampleRate = settingsFile_.read<double>("PCM_SAMPLE_RATE", 44100.0);	
	
	double coeffs2[static_cast<int>(windowSize_/2.0)];
	double coeffs3[static_cast<int>(windowSize_/2.0)];
	double product[static_cast<int>(windowSize_/2.0)];

	for(int j = 0; j < windowSize_/2.0; j++){
		coeffs2[j] = 1;
		coeffs3[j] = 1;
	}

	for(int j = 0; j < (windowSize_/2.0)/2.0; j++){
		coeffs2[j] = (windowVec[2*j] + windowVec[(2*j)+1]) * 0.5;
		if(j < ((windowSize_/2.0) / 3)){
			coeffs3[j] = (windowVec[j*3] + windowVec[(j*3)+1] + windowVec[(j*3)+2]) / 3;
		}
	}
		
	int peak_index = 0;
	double peak = 0.0;
	for(int j = 1; j < windowSize_/2.0; j++){
		product[j] = windowVec[j] * coeffs2[j] * coeffs3[j];
		if(product[j] > peak){
	 		peak_index = j;
	   		peak = product[j];
		}
	}

	int position1_lwr = 0;
	double largest1_lwr = 0;
    	for(int j = 0; j < windowSize_/2.0; j++){
		if(windowVec[j] > largest1_lwr && j != peak_index){
	   		largest1_lwr = windowVec[j];
	    		position1_lwr = j;
		}
    	}

	double ratio1 = windowVec[position1_lwr] / windowVec[peak_index];
	
	if(position1_lwr > (peak_index * 0.4) && position1_lwr < (peak_index * 0.6) && ratio1 > 0.1){
		peak_index = position1_lwr;
	}

	if( peak_index == 0){
		return 0.0;
	}

	return sampleRate / static_cast<double>(peak_index);
}

/** 
*	calcWindowSpectralLoudness
*	@brief calculates the loudness over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre windowSize_ is set
*	@return a double representing the spectral loudness of the parameter vector
*	@par Method Note spectral loudness is calculated over the bark bands as 
*	described below.  As such, we denote the set of bark bands as $b$ and the 
*	number of bark coefficients as $n$.
*
*	\begin{equation}
*		SLo(b) = \sum_{i=0}^{n-1}b_{i}^0.23
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralLoudness(vector<double> windowVec){

	vector<double> barkCoefs = calcWindowBark(windowVec);
	double returnVal = 0;

	for(int i = 0; i < barkCoefs.size(); i++){
		returnVal += powf(barkCoefs[i], 0.23);
	}

	return returnVal;

}

/** 
*	calcWindowSpectralSharpness
*	@brief calculates the sharpness over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre windowSize_ is set, settingsFile is set
*	@return a double representing the spectral sharpness of the parameter vector
*	@par Method Note spectral sharpness is calculated over the bark bands 
*	as described below.  As such, we denote the set of bark bands as $b$ and 
*	the number of bark coefficients as $n$ and the number of bark bands as $N$.
*
*	\begin{equation}
*		SSh(b) = (0.11 * \sum_{i=0}^{n-1} i * b_{i}^{0.23} )* N ~where~i < 15
*	
*		SSh(b) = (0.11 * \sum_{i=0}^{n-1} i * 0.066 *e^{0.171*i} * b_{i}^{0.23} )* N ~where~i \geq 15 
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowSpectralSharpness(vector<double> windowVec){
	vector<double> barkCoefs = calcWindowBark(windowVec);
	double returnVal = 0;

	for(int i = 0; i < barkCoefs.size(); i++){
		double sl = powf(barkCoefs[i], 0.23); //specific loudness
		double g = 1.0;
		if( i >= 15 ){
			g = 0.066 * expf(0.171 * i);
		}
		returnVal += i * g * sl;
	}
	returnVal = 0.11 * returnVal / settingsFile_.read<int>("NUM_BARK_BANDS", 26);

	return returnVal;
}


/** 
*	calcWindowTristimulus1
*	@brief calculates the tristimulus 1 over the parameter vector
*	@param windowVec a vector of doubles representing peak domain values
*	@post
*	@pre windowSize_ is set
*	@return a double representing the tristimulus 1 of the parameter vector
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowTristimulus1(vector<double> windowVec){
	
	int n = static_cast<int>(windowSize_ / 2.0);
	
	double den = 0.0;
	double p1 = 0.0;
	for( int i = 0; i < n; i++){
		if(windowVec[i]){
			den += windowVec[i];
			if(!p1){
				p1 = windowVec[i];
			}
		}
	}

	if(den == 0.0 || p1 == 0.0){
		return 0.0;
	}

	return p1 / den;
}

/** 
*	calcWindowTristimulus2
*	@brief calculates the tristimulus 2 over the parameter vector
*	@param windowVec a vector of doubles representing peak domain values
*	@post
*	@pre windowSize_ is set
*	@return a double representing the tristimulus 2 of the parameter vector
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowTristimulus2(vector<double> windowVec){
	
	int n = static_cast<int>(windowSize_ / 2.0);
	
	double den = 0.0;
	double ps = 0.0;
	double p2 = 0.0;
	double p3 = 0.0;
	double p4 = 0.0;
	for( int i = 0; i < n; i++){
		if(windowVec[i]){
			den += windowVec[i];
			if(!p2){ 
				p2 = windowVec[i];
			} else if(!p3){
				p3 = windowVec[i];
			} else if(!p4){
				p4 = windowVec[i];
			}
		}
	}

	ps = p2 + p3 + p4;

	if(den == 0.0 || ps == 0.0){
		return 0.0;
	}

	return ps / den;
}

/** 
*	calcWindowTristimulus3
*	@brief calculates the tristimulus 3 over the parameter vector
*	@param windowVec a vector of doubles representing peak domain values
*	@post
*	@pre windowSize_ is set
*	@return a double representing the tristimulus 3 of the parameter vector
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowTristimulus3(vector<double> windowVec){
	
	int n = static_cast<int>(windowSize_ / 2.0);
	int count = 0;	

	double den = 0.0;
	double num = 0.0;
	for( int i = 0; i < n; i++){
		if(windowVec[i]){
			den += windowVec[i];
			if(count >= 5){ 
				num += windowVec[i];
			}
			count++;
		}
	}

	if(den == 0.0 || num == 0.0){
		return 0.0;
	}

	return num / den;
}

/** 
*	calcWindowHarmonicOddEvenRatio
*	@brief calculates the odd even ratio over the parameter vector
*	@param windowVec a vector of doubles representing peak domain values
*	@post
*	@pre windowSize_ is set
*	@return a double representing the odd even ratio of the parameter vector
*	@par Method The Harmonic Odd Even Ratio is calculated over the Harmonic 
*	domain and as such we denote the window of input here as $h$ and its length as $n$.
*
*	\begin{equation}
*		HOEr(h) =  \frac{ \mid\mid h_{i} \mod 2 \neq 0 \mid\mid }{ \mid\mid h_{i} \mod 2 = 0 \mid\mid }
*	\end{equation}
*	@note	
*	@remarks
*	@warning
*	@todo
**/
double Feature::calcWindowHarmonicOddEvenRatio(vector<double> windowVec){
	
	int M = static_cast<int>(windowSize_ / 2.0);

    	double odd = 0.0;
	double even = 0.0;

    	for(int i = 0; i < M; i++){
		if(windowVec[i]){
	    		if(i % 2 != 0){
				odd += windowVec[i];
	    		} else {
				even += windowVec[i];
	    		}
		}
    	}

	if(odd == 0.0 || even == 0.0){
        	return 0.0;
    	}

    	return odd / even;
}

/** 
*	calcWindowMFCC
*	@brief calculates the MFCC's over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre settingsFile_ is linked, windowSize_ is set
*	@return a vector of double representing the Mel Scale coeficients over the parameter vector
*	@par Method MFCCs is one of the most popular features used in MIR.  The idea behind MFCC's 
*	is to bin frequencies into groups based on the Mel Scale.  The Mel Scale attempts to simulates 
*	the changes in human perception of sound as frequency changes. A good paper discussing MFCCs 
*	is \cite{logan2000} and our implementation is based on the description in \cite{ETSI2000}.
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> Feature::calcWindowMFCC(vector<double> windowVec){
	
	//extraction vals here
	int numCepstra = settingsFile_.read<int>("NUM_CEPSTRA", 13);
	double PI = settingsFile_.read<double>("PI", 3.14159265358979323846);
	int numMelFilters = settingsFile_.read<int>("NUM_MEL_FILTERS", 23);
	double lowerFilterFreq = settingsFile_.read<double>("LOWER_FILTER_FREQ", 133.3334);
	double upperFilterFreq = settingsFile_.read<double>("UPPER_FILTER_FREQ", 6855.4976);
	double samplingRate = settingsFile_.read<double>("PCM_SAMPLE_RATE", 44100.0);

	//Mel Filtering : center frequencies
	int cbin[numMelFilters + 2];
	cbin[0] = static_cast<int>(round(lowerFilterFreq / samplingRate * windowSize_)); //calculating cbin at 0
        cbin[numMelFilters+1] = static_cast<int>(windowSize_ / 2.0); //calculating cbin at 24
	for (int i = 1; i <= numMelFilters; i++){
		
		double mel[2];
        	mel[0] = 2595 * log10(1 + lowerFilterFreq / 700);
        	mel[1] = 2595 * log10(1 + (samplingRate / 2) / 700);
        
        	// take inverse mel of:
        	double inverse = mel[0] + ((mel[1] - mel[0]) / (numMelFilters + 1)) * i;
		inverse = pow(10, inverse / 2595) - 1;
        	double fcent = 700 * inverse;

            	cbin[i] = static_cast<int>(round(fcent / samplingRate * windowSize_));
        }

	double fbtemp[numMelFilters + 2];
        for (int k = 1; k <= numMelFilters; k++){
            	double num1 = 0, num2 = 0;
            	for (int m = cbin[k - 1]; m <= cbin[k]; m++){
			num1 += ((m - cbin[k - 1] + 1) / (cbin[k] - cbin[k-1] + 1)) * windowVec[m];
            	}
		
            	for (int m = cbin[k] + 1; m <= cbin[k + 1]; m++){
			num2 += (1 - ((m - cbin[k]) / (cbin[k + 1] - cbin[k] + 1))) * windowVec[m];
            	}

            	fbtemp[k] = num1 + num2;
        }

        double fbank[numMelFilters];
        for (int j = 0; j < numMelFilters; j++){
            	fbank[j] = fbtemp[j + 1];
        }
		
        // Non-linear transformation
	double f[numMelFilters];
        for (int j = 0; j < numMelFilters; j++){
            	f[j] = log(fbank[j]);
            	// check if ln() returns a value less than the floor
            	if (f[j] < -50){
			 f[j] = -50;
		}
	}

	// Cepstral coefficients
	double cepc[numCepstra];
	for(int j = 0; j <numCepstra; j++){
		cepc[j] = 0;
	}		

	vector<double> mfccVec(numCepstra);

        for (int j = 0; j < numCepstra; j++){
            	for (int k = 1; k <= numMelFilters; k++){
                	cepc[j] += f[k - 1] * cos(PI * j / numMelFilters * (k - 0.5));
            	}
		mfccVec[j] = cepc[j]; //add MFCC to the segment sum for averaging
        }

	return mfccVec;
}



/** 
*	calcWindowBark
*	@brief calculates the Bark Scale over the parameter vector
*	@param windowVec a vector of doubles representing spectral domain values
*	@post
*	@pre settingsFile_ is linked, windowSize_ is set
*	@return a vector of doubles representing the bark coefficients for the parameter vector
*	@par Method Another highly popular feature function is Bark Bands which separates the 
*	audio source into a series of bands which mimic the psychological bands of hearing. The 
*	Bark scale was proposed by Zwicker in \cite{zwicker1961}.
*	@note results from this function vary amoung extraction libraries based on precision of rounding
*	@remarks
*	@warning
*	@todo
**/
vector<double> Feature::calcWindowBark(vector<double> windowVec){
	
	//extraction vals here
	int numBarkBands = settingsFile_.read<int>("NUM_BARK_BANDS", 26);
	double sampleRate = settingsFile_.read<double>("PCM_SAMPLE_RATE", 44100.0);
	double edges[] = {0, 100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000, 15500, 20500, 27000}; 
	int band_limits[numBarkBands];

	for ( int i = 0; i <= numBarkBands; i++){
		band_limits[i] = static_cast<int>(round(edges[i] / sampleRate * (windowSize_/2.0)));
    	}


	vector<double> barkVec(numBarkBands);
	for(int j= 0; j < numBarkBands - 1; j++){ //only goes to numBarkBands - 1 = 0 for last element
		barkVec[j] = 0;
        	for(int n = band_limits[j]; n < band_limits[j + 1]; n++){
            		barkVec[j] += windowVec[n];
		}
	}

	return barkVec;
}

