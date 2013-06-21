/**
*	@file feature.h
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
*	std::math.h
*	domain.h
*	configFile.h
*
*	Example Usage
*	-----------
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
*	@todo
*	
*
**/

#ifndef __FEATURE_H
#define __FEATURE_H

#include <vector>

#include "configFile.h"

using std::vector;

class Feature{
    public:
	//constructor
	Feature(int windowSize = 1024); ///constructor - inits private member windowSize_

	//accessors
	void setWindowSize(int windowSize);  /// allows for setting of windowSize private member

	/* single return time domain values */
	double calcWindowMean(vector<int> windowVec); ///calculates the mean of a vector of ints
	double calcWindowVariance(vector<int> windowVec); ///calculates the variance of a vector of ints
	double calcWindowStandardDeviation(vector<int> windowVec); ///calculates the standard deviation of a vector of ints
	double calcWindowAverageDeviation(vector<int> windowVec); ///calculates the average deviation of a vector of ints
	double calcWindowSkewness(vector<int> windowVec); ///calculates the skewness of a vector of ints
	double calcWindowKurtosis(vector<int> windowVec); ///calculates the kurtosis of a vector of ints
	double calcWindowZCR(vector<int> windowVec); ///calcualtes the zero crossing rate of a vector of ints
	double calcWindowRMS(vector<int> windowVec); ///calculates the root mean squared (energy) of a vector of ints
	double calcWindowNonZeroCount(vector<int> windowVec); ///calculates the non zero count of a vector of ints

	/* single return spectral domain values */
	double calcWindowSpectralCentroid(vector<double> windowVec); ///calculates the centroid of a vector of doubles
	double calcWindowSpectralVariance(vector<double> windowVec); ///calculates the variance of a vector of doubles
	double calcWindowSpectralStandardDeviation(vector<double> windowVec); ///calculates the standard deviation of a vector of doubles
	double calcWindowSpectralAverageDeviation(vector<double> windowVec); ///calculates the average deviation of a vector of doubles
	double calcWindowSpectralSkewness(vector<double> windowVec); ///calculates the skewness of a vector of doubles
	double calcWindowSpectralKurtosis(vector<double> windowVec); ///calculates the Kurtosis of a vector of doubles
	double calcWindowSpectralIrregularityK(vector<double> windowVec); ///calculates the Irregularity K of a vector of doubles
	double calcWindowSpectralIrregularityJ(vector<double> windowVec); ///calculates the Irregularity J of a vector of doubles
	double calcWindowSpectralFlatness(vector<double> windowVec); ///calculates the Flatness of a vector of doubles
	double calcWindowSpectralTonality(vector<double> windowVec); ///calculates the Tonality of a vector of doubles
	double calcWindowSpectralMin(vector<double> windowVec); ///calculates the min of a vector of doubles
	double calcWindowSpectralMax(vector<double> windowVec); ///calculates the max of a vector of doubles
	double calcWindowSpectralCrest(vector<double> windowVec); ///calculates the crest of a vector of doubles
	double calcWindowSpectralSlope(vector<double> windowVec); ///calculates the slope of a vector of doubles
	double calcWindowSpectralSpread(vector<double> windowVec); ///calculates the spread of a vector of doubles
	double calcWindowSpectralRolloff(vector<double> windowVec); ///calculates the rolloff of a vector of doubles
	double calcWindowSpectralHPS(vector<double> windowVec); ///calculates the harmonic product of a vector of doubles
	double calcWindowSpectralLoudness(vector<double> windowVec); ///calculates the loudness of a vector of doubles
	double calcWindowSpectralSharpness(vector<double> windowVec); ///calculates the sharpness of a vector of doubles
	double calcWindowSpectralPitch(vector<double> windowVec); ///calculates the pitch of a vector of doubles

	/* single return Harmonic domain values */
	double calcWindowTristimulus1(vector<double> windowVec); ///calculates the tristimulus 1 of a vector of doubles
	double calcWindowTristimulus2(vector<double> windowVec); ///calculates the tristimulus 2 of a vector of doubles
	double calcWindowTristimulus3(vector<double> windowVec); ///calculates the tristimulus 3 of a vector of doubles
	double calcWindowHarmonicOddEvenRatio(vector<double> windowVec);///calculates the odd even ratio of a vector of doubles

	/* multi return spectral domain values */
	vector<double> calcWindowMFCC(vector<double> windowVec); ///calculates the MFCC of a vector of doubles
	vector<double> calcWindowBark(vector<double> windowVec); ///calculates the Bark scale of a vector of doubles

    private:
	/* Private Members */
	int windowSize_; ///the size of a single window of values
	ConfigFile settingsFile_; ///the connection object to the settings file


};
#endif
