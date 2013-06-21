/**
*	@file featureExtract.h
*	@class FeatureExtract
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
*	@brief A controller class for the domain and feature classes
*
*	
*	Given a PCM file, a start position, end position, and feature selection, this 
*	class extracts the requested feature from the requested location in the file.
*
*	Required Dependencies
*	-----------
*	std::vector
*	std::math.h
*	featureExtract.h
*	configFile.h
*
*	Example Usage
*	-----------
*	#include "featureExtract.h"
*	#include <vector>
*	using namespace std;
*
*	int main(){	 
*
*		FeatureExtract fe;
*		fe.setFileName("myPCMFile.txt");
*		fe.setWindowSize(1024);
*		fe.setup();
*		fe.getFeature(44100, 88200, 3);
*		vector<double> vecFeatureResults = fe.getValues();
*		cout << vecFeatureResults[0] << endl;
*
*	}
*
*
*	More examples to be found in example folder
*
*	@todo replace namespace std with localised std
*	
*
**/

#include "featureExtract.h"

#include <math.h>
#include <vector>

#include "configFile.h"
#include "feature.h"
#include "domain.h"

using namespace std;

/** 
*	FeatureExtract
*	@brief Constructor - inits private members fileName_ and windowSize_
*	@param fileName	a string representing a pcm file object
*	@param windowSize an int representing the size of a window
*	@post fileName_ and windowSize_ are set
*	@pre fileName names an existing file, and windowSize is a positive integer
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
FeatureExtract::FeatureExtract(string fileName, int windowSize)
	:fileName_(fileName),
	windowSize_(windowSize){}

/** 
*	setup
*	@brief inits the settingsFile_, feature_, and domain_ objects
*	@param
*	@post settingsFile_, feature_, and domain_ are set up
*	@pre
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void FeatureExtract::setup(){
	settingsFile_.setFileName("settings.txt");
	feature_.setWindowSize(windowSize_);
	domain_.setFileName(fileName_);
	domain_.setWindowSize(windowSize_);
	domain_.setup();
}

/** 
*	getValues
*	@brief returns the feature results
*	@param
*	@post
*	@pre featureValues_ is not null
*	@return a vector of double: the feature results
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<double> FeatureExtract::getValues(){
	return featureValues_;
}

/** 
*	setFileName
*	@brief sets the fileName_ private member to the parameter value
*	@param fileName a string representing the new fileName
*	@post fileName_ is set
*	@pre fileName represents a pcm file object
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void FeatureExtract::setFileName(string fileName){
	fileName_ = fileName;
}

/** 
*	setWindowSize
*	@brief sets the windowSize_ private member with the parameter value
*	@param windowSize an int representing the window size
*	@post windowSize_ is set to parameter value
*	@pre windowSize is a positive integer
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void  FeatureExtract::setWindowSize(int windowSize){
	windowSize_ = windowSize;
}

/** 
*	getFeature
*	@brief A switch redirecting to the correct feature extraction method
*	@param startPos an int representing the start location of extraction
*	@param endPos an int representing the end location of extractino
*	@param feature an int representing the feature to extract
*	@post domain_ start and end postitions are set, featureValues_ is filled
*	@pre parameters are all positive integers
*	@return
*	@par Method Each case contains three nested functions.  The first (outer-
*	most) is the setValues funtion for setting our final results into the
*	featureValues_ vector.  The second (mid) is the appropriate averaging 
*	function for the requested feature (Note this means the feature values
*	are averaged over all windows between the start and end position). The
*	third (inner) is a function pointer to the Feature class with the requested
*	feature extraction algorithm which is passed as a parameter to the averaging
*	function.
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void FeatureExtract::getFeature(int startPos, int endPos, int feature){

	domain_.setEndPosition(endPos);
	domain_.setStartPosition(startPos);

	switch (feature){
		case 1:{
			//cout << "            Extracting Mean" << endl;
			setValues( avgOverStatisticalWindow( &Feature::calcWindowMean ) );
			break;
		}
		case 2:{
			//cout << "            Extracting Variance" << endl;
			setValues( avgOverStatisticalWindow( &Feature::calcWindowVariance ) );
			break;
		}
		case 3:{
			//cout << "            Extracting Standard Deviation" << endl;
			setValues( avgOverStatisticalWindow( &Feature::calcWindowStandardDeviation ) );
			break;
		}
		case 4:{
			//cout << "            Extracting Average Deviation" << endl;
			setValues( avgOverStatisticalWindow( &Feature::calcWindowAverageDeviation ) );
			break;
		}
		case 5:{
			//cout << "            Extracting Skewness" << endl;
			setValues( avgOverStatisticalWindow( &Feature::calcWindowSkewness ) );
			break;
		}
		case 6:{
			//cout << "            Extracting Kurtosis" << endl;
			setValues( avgOverStatisticalWindow( &Feature::calcWindowKurtosis ) );
			break;
		}
		case 7:{
			//cout << "            Extracting ZCR" << endl;
			setValues( avgOverStatisticalWindow( &Feature::calcWindowZCR ) );
			break;
		}
		case 9:{
			//cout << "            Extracting RMS" << endl;
			setValues( avgOverStatisticalWindow( &Feature::calcWindowRMS ) );
			break;
		}
		case 10:{
			//cout << "            Extracting Non-Zero Count" << endl;
			setValues( avgOverStatisticalWindow( &Feature::calcWindowNonZeroCount ) );
			break;
		}
		case 11:{
			//cout << "            Extracting Spectral Centroid" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralCentroid ) );
			break;
		}
		case 12:{
			//cout << "            Extracting Spectral Variance" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralVariance ) );
			break;
		}
		case 13:{
			//cout << "            Extracting Spectral Standard Deviation" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralStandardDeviation ) );
			break;
		}
		case 14:{
			//cout << "            Extracting Spectral Average Deviation" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralAverageDeviation ) );
			break;
		}
		case 15:{
			//cout << "            Extracting Spectral Skewness" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralSkewness ) );
			break;
		}
		case 16:{
			//cout << "            Extracting Spectral Kurtosis" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralKurtosis ) );
			break;
		}
		case 17:{
			//cout << "            Extracting Spectral Irregularity K" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralIrregularityK ) );
			break;
		}
		case 18:{
			//cout << "            Extracting Spectral Irregularity J" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralIrregularityJ ) );
			break;
		}
		case 19:{
			//cout << "            Extracting Spectral Flatness" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralFlatness ) );
			break;
		}
		case 20:{
			//cout << "            Extracting Spectral Tonality" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralTonality ) );
			break;
		}
		case 21:{
			//cout << "            Extracting Spectral Min" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralMin ) );
			break;
		}
		case 22:{
			//cout << "            Extracting Spectral Max" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralMax ) );
			break;
		}
		case 23:{
			//cout << "            Extracting Spectral Crest" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralCrest ) );
			break;
		}
		case 24:{
			//cout << "            Extracting Spectral Slope" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralSlope ) );
			break;
		}
		case 25:{
			//cout << "            Extracting Spectral Spread" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralSpread ) );
			break;
		}
		case 26:{
			//cout << "            Extracting Spectral Rolloff" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralRolloff ) );
			break;
		}
		case 27:{
			//cout << "            Extracting Spectral HPS" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralHPS ) );
			break;
		}
		case 28:{
			//cout << "            Extracting Spectral Loudness" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralLoudness ) );
			break;
		}
		case 29:{
			//cout << "            Extracting Spectral Sharpness" << endl;
			setValues( avgOverSpectralWindow( &Feature::calcWindowSpectralSharpness ) );
			break;
		}
		case 31:{
			//cout << "            Extracting Peak Tristimulus 1" << endl;
			setValues( avgOverPeakWindow( &Feature::calcWindowTristimulus1 ) );
			break;
		}
		case 32:{
			//cout << "            Extracting Peak Tristimulus 2" << endl;
			setValues( avgOverPeakWindow( &Feature::calcWindowTristimulus2 ) );
			break;
		}
		case 33:{
			//cout << "            Extracting Peak Tristimulus 3" << endl;
			setValues( avgOverPeakWindow( &Feature::calcWindowTristimulus3 ) );
			break;
		}
		case 40:{
			//cout << "            Extracting MFCC" << endl;
			setValues( avgVectorOverSpectralWindow( &Feature::calcWindowMFCC, settingsFile_.read<int>("NUM_CEPSTRA", 13)));
			break;
		}
		case 41:{
			//cout << "            Extracting Bark" << endl;
			setValues( avgVectorOverSpectralWindow( &Feature::calcWindowBark, settingsFile_.read<int>("NUM_BARK_BANDS", 26)));
			break;
		}
		default:{
			cout << "            Feature set for extraction Does not Exist." << endl;
			break;
		}
	}
}

/** 
*	avgOverStatisticalWindow
*	@brief Averages values over several windows of time domain statistics
*	@param *pt2Func is a pointer to a function in the Feature Class which takes a vector of ints and returns a double
*	@post
*	@pre domain is setup, *pt2Func relates to a function, windowSize_ is set
*	@return a double representing the average value over each of the windows
*	@par Method Using the start position and end position calculate the number of windows
*	required (using windowSize_) to cover the gap. For each of those windows, get the vector
*	of time domain values that make up the window.  For those values extract the appropriate
*	feature from the feature function pointed to by *pt2Func and add it to the avgSum. Finally
*	divide the avgSum by the number of windows to get the average and return this value.
*	@note	
*	@remarks Note using the function pointer saves us having to make a separate averaging function
*	for every feature function, or, even worse, putting the averaging code into each feature functino
*	itself.
*	@warning
*	@todo
**/
//avgOverTimeDomainWindow
double FeatureExtract::avgOverStatisticalWindow(double (Feature::*pt2Func)(vector<int>)){

	int startPosition = domain_.getStartPosition();
	int endPosition = domain_.getEndPosition();
	int numWindows = static_cast<int>(ceil((endPosition - startPosition)/static_cast<double>(windowSize_)));
	double avgSum = 0;

	for(int i = 0; i < numWindows; i++){
		vector<int> windowVec = domain_.getTimeDomainWindow(startPosition+(i*windowSize_));
		avgSum += (feature_.*pt2Func)(windowVec);	
	}

	if(avgSum == 0.0){
		return 0.0;
	}

	return avgSum/numWindows;
}

/** 
*	avgOverSpectralWindow
*	@brief Averages values over several windows of spectral domain statistics
*	@param *pt2Func is a pointer to a function in the Feature Class which takes a vector of doubles and returns a double
*	@post
*	@pre domain is setup, *pt2Func relates to a function, windowSize_ is set
*	@return a double representing the average value over each of the windows
*	@par Method Using the start position and end position calculate the number of windows
*	required (using windowSize_) to cover the gap. For each of those windows, get the vector
*	of spectral domain values that make up the window.  For those values extract the appropriate
*	feature from the feature function pointed to by *pt2Func and add it to the avgSum. Finally
*	divide the avgSum by the number of windows to get the average and return this value.
*	@note	
*	@remarks Note using the function pointer saves us having to make a separate averaging function
*	for every feature function, or, even worse, putting the averaging code into each feature functino
*	itself.
*	@warning
*	@todo
**/
//avgOverFrequencyDomainWindow
double FeatureExtract::avgOverSpectralWindow(double (Feature::*pt2Func)(vector<double>)){

	int startPosition = domain_.getStartPosition();
	int endPosition = domain_.getEndPosition();
	int numWindows = static_cast<int>(ceil((endPosition - startPosition)/static_cast<double>(windowSize_)));
	double avgSum = 0;

	for(int i = 0; i < numWindows; i++){
		vector<double> windowVec = domain_.getFrequencyDomainWindow(startPosition+(i*windowSize_));
		avgSum += (feature_.*pt2Func)(windowVec);
	}

	if(avgSum == 0.0){
		return 0.0;
	}

	return avgSum/numWindows;
}

/** 
*	avgOverPeakWindow
*	@brief Averages values over several windows of peak domain statistics
*	@param *pt2Func is a pointer to a function in the Feature Class which takes a vector of doubles and returns a double
*	@post
*	@pre domain is setup, *pt2Func relates to a function, windowSize_ is set
*	@return a double representing the average value over each of the windows
*	@par Method Using the start position and end position calculate the number of windows
*	required (using windowSize_) to cover the gap. For each of those windows, get the vector
*	of peak domain values that make up the window.  For those values extract the appropriate
*	feature from the feature function pointed to by *pt2Func and add it to the avgSum. Finally
*	divide the avgSum by the number of windows to get the average and return this value.
*	@note	
*	@remarks Note using the function pointer saves us having to make a separate averaging function
*	for every feature function, or, even worse, putting the averaging code into each feature functino
*	itself.
*	@warning
*	@todo
**/
//avgOverPeakDomainWindow
double FeatureExtract::avgOverPeakWindow(double (Feature::*pt2Func)(vector<double>)){

	int startPosition = domain_.getStartPosition();
	int endPosition = domain_.getEndPosition();
	int numWindows = static_cast<int>(ceil((endPosition - startPosition)/static_cast<double>(windowSize_)));
	double avgSum = 0;

	for(int i = 0; i < numWindows; i++){
		vector<double> windowPeakSpectrum = domain_.getPeakDomainWindow(startPosition+(i*windowSize_));
		avgSum += (feature_.*pt2Func)(windowPeakSpectrum);
	}

	if(avgSum == 0.0){
		return 0.0;
	}

	return avgSum/numWindows;
}

/** 
*	avgVectorOverSpectralWindow
*	@brief Averages values over several windows of multi dimensional spectral domain statistics
*	@param *pt2Func is a pointer to a function in the Feature Class which takes a vector of doubles 
*	and returns a  vector  of doubles
*	@post
*	@pre domain is setup, *pt2Func relates to a function, windowSize_ is set
*	@return a double representing the average value over each of the windows for each feature dimension
*	@par Method Using the start position and end position calculate the number of windows
*	required (using windowSize_) to cover the gap. For each of those windows, get the vector
*	of spectral domain values that make up the window.  For those values extract the appropriate
*	feature from the feature function pointed to by *pt2Func and add it to the avgSum. To add the 
*	multidimensional values simply loop through the dimensions and add them in turn. Finally
*	divide the avgSum of each dimension by the number of windows to get the average and return this value.
*	@note	
*	@remarks Note using the function pointer saves us having to make a separate averaging function
*	for every feature function, or, even worse, putting the averaging code into each feature functino
*	itself.
*	@warning
*	@todo
**/
vector<double> FeatureExtract::avgVectorOverSpectralWindow(vector<double> (Feature::*pt2Func)(vector<double>), int numReturn){

	int startPosition = domain_.getStartPosition();
	int endPosition = domain_.getEndPosition();
	int numWindows = static_cast<int>(ceil((endPosition - startPosition)/static_cast<double>(windowSize_)));

	vector<double> avgSum;
	for(int i = 0; i < numReturn; i++){
		avgSum.push_back(0.0);
	}

	for(int i = 0; i < numWindows; i++){
		vector<double> windowVec = domain_.getFrequencyDomainWindow(startPosition+(i*windowSize_));
		vector<double> windowFunc = (feature_.*pt2Func)(windowVec);
		for(int j = 0; j < numReturn; j++){
			avgSum[j] += windowFunc[j];
		}
	}
	
	for(int i = 0; i < numReturn; i++){
		if(avgSum[i] != 0.0){
			avgSum[i] = avgSum[i] / numWindows;
		}
	}

	return avgSum;

}


/** 
*	setValues
*	@brief Sets the results vector with a single value
*	@param value a double representing the value
*	@post featureValues_ is set
*	@pre value parameter is not null
*	@return
*	@par Method first clear the results vector and then set its first value to the parameters
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void FeatureExtract::setValues(double value){
	featureValues_.clear();
	featureValues_.push_back(value);
}

/** 
*	setValues
*	@brief Sets the results vector with a multi-dimensional value
*	@param values a vector of double representing the results values
*	@post featureValues_ is set
*	@pre values parameter is not null
*	@return
*	@par Method first clear the results vector and then set it to equal the parameter vector
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void FeatureExtract::setValues(vector<double> values){
    	featureValues_.clear();
	featureValues_ = values;
}

