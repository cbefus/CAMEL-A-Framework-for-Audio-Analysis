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
*	@todo
*	
*
**/
#ifndef __FEATUREEXTRACT_H
#define __FEATUREEXTRACT_H

#include <string>
#include <vector>

#include "configFile.h"
#include "domain.h"
#include "feature.h"

using namespace std;

class FeatureExtract{
    public:
	//constructor
        FeatureExtract(string fileName = "null.txt", int windowSize = 1024); ///constructor - inits private members
	
	//methods
	void setup(); ///inits dependency classes
	void getFeature(int startPos, int endPos, int feature); ///extracts the requested feature from between start and endPos

	//accessors     
	void setFileName(string fileName); ///sets the member fileName
	void setWindowSize(int windowSize); ///sets the member windowSize
	vector<double> getValues(); ///gets the results vector
             
    private:
	// Private Members
	Feature feature_; ///the Feature Object for the requested feature type
	Domain domain_; ///the Domain Object for the requested feature type
	vector<double> featureValues_; ///the vector of resulting feature values
	ConfigFile settingsFile_;///the confiFile object for the settings file
	string fileName_; ///the fileName for the PCM data file
	int windowSize_; ///the windowSize for a single window of PCM data

	/* functions to avg another functions return values over a specific data vector input */
	double avgOverStatisticalWindow(double (Feature::*pt2Func)(vector<int>)); ///Averages the values over statistics from the time domain
	double avgOverSpectralWindow(double (Feature::*pt2Func)(vector<double>)); ///Averages the values over statistics from the spectral domain
	double avgOverPeakWindow(double (Feature::*pt2Func)(vector<double>)); ///Averages the values over statistics from the Peak domain
	vector<double> avgVectorOverSpectralWindow(vector<double> (Feature::*pt2Func)(vector<double>), int numReturn); ///Averages multi-dimensional statistics over the spectral domain
			
	/* private accessor functions */
        void setValues(double value); ///sets the results vector with a single dimension result
        void setValues(vector<double> values); ///sets the results vector with a multi dimensional result
           
};
#endif
