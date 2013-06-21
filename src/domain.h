/**
*	@file domain.h
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
*		domain.setup();	
*		vector<double> domainVec = domain.getFrequencyDomainWindow(44100);
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
*
**/


#ifndef __DOMAIN_H
#define __DOMAIN_H

#include <vector>
#include <string>

#include "fileVector.h"
#include "configFile.h"

using std::vector;
using std::string;

class Domain{
    public:
	//constructor
	Domain(string fileName = "null.txt", int windowSize = 1024); ///Constructor - inits private members

	//methods
	void setup(); ///sets up the FileVector (a static dependency object)

	//accessors
	void setFileName(string fileName); ///Accessor to set private member fileName_
	void setWindowSize(int windowSize); ///Accessor to set private member windowSize_
	void setStartPosition(int begin); ///Accessor to set private member begin_
	void setEndPosition(int end); ///Accessor to set private member end_
	int getStartPosition(); ///Accessor to recover private member begin_
	int getEndPosition(); ///Accessor to recover private member end_

	/* functions to get different domains */
	vector<int> getTimeDomainWindow(int startingAt); ///returns the time domain values from a position		
	vector<double> getFrequencyDomainWindow(int startingAt); ///returns the frequency domain values from a position
	vector<double> getPeakDomainWindow(int startingAt); ///returns the peak domain values from a position

    private:			
	//members
	int begin_; ///the start location of the domain values to calculate
	int end_; ///the end location of the domain values to calculate
	int windowSize_; ///the size of each window over which to calculate domain values
	string fileName_; /// the filename of the pcm file on which to operate

	FileVector<int> pcmVec_; ///The fileVector of values from the pcm file
	ConfigFile settingsFile_; ///The configFile for the settings file

	vector<double> calcFrequencyDomainWindow(vector<int> pcmWindow); ///Calculates the freq domain over a pcm window
	vector<double> calcPeakDomainWindow(vector<double> windowVec); ///Calculates the Peak domain over a freq window
	vector<double> calcFFT(float pcmVals[]); ///Calculates the FFT values over an array of pcmVals
	vector<double> applySpectrumFunc(vector<double> fftVec); ///Applies the appropriate Spectrum function for freq domain
	vector<double> applyWindowFunc(vector<double> windowVec); ///Applies the appropriate window function for freq domain

	/* different fft spectrum functions */
	vector<double> applyLogMagnitude(vector<double> fftVec); ///Calculates the Log Magnitude spectrum over fft values
	vector<double> applyPower(vector<double> fftVec); ///Calculates the Power spectrum over fft values
	vector<double> applyLogPower(vector<double> fftVec); ///Calculates the Log Power spectrum over fft values
	vector<double> applyMagnitude(vector<double> fftVec); ///Calculates the Magnitude spectrum over fft values

	/* windowing functions */
	vector<double> applyRectangular(vector<double> windowVec); ///Calculates the Rectangular Window function over a window
	vector<double> applyHamming(vector<double> windowVec);///Calculates the Hamming Window function over a window
	vector<double> applyHann(vector<double> windowVec);///Calculates the Hann Window function over a window
	vector<double> applyBartlett(vector<double> windowVec);///Calculates the Bartlett Window function over a window
	vector<double> applyTriangular(vector<double> windowVec);///Calculates the Triangular Window function over a window
	vector<double> applyBartlett_hann(vector<double> windowVec);///Calculates the Bartlett-Hann Window function over a window
	vector<double> applyBlackman(vector<double> windowVec);///Calculates the Blackman Window function over a window
	vector<double> applyBlackman_harris(vector<double> windowVec);///Calculates the Blackman-Harris Window function over a window

};
#endif
