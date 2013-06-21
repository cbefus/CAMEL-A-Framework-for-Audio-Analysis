/**
*	@file segmenter.h
*	@class Segmenter
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
*	@brief Segments and extracts features from a PCM file
*
*	
*	Given a PCM file, a segmentation method index, feature method index, and a 
*	number of segments to extract this class segments the file according to the
*	method decided on and extracts the resulting features from each segment.	
*
*	Required Dependencies
*	-----------
*	std::vector
*	std::string
*	std::algorithm
*	std::math.h
*	configFile.h
*	segmenter.h
*
*	Example Usage
*	-----------
*	#include "segmenter.h"
*	#include <vector>
*	using namespace std;
*
*	int main(){	 
*
*		Segmenter s("myPCMfile.txt");
*		s.setup();
*		s.run(3, 2, 10);
*		vector< vector<double> > vvdFeatResults = s.getFeatureResults();
*		vector<int> vviPosResults = s.getPositionResults();	
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

#ifndef __SEGMENTER_H
#define __SEGMENTER_H

#include <string>
#include <vector>


#include "configFile.h"

#include "featureExtract.h"

using namespace std;

class Segmenter{
	public:

		//constructor
		Segmenter(string strPCMFileName = "null.txt"); ///constructor - inits private members

		//methods
		void setup(); ///inits dependencies classes
		void run(int intfeatureType, int intSegmentType = 0, int intNumSegments=20); ///runs specified actions

		//accessors
		vector< vector<double> > getFeatureResults(); ///returns the feature values for each segment
		vector<int> getPositionResults(); ///returns the position values for each segment
		void setFileName(string strFileName); ///sets the filename for the pcm file

	private:
		
		//members
		ConfigFile settingsFile_; ///the ConfigFile object holding settings information

		string strPCMFileName_; ///the file name of the PCM file object

		int intNumSamples_; /// the number of PCM samples in the file
		int intBitsPerChannel_; /// the bits per channel
		int intChannels_; /// the number of channels
		int intSampleRate_; /// the sampleRate
		bool blnNormalized_; /// whether or not the samples are normalised
	
		FeatureExtract fexExtractor_; ///FeatureExtract Object hooked into PCM file

		vector<int> vecPositionResults_; ///results vector for segment positions
		vector< vector<double> > vecFeatureResults_; ///results vector for segment feature values

		//member functions
		/* segmentation methods */
		void staticBySizeSegmentation(int intFeatureType); ///segment into groups of a given size
		void staticByNumSegmentation(int intNumSegments, int intFeatureType); ///segment into a given number of segments
		void topkFarthestNeighborSegmentation(int intNumSegments, int intFeatureType);///return segments based on algo
		void topkLargestMatrixScoreSegmentation(int intNumSegments, int intFeatureType);///return segments based on algo
		void greedyMergeSegmentation(int intNumSegments, int intFeatureType);///return segments based on algo

		/* segmentation helpers */
		void segmentGivenPositionsVector(vector<int> segmentPositions, int intFeatureType); ///perform segmentation given position vector
		double euclideanDist(vector<double> featureA, vector<double> featureB); ///calculate the euclidean distance between two feature vectors
		vector< vector<double> > multiplyVectors(vector< vector<double> > distVec, vector< vector<int> > kernelVec); /// multiply a distance vector against a kernel vector
		vector<int> getTopkWithHueristic(vector<double> valueVector, int intNumSegments, int intNumHueristicFrames, int intFrameSize); ///given a series of segment distances return the topK that are not within a hueristic distance of eachother
		void setPositionResultsVec(vector<int> vecPositionResults);///set the results vector for position
		void setFeatureResultsVec(vector< vector<double> > vecFeatureResults);///set the results vector for value
		void readCheckPCMHeader(); ///read the PCM file header and grab values

};
#endif
