/**
*	@file segmenter.cpp
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
*	@todo replace namespace with localised std
*	
*
**/
#include "segmenter.h"

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <math.h>
#include "domain.h"

using namespace std;

/** 
*	Segmenter
*	@brief Constructor method sets the PCM filename
*	@param strPCMFileName a string representing the PCM file object
*	@post strPCMFileName_ is set to parameter
*	@pre PCM file actually exists
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
Segmenter::Segmenter(string strPCMFileName){
	strPCMFileName_ = strPCMFileName;
}

/** 
*	setup
*	@brief sets up dependency objects for the class
*	@param 
*	@post fexExtractor_ is setup, settingsFile_ is setup
*	@pre strPCMFileName_ is set
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void Segmenter::setup(){
	settingsFile_.setFileName("settings.txt");
	readCheckPCMHeader();
	fexExtractor_.setFileName(strPCMFileName_);
	fexExtractor_.setWindowSize(settingsFile_.read<int>("WINDOW_SIZE", 1024));
	fexExtractor_.setup();
}

/** 
*	run
*	@brief Switch to the appropriate segmentation method funtion
*	@param intFeatureType an int representing what feature (from featureExtract) to use
*	@param intSegmentType an int representing what Segmentation algo to use
*	@param intNumSegments an int representing the number of segments desired
*	@post
*	@pre
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void Segmenter::run(int intFeatureType, int intSegmentType, int intNumSegments){

	switch(intSegmentType){
		case 1:{ //do static
			staticBySizeSegmentation(intFeatureType);
			break;
		}
		case 2:{
			staticByNumSegmentation(intNumSegments, intFeatureType);
			break;
		}
		case 3:{ //do top-k farthest neighbors
			topkFarthestNeighborSegmentation(intNumSegments, intFeatureType);
			break;
		}
		case 4:{ //do top-k biggest matrix scores
			topkLargestMatrixScoreSegmentation(intNumSegments, intFeatureType);
			break;
		}
		case 5:{ //do greedy merge based
			greedyMergeSegmentation(intNumSegments, intFeatureType);			
			break;
		}
		default:{
			break;
		}
	}

}

/** 
*	staticBySizeSegmentation
*	@brief static by size segmentation just separates the values into groups of even length (l)
*	@param intFeatureType an int representing the feature type for extraction
*	@post
*	@pre settingsFile_ is setup, intNumSamples_ is set
*	@return
*	@par Method Just discover how many segments would be created by splitting the PCM file into
*	groups of length l.  Then call staticByNumSegmentation with this value.
*	@note size l is set in the settings file
*	@remarks
*	@warning
*	@todo
**/
void Segmenter::staticBySizeSegmentation(int intFeatureType){
	
	vector< vector<double> > vecCurrFeatureResults;
	int intSegmentSize = settingsFile_.read<int>("SEGMENT_SIZE", 102400);
	int intNumSegments = static_cast<int>(ceil((intNumSamples_)/static_cast<double>(intSegmentSize)));
	
	staticByNumSegmentation(intNumSegments, intFeatureType);

}


/** 
*	staticByNumSegmentation
*	@brief static by num segmentation just separates the values into k groups of even length
*	@param intNumSegments an int representing the number of resulting segments (k)
*	@param intFeatureType an int representing feature type for extraction
*	@post
*	@pre
*	@return
*	@par Method Calculate the length of each of the segments given the PCM file size.  For
*	each segment simply calculate the start and end positions using the length of segments
*	and the current segment iterator.  From this start/end position calculate features and store.
*	@note	
*	@remarks
*	@warning
*	@todo
**/

void Segmenter::staticByNumSegmentation(int intNumSegments, int intFeatureType){

	vector<int> vecPositionResults;
	vecPositionResults.push_back(0);
	vector< vector<double> > vecFeatureResults;
	int intSegmentSize = static_cast<int>(ceil((intNumSamples_)/static_cast<double>(intNumSegments)));
	
	int intSegStartPos = 0;
	int intSegEndPos = 0;
	
	for(int i = 0; i < intNumSegments; i++){
		intSegStartPos = i * intSegmentSize;
		intSegEndPos = intSegStartPos + intSegmentSize;
		
		if(intSegEndPos > intNumSamples_){
			intSegEndPos = intNumSamples_;
		}
		
		fexExtractor_.getFeature(intSegStartPos, intSegEndPos, intFeatureType);
		vecPositionResults.push_back(intSegEndPos);
		vecFeatureResults.push_back(fexExtractor_.getValues());
		
	}

	setPositionResultsVec(vecPositionResults);
	setFeatureResultsVec(vecFeatureResults);

}


/** 
*	topkFarthestNeighborSegmentation
*	@brief Find the top-K farthest neighbors in terms of a given feature and segment there
*	@param intNumSegments an int representing the requested number of segments output (K)
*	@param intFeatureType an int representing the feature type for evaluation/extraction
*	@post
*	@pre
*	@return
*	@par Method top-k farthest neighbors segmentation measures the euclidean distance between 
*	each frame (size f) and picks the top-k as the segmentation locations such that they are 
*	at least b distance apart (simple hueristic)
*	@note f and b are set in the settings file
*	@remarks loosely based on several papers which mimic George Tzanetakis of 1999
*	@warning
*	@todo
**/
void Segmenter::topkFarthestNeighborSegmentation(int intNumSegments, int intFeatureType){

	//first frame and get feature vals
	int intFrameSize = settingsFile_.read<int>("TKFN_FRAME_SIZE", 1024);
	int intNumFrames = static_cast<int>(ceil((intNumSamples_)/static_cast<double>(intFrameSize)));
	
	staticByNumSegmentation(intNumFrames, intFeatureType);
	vector< vector<double> > vvdCurrResults = getFeatureResults();
	
	//then calc neighboring distances
	vector<double> vdblDistVals;
	for(int i = 0; i < vvdCurrResults.size()-1; ++i){
		vdblDistVals.push_back(euclideanDist(vvdCurrResults[i], vvdCurrResults[i+1]));
		//cout << vdblDistVals[i] << endl;
	}
	
	//get top-k max but not any within hueristic percentage of eachother
	int intNumHueristicFrames = static_cast<int>(intNumFrames * settingsFile_.read<double>("TKFN_HUERISTIC_PERCENT", 0.05));

	vector<int> vintSegmentLocations = getTopkWithHueristic(vdblDistVals, intNumSegments, intNumHueristicFrames, intFrameSize);

	segmentGivenPositionsVector(vintSegmentLocations, intFeatureType);

}


/** 
*	topkLargestMatrixScoreSegmentation
*	@brief find top-K largest matrix scores and segment there
*	@param intNumSegments an int representing the number of required output segments (K)
*	@param intFeatureType an int representing the feature type for extraction
*	@post
*	@pre
*	@return
*	@par Method top-k largest matrix score segmentation makes a self similarity and 
*	cross similarity matrix of each frame in some spectrum (S) and then measures the 
*	difference between them as a score. the top-k largest of these are selected as 
*	segmentation locations (with hueristics to prevent small segments using b). Note 
*	that the extracted feature from this segmentation is not the same as the spectrum 
*	used for the matricies. Also note the spectrum type and window function for S 
*	are both set under feature extraction settings in the settings file.  Furthermore, 
*	this algorithm directly accesses the domain class
*	@note S and b are set in the settings file
*	@remarks
*	@warning
*	@todo
**/

void Segmenter::topkLargestMatrixScoreSegmentation(int intNumSegments, int intFeatureType){

	//first frame and get domain vectors -- note for this algo frames = windows
	int intFrameSize = settingsFile_.read<int>("WINDOW_SIZE", 1024);
	int intNumFrames = static_cast<int>(ceil((intNumSamples_)/static_cast<double>(intFrameSize)));
	int intKernelSize = settingsFile_.read<int>("TKLM_KERNEL_SIZE", 44);	

	if(intKernelSize % 2 != 0){
		intKernelSize -= 1;
	}

	Domain domExtractor(strPCMFileName_ , intFrameSize);
	domExtractor.setup();
	domExtractor.setEndPosition(intNumSamples_);
	domExtractor.setStartPosition(0);


	/**** precalculate the domain values for all frames ****/
	vector< vector<double> > freqDomWindows;
	for(int i = 0; i < intNumFrames; i++){
		vector<double> domWind = domExtractor.getFrequencyDomainWindow(i*intFrameSize);
		freqDomWindows.push_back(domWind);
	}


	//create a entire multidimensional distance vector as a series of 2 dimensional vectors with repetition 
	//(where we center kernel only moves by 1 frame at a time) 

	vector< vector< vector<double> > > distanceVec;
	for(int i  = intKernelSize-1; i < intNumFrames; ++i){
		vector< vector<double> > tempVecOut;
		for(int j = i; j > i-intKernelSize; --j){
			vector<double> tempVecIn;
			for(int k = i-intKernelSize+1; k < i+1; ++k){				
				//vector<double> firstValVec = domExtractor.getFrequencyDomainWindow(k*intFrameSize);
				//vector<double> secondValVec = domExtractor.getFrequencyDomainWindow(j*intFrameSize);
				tempVecIn.push_back( euclideanDist(freqDomWindows[k], freqDomWindows[j]));
			}
			tempVecOut.push_back(tempVecIn);
		} 
		distanceVec.push_back(tempVecOut);
	}

	//construct entire multidimensional checkerboard kernel vector as 2 dimensional  vector 
	vector< vector<int> > kernelVec;	
	//first half of kernal is square of 1's followed by -1's
	for(int i = 0; i < intKernelSize/2; ++i){
		vector<int> tempVecA;
		for(int j = 0; j < intKernelSize/2; ++j){
			tempVecA.push_back(1);
		}
		kernelVec.push_back(tempVecA);
		vector<int> tempVecB;
		for(int j = 0; j < intKernelSize/2; ++j){
			tempVecB.push_back(-1);
		}
		kernelVec.push_back(tempVecB);
	}
	//second half of kernal is square of -1's followed by 1's
	for(int i = intKernelSize/2; i < intKernelSize; ++i){
		vector<int> tempVecA;
		for(int j = 0; j < intKernelSize/2; ++j){
			tempVecA.push_back(-1);
		}
		kernelVec.push_back(tempVecA);
		vector<int> tempVecB;
		for(int j = 0; j < intKernelSize/2; ++j){
			tempVecB.push_back(1);
		}
		kernelVec.push_back(tempVecB);
	}

	//multiply kernel (using matrix multiplication) values onto distance values
	//note that if all distances = 0 then we get a score of 0... big distances = big or little score
	for(int i = 0; i < distanceVec.size(); ++i){
		distanceVec[i] = multiplyVectors(distanceVec[i], kernelVec);
	}
	
	//sum each kernel^2 set of mixed distance values to create novelty score
	vector<double> noveltyVec;
	//first pad with 0's for values skipped by centering
	for(int i = 0; i < intKernelSize; ++i){
		noveltyVec.push_back(0.0);
	}

	for(int i = 0; i < distanceVec.size(); ++i){
		double noveltyScore = 0.0;
		for(int j = 0; j < distanceVec[i].size(); ++j){
			for(int k = 0; k < distanceVec[i][j].size(); ++k){
				noveltyScore += distanceVec[i][j][k];
			}
		}
		noveltyVec.push_back(noveltyScore);
	}

	//find top-k extrema of these scores (considering hueristic) and mark as segmentation locations
	int intNumHueristicFrames = static_cast<int>(intNumFrames * settingsFile_.read<double>("TKLM_HUERISTIC_PERCENT", 0.05));

	vector<int> vintSegmentLocations = getTopkWithHueristic(noveltyVec, intNumSegments, intNumHueristicFrames, intFrameSize);

	segmentGivenPositionsVector(vintSegmentLocations, intFeatureType);
}


/** 
*	greedyMergeSegmentation
*	@brief Merge similar frames until only K are left unmerged... these are the segments
*	@param intNumSegments an int representing the number of required output segments
*	@param intFeatureType an int representing the feature used for evaluation/extraction
*	@post
*	@pre
*	@return
*	@par Method greedy merge based segmentation just merges frames which are most similar 
*	until only k frames are left
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void Segmenter::greedyMergeSegmentation(int intNumSegments, int intFeatureType){

	//first frame and get feature vals
	int intFrameSize = settingsFile_.read<int>("GM_FRAME_SIZE", 1024);
	int intNumFrames = static_cast<int>(ceil((intNumSamples_)/static_cast<double>(intFrameSize)));
	
	staticByNumSegmentation(intNumFrames, intFeatureType);
	vector< vector<double> > vvdCurrResults = getFeatureResults();
	vector<int> vviCurrPositionResults = getPositionResults();
	
	//then calc neighboring distances
	vector<double> vdblDistVals;
	for(int i = 0; i < vvdCurrResults.size()-1; ++i){
		vdblDistVals.push_back(euclideanDist(vvdCurrResults[i], vvdCurrResults[i+1]));
	}

	//now find min distance and merge out -- rextracting new value and calculating new distances
	//do this until only k segments left (or k-1 distances)

	double dblMaxDist = 0.0;
	for(int i = 0; i < vdblDistVals.size(); ++i){
		if( vdblDistVals[i] > dblMaxDist ){
			dblMaxDist = vdblDistVals[i];
		}
	}
	
	int intCurrNumDistances = vdblDistVals.size();
	while(intCurrNumDistances > intNumSegments-1){

		double dblMinDist = dblMaxDist; 
		int intMinPos = 0;
		for(int i = 0; i < vdblDistVals.size(); ++i){
			if(vdblDistVals[i] < dblMinDist){
				dblMinDist = vdblDistVals[i];
				intMinPos = i;
			}
		}

		
		//get new segment feature value
		fexExtractor_.getFeature(vviCurrPositionResults[intMinPos], vviCurrPositionResults[intMinPos+2], intFeatureType);

		vvdCurrResults[intMinPos+1] = fexExtractor_.getValues(); //set new value in right segment value area
		vvdCurrResults.erase(vvdCurrResults.begin()+intMinPos); //remove left value area
		
		//calc new left and right distance -- if needed
		if(intMinPos-1 > -1){
			vdblDistVals[intMinPos-1] = euclideanDist(vvdCurrResults[intMinPos-1], vvdCurrResults[intMinPos]); 
		}
		if(intMinPos+1 < intCurrNumDistances){		
			vdblDistVals[intMinPos+1] = euclideanDist(vvdCurrResults[intMinPos], vvdCurrResults[intMinPos+1]);
		}

		vdblDistVals.erase(vdblDistVals.begin()+intMinPos); //remove extra distance
		vviCurrPositionResults.erase(vviCurrPositionResults.begin()+intMinPos+1); //remove extra position
				
		intCurrNumDistances--;
	}

	vviCurrPositionResults.erase(vviCurrPositionResults.begin());
	vviCurrPositionResults.erase(vviCurrPositionResults.end()-1);
	segmentGivenPositionsVector(vviCurrPositionResults, intFeatureType);


}

/** 
*	segmentGivenPositionsVector
*	@brief Method given a vector of pcm positions (not including start and end) and a feature type stores the results vecs
*	@param vinSegmentPositions a vector of int representing the positions in the file for segmentation
*	@param intFeatureType an int representing the feature for extraction
*	@post
*	@pre
*	@return
*	@par Method Add 0 to the start and intNumSamples_ to the end of the segment positions. For each of the segments
*	extract the feature from that segment and add it to the values vector. Set the final results vectors.
*	@note	
*	@remarks
*	@warning
*	@todo clear results vectors before setting
**/
void Segmenter::segmentGivenPositionsVector(vector<int> vintSegmentPositions, int intFeatureType){
	
	vector<int> vintSegmentPositionsFull;
	vintSegmentPositionsFull.push_back(0);
	for(int i = 0; i < vintSegmentPositions.size(); ++i){
		vintSegmentPositionsFull.push_back(vintSegmentPositions[i]);
	}
	vintSegmentPositionsFull.push_back(intNumSamples_);

	vector< vector<double> > vecFeatureResults;	

	for(int i = 0; i < vintSegmentPositionsFull.size()-1; ++i){

		if(vintSegmentPositionsFull[i+1] > intNumSamples_){
			vintSegmentPositionsFull[i+1] = intNumSamples_;
		}

		fexExtractor_.getFeature(vintSegmentPositionsFull[i], vintSegmentPositionsFull[i+1], intFeatureType);
		vecFeatureResults.push_back(fexExtractor_.getValues());
	}

	setPositionResultsVec(vintSegmentPositionsFull);
	setFeatureResultsVec(vecFeatureResults);
}


/** 
*	getTopkWithHueristic
*	@brief Given a vector of distances calculate the top-k positions taking account for the hueristic, and sorted
*	@param valueVector a vector of double with the distances for all frames
*	@param intNumSegments an int representing the number of segments (K)
*	@param intNumHueristicFrames the number of frames between any selection (min)
*	@param intFrameSize the size of each frame
*	@post
*	@pre
*	@return a vector of ints representing the location of the top-k segment locations
*	@par Method Find the largest location, remove any locations within the hueristic. repeat. return locations.
*	@note	
*	@remarks
*	@warning
*	@todo
**/

vector<int> Segmenter::getTopkWithHueristic(vector<double> valueVector, int intNumSegments, int intNumHueristicFrames, int intFrameSize){
	
	if(intNumSegments * (2*intNumHueristicFrames) + intNumSegments > valueVector.size()){
		cout << "Err:\\> Segmenter::topkFarthestNeighborSegmentation --> not enough frames for parameters" << endl;
		cout << "    \\>   " <<  intNumSegments * (2*intNumHueristicFrames) + intNumSegments << " > " << valueVector.size() << endl;
	}

	vector<int> vintSegmentLocations;
	vector<int> vintDontLookIndicies;
	for(int i = 0; i < intNumSegments-1; ++i){ //selecting 5 locations gives 6 segments...
		double dblMaxDist = -1.0;
		int intMaxIndex = -1;
		for(int j = 0; j < valueVector.size(); ++j){
			bool blnSkip = false;
			for(int k = 0; k < vintDontLookIndicies.size(); ++k){
				if(vintDontLookIndicies[k] == j){
					blnSkip = true;
					break;
				}
			}

			if(!blnSkip && valueVector[j] > dblMaxDist){
				dblMaxDist = valueVector[j];
				intMaxIndex = j;
			}
		}

		vintSegmentLocations.push_back(intMaxIndex);

		int intDontLookStart = intMaxIndex-intNumHueristicFrames;
		int intDontLookEnd = intDontLookStart + (2 * intNumHueristicFrames) + 1;
		for(int j = intDontLookStart; j < intDontLookEnd; ++j){		
			vintDontLookIndicies.push_back(j);
		}
	}

	sort(vintSegmentLocations.begin(), vintSegmentLocations.end());

	for(int i = 0; i < intNumSegments-1; ++i){
		vintSegmentLocations[i] *= intFrameSize;
	}

	return vintSegmentLocations;
}


/** 
*	multiplyVectors
*	@brief multiplies a distance vector against a kernel vector
*	@param distVec is a 2-dim vector of doubles representing distances
*	@param kernelVec is a 2-dim vector of ints representing a kernel
*	@post
*	@pre assumes both vectors are in same dimensions and that the vectors are square
*	@return a 2-dim vector of double representing the outcome of matrix multiplication
*	@par Method iterate through the distance vec applying matrix multiplication
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector< vector<double> > Segmenter::multiplyVectors(vector< vector<double> > distVec, vector< vector<int> > kernelVec){

	vector< vector<double> > returnVec;

	for(int i = 0; i < distVec.size(); ++i){
		for(int j = 0; j < distVec.size(); ++j){
			vector<double> tempVec;
			double sumVals = 0.0;
			for( int k = 0; k < distVec.size(); ++k){
				sumVals += distVec[i][k]*kernelVec[k][j];
				tempVec.push_back(sumVals);
			}
			returnVec.push_back(tempVec);
		}
	}

	return returnVec;
}


/** 
*	euclideanDist
*	@brief calculates the euclidean distance between two vectors
*	@param featureA a vector of double representing the values of a feature
*	@param featureB a vector of double representing the values of a feature
*	@post
*	@pre featureA and featureB are of ths same length
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning if param vectors are not the same length we warn and return a dist of 0.0
*	@todo
**/
double Segmenter::euclideanDist(vector<double> featureA, vector<double> featureB){
	double dblEuclidDist = 0.0;	
	if(featureA.size() != featureB.size()){
		cout << "Err:\\> Segmenter:: euclideanDist --> unmatched feature vectors" << endl;
		return dblEuclidDist;
	}

	
	for( int i = 0; i < featureA.size(); ++i){
		dblEuclidDist += pow(featureA[i] - featureB[i], 2);
	}

	return sqrt(dblEuclidDist);

}

/** 
*	getPositionResults
*	@brief returns the position result vector
*	@param
*	@post
*	@pre
*	@return a vector of ints representing segmentation locations
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector<int> Segmenter::getPositionResults(){
	
	return vecPositionResults_;

}

/** 
*	getFeatureResults
*	@brief returns the feature results vector
*	@param name descr
*	@post
*	@pre
*	@return a 2-dim vector of double representing features for each segment
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
vector< vector<double> > Segmenter::getFeatureResults(){
	
	return vecFeatureResults_;

}

/** 
*	setFileName
*	@brief sets the PCM file name
*	@param strFileName a string representing the file object
*	@post
*	@pre
*	@return
*	@par Method just sets the member
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void Segmenter::setFileName(string strFileName){
	strPCMFileName_ = strFileName;
}

/** 
*	setPositionResultsVec
*	@brief sets the position results vector
*	@param vecPositionResults a vector of ints representing segmentation positions
*	@post
*	@pre
*	@return
*	@par Method clear the results vector and then set it
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void Segmenter::setPositionResultsVec(vector<int> vecPositionResults){
	vecPositionResults_.clear();
	vecPositionResults_ = vecPositionResults;
}

/** 
*	setFeatureResultsVec
*	@brief sets the feature results vector
*	@param vecFeatureResults a 2-dim vector of feature values for each segment
*	@post
*	@pre
*	@return
*	@par Method clear the vector and set with new values
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void Segmenter::setFeatureResultsVec(vector< vector<double> > vecFeatureResults){
	vecFeatureResults_.clear();
	vecFeatureResults_ = vecFeatureResults;
}

/** 
*	readCheckPCMHeader
*	@brief Checks that the PCM header is correct and grabs values
*	@param
*	@post intNumSamples_, intBitsPerChannel_, intChannels_, intSampleRate_, blnNormalized_ are set
*	@pre strPCMFileName_ is set.
*	@return
*	@par Method 
*	@note	
*	@remarks
*	@warning
*	@todo
**/
void Segmenter::readCheckPCMHeader(){

	//cout << "        Reading and Checking PCM file Header Properties." << endl;
	std::ifstream fin(strPCMFileName_.c_str());
	string junk = "";
	fin >> junk >> intNumSamples_;
	fin >> junk >> intBitsPerChannel_;
	if (intBitsPerChannel_ != 16) {
		cout << "Wrn:\\> Segmenter::readCheckPCMHeader-->: BitsPerChannel not 16: "<< intBitsPerChannel_ << endl;
	}
	fin >> junk >> intChannels_;
	if (intChannels_ != 1) {
		cout << "Wrn:\\> Segmenter::readCheckPCMHeader-->: intChannels not 1:  "<< intChannels_ << endl;
	}
	fin >> junk >> intSampleRate_;
	if (intSampleRate_ != 44100) {
		cout << "Wrn:\\> Segmenter::readCheckPCMHeader-->: intSampleRate not 44100: " << intSampleRate_ << endl;
	}
	fin >> junk >> junk;
	blnNormalized_ = (junk == "TRUE")? true : false;
	if (blnNormalized_ != false) {
		cout << "Wrn:\\> Segmenter::readCheckPCMHeader-->: blnNormalized not false: "<< junk << endl;
	}
	fin.close();
}
