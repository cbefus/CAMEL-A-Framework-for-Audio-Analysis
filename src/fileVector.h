/**
*	@file fileVector.h
*	@class FileVector
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
*	@brief reads in a file of values and puts them into a vector for easy access
*	
*
*	Required Dependencies
*	-----------
*	std::vector
*	std::string
*	std::fstream
*	configFile.h
*
*	Example Usage
*	-----------
*	#include "fileVector.h"
*	#include <vector>
*	using namespace std;
*
*	int main(){
*		int vecStart = 44100;
*		int vecEnd = 88200;
*		FileVector<int> p;
*		p.setFileName("myPCMfile.txt");
*		p.setup();
*		for(int i = vecStart; i < p.size() && i < vecEnd; ++i){
*			cout << pcmVec[i] << endl;
*		}	
*	}
*
*
*	More examples to be found in example folder
*
*	@todo
*	
*
**/
#ifndef FILEVECTOR_H
#define FILEVECTOR_H

#include <fstream>
#include <string>
#include <vector>

#include "configFile.h"

//precondition of file is that it is of only one type (T) exclusive of the header separated by whitespace or newline

template <class T>
class FileVector{

	public:
		FileVector(std::string fileName); ///Constructor - intializes the private member fileName_
		void setup(); ///copies over the content of the file into private vector dataVector_
		T operator[](int position) {return getValueAt(position);}; ///assigns the operator [] to function getValueAt
		int size() { return dataVector_.size(); }; ///returns the size of the dataVector
		void setFileName(std::string fileName) {fileName_ = fileName;}; ///public access to private member... not needed

	private:
		std::string fileName_; ///stores the string name of the file to be copied to vector format
		std::vector<T> dataVector_; ///the vector to hold the file contents
		ConfigFile settingsFile_; ///the configFile object for getting file header size

		T getValueAt(int position) { return dataVector_[position]; }; ///returns the value of the vector at a position
		bool fileExists(); ///returns true or false to the existance of the private member file
};

//public
/** 
*	FileVector
*	@brief Constructor for FileVector Object
*	@param fileName the file name as a string to attach the vector to
*	@post FileVector object exists for access, though perhaps attached to a non existant file
*	@pre The requested file exists
*	@par Method Simply sets the private member fileName_ to the input parameter
*	@note 
*	@remarks This probably doesn't need to be a template class since we restrict the file format elsewhere
*	@warning
*	@todo add a check for file existance
**/
template <class T>
FileVector<T>::FileVector(std::string fileName = "null.txt"){
	fileName_ = fileName;
}


/** 
*	setup
*	@brief Copies all but the file header into a private vector object (dataVector_)
*	
*	@post dataVector_ contains the content of the file (minus the header)
*	@pre The requested file exists and has more content than the header lines
*	@par Method 1) Open the input file, 2) Read in number of header lines from settings file, 
*	3) grab the header lines and do nothing with them, 4) grab all following values and push onto vector 
*	@note 
*	@remarks This probably doesn't need to be a template class since we restrict the file format elsewhere
*	@warning
*	@todo see if vector size changes -> implies success
**/
template <class T>
void FileVector<T>::setup(){
	
	if(fileExists()){
		std::ifstream fin(fileName_.c_str());
		std::string headerValue;	

		int headerSize = settingsFile_.read<int>("PCM_HEADER_LINES", 5);
		for(int i = 0; i < headerSize; ++i){
			std::getline(fin, headerValue);
		}

		T fileValue;
		while(fin >> fileValue){
			dataVector_.push_back(fileValue);
		}
	
		fin.close();
	}

}


/** 
*	fileExists
*	@brief Checks to see if the private member file exists
*	
*	
*	@pre the private member fileName_ has been initialised
*	@return boolean indicating if the file exists (true) or not (false)
*	@par Method Create a file stream to the object, if we can open it (close stream) and return true.
*	Otherwise close the stream and return false.
*	@note	
*	@remarks
*	@warning
*	@todo Actually use this above
**/

//return true if a file exists given a name
template <class T>
bool FileVector<T>::fileExists(){
	
	std::ifstream fin(fileName_.c_str());
	if( fin.is_open()) {	
		fin.close();
		return true;
	}
	fin.close();
	return false;
}





#endif
