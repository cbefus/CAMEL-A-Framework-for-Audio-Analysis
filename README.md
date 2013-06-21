

Welcome to Content-based Audio and Music Extraction Library.  The main purpose of this library is to give simple access to the common MIR extraction functions without a massive learning curve or the restrictions of GUIs. For more information on the reasoning and design considerations behind CAMEL see C. Sanden, C. Befus, and J. Zhang, CAMEL: A lightweight Framework for Content-based Audio and Music Analysis, In Proceedings of Audio Mostly, Piteå, Sweden, 2010. 

This readme will outline the following

0. License

1. Using CAMEL
1.1 Building the library and Dependencies
1.2 Creating the proper input file format (PCM ASCII text)
1.3 Extracting a Domain from a file using CAMEL
1.4 List of Domains in CAMEL (with identifiers for switch/settings)
1.5 Extracting a Feature from a file using CAMEL
1.6 List of Features in CAMEL (with identifiers for switch/settings)
1.7 Extracting Segmentations from a file using CAMEL
1.8 List of Segmentations in CAMEL (with identifiers for switch/settings)

2. Extending CAMEL
2.1 Adding a Domain to CAMEL
2.2 Adding a Feature to CAMEL
2.3 Adding a Segmentation to CAMEL


0. License
-----------

Copyright (c) 2010 Chad R. Befus, Chris Sanden, Cody Rioux
	
CAMEL (Content-based Audio and Music Extraction Library) is unrestricted, open 
source "software", with respect, but not limited, to modification, use, 
publishing, and distribution, subject to the following conditions:

1. This copyright notice is maintained in its current form throughout CAMEL.

2. Any publication or distribution credits the usage of CAMEL appropriately.
	
This software is provided "as is" without warranty. The authors are in no way
liable for any misuse or damages arising from the use, modification or 
distribution of it or any part thereof.

-----------

1. Using CAMEL
###########
-----------

Every file in CAMEL begins with a header which includes License information, 
file dependencies, and an example usage of the file.  Every function in CAMEL 
begins with commenting which outlines the parameters, conditions, and methods 
used throughout the function.  Reading these comments can help with the use of 
CAMEL.  Many of the comments for features and domains has the math behind the 
methodology writting in latex form.  All the code in CAMEL has been evaluated 
against the output of several other mainstream extraction libraries and is, 
to that end, correct (though several errors were discoved in other libraries 
during this testing).

Note: this readme is written assuming the user is using a linux machine.  
However CAMEL has been tested and used successfully on both a windows and 
Mac machine.  The content of this readme is still applicable (in most cases) 
to any operating system.

-----------

1.1 Building the library and Dependencies
-----------
First off, at this point (version 1.0) the only dependency for building CAMEL
is that the FFTW library be installed on your machine. FFTW is freely available
and can be found at http://www.fftw.org/.  For the purposes of creating CAMEL 
libfftw3-dev was used and is recommended for future use.  On a linux terminal 
this can be installed using: 

	~$ sudo apt-get install libfftw3-dev

Building CAMEL in its current form is extremely simple. Assuming your main file is
named main.cpp and you are preforming a segmentation/extraction task the makefile 
provided should handle the full build process. simply typing:

	~$ make

should handle the process.

If you wish to skip the segmentation task for feature extraction then the lines in the makefile: 

	segmenter.o : src/segmenter.h src/featureExtract.h
	main.o : src/segmenter.h

can be replaced with the lines:

	main.o : src/featureExtract.h

and the line:

	featExtract : main.o src/segmenter.o src/featureExtract.o src/configFile.o src/domain.o src/feature.o

can be replaced with:

	featExtract : main.o src/featureExtract.o src/configFile.o src/domain.o src/feature.o

similarily if you wish to skip the featureExtract step (only use the domain class) you could 
do the same steps but removing the featureExtract.o content as well.

If you wish for your executable to be named something other than featExtract simply rename 
the word featExtract in the makefile with the one you prefer.

If you wish to use a file other than main.cpp as your main code file simply replace the word 
main.o in the makefile with the one you are using.

After building CAMEL use: 

	~$ make clean 

to remove the temporary files from the build process.  Also:

	~$ make clean-all

removes the executable along with the temporary files if you want to do a clean build.

-----------

1.2 Creating the proper input file format
-----------
The input files for CAMEL are in a particular format.  We are working on making a .mp3 decoder
however for the time being we use a PCM  ASCII text print out.  An example file is provided called 
samplePCM.txt.  The format of these files is a 5 line header explaining the content of the rest 
of the file:

SAMPLES:	221184
BITSPERSAMPLE:	16
CHANNELS:	1
SAMPLERATE:	44100
NORMALIZED:	FALSE

The samples line gives the number of lines in the file after the header (PCM samples). For CAMEL 
this line is extremely important as this value is used for much of the processing.  The next 
four lines are not actually used by CAMEL but are checked and warned about if not matching the 
above parameters.  CAMEL expects the audio to be in mono at 44100 without normalization.  However 
everything in CAMEL should still work if these things are not true (though some of the maths might
be off). If these are not the values you wish to use then it is HIGHLY RECOMMENDED that you adjust 
the settings in the settings.txt file accordingly.

Following the header the PCM format has a single integer value on each new line for the number of
samples listed in the header.

We converted .mp3 files to the format above using Adobe Audition 3.0.  Audition should be able to 
convert any format of file your data set is in into this format.  (Note we ran Audition on a linux 
system using Wine 1.3 -- so no complaining that it is windows only!!).

To use Adobe Audition to convert to our input format simply open Audition, and under the File menu 
select to batch processing 

	Under tab 1: add the data file you intend to convert. 
	under tab 3: check off the conversion settings and for destination format select 44100, mono, and 16 bit, 
	under tab 4: select ASCII text data (*.txt) as the output format 
	under tab 5: set up your destination ect. 

Finally hitting Run Batch should make the conversion correctly.  These files are now ready for processing in CAMEL


-----------

1.3 Extracting a Domain from a file using CAMEL
-----------
In this section we are assuming you have adjusted the makefile and preprocessed the data as 
explained above.  Here we give a general explanation and some example code on extracting some 
domain values from a PCM ASCII text file.

The simplest example of extracting a domian from CAMEL is to print to screen a single window
of spectral domain values.  To do this we first open the settings.txt file and adjust the 
appropriate settings for the spectral domain. This includes the settings under headings 
'#PCM Settings' and '#  FFT SPECTRUM'

Next we write a simple main.cpp file.  In it we create an empty domain object, we set the
file name for our PCM file in the domain object, we set the window size for our extraction
to 1024 arbitrarily.  We set the start position of the file to 0. 
We set the end position to the end of the file (or 10000 arbitrarily).
We then call the setup funtion for the domain object allowing it to organize itself.  
Finally we calculate the spectral domain requested (in the settings file) via the 
getFrequencyDomainWindow function, passing in the location we wish to start extraction.  
Then a simple for loop prints the contents of the resulting vector to standard out. Build
and run.


#include "src/domain.h"
#include <vector>
using namespace std;

int main(){	 
	Domain domain;
	domain.setFileName("samplePCM.txt");
	domain.setWindowSize(1024);
	domain_.setStartPosition(0);
	domain_.setEndPosition(10000);
	domain.setup();	
	vector<double> domainVec = domain.getFrequencyDomainWindow(0);
	for(int i = 0; i < 1024; ++i){
		cout << domainVec[i] << endl;
	}
	return 0;
}

-----------


1.4 List of Domains in CAMEL (with identifiers for switch/settings)
-----------

Domains			Function
1 Time			getTimeDomainWindow()
2 Spectral		getFrequencyDomainWindow()
3 Peak			getPeakDomainWindow()

Spectrums		Switch Value		Function
1 Log Magnitude		1			applyLogMagnitude()
2 Power			2			applyPower()
3 Log Power		3			applyLogPower()
4 Magnitude		4			applyMagnitude()
			

Window Functions	Switch Value		Function
1 Hann			1			applyHann()
2 Hamming		2			applyHamming()
3 Triangular		3			applyTriangular()
4 Bartlett		4			applyBartlett()
5 Bartlett Hann		5			applyBartlett_hann()
6 Blackman		6			applyBlackman()
7 Blackman Harris	7			applyBlackman_harris()
8 Rectangular		8			applyRectangular()
			

-----------

1.5 Extracting a Feature from a file using CAMEL
-----------
In this section we are assuming you have adjusted the makefile and preprocessed the data as 
explained above. Here we give a general explanation and some example code on extracting some 
feature values from a PCM ASCII text file.

The simplest example of extracting a feature from CAMEL is to extract the mean from a window
of the time domain.  To do this we first begin by making appropriate alterations in the settings 
file. For feature extraction settings manipulate, as needed, values under headings '#PCM Settings',
and '#Settings for featureExtract class' including '#  FFT SPECTRUM'

Next we write a simple main.cpp.  In our main file we create a feature extraction object.  We 
point it at our PCM ASCII text file. We set our window size to 1024. We call the setup function 
for the feature extraction object. We then use the getFeature function to calculate feature 1 (mean)
from between positions 0 and 7000.  we collect the results using the getValues function.  We can
then print these results to standard. Note that all results are returned as a vector, clearly for 
1-dimensional results such as the mean we only need to output the first position.

#include "src/featureExtract.h"
#include <vector>
using namespace std;

int main(){	 

	FeatureExtract fe;
	fe.setFileName("samplePCM.txt");
	fe.setWindowSize(1024);
	fe.setup();
	fe.getFeature(0, 7000, 1);
	vector<double> vecFeatureResults = fe.getValues();
	cout << vecFeatureResults[0] << endl;

}

-----------


1.6 List of Features in CAMEL (with identifiers for switch/settings)
-----------


Feature			Identifier	Function
1 Mean				1	calcWindowMean
2 Variance			2	calcWindowVariance
3 Standard Deviation		3	calcWindowStandardDeviation
4 Average Deviation		4	calcWindowAverageDeviation
5 Skewness 			5	calcWindowSkewness 
6 Kurtosis 			6	calcWindowKurtosis 
7 ZCR 				7	calcWindowZCR 
8 RMS 				9	calcWindowRMS 
9 Non-Zero Count 		10	calcWindowNonZeroCount 
10 Spectral Centroid 		11	calcWindowSpectralCentroid 
11 Spectral Variance 		12	calcWindowSpectralVariance 
12 Spectral Standard Deviation 	13	calcWindowSpectralStandardDeviation 
13 Spectral Average Deviation 	14	calcWindowSpectralAverageDeviation 
14 Spectral Skewness 		15	calcWindowSpectralSkewness 
15 Spectral Kurtosis 		16	calcWindowSpectralKurtosis 
16 Spectral Irregularity K 	17	calcWindowSpectralIrregularityK 
17 Spectral Irregularity J 	18	calcWindowSpectralIrregularityJ 
18 Spectral Flatness 		19	calcWindowSpectralFlatness 
19 Spectral Tonality 		20	calcWindowSpectralTonality 
20 Spectral Min 		21	calcWindowSpectralMin 
21 Spectral Max 		22	calcWindowSpectralMax 
22 Spectral Crest 		23	calcWindowSpectralCrest 
23 Spectral Slope 		24	calcWindowSpectralSlope 
24 Spectral Spread 		25	calcWindowSpectralSpread 
25 Spectral Rolloff 		26	calcWindowSpectralRolloff 
26 Spectral HPS 		27	calcWindowSpectralHPS 
27 Spectral Loudness 		28	calcWindowSpectralLoudness 
28 Spectral Sharpness 		29	calcWindowSpectralSharpness 
29 Peak Tristimulus 1 		31	calcWindowTristimulus1 
30 Peak Tristimulus 2 		32	calcWindowTristimulus2 
31 Peak Tristimulus 3 		33	calcWindowTristimulus3 
32 MFCC 			40	calcWindowMFCC 
33 Bark 			41	calcWindowBark


-----------



1.7 Extracting Segmentations from a file using CAMEL
-----------
In this section we are assuming you have adjusted the makefile and preprocessed the data as 
explained above. Here we give a general explanation and some example code on segmenting and 
extracting some domain values from a PCM ASCII text file.

The simplest example of segmentation is to create a static segmentation of the file and
extract a feature value of mean from each segment.

First, as always, we make changes to the appropriate locations in the settings.txt file for
segmentation.  In the case of segmentation this generally means all the settings are fair
game depending on your needs.

Next we create a simple main.cpp file.  In our main file we create a segmenter object and
initialise it with the PCM ASCII text file on which we are going to work with.  Next we
run the setup method for the segmentation object.  Then we call the run function on the 
segmentation object telling it to extract feature 1 (mean) from segmentation method 2 (static 
by number of segments) and to make 10 segments.

We then get the results (a vector of 10 vectors, each containing the mean value for a segment)
for the extraction process using the getFeatureResults method and we get the position result
using the getPositionResults method.


#include "src/segmenter.h"
#include <vector>
using namespace std;

int main(){	 

	Segmenter s("samplePCM.txt");
	s.setup();
	s.run(1, 2, 10);
	vector< vector<double> > vvdFeatResults = s.getFeatureResults();
	vector<int> vviPosResults = s.getPositionResults();	

}

-----------


1.8 List of Segmentations in CAMEL (with identifiers for switch/settings)
-----------

Segmentation		Identifier	Function
staticBySize		1		staticBySizeSegmentation
staticByNum		2		staticByNumSegmentation
topkFarthestNeighbor	3		topkFarthestNeighborSegmentation
topkLargestMatrixScore	4		topkLargestMatrixScoreSegmentation
greedyMerge		5		greedyMergeSegmentation	


-----------



2. Extending CAMEL
###########
-----------

CAMEL was also designed with simplified extension in mind.  It is relatively easy to 
add new features, domains, and segmentation algorithms to CAMEL.  If you wish your 
extention to be added to the main version of CAMEL please contact the authors.

-----------

2.1 Adding a Domain to CAMEL
-----------

There are of course several different things one might want to add into the domain class. This 
includes a new windowing function, a new spectrum function, or a new domain all together.  We will
discuss each of these here.

2.1.1 Adding a window function

Adding a new windowing function is extremely simple.  First create the function itself in the 
domain.cpp. The return type should be a vector of doubles. This function should be written to 
evaluate over a single window of data passed in  as a vector of doubles (the only parameter). 
Any settings or parameters for this function beyond the vector can be added to the settings 
file and be sourced via 

	mySetting = settingsFile_.read<mySettingType>("MYSETTINGNAME", mySettingDefaultValue);

Next add the function prototype to the domain.h header file.  Next add a new switch to the 
applyWindowFunc() function in the domain.cpp class which references your new function. Finally
use the setting SPECTRUM_WINDOW in the settings file to select your new method by its switch
value

2.1.2 Adding a spectrum function

Adding a new spectrum function is very similar to adding a new window function.  First create
a new function in the domain.cpp class for your spectrum where the only parameter is a vector
of the type you expect to calculate over (FFT values would be a vector of double).  As a return 
value set it to a vector of doubles. For any settings which are required beyond the parameter
vector use the settings file as explained abive under windowing functions.

Next add the function prototype to the domain.h header file.  Then add a new switch to the 
applySpectrumFunc() function in the domain.cpp class which references your new function. Finally
use the setting SPECTRUM_TYPE in the settings file to select your new method by its swith value.


2.1.3 Adding a new domain

Adding a new domain is also simple (assuming it is calculated of the time domain, spectral domain,
FFT, or the peak domain).

As with windows and spectrums assume that you are writing a function only concerned with calculating
your domain over a single vector of values (ints if its off the time domain, double for all else).
Your only parameter should be that vector. All other settings should be dealt with as above. Your 
function should return a vector of ints or doubles.

Next add the function prototype to the domain.h header file.  Source your function directly as seen
above in the examples on usage or add it to the featureExtract.cpp class for calculating/extracting
features over it (create a avgOverYOURDOMAINWindow function following the same methods as the others
in the featureExtract class and then add it to the featureExtract switch).


-----------

2.2 Adding a Feature to CAMEL
-----------

Adding a feature to CAMEL is extremely simple.  For the sake of simplicity we will assume you are
adding a feature over one of the domains already built into the system, otherwise simply see the
adding a new domain section first.

The first step is to add the statistical feature calculation to the feature.cpp file.  Just 
write the function to calculate the feature over a single vector of the correct domain values.  Your 
return type should be either a single value (cast it to a double if its not already) or an array of 
doubles (if it returns multi dimensionals like MFCC).  Your only parameter should be the vector 
of domain values.  Any other settings should use the settingsfile method outlined above in the 
adding a domain section.

Next add the prototype to the feature.h header file.  The add a switch to your feature in the 
featureExtract.cpp getFeature() function.  Note use the correct avgOverDomainWindow function and 
pass your feature to it as a pointer... see the other switch cases for examples.  You can now source
your feature directly as in the extracting a feature example above, or through segmentation.


-----------

2.3 Adding a Segmentation to CAMEL
-----------

Adding a segmentation algorithm to CAMEL simply involves creating a function in the 
segmenter.cpp file.  Note your parameters can fluctuate here but by standard they should be
two values: the number of segments, and the feature extraction method to process on. All other 
settings should be collected from the settings.txt system described above.

Once you have added your segmentation algorithms function, add its prototype to the segmentation.h 
header file and add a switch to it in the run() function of the segmentation class.  Most of the 
segmentation functions in segmenter.cpp make use of several of the other functions (such as the
segmentGivenPositionsVector() function which given a vector of integer positions returns fills the
results vectors for both positions and feature values).

You can now source your segmentation function as explained above.



-----------



