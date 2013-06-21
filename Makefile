# definitions
CC = g++
CCFLAGS = -Wall

## targets and dependencies
all : featExtract

featExtract : main.o src/segmenter.o src/featureExtract.o src/configFile.o src/domain.o src/feature.o
	$(CC) $(CCFLAGS) $^ -o $@ -lfftw3f -ldl

# default rule for compiling .cc to .o
%.o: %.cc                                ## next line must begin with a TAB
	$(CC) -c $(CCFLAGS) $< 



domain.o : src/fileVector.h src/configFile.h src/feature.h
feature.o : src/fileVector.h src/configFile.h src/domain.h
featureExtract.o : src/fileVector.h src/configFile.h src/feature.h src/domain.h
segmenter.o : src/segmenter.h src/featureExtract.h
main.o : src/segmenter.h

clean:                 ## next lines must begin with a TAB
	rm -f src/*.o
	rm -f *.o
	rm -f *~ *# .#*

clean-all : clean      ## next line must begin with a TAB
	rm -f featExtract

