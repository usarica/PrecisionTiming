{
	gSystem->Load("libFWCoreFWLite.so"); 
	AutoLibraryLoader::enable();
	gSystem->Load("libDataFormatsFWLite.so");
	gSystem->Load("libDataFormatsPatCandidates.so");
	
#include "DataFormats/FWLite/interface/Handle.h"
}
