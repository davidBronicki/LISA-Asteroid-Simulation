#include <vector>
#include <string>

#define INT 0
#define FLOAT 1
#define STRING 2
#define BOOL 3

//data size: if -1 will change to whatever length is passed
//if 0, must pass 0 arguments and must be bool
//other must pass <dataSize> aruments
namespace arguments{
	class argumentsHandler{
		std::vector<bool> dataCalled;
		std::vector<void*> data;
		std::vector<int> dataTypes;
		std::vector<int> dataSizes;
		std::vector<std::vector<std::string>> flags;
		std::string helpString;
		bool helpCalled;
		int getFlagIndex(std::string flag);
		bool setFlagValue(std::string flag, std::vector<std::string> rawDatum);
	public:
		argumentsHandler();
		argumentsHandler& addFlag(std::vector<std::string> flagNames, int dataType, int dataSize);
		argumentsHandler& setHelp(std::string inputHelpString);
		void* getFlagValue(std::string);
		void passArgs(int argc, char** argv);
		bool isFlagCalled(std::string flag);
		bool isHelpCalled();
	};
}