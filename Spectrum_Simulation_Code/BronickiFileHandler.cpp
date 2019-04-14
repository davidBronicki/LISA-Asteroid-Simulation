#include <fstream>
#include <vector>
#include <memory>
#include <string>
#include <exception>

namespace bfh{
	template<typename T>
	class DataHandler;
	class DataHandler<T>{
		vector<T> dataList;
		char typeName;
	public:
		DataHandler(T inputData, char inputName){
			dataList = inputData;
			typeName = inputName;
			if (inputName != 'c' && inputName != 'f' && inputName != 's'){
				throw
			}
		}
	}
	class FileHandler{

	};
}
