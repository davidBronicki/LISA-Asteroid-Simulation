#include "cuda_runtime.h"

#include <vector>

namespace cudaUtil{
	template <class T>
	class hdmem
	{
		unsigned int* references;
		T* hostMemT;
		T* deviceMemT;
		int numElements;
		int dataSize;
		cudaError_t cudaStatus;
		bool* deviceMemLoaded;
	public:
		hdmem();
		hdmem(int numberOfEntries);
		hdmem(T* instantiationData, int numberOfEntries);
		hdmem(T instantiationData, int blank);
		hdmem(std::vector<T> instantiationData);
		hdmem(const hdmem<T> &other);
		~hdmem();
		hdmem<T>& operator=(hdmem<T> other);
		hdmem<T>& loadMemory();
		hdmem<T>& unloadMemory();
		hdmem<T>& updateHost();
		hdmem<T>& updateDevice();
		T get(int index);
		std::vector<T> getAll();
		hdmem<T>& set(T value, int index);
		hdmem<T>& setAll(T* refValue);
		hdmem<T>& setAll(std::vector<T> refValue);
		operator T* ();
		int size();
		bool isGood();
		int deviceMemoryAllocated();
	};

	template <class T>
	hdmem<T>::hdmem()
	{
		references = new unsigned int(1);
		deviceMemLoaded = new bool(false);
		hostMemT = new T[1];
	}

	template <class T>
	hdmem<T>::hdmem(int numberOfEntries)
	{
		references = new unsigned int(1);
		numElements = numberOfEntries;
		dataSize = numberOfEntries * sizeof(T);
		hostMemT = new T[numberOfEntries];
		deviceMemLoaded = new bool(false);
	}

	template <class T>
	hdmem<T>::hdmem(T* instantiationData, int numberOfEntries)
	{
		references = new unsigned int(1);
		numElements = numberOfEntries;
		dataSize = numberOfEntries * sizeof(T);
		hostMemT = new T[numberOfEntries];
		for (int i = 0; i < numberOfEntries; i++)
		{
			hostMemT[i] = instantiationData[i];
		}
		deviceMemLoaded = new bool(false);
	}

	template <class T>
	hdmem<T>::hdmem(T instantiationData, int numberOfEntries)
	{
		references = new unsigned int(1);
		numElements = 1;
		dataSize = sizeof(T);
		hostMemT = new T[1];
		*hostMemT = instantiationData;
		deviceMemLoaded = new bool(false);
	}

	template <class T>
	hdmem<T>::hdmem(std::vector<T> instantiationData)
	{
		references = new unsigned int(1);
		numElements = instantiationData.size();
		dataSize = numElements * sizeof(T);
		hostMemT = new T[numElements];
		for (int i = 0; i < numElements; i++)
		{
			hostMemT[i] = instantiationData[i];
		}
		deviceMemLoaded = new bool(false);
	}

	template <class T>
	hdmem<T>::hdmem(const hdmem<T> &other)
	{
		dataSize = other.dataSize;
		numElements = other.numElements;
		deviceMemLoaded = other.deviceMemLoaded;
		cudaStatus = other.cudaStatus;
		deviceMemT = other.deviceMemT;
		hostMemT = other.hostMemT;
		references = other.references;
		(*references)++;
	}

	template <class T>
	hdmem<T>::~hdmem()
	{
		(*references)--;
		if (*references == 0)
		{
			delete[] hostMemT;
			if (deviceMemLoaded) cudaFree(deviceMemT);
			delete references;
			delete deviceMemLoaded;
		}
	}

	template <class T>
	hdmem<T>& hdmem<T>::operator=(hdmem<T> other)
	{
		(*this).~hdmem();
		dataSize = other.dataSize;
		numElements = other.numElements;
		deviceMemLoaded = other.deviceMemLoaded;
		cudaStatus = other.cudaStatus;
		deviceMemT = other.deviceMemT;
		hostMemT = other.hostMemT;
		references = other.references;
		(*references)++;
		return *this;
	}

	template <class T>
	hdmem<T>& hdmem<T>::loadMemory()
	{
		if (!*deviceMemLoaded)
		{
			cudaStatus = cudaMalloc((void**)&(deviceMemT), dataSize);
			if (cudaStatus != cudaSuccess)
			{
				fprintf(stderr, "Cuda Memory Allocation Failed!\n");
				return *this;
			}
			*deviceMemLoaded = true;
		}
		return *this;
	}

	template <class T>
	hdmem<T>& hdmem<T>::unloadMemory()
	{
		if (*deviceMemLoaded)
		{
			cudaFree(deviceMemT);
			*deviceMemLoaded = false;
		}
		return *this;
	}

	template <class T>
	hdmem<T>& hdmem<T>::updateDevice()
	{
		cudaStatus = cudaMemcpy(deviceMemT, hostMemT, dataSize, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Cuda Host to Device Memory Copy Failed!\n");
			return *this;
		}
		return *this;
	}

	template <class T>
	hdmem<T>& hdmem<T>::updateHost()
	{
		cudaStatus = cudaMemcpy(hostMemT, deviceMemT, dataSize, cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "Cuda Device to Host Memory Copy Failed!\n");
			return *this;
		}
		return *this;
	}

	template <class T>
	T hdmem<T>::get(int index)
	{
		return hostMemT[index];
	}

	template <class T>
	std::vector<T> hdmem<T>::getAll()
	{
		std::vector<T> output = std::vector<T>();
		for (int i = 0; i < numElements; i++)
		{
			output.push_back(hostMemT[i]);
		}
		return output;
	}

	template <class T>
	hdmem<T>& hdmem<T>::set(T value, int index)
	{
		if (index < numElements)
		{
			hostMemT[index] = value;
			return *this;
		}
		else return *this;
	}

	template <class T>
	hdmem<T>& hdmem<T>::setAll(T* refValue)
	{
		for (int i = 0; i < numElements; i++)
		{
			try
			{
				hostMemT[i] = refValue[i];
			}
			catch (...)
			{
				return *this;
			}
		}
		if (*deviceMemLoaded) return updateDevice();
		else return *this;
	}

	template <class T>
	hdmem<T>& hdmem<T>::setAll(std::vector<T> refValue)
	{
		if (numElements == refValue.size())
		{
			for (int i = 0; i < numElements; i++)
			{
				hostMemT[i] = refValue[i];
			}
			if (*deviceMemLoaded) return updateDevice();
			else return *this;
		}
		else return *this;
	}

	template <class T>
	hdmem<T>::operator T* ()
	{
		return deviceMemT;
	}

	template <class T>
	int hdmem<T>::size()
	{
		return numElements;
	}

	template <class T>
	bool hdmem<T>::isGood()
	{
		return cudaStatus == cudaSuccess;
	}

	template <class T>
	int hdmem<T>::deviceMemoryAllocated()
	{
		if (*deviceMemLoaded) return (long long)dataSize;
		else return 0;
	}
}