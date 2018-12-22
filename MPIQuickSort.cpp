#include "pch.h"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;


void lessThan(int *array, int &arraySize, int *arrayNew, int &arrayNewSize, int *arraySend, int &arraySendSize, int pivot) {

	for (int t = 0; t < arraySize; ++t) {
		if (array[t] <= pivot) {
			arrayNew[arrayNewSize] = array[t];
			++arrayNewSize;
		}
		else {
			arraySend[arraySendSize] = array[t];
			++arraySendSize;
		}
	}
}


void greaterThan(int *array, int &arraySize, int *arrayNew, int &arrayNewSize, int *arraySend, int &arraySendSize, int pivot) {

	for (int t = 0; t < arraySize; ++t) {
		if (array[t] > pivot) {
			arrayNew[arrayNewSize] = array[t];
			++arrayNewSize;
		}
		else {
			arraySend[arraySendSize] = array[t];
			++arraySendSize;
		}
	}
}


void QuickSort(int *array, int arraySize, int rootProc) {
	const int idMsgPivot = 1;
	const int idMsgArray = 2;
	const int idMsgArraySize = 3;

	int procSize, procRank, arrayNum, size, procArraySize, procArrayNewSize, procArraySendSize, procArrayRecvSize;
	int *procArray, *procArrayNew, *procArraySend;

	MPI_Comm_size(MPI_COMM_WORLD, &procSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

	if (procRank == 0) {
		size = arraySize;
		MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
		arrayNum = size / procSize;
	}
    // сообщение размера массива всем процессам
	MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	arrayNum = size / procSize;

	int diff = size - arrayNum * procSize;

	int *scounts = new int[procSize];
	for (int i = 0; i < procSize - 1; ++i)
		scounts[i] = arrayNum;
	scounts[procSize - 1] = arrayNum + diff;

	int *displs = new int[procSize];
	displs[0] = 0;
	for (int i = 1; i < procSize; ++i) {
		displs[i] = displs[i - 1] + scounts[i - 1];
	}
	procArraySize = scounts[procRank];
	procArray = new int[procArraySize];
	MPI_Scatterv(array, scounts, displs, MPI_INT, procArray, scounts[procRank], MPI_INT, 0, MPI_COMM_WORLD);
	delete[] scounts;
	delete[] displs;

	// количество итераций деления гиперкуба, фактически это логарифм по основанию 2 от количества процессов
	int iterationsNum = 0;
	while ((1 << iterationsNum) != procSize && (1 << iterationsNum) < procSize) {
		++iterationsNum;
	}

	int currentBit, stepMainOfCube, pairProcess, pivot;
	for (int i = iterationsNum - 1; i >= 0; --i) {
		currentBit = 1 << i;
		stepMainOfCube = 1 << (i + 1); // расстояние между главными процессами гиперкубов
		bool isMainOfCube = false;

		// определение парного процесса (процесс парный - если разница в одном бите)
		pairProcess = procRank | currentBit;
		if (procRank == pairProcess)
			pairProcess -= currentBit;

		// поиск главных процессов каждого гиперкуба
		for (int j = 0; j < procSize; j += stepMainOfCube) {
			if (procRank == j) {
				isMainOfCube = true;
				// подсчет опорного элемента
				if (procArraySize != 0) {
					pivot = procArray[procArraySize / 2];
				}
				else {
					pivot = 0;
				}
					
				// сообщение опорного элемента другим процессам текущего гиперкуба
				for (int t = procRank + 1; t < procRank + stepMainOfCube; ++t) {
					MPI_Send(&pivot, 1, MPI_INT, t, idMsgPivot, MPI_COMM_WORLD);
				}
				break;
			}
		}

		// получение опорного элемента от главного гиперкуба
		if (!isMainOfCube) {
			MPI_Status status;
			MPI_Recv(&pivot, 1, MPI_INT, MPI_ANY_SOURCE, idMsgPivot, MPI_COMM_WORLD, &status);
		}

		procArrayNew = new int[procArraySize];
		procArraySend = new int[procArraySize];
		procArrayNewSize = 0;
		procArraySendSize = 0;

		// обмен блоками данных между парными процессами
		// отправка данных парному процессу с помощью MPI_Send
		// получение данных от парного процесса с помощью MPI_Recv
		if (procRank < pairProcess) {
			lessThan(procArray, procArraySize, procArrayNew, procArrayNewSize, procArraySend, procArraySendSize, pivot);
			delete[] procArray;

			MPI_Send(&procArraySendSize, 1, MPI_INT, pairProcess, idMsgArraySize, MPI_COMM_WORLD);
			MPI_Send(procArraySend, procArraySendSize, MPI_INT, pairProcess, idMsgArray, MPI_COMM_WORLD);
			delete[] procArraySend;

			MPI_Status status;
			MPI_Recv(&procArrayRecvSize, 1, MPI_INT, pairProcess, idMsgArraySize, MPI_COMM_WORLD, &status);

			procArraySize = procArrayNewSize + procArrayRecvSize;
			procArray = new int[procArraySize];
			memcpy(procArray, procArrayNew, procArrayNewSize * sizeof(int));
			delete[] procArrayNew;

			MPI_Recv(procArray + procArrayNewSize, procArrayRecvSize, MPI_INT, pairProcess, idMsgArray, MPI_COMM_WORLD, &status);
		}
		else {
			greaterThan(procArray, procArraySize, procArrayNew, procArrayNewSize, procArraySend, procArraySendSize, pivot);
			delete[] procArray;

			MPI_Status status;
			MPI_Recv(&procArrayRecvSize, 1, MPI_INT, pairProcess, idMsgArraySize, MPI_COMM_WORLD, &status);
			MPI_Send(&procArraySendSize, 1, MPI_INT, pairProcess, idMsgArraySize, MPI_COMM_WORLD);

			procArraySize = procArrayNewSize + procArrayRecvSize;
			procArray = new int[procArraySize];
			memcpy(procArray, procArrayNew, procArrayNewSize * sizeof(int));
			delete[] procArrayNew;

			MPI_Recv(procArray + procArrayNewSize, procArrayRecvSize, MPI_INT, pairProcess, idMsgArray, MPI_COMM_WORLD, &status);

			MPI_Send(procArraySend, procArraySendSize, MPI_INT, pairProcess, idMsgArray, MPI_COMM_WORLD);
			delete[] procArraySend;
		}
	}

	// сортировка части массива
	sort(procArray, procArray + procArraySize);

	int *rcounts = new int[procSize];
	displs = new int[procSize];
	// сбор размеров блоков данных всех процессов
	MPI_Gather(&procArraySize, 1, MPI_INT, rcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
	displs[0] = 0;
	for (int i = 1; i < procSize; ++i) {
		displs[i] = displs[i - 1] + rcounts[i - 1];
	}
	// сбор блоков данных всех процессов (MPI_Gatherv, т.к. умеет собирать блоки разной длины)
	MPI_Gatherv(procArray, procArraySize, MPI_INT, array, rcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

	delete[] procArray, rcounts, displs;
}


int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	int procRank, procSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procSize);

	int *array;
	int arraySize = 1000000;
	array = new int[arraySize];
	double startTime, finishTime;

	if (procRank == 0) {
		ifstream input("D:\input1000000.txt");
		for (int i = 0; i < arraySize; ++i) {
			input >> array[i];
		}
		startTime = MPI_Wtime();
	}
        
	QuickSort(array, arraySize, 0);

	if (procRank == 0) {
		finishTime = MPI_Wtime();

		ofstream out("D:\work_time.txt");
		out << finishTime - startTime << endl;
		out.close();

		ofstream output("D:\output1000000.txt");
		for (int i = 0; i < arraySize; ++i) {
			output << array[i] << " ";
		}
		output.close();
	}

	delete[] array;
	MPI_Finalize();

	return 0;
}

