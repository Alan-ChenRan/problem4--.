#include<stdio.h>
#include<iostream>
#include<mpi.h>
#include <stdlib.h>
#include <time.h> 
#include<algorithm>
#include<math.h>
#include<assert.h>
#include<malloc.h>
#include<sys/time.h>
#define For(i,j,n) for(int i=j;i<n;++i)
using namespace std;
int ProcNum, ProcRank;
int num_add, num_Multiply, num_exchange;
int Sum_add, Sum_Multiply, Sum_exchange;

int flag=0,Na=0,Nb=0,N=0,Type=0;
int pos(int x,int y){return y+N*x;}
template<class T>
T get_num(FILE* fd,T a)
{
    char t[1];
    T sum=0;a=0;
    fread(t,1,1,fd);
    while(t[0]<'0'||t[0]>'9')fread(t,1,1,fd);
    while(t[0]>='0'&&t[0]<='9')sum=sum*10+t[0]-'0',fread(t,1,1,fd);
    
    return sum;
}
template<class T>
void Reada(T*a,FILE *fa)
{ 
   T k;k=0;
   For(i,0,N)For(j,0,N)a[pos(i,j)]=rand()%100;//get_num(fa,k);
}

template<class T>
void Readb(T*a,FILE *fa)
{ 
   T k;k=0;
   For(i,0,N)a[i]=rand()%100;//get_num(fa,k);
}
char buf[20];
int L=0;
template<class T>
void Writebuf(T x)
{
    L=0;
    if(x==0){buf[0]='0';buf[1]='\0';L=1;return;}
    while(x>0)
    {
        buf[L]=x%10+'0';x/=10;++L;
    }
    For(i,0,L/2){char t=buf[i];buf[i]=buf[L-1-i];buf[L-1-i]=t;}
    buf[L]=' ';++L;
    buf[L]='\0';
}
template<class T>
void Printing(T*c,FILE*fc)
{
    For(i,0,N)
    {
        For(j,0,N)
        {
            Writebuf(c[pos(i,j)]);
            fwrite(buf,1,L,fc);
        }
        fwrite("\n",1,1,fc);
    }
}

template<class T>
void ParallelResultCalculation(T* pProcRows, T* pVector, T*pProcResult, int Size, int RowNum)
{
	int i, j; // Loop variables 
	for (i = 0; i < RowNum; i++) {
		pProcResult[i] = 0;
		for (j = 0; j < Size; j++)
			pProcResult[i] += pProcRows[i * Size + j] * pVector[j], ++num_add, ++num_Multiply;
	}

}

template<class T>
void ResultReplication(T* pProcResult, T* pResult, int Size,int RowNum)
{
	int i; // Loop variable 
	int* pReceiveNum; // Number of elements, that current process sends 
	int* pReceiveInd; /* Index of the first element from current process
	in result vector */
	int RestRows = Size; // Number of rows, that haven’t been distributed yet 
	//Alloc memory for temporary objects 
	pReceiveNum = new int[ProcNum];
	pReceiveInd = new int[ProcNum];
	//Define the disposition of the result vector block of current processor 
	pReceiveInd[0] = 0;
	pReceiveNum[0] = Size / ProcNum;
	for (i = 1; i < ProcNum; i++) {
		RestRows -= pReceiveNum[i - 1];
		pReceiveNum[i] = RestRows / (ProcNum - i);
		pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
	}
	//Gather the whole result vector on every processor 
        if(flag)MPI_Allgatherv(pProcResult, pReceiveNum[ProcRank], MPI_LONG_LONG_INT , pResult,pReceiveNum, pReceiveInd, MPI_LONG_LONG_INT , MPI_COMM_WORLD);
        else MPI_Allgatherv(pProcResult, pReceiveNum[ProcRank], MPI_INT , pResult,pReceiveNum, pReceiveInd, MPI_INT , MPI_COMM_WORLD);
	//Free the memory 
	delete[] pReceiveNum;
	delete[] pReceiveInd;
}

template<class T>
void DataDistribution(T* pMatrix, T* pProcRows, T* pVector,int Size, int RowNum)
{
	int* pSendNum; // the number of elements sent to the process 
	int* pSendInd; // the index of the first data element sent to the process 
	int RestRows = Size; // Number of rows, that haven’t been distributed yet 
	if(flag)MPI_Bcast(pVector, Size, MPI_LONG_LONG_INT , 0, MPI_COMM_WORLD);
        else MPI_Bcast(pVector, Size, MPI_INT , 0, MPI_COMM_WORLD);
	// Alloc memory for temporary objects 
	pSendInd = new int[ProcNum];
	pSendNum = new int[ProcNum];
	// Define the disposition of the matrix rows for current process 
	RowNum = (Size / ProcNum);
	pSendNum[0] = RowNum * Size;
	pSendInd[0] = 0;
	for (int i = 1; i < ProcNum; i++) {
		RestRows -= RowNum;
		RowNum = RestRows / (ProcNum - i);
		pSendNum[i] = RowNum * Size;
		pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
	}
	// Scatter the rows 
	if(flag)MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_LONG_LONG_INT , pProcRows,pSendNum[ProcRank], MPI_LONG_LONG_INT , 0, MPI_COMM_WORLD);
        else MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_INT , pProcRows,pSendNum[ProcRank], MPI_INT , 0, MPI_COMM_WORLD);
	// Free the memory 
	delete[] pSendNum;
	delete[] pSendInd;
}

template<class T>
void RandomDataInitialization(T *pMatrix, T *pVector, int Size)
{
	for (int i = 0; i < Size; ++i)pVector[i] = rand()%10;
	for (int i = 0; i < Size* Size; ++i)pMatrix[i] = rand()%10;
}

template<class T>
void ProcessInitialization(T*& pProcRows, T*& pProcResult,int& Size, int& RowNum)
{
	int RestRows; // Number of rows, that haven’t been distributed yet 
	int i; // Loop variable 
	MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	RestRows = Size;
	for (i = 0; i < ProcRank; i++)RestRows = RestRows - RestRows / (ProcNum - i);
	RowNum = RestRows / (ProcNum - ProcRank);
	pProcRows = new T[RowNum * Size];
	pProcResult = new T[RowNum];
}

template<class T>
void ProcessTermination(T* pMatrix, T* pVector, T* pResult, T* pProcRows, T* pProcResult,int Size)
{
	if (ProcRank == 0)
	{
		//cout <<"0---------" << endl;
		cout <<"pMatrix：--------" << endl;
		for (int i = 0; i < Size; ++i)
		{
			for (int j = 0; j < Size; ++j)cout << pMatrix[i * Size + j]<<" "; cout << endl;
		}
		cout << "pVector：--------" << endl;
		for (int i = 0; i < Size; ++i)cout << pVector[i]<<" ";; cout << endl;
		cout << "pResult：--------" << endl;
		for (int i = 0; i < Size; ++i)cout << pResult[i]<<" ";; cout << endl;
	}
}
int main(int argc, char* argv[]) 
{
    srand((unsigned)time(NULL));
    int Size; // Sizes of initial matrix and vector 
    int RowNum;
	
    FILE *fa=fopen(argv[1],"r");
    if(fa==NULL){printf("file a error\n");return 1;}
    FILE *fb=fopen(argv[2],"r");
    if(fb==NULL){printf("file b error\n");return 1;}
    FILE *fc=fopen(argv[3],"w");
    if(fc==NULL){printf("file c error\n");return 1;}

    char type_a[1],type_b[1];
    fread(type_a,1,1,fa);
    fread(type_b,1,1,fb);
    Na=get_num(fa,Na);
    Nb=get_num(fb,Nb);
    N=max(Na,Nb);Size=N;
    if(type_a[0]=='l'&&type_b[0]=='l')flag=1;
 
    flag=0;
    //flag=1;
    N=1000;Size=N;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    
    if(flag==0)
    {
        int*pMatrix,*pVector,*pResult,*pProcRows,*pProcResult;
        pMatrix=(int*)malloc(N*N*sizeof(int));
        pVector=(int*)malloc(N*sizeof(int));
        pResult=(int*)malloc(N*N*sizeof(int));
        For(i,0,N)For(j,0,N)pResult[pos(i,j)]=0;

        if(ProcRank==0)Reada(pMatrix,fa),Readb(pVector,fb);

    	ProcessInitialization(pProcRows, pProcResult,Size, RowNum);
    	DataDistribution(pMatrix, pProcRows, pVector, Size, RowNum);
    	ParallelResultCalculation(pProcRows, pVector, pProcResult, Size, RowNum);
    	ResultReplication(pProcResult, pResult, Size, RowNum);
    	ProcessTermination(pMatrix, pVector, pResult, pProcRows, pProcResult, Size);
    	if(ProcRank==0)
        {
           Printing(pResult,fc);
           fflush(fc);fclose(fa);fclose(fb);fclose(fc);
        }
        MPI_Finalize();
    }
    else
    {
        long long*pMatrix,*pVector,*pResult,*pProcRows,*pProcResult;
        pMatrix=(long long*)malloc(N*N*sizeof(long long));
        pVector=(long long*)malloc(N*sizeof(long long));
        pResult=(long long*)malloc(N*N*sizeof(long long));
        For(i,0,N)For(j,0,N)pResult[pos(i,j)]=0;

        if(ProcRank==0)Reada(pMatrix,fa),Readb(pVector,fb);

    	ProcessInitialization(pProcRows, pProcResult,Size, RowNum);
    	DataDistribution(pMatrix, pProcRows, pVector, Size, RowNum);
    	ParallelResultCalculation(pProcRows, pVector, pProcResult, Size, RowNum);

    	ResultReplication(pProcResult, pResult, Size, RowNum);

    	ProcessTermination(pMatrix, pVector, pResult, pProcRows, pProcResult, Size);

    	if(ProcRank==0)
        {
           Printing(pResult,fc);
           fflush(fc);fclose(fa);fclose(fb);fclose(fc);
        }
        MPI_Finalize();
    }
    return 0;
}
