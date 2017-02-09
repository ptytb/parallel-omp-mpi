#include <mpi.h>
#include <iostream>
#include <algorithm>
#include <math.h>

#define MEASURE_BEGIN startTime = MPI_Wtime();
#define MEASURE_END(M) \
    deltaTime = MPI_Wtime() - startTime; \
    log(#M << ' ' << result << \
        " за время "<<  deltaTime << std::endl); 

#define READ_USER(W) \
std::cout << "Введите параметр " << #W << "=";\
std::cin >> W;

#define LIMIT_VALUE_LEFT 1.
#define LIMIT_VALUE_RIGHT 0.

#define REGION_LEFT 0.
#define REGION_RIGHT 1.

#define ROOT 0
#define TAIL nodeCount - 1

#define TAG_VALUE 0

#define CONVERGENCE_CHECK_PERIOD 200
#define CONVERGENCE_EPSILON 1E-3

#define log(M) \
{\
    int _nodeId;\
    MPI_Comm_rank(MPI_COMM_WORLD, &_nodeId);\
    std::cout << '#' << _nodeId << " " <<  M;\
    std::cout.flush();\
}

inline double step(double &t, double &a, double &hh, double &left,
                   double &mid, double &right);
inline void sendDouble(double &value, int id);
inline double recvDouble(int id);
//inline void sendInt(int &value, int id);
//inline int recvInt(int id);
void dump(double what[], int size, int n);
inline bool converges(double a[], double b[], int size);

int main(int argc, char *argv[])
{
    int nodeId, nodeCount;
    double startTime, deltaTime, result;
    double a; // Теплопроводность
    int N, arraySize, timeLimit;
    double t, h, hh;
    double *uPrev, *uNext, *buffer;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nodeCount);
    MPI_Comm_rank(MPI_COMM_WORLD, &nodeId);

    bool isHead = (nodeId == ROOT);
    bool isTail = (nodeId == TAIL);

    if (isHead)
    {
        READ_USER(a);
        READ_USER(N);
        READ_USER(timeLimit);
    }

    MPI_Bcast(&a, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 
    MPI_Bcast(&timeLimit, 1, MPI_INT, ROOT, MPI_COMM_WORLD); 

    arraySize = N;

    std::cout << std::fixed;
    std::cout.precision(2);

    log("Начало вычислений" << std::endl);
    MPI_Barrier(MPI_COMM_WORLD);

    a *= a;
    h = (REGION_RIGHT - REGION_LEFT) / N;

    int tail = (N - (N / nodeCount) * nodeCount);
    N /= nodeCount;
    N += (isTail ? tail : 0);
    buffer = new double[2 * N];

    //log("N= " << N << " nc=" << nodeCount);

    hh = h * h;

    // Из условия устойчивости
    t = 0.5 * hh / a;

    uPrev = buffer;
    uNext = buffer + N;

    // Краевые условия первого порядка

    if (isHead)
    {
        uNext[0] = uPrev[0] = LIMIT_VALUE_LEFT;
    }

    if (isTail)
    {
        uNext[N - 1] = uPrev[N - 1] = LIMIT_VALUE_RIGHT;
    }

    memset(uPrev + isHead, '\0', (N - 1 - isTail) * sizeof(double));
    //dump(uPrev, N, 0);

    MEASURE_BEGIN;

    int j;

    for (j = 0;
         j < timeLimit && !converges(uPrev, uNext, N);
         ++j, std::swap(uPrev, uNext))
    {
        for (int i = 1; i < N - 1; ++i)
        {
            uNext[i] = step(t, a, hh, uPrev[i - 1], uPrev[i], uPrev[i + 1]);
        }

        double leftEdge, rightEdge;

        if (nodeId % 2 == 0)
        {
            if (!isTail)
            {
                sendDouble(uPrev[N - 1], nodeId + 1);
                rightEdge = recvDouble(nodeId + 1);
            }

            if (!isHead)
            {
                leftEdge = recvDouble(nodeId - 1);
                sendDouble(uPrev[0], nodeId - 1);
            }
        }
        else
        {
            if (!isHead)
            {
                leftEdge = recvDouble(nodeId - 1);
                sendDouble(uPrev[0], nodeId - 1);
            }

            if (!isTail)
            {
                sendDouble(uPrev[N - 1], nodeId + 1);
                rightEdge = recvDouble(nodeId + 1);
            }
        }

        if (!isHead)
        {
            uNext[0] = step(t, a, hh, leftEdge, uPrev[0], uPrev[1]);
        }

        if (!isTail)
        {
            uNext[N - 1] = step(t, a, hh, uPrev[N - 2], uPrev[N - 1], rightEdge);
        }

        //dump(uNext, N, j);
    }

    MEASURE_END(Конец вычислений);

    // Результат в uPrev
    
    //dump(uPrev, N, timeLimit);

    /**************************************************/
    // Собираем фрагменты на 0
    /**************************************************/
    
    double *resultArray;

    if (isHead)
    {
        resultArray = new double[arraySize];

        for (int i = 0; i < arraySize; ++i)
        {
            resultArray[i] = 666.;
        }
    }

    const int slice = arraySize / nodeCount;

    MPI_Gather(uPrev, slice, MPI_DOUBLE,
               resultArray, slice, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    if (isTail && tail)
    {
        MPI_Send(uPrev + slice, tail, MPI_DOUBLE, ROOT, TAG_VALUE, MPI_COMM_WORLD);
    }

    if (isHead)
    {
        if (tail)
        {
            MPI_Recv(resultArray + arraySize - tail, tail, MPI_DOUBLE,
                     TAIL, TAG_VALUE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        dump(resultArray, arraySize, j);
        log("Результат: t=" << t * j);
        delete resultArray;
    }

    delete buffer;

    MPI_Finalize();
    return 0;
}

inline double step(double &t, double &a, double &hh, double &left,
                   double &mid, double &right)
{
    return mid + t * a * (left - 2. * mid + right) / hh;
}

void sendDouble(double &value, int id)
{
    //log("send..." << std::endl);
    MPI_Send(&value, 1, MPI_DOUBLE, id, TAG_VALUE, MPI_COMM_WORLD);
    //log("sent " << value <<  " to " << id << std::endl);
}

double recvDouble(int id)
{
    double value;

    //log("recv..." << std::endl);
    MPI_Recv(&value, 1, MPI_DOUBLE, id, TAG_VALUE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //log("recv " << value << " from " << id << std::endl);
    return value;
}

//void sendInt(int &value, int id)
//{
    //MPI_Send(&value, 1, MPI_INT, id, TAG_VALUE, MPI_COMM_WORLD);
//}

//int recvInt(int id)
//{
    //MPI_Status stat;
    //int value;

    //MPI_Recv(&value, 1, MPI_INT, id, TAG_VALUE, MPI_COMM_WORLD, &stat);
    //return value;
//}

void dump(double what[], int size, int n)
{
    log("Итерация " << n << ": ");
    
    for (int i = 0; i < size; ++i)
    {
        std::cout << what[i] << ' ';
    }

    std::cout << std::endl;
}

bool converges(double a[], double b[], int size)
{
    static unsigned int n = 0;
    n++;

    if (n % CONVERGENCE_CHECK_PERIOD == 0)
    {
        double sum = 0.;

        for (int i = 0; i < size; ++i)
        { 
            double d = a[i] - b[i];
            sum += d * d;
        }

        return sqrt(sum) <= CONVERGENCE_EPSILON;
    }

    return false;
}

