#include "Taylor.hpp"
#include "fobj.hpp"
//#include "stack.hpp"

#include <string>
#include <iostream>
#include <cstdlib>
#include <functional>
#include <cmath>
#include <omp.h>
#include <stack>
#include <vector>

typedef std::binder1st<std::mem_fun1_t<double, fobj_t, double> > fun_t;

void readUserF(std::string &userF);
Taylor::Computable *parseUserF(std::string &userF);
double integrate(double a, double b, double e, std::string &userF);
double integrate_parallel(double a, double b, double e, std::string &userF);

struct stack_item_t
{
    stack_item_t(double a, double b, double fa, double fb, double s)
    {
        this->a = a;
        this->b = b;
        this->fa = fa;
        this->fb = fb;
        this->s = s;
    }

    double a, b, fa, fb, s;
};

#define MEASURE_BEGIN startTime = omp_get_wtime();
#define MEASURE_END(M) \
    deltaTime = omp_get_wtime() - startTime; \
    std::cout << #M << ' ' << result << \
        " за время "<<  deltaTime << std::endl; \
        std::cout.flush();

int main(int argc, char *argv[])
{
    std::string userF;
    double startTime, deltaTime, singleTime, multiTime, result, LEFT, RIGHT,
           EPSILON;

    std::cout << "Введите левую границу, правую, точность: ";
    std::cin >> LEFT;
    std::cin >> RIGHT;
    std::cin >> EPSILON;

    readUserF(userF);

    std::cout << "Начало расчета" << std::endl;

    MEASURE_BEGIN
    result = integrate(LEFT, RIGHT, EPSILON, userF);
    MEASURE_END(Последовательный расчет)
    singleTime = deltaTime;

    MEASURE_BEGIN
    result = integrate_parallel(LEFT, RIGHT, EPSILON, userF);
    MEASURE_END(Параллельный расчет)
    multiTime = deltaTime;

    std::cout << "Конец расчета" << std::endl;

    std::cout << "Ускорение: single/multi = " << singleTime / multiTime <<
        std::endl;

    std::cout << "Эффективность, %: (single/multi)/p = " <<
        100. * (singleTime / multiTime) / omp_get_num_procs() << '%' <<
             std::endl;

    return EXIT_SUCCESS;
}

void readUserF(std::string &userF)
{
    std::cout << "Введите функцию: ";
    std::cin >> userF;
}

double integrate(double a, double b, double e, std::string &userF)
{
    fobj_t fobj(userF);
    fun_t f = std::bind1st(std::mem_fun(&fobj_t::y), &fobj);

    std::stack<stack_item_t> stk;

    double j = 0.;

    double fa = f(a);
    double fb = f(b);

    double sab = (fa + fb) * (b - a) / 2.;

    for ( ; ; )
    {
        double c = (a + b) / 2.;

        double fc = f(c);

        double sac = (fa + fc) * (c - a) / 2.;
        double scb = (fc + fb) * (b - c) / 2.;
        double sacb = sac + scb;
        double a_sacb = fabs(sacb);

        // Относительное отклонение
        if (fabs(sab - sacb) >= e * a_sacb && a_sacb >= e)
        {
            // Левую трапецию в стек
            stk.push(stack_item_t(a, c, fa, fc, sac));

            // Далее разбиваем правую
            a = c;
            fa = fc;
            sab = scb;
        }
        else
        {
            j += sacb;

            if (stk.empty())
            {
                break;
            }

            stack_item_t &i = stk.top(); 
            a = i.a; 
            b = i.b; 
            fa = i.fa; 
            fb = i.fb; 
            sab = i.s; 
            stk.pop(); 
        }
    }

    return j;
}

double integrate_parallel(double a, double b, double e, std::string &userF)
{
    double j = 0.;
    int p = omp_get_num_procs();
    double slice = (b - a) / static_cast<double> (p);

    std::cout << p << " кусков по " << slice << std::endl;

    omp_set_dynamic(false);

#pragma omp parallel num_threads(p)
    {
        double tid = omp_get_thread_num();

#pragma omp critical
        {
            std::cout << tid << "# " << a + slice * tid << ".." << a + slice * (tid + 1.) << std::endl;
        }

        double jl = integrate(a + slice * tid, a + slice * (tid + 1.), e, userF);

#pragma omp atomic
        j += jl; 

        std::cout << tid << " закончил" << std::endl;
    }

    return j;
}

