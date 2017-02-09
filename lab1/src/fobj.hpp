#ifndef FOBJ_H_1NRBK0OY
#define FOBJ_H_1NRBK0OY

#include "Taylor.hpp"

class fobj_t
{
    public:
        fobj_t(std::string &userF);
        virtual ~fobj_t();
        double y(double x);
        bool isOkay() const;

    private:
        Taylor::Computable *parseUserF(std::string &userF);
        Taylor::Computable *f;
        Taylor::Symbol *x;
        Taylor::Link *linkX;
        bool _okay;
};

#endif /* end of include guard: FOBJ_H_1NRBK0OY */
