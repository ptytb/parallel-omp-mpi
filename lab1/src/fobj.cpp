#include "fobj.hpp"

#include "iostream"

fobj_t::fobj_t(std::string &userF)
{
    _okay = (f = parseUserF(userF));
}

fobj_t::~fobj_t()
{
    delete f;
    delete linkX;
    delete x;
}

bool fobj_t::isOkay() const
{
    return _okay;
}

Taylor::Computable *fobj_t::parseUserF(std::string &userF)
{
    Taylor::Taylor taylor;
    x = new Taylor::Symbol("x", 0.0, true);
    linkX = new Taylor::Link(x);
    taylor.bindSymbol(linkX);

    Taylor::Computable *rawf = NULL;
    Taylor::Error *e = taylor.parse(userF.c_str(), rawf);

    if (e)
    {
        std::cerr << "Ошибка в позиции " << e->getBegin() << ".." << e->getEnd()
            << " : " << e->message() << " : " << e->information() << std::endl;
        return NULL;
        delete e;
    }

    Taylor::Computable *f = rawf->getSimplified();
    delete rawf;
    return f;
}

double fobj_t::y(double x)
{
    //linkX->setLink(this->x);
    this->x->setValue(x);	
    double value = f->getValue();

    if (!isnormal(value) && fpclassify(value) != FP_ZERO)
    {
        std::cerr << "Значение функции не является числом: f(" << x
            << ") = " << value << std::endl;
    }

    return value;
}

