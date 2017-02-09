#ifndef STACK_HPP_OQZ1EBAW
#define STACK_HPP_OQZ1EBAW

// данные, описывающие стек
extern int sp;

struct stack_item_t
{
    double a, b, fa, fb, s;
};

#define STK_SIZE 10000000
extern stack_item_t stk[STK_SIZE];

// макроопределения доступа к данным стека

#define STACK_IS_FREE (sp <= 0)

#define PUT_INTO_STACK(A,B,FA,FB,S) \
{ \
    if (sp == STK_SIZE) \
    { \
	std::cerr << "Переполнение локального стека" << std::endl; \
	exit(EXIT_FAILURE); \
    } \
    stk[sp].a = A; \
    stk[sp].b = B; \
    stk[sp].fa = FA; \
    stk[sp].fb = FB; \
    stk[sp].s = S; \
    sp++; \
}

#define GET_FROM_STACK(A,B,FA,FB,S) \
{ \
    sp--; \
    A = stk[sp].a; \
    B = stk[sp].b; \
    FA = stk[sp].fa; \
    FB = stk[sp].fb; \
    S = stk[sp].s; \
}

#endif /* end of include guard: STACK_HPP_OQZ1EBAW */
