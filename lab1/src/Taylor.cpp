#include "Taylor.hpp"
#include <stdlib.h>
#include <string.h>
#include <stack>

#define STATE_CLEAR 0
#define STATE_NUMBER 1
#define STATE_SYMBOL 2
#define STATE_SYM_NW 3

#define iszero Taylor::Taylor::iszero
#define isone Taylor::Taylor::isone

namespace Taylor {

	int Taylor::accuracy = 20;
	double Taylor::epsilon = 10.0E-20;
	int Taylor::displayAccuracy = 2;
	bool Taylor::expoPres = false;

	void Taylor::setAccuracy(int acc)
	{
		accuracy = acc;
		epsilon = pow(10.0, -acc);
	}
	
	int Taylor::getAccuracy() { return accuracy; }
	void Taylor::setDisplayAccuracy(int acc) { displayAccuracy = acc; }
	int Taylor::getDisplayAccuracy() { return displayAccuracy; }
	bool Taylor::iszero(double x) { return fabs(x) < epsilon; }
	bool Taylor::isone(double x) { return fabs(x - 1.0) < epsilon; }
	void Taylor::setExpoPres(bool ep) { expoPres = ep; }	
	bool Taylor::getExpoPres() { return expoPres; }

	int syntaxMatrix[10][12] = {
		// Может ли после A следовать B
		// A в позиции:				       ^  $
		// A:  / B:	fu  va nu +  -  *  /  ^  (  )
		/* fu */	{0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0},
		/* va */	{0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1},
		/* nu */	{0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1},
		/* +  */	{1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0},
		/* -  */	{1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0},
		/* *  */	{1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0},
		/* /  */	{1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0},
		/* ^  */	{1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0},
		/* (  */	{1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0},
		/* )  */	{0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1}
		// Смещение (строка и столбец массива) у объектов получается методом getCode()
	};

	Error *Taylor::checkSyntax(std::list<Computable*> &infix)
	{
		std::list<Computable*>::iterator i;
		Computable *a = NULL;
		Computable *b = NULL;

		for (i = infix.begin(); i != infix.end(); i++) {
			a = b;
			b = *i;
			if (!a) {
				// В начале выражения
				if (!syntaxMatrix[b->getCode()][10]) {
					return new ESyntax(*(b->getPosition()));
				}
				continue;
			}
			if (!syntaxMatrix[a->getCode()][b->getCode()]) {
				return new ESyntax(*(b->getPosition()));
			}
		}
		// В конце выражения
		if (!syntaxMatrix[b->getCode()][11]) {
			return new ESyntax(*(b->getPosition()));
		}
		return NULL;
	}

	bool isword(char c)
	{
		return c >= 'a' && c <= 'z' ||
			c >= 'A' && c <= 'Z' ||
			c >= '0' && c <= '9' ||
			c == '_';
	}

	bool isnumber(char c)
	{
		return c >= '0' && c <= '9' ||
			c == '.' || c == ',';
	}

	bool isnwsym(char c)
	{
		return c == '+' || c == '-' || c == '*' || c == '/' ||
			c == '^' || c == '(' || c == ')';
	}

	Computable *Taylor::getObject(const char *expr, int argCount)
	{
		Computable *symbol;
		if (strcmp(expr, "Pi") == 0) {
			return new Pi();        
		} else if (strcmp(expr, "Exp") == 0) {
			return new Exponent();        
		} else if (strcmp(expr, "+") == 0) {
			return new OpAdd();
		} else if (strcmp(expr, "-") == 0) {
			if (argCount == 1) {
				return new OpNeg();
			} else {
				return new OpSub();
			}
		} else if (strcmp(expr, "*") == 0) {
			return new OpMul();
		} else if (strcmp(expr, "/") == 0) {
			return new OpDiv();
		} else if (strcmp(expr, "^") == 0) {
			return new OpPow();
		} else if (strcmp(expr, "(") == 0) {
			return new OpOpenRoundBracket();
		} else if (strcmp(expr, ")") == 0) {
			return new OpCloseRoundBracket();
		} else if (strcmp(expr, "sin") == 0) {
			return new FunSin();
		} else if (strcmp(expr, "cos") == 0) {
			return new FunCos();
		} else if (strcmp(expr, "ln") == 0) {
			return new FunLn();
		} else if (strcmp(expr, "tg") == 0) {
			return new FunTg();
		} else if (strcmp(expr, "ctg") == 0) {
			return new FunCtg();
		} else if (strcmp(expr, "sh") == 0) {
			return new FunSinHyp();
		} else if (strcmp(expr, "ch") == 0) {
			return new FunCosHyp();
		} else if (strcmp(expr, "th") == 0) {
			return new FunTgHyp();
		} else if (strcmp(expr, "cth") == 0) {
			return new FunCtgHyp();
		} else if (strcmp(expr, "arcsin") == 0) {
			return new FunArcSin();
		} else if (strcmp(expr, "arccos") == 0) {
			return new FunArcCos();
		} else if (strcmp(expr, "arctg") == 0) {
			return new FunArcTg();
		} else if (strcmp(expr, "arcctg") == 0) {
			return new FunArcCtg();
		} else if (strcmp(expr, "Arsh") == 0) {
			return new FunArSh();
		} else if (strcmp(expr, "Arch") == 0) {
			return new FunArCh();
		} else if (strcmp(expr, "Arth") == 0) {
			return new FunArTh();
		} else if (strcmp(expr, "Arcth") == 0) {
			return new FunArCth();
		} else if (symbol = getSymbol(expr)) {
			return new Link(symbol);
		}
		return NULL;
	}

	Taylor::Taylor() { }

	Taylor::~Taylor() { } 

	Error *Taylor::parse(const char *expr, Computable *&to)
	{
		int len = strlen(expr);
		int beg = 0;
		int end = 0;
		int state = STATE_CLEAR;
		bool unaryModifier = true;

		if (!len) {
			return new EEmptyExpression();
		}

		std::list<Computable*> infix;
		// Входной список лексем infix в инфиксной нотации
		for (int i = 0; i <= len; ++i) {
			switch (state) {
NEW_STATE:
				case STATE_CLEAR:
					if (!isnumber(expr[i]) && isword(expr[i])) {
						state = STATE_SYMBOL;
						unaryModifier = false;
					} else if (isnumber(expr[i])) {
						state = STATE_NUMBER;
						unaryModifier = false;
					} else if (isnwsym(expr[i])) {
						state = STATE_SYM_NW;
					} else if (i == len) {
						continue;
					} else {
						return new EInvalidChar(Position(i, i + 1));
					}
					beg = i;
					break;

				case STATE_SYMBOL:
					end = i;
					if (!isword(expr[i]) || i == len) {
						int symlen = end - beg;
						char *symbol = (char*) malloc(symlen + 1);
						memcpy(symbol, expr + beg, symlen);
						symbol[symlen] = '\0';

						Computable *object = getObject(symbol);
						if (!object) {
							// Освободить память
							std::list<Computable*>::iterator i;
							for (i = infix.begin(); i != infix.end(); i++) {
								delete (*i);
							}
							// Символ не определен
							return new EUndefinedSymbol(Position(beg, end), symbol); 
						}
						object->setPosition(new Position(beg, end));
						infix.push_back(object);

						free(symbol);
						state = STATE_CLEAR;
						goto NEW_STATE;
					}
					break;

				case STATE_NUMBER:
					end = i;
					if (!isnumber(expr[i]) || i == len) {
						int numlen = end - beg;
						char *number = (char*) malloc(numlen + 1);
						memcpy(number, expr + beg, numlen);
						number[numlen] = '\0';

						double value;
						sscanf(number, "%lf", &value);
						Number *object = new Number(value);
						object->setPosition(new Position(beg, end));
						infix.push_back(object);

						free(number);
						state = STATE_CLEAR;
						goto NEW_STATE;
					}
					break;

				case STATE_SYM_NW:
					end = i;

					char op[2];
					op[0] = expr[beg];
					op[1] = '\0';
					Computable *object = getObject(op,
							unaryModifier ? 1 : -1);
					object->setPosition(new Position(beg, end));
					infix.push_back(object);
					if (op[0] == '(') {
						unaryModifier = true;
					} else {
						unaryModifier = false;
					}

					state = STATE_CLEAR;
					goto NEW_STATE;
			}
		}

		// Проверить синтаксис, проверяет только совместимость соседствующих
		// лексем
		Error *e = checkSyntax(infix);
		if (e) {
			// Освободить память
			std::list<Computable*>::iterator i;
			for (i = infix.begin(); i != infix.end(); i++) {
				delete (*i);
			}
			return e;
		}

		std::list<Computable*> outrpn;
		std::stack<Computable*> operators;
		// Выходной список outrpn в обратной польской нотации
		std::list<Computable*>::iterator i;
		for (i = infix.begin(); i != infix.end(); i++) {
			Computable *object = *i;
			if (dynamic_cast<Number*> (object) || dynamic_cast<Link*> (object)) {
				outrpn.push_back(object);
			} else if (dynamic_cast<OpOpenRoundBracket*> (object)) {
				operators.push(object);
			} else if (dynamic_cast<OpCloseRoundBracket*> (object)) {
				bool bracketsMatches = false;
				while (!operators.empty()) {
					Computable *objectInBrackets;
					objectInBrackets = operators.top();
					operators.pop();
					if (dynamic_cast<OpOpenRoundBracket*> (objectInBrackets)) {
						bracketsMatches = true;
						break;
					}
					outrpn.push_back(objectInBrackets);
				}
				if (!bracketsMatches) {
					// Несоответствие скобок, не хватает открывающей
					Position pos = *(object->getPosition());
					// Освободить память
					std::list<Computable*>::iterator i;
					for (i = infix.begin(); i != infix.end(); i++) {
						delete (*i);
					} 
					return new EParenthesesMismatchOpen(pos);
				}
				if (!operators.empty() && dynamic_cast<Function*> (operators.top())) {
					outrpn.push_back(operators.top());
					operators.pop();
				}
			} else if (dynamic_cast<Operator*> (object)) {
				Operator *op1 = (Operator*) object;
				Operator *op2;
				while (!operators.empty() &&
						(op2 = dynamic_cast<Operator*> (operators.top())) &&
						(op1->getAssociativity() == ALEFT &&
						 op1->getPriority() <= op2->getPriority() ||
						 op1->getAssociativity() == ARIGHT &&
						 op1->getPriority() < op2->getPriority())) {
					outrpn.push_back(op2);
					operators.pop();
				}
				operators.push(op1);
			} else {
				// Функция
				operators.push(object);
			}
		}
		while (!operators.empty()) {
			Computable *object = operators.top();
			operators.pop();
			outrpn.push_back(object);

			if (dynamic_cast<OpOpenRoundBracket*> (object)) {
				// Несоответствие скобок, не хватает закрывающей
				Position pos = *(object->getPosition());
				// Освободить память
				std::list<Computable*>::iterator i;
				for (i = infix.begin(); i != infix.end(); i++) {
					delete (*i);
				} 
				return new EParenthesesMismatchClose(pos);
			}
		}

		// Создаем программу Канторовича (дерево) из промежуточной ПОЛИЗ
		std::list<Computable*> args;
		std::list<Computable*>::iterator k;
		for (k = outrpn.begin(); k != outrpn.end(); k++) {
			Computable *object = *k;
			ComputableWithArguments *op = dynamic_cast<ComputableWithArguments*> (object);
			if (op) {
				int arity = op->getArity();
				for (int i = 0; i < arity; ++i) {
					if (*k == outrpn.front()) break;
					k--;
					args.push_front(*k);
				}
				if (args.size() < arity) {
					// Не хватает аргументов
					Position pos = *(object->getPosition());
					// Освободить память
					std::list<Computable*>::iterator i;
					// До текущего объекта все объекты уже в дереве, их нельзя
					// удалять через infix
					bool currentPassed = false;
					for (i = infix.begin(); i != infix.end(); i++) {
						if (currentPassed) {
							delete (*i);
						} else if (*i == object) {
							currentPassed = true;
						}
					} 
					delete object;
					return new ENotEnoughArgs(pos, arity, args.size());
				}
				op->setArguments(args);
				for (int i = 0; i < arity; ++i) {
					k = outrpn.erase(k);
				}
				args.clear();
			}
		}

		to = outrpn.front();
		return NULL;
	}

	Computable *Taylor::getSeries(Computable *source, int degree, Link *x, Symbol *arg, Symbol *a)
	{
		// Считаем значение протзводной dnf в точке a, dnf(arg), arg --> a. a - константа
		x->setLink(a);
		Computable *series = new Number(source->getValue());

		Computable *dnf = source->clone();
		double factorial = 1.0;

		for (int i = 1; i <= degree; ++i) {
			// Берем производную dnf, dnf(arg), arg --> x. x - переменная
			x->setLink(arg);

			Computable *rawdnf = dnf->getDerivative();
			delete dnf;
			dnf = rawdnf->getSimplified();
			delete rawdnf;

			factorial *= i;

			OpAdd *next = new OpAdd();
			OpMul *kan = new OpMul();
			OpPow *pow = new OpPow();
			OpSub *sub = new OpSub();

			// Считаем значение производной dnf в точке a, dnf(arg), arg --> a. a - константа
			x->setLink(a);
			Number *an = new Number(dnf->getValue() / factorial);

			std::list<Computable*> args;
			args.push_back(new Link(arg));
			args.push_back(new Link(a));
			sub->setArguments(args);

			args.clear();
			args.push_back(sub);
			args.push_back(new Number((double) i));
			pow->setArguments(args);

			args.clear();
			args.push_back(an);
			args.push_back(pow);
			kan->setArguments(args);

			args.clear();
			args.push_back(series);
			args.push_back(kan);
			next->setArguments(args);

			series = next;
		}
		delete dnf;
		Computable *seriesSimp = series->getSimplified();
		delete series;
		return seriesSimp;
	}

	void Taylor::getSeriesList(std::list<Computable*> &tsList, Computable *source,
			int degree, Link *x, Symbol *arg, Symbol *a)
	{
		// Считаем значение функции в точке a, source(arg), arg --> a. a - константа
		x->setLink(a);
		Computable *series = new Number(source->getValue());
		tsList.push_back(series->clone());

		if (iszero(series->getValue())) {
			delete series;
			series = NULL;
		}

		Computable *dnf = source->clone();
		double factorial = 1.0;

		for (int i = 1; i <= degree; ++i) {
			// Берем производную dnf, dnf(arg), arg --> x. x - переменная
			x->setLink(arg);

			Computable *rawdnf = dnf->getDerivative();
			delete dnf;
			dnf = rawdnf->getSimplified();
			delete rawdnf;

			factorial *= i;

			OpMul *kan = new OpMul();
			OpPow *pow = new OpPow();
			OpSub *sub = new OpSub();

			// Считаем значение производной dnf в точке a, dnf(arg), arg --> a. a - константа
			x->setLink(a);
			Number *an = new Number(dnf->getValue() / factorial);

			std::list<Computable*> args;
			args.push_back(new Link(arg));
			args.push_back(new Link(a));
			sub->setArguments(args);

			args.clear();
			args.push_back(sub);
			args.push_back(new Number((double) i));
			pow->setArguments(args);

			args.clear();
			args.push_back(an);
			args.push_back(pow);
			kan->setArguments(args);

			Computable *kansimp = kan->getSimplified();
			delete kan;

			if (!(dynamic_cast<Number*> (kansimp) && iszero(kansimp->getValue()))) {
				if (!series) {
					series = kansimp;
				} else {
					OpAdd *next = new OpAdd();
					args.clear();
					args.push_back(series);
					args.push_back(kansimp);
					next->setArguments(args);
					series = next;
				}
				tsList.push_back(i != degree ? series->clone() : series);
			} else {
				delete kansimp;
				if (series) {
					tsList.push_back(series->clone());
				} else {
					tsList.push_back(new Number(0.0));
				}
			}

		}
		delete dnf;
		if (degree == 0 && series) delete series;
	}

	//////////////////////////////////////////////////////

	double OpNeg::getValue()
	{
		return -arguments.front()->getValue();
	}

	Computable *OpNeg::getDerivative()
	{
		OpNeg *neg = new OpNeg();
		std::list<Computable*> args;
		args.push_back(arguments.front()->getDerivative());
		neg->setArguments(args);
		return neg;
	}

	Computable *OpNeg::getSimplified() {
		OpNeg *opNeg = new OpNeg();
		copySimplifiedArgumentsTo(opNeg);
		return opNeg;
	}

	double OpAdd::getValue() {
		Computable *a = arguments.front();
		Computable *b = arguments.back();
		return a->getValue() + b->getValue();
	}

	Computable *OpAdd::getDerivative() {
		Computable *a = arguments.front();
		Computable *b = arguments.back();
		OpAdd *addDeriv = new OpAdd();
		std::list<Computable*> args;
		args.push_back(a->getDerivative());
		args.push_back(b->getDerivative());
		addDeriv->setArguments(args);
		return addDeriv;
	}

	double OpSub::getValue() {
		Computable *a = arguments.front();
		Computable *b = arguments.back();
		return a->getValue() - b->getValue();
	}

	Computable *OpSub::getDerivative() {
		Computable *a = arguments.front();
		Computable *b = arguments.back();
		OpSub *subDeriv = new OpSub();
		std::list<Computable*> args;
		args.push_back(a->getDerivative());
		args.push_back(b->getDerivative());
		subDeriv->setArguments(args);
		return subDeriv;
	}

	double OpMul::getValue() {
		Computable *a = arguments.front();
		Computable *b = arguments.back();
		return a->getValue() * b->getValue();
	}

	Computable *OpMul::getDerivative() {
		Computable *a = arguments.front();
		Computable *b = arguments.back();

		OpAdd *mulDeriv = new OpAdd();
		OpMul *da = new OpMul();
		OpMul *db = new OpMul();

		std::list<Computable*> args;
		args.push_back(a->getDerivative());
		args.push_back(b->clone());
		da->setArguments(args);

		args.clear();
		args.push_back(a->clone());
		args.push_back(b->getDerivative());
		db->setArguments(args);

		args.clear();
		args.push_back(da);
		args.push_back(db);
		mulDeriv->setArguments(args);
		return mulDeriv;
	}

	double OpDiv::getValue() {
		Computable *a = arguments.front();
		Computable *b = arguments.back();
		return a->getValue() / b->getValue();
	}

	Computable *OpDiv::getDerivative() {
		Computable *a = arguments.front();
		Computable *b = arguments.back();

		OpDiv *divDeriv = new OpDiv();
		OpSub *numerator = new OpSub();
		OpMul *da = new OpMul();
		OpMul *db = new OpMul();
		OpPow *denumerator = new OpPow();

		std::list<Computable*> args;
		args.push_back(a->getDerivative());
		args.push_back(b->clone());
		da->setArguments(args);

		args.clear();
		args.push_back(a->clone());
		args.push_back(b->getDerivative());
		db->setArguments(args);

		args.clear();
		args.push_back(da);
		args.push_back(db);
		numerator->setArguments(args);

		args.clear();
		args.push_back(b->clone());
		args.push_back(new Number(2.0));
		denumerator->setArguments(args);

		args.clear();
		args.push_back(numerator);
		args.push_back(denumerator);
		divDeriv->setArguments(args);
		return divDeriv;
	}

	double FunSin::getValue() {
		Computable *arg = arguments.front();
		return sin(arg->getValue());
	}

	Computable *FunSin::getDerivative() {
		OpMul *der = new OpMul();
		FunCos *cos = new FunCos();
		copyArgumentsTo(cos);

		std::list<Computable*> args;
		args.push_back(arguments.front()->getDerivative());
		args.push_back(cos);
		der->setArguments(args);
		return der;
	}

	double FunCos::getValue() {
		Computable *arg = arguments.front();
		return cos(arg->getValue());
	}

	Computable *FunCos::getDerivative() {
		OpMul *dcos = new OpMul();
		OpMul *neg = new OpMul();
		FunSin *sin = new FunSin();

		copyArgumentsTo(sin);

		std::list<Computable*> args;
		args.push_back(new Number(-1.0));
		args.push_back(sin);
		neg->setArguments(args);

		args.clear();
		args.push_back(arguments.front()->getDerivative());
		args.push_back(neg);
		dcos->setArguments(args);
		return dcos;
	}

	double FunLn::getValue() {
		Computable *arg = arguments.front();
		return log(arg->getValue());
	}

	Computable *FunLn::getDerivative() {
		Computable *arg = arguments.front();
		OpDiv *dln = new OpDiv();

		std::list<Computable*> args;
		args.push_back(arg->getDerivative());
		args.push_back(arg->clone());
		dln->setArguments(args);
		return dln;
	}

	double FunTg::getValue() {
		Computable *arg = arguments.front();
		return tan(arg->getValue());
	}

	Computable *FunTg::getDerivative() {
		OpDiv *div = new OpDiv();
		OpPow *pow = new OpPow();
		FunCos *cos = new FunCos();
		copyArgumentsTo(cos);

		std::list<Computable*> args;
		args.push_back(cos);
		args.push_back(new Number(2.0));
		pow->setArguments(args);

		args.clear();
		args.push_back(arguments.front()->getDerivative());
		args.push_back(pow);
		div->setArguments(args);
		return div;
	}

	double FunCtg::getValue() {
		Computable *arg = arguments.front();
		return 1.0/tan(arg->getValue());
	}

	Computable *FunCtg::getDerivative() {
		OpDiv *div = new OpDiv();
		OpPow *pow = new OpPow();
		FunSin *sin = new FunSin();
		OpNeg *neg = new OpNeg();
		copyArgumentsTo(sin);

		std::list<Computable*> args;
		args.push_back(sin);
		args.push_back(new Number(2.0));
		pow->setArguments(args);

		args.clear();
		args.push_back(arguments.front()->getDerivative());
		args.push_back(pow);
		div->setArguments(args);

		args.clear();
		args.push_back(div);
		neg->setArguments(args);
		return neg;
	}

	double FunTgHyp::getValue() {
		Computable *arg = arguments.front();
		return tanh(arg->getValue());
	}

	Computable *FunTgHyp::getDerivative() {
		OpDiv *div = new OpDiv();
		OpPow *pow = new OpPow();
		FunCosHyp *cosh = new FunCosHyp();
		copyArgumentsTo(cosh);

		std::list<Computable*> args;
		args.push_back(cosh);
		args.push_back(new Number(2.0));
		pow->setArguments(args);

		args.clear();
		args.push_back(arguments.front()->getDerivative());
		args.push_back(pow);
		div->setArguments(args);
		return div;
	}

	double FunCtgHyp::getValue() {
		Computable *arg = arguments.front();
		return 1.0/tanh(arg->getValue());
	}

	Computable *FunCtgHyp::getDerivative() {
		OpDiv *div = new OpDiv();
		OpPow *pow = new OpPow();
		FunSinHyp *sinh = new FunSinHyp();
		OpNeg *neg = new OpNeg();
		copyArgumentsTo(sinh);

		std::list<Computable*> args;
		args.push_back(sinh);
		args.push_back(new Number(2.0));
		pow->setArguments(args);

		args.clear();
		args.push_back(arguments.front()->getDerivative());
		args.push_back(pow);
		div->setArguments(args);

		args.clear();
		args.push_back(div);
		neg->setArguments(args);
		return neg;
	}

	double FunSinHyp::getValue()
	{
		return sinh(arguments.front()->getValue());
	}

	Computable *FunSinHyp::getDerivative()
	{
		OpMul *dsh = new OpMul();
		FunCosHyp *ch = new FunCosHyp();
		
		copyArgumentsTo(ch);

		std::list<Computable*> args;
		args.push_back(arguments.front()->getDerivative());
		args.push_back(ch);
		dsh->setArguments(args);		
		return dsh;
	}

	double FunCosHyp::getValue()
	{
		return cosh(arguments.front()->getValue());
	}

	Computable *FunCosHyp::getDerivative()
	{
		OpMul *dch = new OpMul();
		FunSinHyp *sh = new FunSinHyp();
		
		copyArgumentsTo(sh);

		std::list<Computable*> args;
		args.push_back(arguments.front()->getDerivative());
		args.push_back(sh);
		dch->setArguments(args);		
		return dch;
	}

	double FunArcSin::getValue()
	{
		return asin(arguments.front()->getValue());
	}

	Computable *FunArcSin::getDerivative()
	{
		Computable *arg = arguments.front();
		OpDiv *darcsin = new OpDiv();
		OpPow *powDenum = new OpPow();
		OpSub *sub = new OpSub();
		OpPow *powArg = new OpPow();
		
		std::list<Computable*> args;
		args.push_back(arg->clone());
		args.push_back(new Number(2.0));
		powArg->setArguments(args);

		args.clear();
		args.push_back(new Number(1.0));
		args.push_back(powArg);
		sub->setArguments(args);

		args.clear();
		args.push_back(sub);
		args.push_back(new Number(0.5));
		powDenum->setArguments(args);

		args.clear();
		args.push_back(arg->getDerivative());
		args.push_back(powDenum);
		darcsin->setArguments(args);
		return darcsin;
	}

	double FunArcCos::getValue()
	{
		return acos(arguments.front()->getValue());
	}

	Computable *FunArcCos::getDerivative()
	{
		Computable *arg = arguments.front();
		OpNeg *negfrac = new OpNeg();
		OpDiv *frac = new OpDiv();
		OpPow *powDenum = new OpPow();
		OpSub *sub = new OpSub();
		OpPow *powArg = new OpPow();
		
		std::list<Computable*> args;
		args.push_back(arg->clone());
		args.push_back(new Number(2.0));
		powArg->setArguments(args);

		args.clear();
		args.push_back(new Number(1.0));
		args.push_back(powArg);
		sub->setArguments(args);

		args.clear();
		args.push_back(sub);
		args.push_back(new Number(0.5));
		powDenum->setArguments(args);

		args.clear();
		args.push_back(arg->getDerivative());
		args.push_back(powDenum);
		frac->setArguments(args);

		args.clear();
		args.push_back(frac);
		negfrac->setArguments(args);
		return negfrac;
	}

	double FunArcTg::getValue()
	{
		return atan(arguments.front()->getValue());
	}
	
	Computable *FunArcTg::getDerivative()
	{
		Computable *arg = arguments.front();
		OpDiv *datan = new OpDiv();
		OpAdd *denum = new OpAdd();
		OpPow *pow = new OpPow();

		std::list<Computable*> args;
		args.push_back(arg->clone());
		args.push_back(new Number(2.0));
		pow->setArguments(args);

		args.clear();
		args.push_back(new Number(1.0));
		args.push_back(pow);
		denum->setArguments(args);

		args.clear();
		args.push_back(arg->getDerivative());
		args.push_back(denum);
		datan->setArguments(args);

		return datan;
	}

	double FunArcCtg::getValue()
	{
		double val = arguments.front()->getValue();
		return val >= 0.0 ?
			asin(1.0/sqrt(1.0+val*val)) :
			PI_VALUE - asin(1.0/sqrt(1.0+val*val));
	}

	Computable *FunArcCtg::getDerivative()
	{
		Computable *arg = arguments.front();
		OpNeg *dactan = new OpNeg();
		OpDiv *frac = new OpDiv();
		OpAdd *denum = new OpAdd();
		OpPow *pow = new OpPow();

		std::list<Computable*> args;
		args.push_back(arg->clone());
		args.push_back(new Number(2.0));
		pow->setArguments(args);

		args.clear();
		args.push_back(new Number(1.0));
		args.push_back(pow);
		denum->setArguments(args);

		args.clear();
		args.push_back(arg->getDerivative());
		args.push_back(denum);
		frac->setArguments(args);

		args.clear();
		args.push_back(frac);
		dactan->setArguments(args);
		return dactan;
	}

	double FunArSh::getValue()
	{
		return asinh(arguments.front()->getValue());
	}

	Computable *FunArSh::getDerivative()
	{
		Computable *arg = arguments.front();
		OpDiv *darsh = new OpDiv();
		OpPow *powDenum = new OpPow();
		OpAdd *add = new OpAdd();
		OpPow *powArg = new OpPow();
		
		std::list<Computable*> args;
		args.push_back(arg->clone());
		args.push_back(new Number(2.0));
		powArg->setArguments(args);

		args.clear();
		args.push_back(new Number(1.0));
		args.push_back(powArg);
		add->setArguments(args);

		args.clear();
		args.push_back(add);
		args.push_back(new Number(0.5));
		powDenum->setArguments(args);

		args.clear();
		args.push_back(arg->getDerivative());
		args.push_back(powDenum);
		darsh->setArguments(args);
		return darsh;
	}

	double FunArCh::getValue()
	{
		return acosh(arguments.front()->getValue());
	}

	Computable *FunArCh::getDerivative()
	{
		Computable *arg = arguments.front();
		OpDiv *darch = new OpDiv();
		OpPow *powDenum = new OpPow();
		OpSub *sub = new OpSub();
		OpPow *powArg = new OpPow();
		
		std::list<Computable*> args;
		args.push_back(arg->clone());
		args.push_back(new Number(2.0));
		powArg->setArguments(args);

		args.clear();
		args.push_back(powArg);
		args.push_back(new Number(1.0));
		sub->setArguments(args);

		args.clear();
		args.push_back(sub);
		args.push_back(new Number(0.5));
		powDenum->setArguments(args);

		args.clear();
		args.push_back(arg->getDerivative());
		args.push_back(powDenum);
		darch->setArguments(args);
		return darch;
	}

	double FunArTh::getValue()
	{
		return atanh(arguments.front()->getValue());
	}

	Computable *FunArTh::getDerivative()
	{
		Computable *arg = arguments.front();
		OpDiv *darth = new OpDiv();
		OpSub *denum = new OpSub();
		OpPow *pow = new OpPow();

		std::list<Computable*> args;
		args.push_back(arg->clone());
		args.push_back(new Number(2.0));
		pow->setArguments(args);

		args.clear();
		args.push_back(new Number(1.0));
		args.push_back(pow);
		denum->setArguments(args);

		args.clear();
		args.push_back(arg->getDerivative());
		args.push_back(denum);
		darth->setArguments(args);

		return darth;
	}

	double FunArCth::getValue()
	{
		double val = arguments.front()->getValue();
		return 0.5 * log((val + 1.0)/(val - 1.0));
	}

	Computable *FunArCth::getDerivative()
	{
		Computable *arg = arguments.front();
		OpNeg *darcth = new OpNeg();
		OpDiv *frac = new OpDiv();
		OpSub *denum = new OpSub();
		OpPow *pow = new OpPow();

		std::list<Computable*> args;
		args.push_back(arg->clone());
		args.push_back(new Number(2.0));
		pow->setArguments(args);

		args.clear();
		args.push_back(pow);
		args.push_back(new Number(1.0));
		denum->setArguments(args);

		args.clear();
		args.push_back(arg->getDerivative());
		args.push_back(denum);
		frac->setArguments(args);

		args.clear();
		args.push_back(frac);
		darcth->setArguments(args);
		return darcth;
	}

	double OpPow::getValue() {
		Computable *a = arguments.front();
		Computable *b = arguments.back();
		return pow(a->getValue(), b->getValue());
	}

	Computable *OpPow::getDerivative() {
		Computable *powDerivRes;
		Computable *a = arguments.front();
		Computable *b = arguments.back();
		Link *al = dynamic_cast<Link*> (a);
		Link *bl = dynamic_cast<Link*> (b);
		Number *an = NULL;
		Number *bn = NULL;

		if (al) {
			an = dynamic_cast<Number*> (al->getSource());
		} else {
			an = dynamic_cast<Number*> (a);
		}
		if (bl) {
			bn = dynamic_cast<Number*> (bl->getSource());
		} else {
			bn = dynamic_cast<Number*> (b);
		}

		bool aIsVar = an && an->isVariable() || !an;
		bool bIsVar = bn && bn->isVariable() || !bn;

		if (aIsVar && !bIsVar) {
			OpMul *powDerivFull = new OpMul();
			OpMul *powDeriv = new OpMul();
			OpPow *bpow = new OpPow();
			OpSub *bsub = new OpSub();

			std::list<Computable*> args;
			args.push_back(b->clone());
			args.push_back(new Number(1.0));
			bsub->setArguments(args);

			args.clear();
			args.push_back(a->clone());
			args.push_back(bsub);
			bpow->setArguments(args);

			args.clear();
			args.push_back(b->clone());
			args.push_back(bpow);
			powDeriv->setArguments(args);

			args.clear();
			args.push_back(a->getDerivative());
			args.push_back(powDeriv);
			powDerivFull->setArguments(args);

			powDerivRes = powDerivFull;
		} else if (bIsVar && !aIsVar) {
			OpMul *powDerivFull = new OpMul();
			OpMul *powDeriv = new OpMul();
			OpPow *powab = new OpPow();
			FunLn *lna = new FunLn();

			std::list<Computable*> args;
			args.push_back(a->clone());
			args.push_back(b->clone());
			powab->setArguments(args);

			args.clear();
			args.push_back(a->clone());
			lna->setArguments(args);

			args.clear();
			args.push_back(powab);
			args.push_back(lna);
			powDeriv->setArguments(args);

			args.clear();
			args.push_back(b->getDerivative());
			args.push_back(powDeriv);
			powDerivFull->setArguments(args);

			powDerivRes = powDerivFull;
		} else if (an && bn && aIsVar && bIsVar) {
			OpMul *powDeriv = new OpMul();
			OpPow *pow = new OpPow();
			OpAdd *sum = new OpAdd();
			FunLn *ln = new FunLn();

			std::list<Computable*> args;
			args.push_back(a->clone());
			ln->setArguments(args);

			args.clear();
			args.push_back(ln);
			args.push_back(new Number(1.0));
			sum->setArguments(args);

			args.clear();
			args.push_back(a->clone());
			args.push_back(b->clone());
			pow->setArguments(args);

			args.clear();
			args.push_back(pow);
			args.push_back(sum);
			powDeriv->setArguments(args);

			powDerivRes = powDeriv;
		} else if (an && bn && !aIsVar && !bIsVar) {
			return new Number(0.0);
		} else {
			OpAdd *powDeriv = new OpAdd();
			OpPow *leftPow = new OpPow();
			OpMul *leftMul1 = new OpMul();
			OpMul *leftMul2 = new OpMul();
			FunLn *ln = new FunLn();
			OpSub *sub = new OpSub();
			OpPow *rightPow = new OpPow();
			OpMul *rightMul1 = new OpMul();
			OpMul *rightMul2 = new OpMul();

			std::list<Computable*> args;
			args.push_back(a->clone());
			args.push_back(b->clone());
			leftPow->setArguments(args);

			args.clear();
			args.push_back(leftPow);
			args.push_back(b->getDerivative());
			leftMul1->setArguments(args);

			args.clear();
			args.push_back(a->clone());
			ln->setArguments(args);

			args.clear();
			args.push_back(leftMul1);
			args.push_back(ln);
			leftMul2->setArguments(args);

			args.clear();
			args.push_back(b->clone());
			args.push_back(new Number(1.0));
			sub->setArguments(args);

			args.clear();
			args.push_back(a->clone());
			args.push_back(sub);
			rightPow->setArguments(args);

			args.clear();
			args.push_back(rightPow);
			args.push_back(b->clone());
			rightMul1->setArguments(args);

			args.clear();
			args.push_back(rightMul1);
			args.push_back(a->getDerivative());
			rightMul2->setArguments(args);

			args.clear();
			args.push_back(leftMul2);
			args.push_back(rightMul2);
			powDeriv->setArguments(args);

			powDerivRes = powDeriv;
		}

		return powDerivRes;
	}

	void Scope::bindSymbol(Computable *symbol) {
		symbols[symbol->toString()] = symbol;
	}

	void Scope::releaseSymbol(const char *name) {
		symbols.erase(name);
	}

	Computable *Scope::getSymbol(const char *name) {
		return (symbols.count(name) != 0) ?
			symbols[name] : NULL;
	}

	char *ComputableWithArguments::toString() {
		const char *mask = getMask();
		int resLen = 0;
		std::list<char*> out;
		std::list<Computable*>::iterator arg = arguments.begin();
		for (int i = 0; mask[i] != '\0'; ++i) {
			if (mask[i] == '@') {
				char *argStr = (*arg)->toString();
				resLen += strlen(argStr);
				out.push_back(argStr);
				arg++;					
			} else {
				char *ch = (char*) malloc(sizeof(char) * 2);
				ch[0] = mask[i];
				ch[1] = '\0';
				out.push_back(ch);
				resLen++;
			}
		}
		char *resStr = (char*) malloc(sizeof(char) * resLen + 1);
		resStr[0] = '\0';
		std::list<char*>::iterator str;
		for (str = out.begin(); str != out.end(); str++) {
			strcat(resStr, *str);
			free(*str);                
		}
		return resStr;
	}

	char *Number::toString() {
		double val = getValue();
		if (Taylor::getExpoPres()) {
			int len = snprintf(NULL, 0, "%.*le", Taylor::getDisplayAccuracy(), val);
			char *str = (char*) malloc(sizeof(char) * len + 1);
			sprintf(str, "%.*le", Taylor::getDisplayAccuracy(), val); 
			return str;
		} else {
			int len = snprintf(NULL, 0, "%.*lf", Taylor::getDisplayAccuracy(), val);
			char *str = (char*) malloc(sizeof(char) * len + 1);
			sprintf(str, "%.*lf", Taylor::getDisplayAccuracy(), val); 
			return str;
		}
	}

	Computable *OpAdd::getSimplified() {
		Computable *sa = arguments.front()->getSimplified();
		Computable *sb = arguments.back()->getSimplified();
		Number *na = dynamic_cast<Number*> (sa);
		Number *nb = dynamic_cast<Number*> (sb);
		if (na && nb) {
			Number *res = new Number(na->getValue() + nb->getValue());
			delete na;
			delete nb;
			return res;
		} else if (na && iszero(na->getValue())) {
			delete na;
			return sb;
		} else if (nb && iszero(nb->getValue())) {
			delete nb;
			return sa;
		}
		OpAdd *res = new OpAdd();
		std::list<Computable*> args;
		args.push_back(sa);
		args.push_back(sb);
		res->setArguments(args);
		return res;
	}

	Computable *OpSub::getSimplified() {
		Computable *sa = arguments.front()->getSimplified();
		Computable *sb = arguments.back()->getSimplified();
		Number *na = dynamic_cast<Number*> (sa);
		Number *nb = dynamic_cast<Number*> (sb);
		if (na && nb) {
			Number *res = new Number(na->getValue() - nb->getValue());
			delete na;
			delete nb;
			return res;
		} else if (na && iszero(na->getValue())) {
			OpMul *res = new OpMul();

			std::list<Computable*> args;
			args.push_back(new Number(-1.0));
			args.push_back(sb);
			res->setArguments(args);

			delete na;
			return res;
		} else if (nb && iszero(nb->getValue())) {
			delete nb;
			return sa;
		}
		OpSub *res = new OpSub();
		std::list<Computable*> args;
		args.push_back(sa);
		args.push_back(sb);
		res->setArguments(args);
		return res;
	}

	Computable *OpMul::getSimplified() {
		Computable *sa = arguments.front()->getSimplified();
		Computable *sb = arguments.back()->getSimplified();
		Number *na = dynamic_cast<Number*> (sa);
		Number *nb = dynamic_cast<Number*> (sb);

		if (na && iszero(na->getValue())) {
			delete na;
			delete sb;
			return new Number(0.0);
		} else if (nb && iszero(nb->getValue())) {
			delete nb;
			delete sa;
			return new Number(0.0);
		} else if (na && nb) {
			Number *res = new Number(na->getValue() * nb->getValue());
			delete na;
			delete nb;
			return res;
		} else if (na && isone(na->getValue())) {
			delete na;
			return sb;
		} else if (nb && isone(nb->getValue())) {
			delete nb;
			return sa;
		}
		OpMul *res = new OpMul();
		std::list<Computable*> args;
		args.push_back(sa);
		args.push_back(sb);
		res->setArguments(args);
		return res;
	}

	Computable *OpDiv::getSimplified() {
		Computable *sa = arguments.front()->getSimplified();
		Computable *sb = arguments.back()->getSimplified();
		Number *na = dynamic_cast<Number*> (sa);
		Number *nb = dynamic_cast<Number*> (sb);

		// В числителе 0
		if (na && iszero(na->getValue())) {
			delete na;
			delete sb;
			return new Number(0.0);
			// В знаменателе 1
		} else if (nb && isone(nb->getValue())) {
			delete nb;
			return sa;
			// Числитель и знаменатель числа
		} else if (na && nb) {
			Number *res = new Number(na->getValue() / nb->getValue());
			delete na;
			delete nb;
			return res;
		}

		OpDiv *res = new OpDiv();
		std::list<Computable*> args;
		args.push_back(sa);
		args.push_back(sb);
		res->setArguments(args);
		return res;
	}

	Computable *OpPow::getSimplified()
	{
		Computable *sa = arguments.front()->getSimplified();
		Computable *sb = arguments.back()->getSimplified();
		Number *na = dynamic_cast<Number*> (sa);
		Number *nb = dynamic_cast<Number*> (sb);

		if (nb && iszero(nb->getValue())) {
			delete nb;
			delete sa;
			return new Number(1.0);
		} else if (nb && isone(nb->getValue())) {
			delete nb;
			return sa;
		}

		OpPow *res = new OpPow();
		std::list<Computable*> args;
		args.push_back(sa);
		args.push_back(sb);
		res->setArguments(args);
		return res;
	}

} // namespace Taylor


