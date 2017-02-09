#pragma once

/** Taylor v0.1
 * Разложение функции в ряд Тейлора
 * Автор: Пронин Илья
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <list>
#include <map>
#include <string>
#include <stdarg.h>

#define ALEFT 1
#define ARIGHT 2

#define PI_VALUE	3.14159265358979323846
#define EXP_VALUE	2.71828182845904523536

namespace Taylor {

	class Position
	{
		public:
			Position(int beg, int end) {
				this->beg = beg;
				this->end = end;
			}

			int getBegin() { return beg; }
			int getEnd() { return end; }

		private:
			int beg;
			int end;
	};

	class SymbolInfo
	{
		public:
			SymbolInfo(const char *name) {
				this->name = strdup(name);
			}
			
			SymbolInfo(SymbolInfo &si) {
				name = strdup(si.name);
			}

			virtual ~SymbolInfo() { free(name); }
			virtual char *toString() {
				return strdup(name);
			}

		protected:
			char *name;
	};

	class Computable
	{
		public:
			Computable() { pos = NULL; }
			Computable(Computable &from) { pos = NULL; }
			virtual ~Computable() { if (pos) delete pos; }

			virtual Computable& operator=(Computable &from) {
				pos = NULL;
				return *this;
			}

			virtual double getValue() = 0;
			virtual Computable *getDerivative() = 0;
			virtual char *toString() = 0;
			virtual int getCode() = 0;
			virtual Computable *clone() = 0;
			virtual Computable *getSimplified() = 0;

			virtual void setPosition(Position *pos) { this->pos = pos; }
			virtual Position *getPosition() { return pos; }

		private:
			// Позиция объекта в исходной строке, для вывода сообщений
			// об ошибке
			Position *pos;
	};

	class Number : public Computable
	{
		public:
			Number(double value, bool variable = false) {
				this->value = value;
				this->variable = variable;
			}
			virtual ~Number() { }
			virtual double setValue(double value) { this->value = value; }
			virtual double getValue() { return value; }
			virtual Computable *getDerivative() {
				return isVariable() ?
					new Number(1.0) : new Number(0.0);
			}
			virtual char *toString();
			virtual bool isVariable() { return variable; }
			virtual void setVariable(bool variable) { this->variable = variable; }
			virtual int getCode() { return 2; }

			virtual Computable *clone() { return new Number(*this); }
			virtual Computable *getSimplified() { return clone(); }

		protected:
			double value;
			bool variable;
	};

	class Symbol : public SymbolInfo, public Number
	{
		public:
			Symbol(const char *name, double value, bool variable = false,
					bool displayAsNumber = false) :
				SymbolInfo(name), Number(value, variable) {
					this->displayAsNumber = displayAsNumber;
				}
			virtual ~Symbol() { }
			virtual char *toString() {
				if (displayAsNumber) {
					return Number::toString();
				} else {
					return SymbolInfo::toString();
				}
			}
			virtual int getCode() { return 1; }
			virtual Computable *clone() { return new Symbol(*this); }
			virtual Computable *getSimplified() { return clone(); }
			bool getDisplayAsNumber() { return displayAsNumber; }

		private:
			bool displayAsNumber;
	};

	class Link : public Computable
	{
		public:
			Link(Computable *object) { this->object = object; }
			virtual ~Link() { }
			virtual double getValue() { return getSource()->getValue(); }
			virtual Computable *getDerivative() { return getSource()->getDerivative(); }
			virtual char *toString() { return getSource()->toString(); }
			virtual int getCode() { return getSource()->getCode(); }
			virtual Computable *clone() { return new Link(*this); }
			virtual Computable *getSimplified() { return clone(); }

			Computable *getLink() { return object; }
			void setLink(Computable *object) { this->object = object; }

			Computable *getSource() {
				Computable *src = object;
				Link *ln;
				while (ln = dynamic_cast<Link*> (src)) {
					src = ln->getLink();
				}
				return src;
			}

		protected:
			Computable *object;
	};

	class Scope
	{
		public:
			void bindSymbol(Computable *symbol);
			void releaseSymbol(const char *name);
			Computable *getSymbol(const char *name);

		protected:
			std::map<std::string, Computable*> symbols;
	};

	class ComputableWithArguments : public Computable
	{
		public:
		/*
			ComputableWithArguments(...) {
				va_list ap;
				va_start(ap, NULL);
				for (int i = 0; i < 2; ++i) {
					Computable *object = va_arg(ap, Computable*);
				}
				va_end(ap);
			}
		*/
			virtual ~ComputableWithArguments() {
				std::list<Computable*>::iterator i;
				for (i = arguments.begin(); i != arguments.end(); i++) {
					delete (*i);
				}
			}

			virtual char *toString();

			virtual void setArguments(std::list<Computable*> arguments)
			{
				this->arguments = arguments;
			}

			virtual const char *getMask() = 0;
			virtual int getArity() = 0;

		protected:
			std::list<Computable*> arguments;

			virtual void copyArgumentsTo(ComputableWithArguments *to) {
				std::list<Computable*>::iterator i;
				for (i = arguments.begin(); i != arguments.end(); i++) {
					to->arguments.push_back((*i)->clone());
				}
			}

			virtual void copySimplifiedArgumentsTo(ComputableWithArguments *to) {
				std::list<Computable*>::iterator i;
				for (i = arguments.begin(); i != arguments.end(); i++) {
					to->arguments.push_back((*i)->getSimplified());
				}
			}
	};

	class Operator : public ComputableWithArguments
	{
		public:
			virtual ~Operator() { }
			virtual int getPriority() = 0;
			virtual int getAssociativity() = 0;

	};

	class Function : public ComputableWithArguments
	{
		public:
			virtual ~Function() { }
			virtual int getCode() { return 0; }
	};

	//////////////////////////////////////////////////////

	class Pi : public Symbol
	{
		public:
			Pi() : Symbol("Pi", PI_VALUE) { }
			virtual ~Pi() { }
	};

	class Exponent : public Symbol
	{
		public:
			Exponent() : Symbol("Exp", EXP_VALUE) { }
			virtual ~Exponent() { }
	};

	//////////////////////////////////////////////////////

	class OpNeg : public Operator
	{
		public:
			virtual ~OpNeg() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getPriority() { return 5; }
			virtual int getArity() { return 1; }
			virtual int getAssociativity() { return ARIGHT; }
			virtual const char *getMask() { return "-(@)"; }
			virtual int getCode() { return 4; }
			
			virtual Computable *clone() {
				OpNeg *opNeg = new OpNeg();
				copyArgumentsTo(opNeg);
				return opNeg;
			}
			virtual Computable *getSimplified();
	};

	class OpAdd : public Operator
	{
		public:
			virtual ~OpAdd() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getPriority() { return 2; }
			virtual int getArity() { return 2; }
			virtual int getAssociativity() { return ALEFT; }
			virtual const char *getMask() { return "@+@"; }
			virtual int getCode() { return 3; }
			
			virtual Computable *clone() {
				OpAdd *opAdd = new OpAdd();
				copyArgumentsTo(opAdd);
				return opAdd;
			}
			virtual Computable *getSimplified();
	};

	class OpSub : public Operator
	{
		public:
			virtual ~OpSub() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getPriority() { return 2; }
			virtual int getArity() { return 2; }
			virtual int getAssociativity() { return ALEFT; }
			virtual const char *getMask() { return "@-@"; }
			virtual int getCode() { return 4; }
			
			virtual Computable *clone() {
				OpSub *opSub = new OpSub();
				copyArgumentsTo(opSub);
				return opSub;
			}
			virtual Computable *getSimplified();
	};

	class OpMul : public Operator
	{
		public:
			virtual ~OpMul() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getPriority() { return 3; }
			virtual int getArity() { return 2; }
			virtual int getAssociativity() { return ALEFT; }
			virtual const char *getMask() { return "(@)*(@)"; }
			virtual int getCode() { return 5; }
			
			virtual Computable *clone() {
				OpMul *opMul = new OpMul();
				copyArgumentsTo(opMul);
				return opMul;
			}
			virtual Computable *getSimplified();
	};

	class OpPow : public Operator
	{
		public:
			virtual ~OpPow() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getPriority() { return 4; }
			virtual int getArity() { return 2; }
			virtual int getAssociativity() { return ARIGHT; }
			virtual const char *getMask() { return "(@)^(@)"; }
			virtual int getCode() { return 7; }
			
			virtual Computable *clone() {
				OpPow *opPow = new OpPow();
				copyArgumentsTo(opPow);
				return opPow;
			}
			virtual Computable *getSimplified();
	};

	class OpDiv : public Operator
	{
		public:
			virtual ~OpDiv() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getPriority() { return 3; }
			virtual int getArity() { return 2; }
			virtual int getAssociativity() { return ALEFT; }
			virtual const char *getMask() { return "(@)/(@)"; }
			virtual int getCode() { return 6; }
			
			virtual Computable *clone() {
				OpDiv *opDiv = new OpDiv();
				copyArgumentsTo(opDiv);
				return opDiv;
			}
			virtual Computable *getSimplified();
	};

	class OpOpenRoundBracket : public Operator
	{
		public:
			virtual ~OpOpenRoundBracket() { }
			virtual double getValue() { return 0; }
			virtual Computable *getDerivative() { return NULL; }
			virtual int getPriority() { return 0; }
			virtual int getArity() { return 0; }
			virtual int getAssociativity() { return ALEFT; }
			virtual const char *getMask() { return NULL; }
			virtual int getCode() { return 8; }
			virtual Computable *clone() { return new OpOpenRoundBracket(); }
			virtual Computable *getSimplified() { return clone(); }
	};

	class OpCloseRoundBracket : public Operator
	{
		public:
			virtual ~OpCloseRoundBracket() { }
			virtual double getValue() { return 0; }
			virtual Computable *getDerivative() { return NULL; }
			virtual int getPriority() { return 1; }
			virtual int getArity() { return 0; }
			virtual int getAssociativity() { return ALEFT; }
			virtual const char *getMask() { return NULL; }
			virtual int getCode() { return 9; }
			virtual Computable *clone() { return new OpCloseRoundBracket(); }
			virtual Computable *getSimplified() { return clone(); }
	};

	class FunSin : public Function
	{
		public:
			virtual ~FunSin() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "sin(@)"; }
			
			virtual Computable *clone() {
				FunSin *funSin = new FunSin();
				copyArgumentsTo(funSin);
				return funSin;
			}

			virtual Computable *getSimplified() {
				FunSin *funSin = new FunSin();
				copySimplifiedArgumentsTo(funSin);	
				return funSin;
			}

	};

	class FunCos : public Function
	{
		public:
			virtual ~FunCos() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "cos(@)"; }
			
			virtual Computable *clone() {
				FunCos *funCos = new FunCos();
				copyArgumentsTo(funCos);
				return funCos;
			}
			virtual Computable *getSimplified() {
				FunCos *funCos = new FunCos();
				copySimplifiedArgumentsTo(funCos);
				return funCos;
			}
	};

	class FunLn : public Function
	{
		public:
			virtual ~FunLn() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "ln(@)"; }
			
			virtual Computable *clone() {
				FunLn *funLn = new FunLn();
				copyArgumentsTo(funLn);
				return funLn;
			}
			virtual Computable *getSimplified() {
				FunLn *funLn = new FunLn();
				copySimplifiedArgumentsTo(funLn);
				return funLn;
			}
	};

	class FunTg : public Function
	{
		public:
			virtual ~FunTg() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "tg(@)"; }
			
			virtual Computable *clone() {
				FunTg *funTg = new FunTg();
				copyArgumentsTo(funTg);
				return funTg;
			}
			virtual Computable *getSimplified() {
				FunTg *funTg = new FunTg();
				copySimplifiedArgumentsTo(funTg);
				return funTg;
			}
	};

	class FunCtg : public Function
	{
		public:
			virtual ~FunCtg() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "ctg(@)"; }
			
			virtual Computable *clone() {
				FunCtg *funTg = new FunCtg();
				copyArgumentsTo(funTg);
				return funTg;
			}
			virtual Computable *getSimplified() {
				FunCtg *funCtg = new FunCtg();
				copySimplifiedArgumentsTo(funCtg);
				return funCtg;
			}
	};

	class FunSinHyp : public Function
	{
		public:
			virtual ~FunSinHyp() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "sh(@)"; }
			
			virtual Computable *clone() {
				FunSinHyp *funSin = new FunSinHyp();
				copyArgumentsTo(funSin);
				return funSin;
			}

			virtual Computable *getSimplified() {
				FunSinHyp *funSinHyp = new FunSinHyp();
				copySimplifiedArgumentsTo(funSinHyp);
				return funSinHyp;
			}
	};

	class FunCosHyp : public Function
	{
		public:
			virtual ~FunCosHyp() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "ch(@)"; }
			
			virtual Computable *clone() {
				FunCosHyp *funCosHyp = new FunCosHyp();
				copyArgumentsTo(funCosHyp);
				return funCosHyp;
			}

			virtual Computable *getSimplified() {
				FunCosHyp *funCosHyp = new FunCosHyp();
				copySimplifiedArgumentsTo(funCosHyp);
				return funCosHyp;
			}
	};

	class FunTgHyp : public Function
	{
		public:
			virtual ~FunTgHyp() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "th(@)"; }
			
			virtual Computable *clone() {
				FunTgHyp *funCosHyp = new FunTgHyp();
				copyArgumentsTo(funCosHyp);
				return funCosHyp;
			}

			virtual Computable *getSimplified() {
				FunTgHyp *funTgHyp  = new FunTgHyp();
				copySimplifiedArgumentsTo(funTgHyp );
				return funTgHyp ;
			}
	};

	class FunCtgHyp : public Function
	{
		public:
			virtual ~FunCtgHyp() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "cth(@)"; }
			
			virtual Computable *clone() {
				FunCtgHyp *funCosHyp = new FunCtgHyp();
				copyArgumentsTo(funCosHyp);
				return funCosHyp;
			}

			virtual Computable *getSimplified() {
				FunCtgHyp *funCtgHyp  = new FunCtgHyp();
				copySimplifiedArgumentsTo(funCtgHyp );
				return funCtgHyp ;
			}
	};

	class FunArcSin : public Function
	{
		public:
			virtual ~FunArcSin() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "arcsin(@)"; }
			
			virtual Computable *clone() {
				FunArcSin *funArcSin = new FunArcSin();
				copyArgumentsTo(funArcSin);
				return funArcSin;
			}

			virtual Computable *getSimplified() {
				FunArcSin *funArcSin = new FunArcSin();
				copySimplifiedArgumentsTo(funArcSin );
				return funArcSin ;
			}
	};

	class FunArcCos : public Function
	{
		public:
			virtual ~FunArcCos() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "arccos(@)"; }
			
			virtual Computable *clone() {
				FunArcCos *funArcCos = new FunArcCos();
				copyArgumentsTo(funArcCos);
				return funArcCos;
			}

			virtual Computable *getSimplified() {
				FunArcCos *funArcCos = new FunArcCos();
				copySimplifiedArgumentsTo(funArcCos);
				return funArcCos;
			}
	};

	class FunArcTg : public Function
	{
		public:
			virtual ~FunArcTg() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "arctg(@)"; }
			
			virtual Computable *clone() {
				FunArcTg *funArcTg = new FunArcTg();
				copyArgumentsTo(funArcTg);
				return funArcTg;
			}

			virtual Computable *getSimplified() {
				FunArcTg *funArcTg = new FunArcTg();
				copySimplifiedArgumentsTo(funArcTg);
				return funArcTg;
			}
	};

	class FunArcCtg : public Function
	{
		public:
			virtual ~FunArcCtg() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "arcctg(@)"; }
			
			virtual Computable *clone() {
				FunArcCtg *funArcCtg = new FunArcCtg();
				copyArgumentsTo(funArcCtg);
				return funArcCtg;
			}

			virtual Computable *getSimplified() {
				FunArcCtg *funArcCtg = new FunArcCtg();
				copySimplifiedArgumentsTo(funArcCtg);
				return funArcCtg;
			}
	};

	class FunArSh : public Function
	{
		public:
			virtual ~FunArSh() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "Arsh(@)"; }
			
			virtual Computable *clone() {
				FunArSh *funArSh = new FunArSh();
				copyArgumentsTo(funArSh);
				return funArSh;
			}

			virtual Computable *getSimplified() {
				FunArSh *funArSh = new FunArSh();
				copySimplifiedArgumentsTo(funArSh);
				return funArSh;
			}
	};

	class FunArCh : public Function
	{
		public:
			virtual ~FunArCh() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "Arch(@)"; }
			
			virtual Computable *clone() {
				FunArCh *funArCh = new FunArCh();
				copyArgumentsTo(funArCh);
				return funArCh;
			}

			virtual Computable *getSimplified() {
				FunArCh *funArCh = new FunArCh();
				copySimplifiedArgumentsTo(funArCh);
				return funArCh;
			}
	};

	class FunArTh : public Function
	{
		public:
			virtual ~FunArTh() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "Arth(@)"; }
			
			virtual Computable *clone() {
				FunArTh *funArTh = new FunArTh();
				copyArgumentsTo(funArTh);
				return funArTh;
			}

			virtual Computable *getSimplified() {
				FunArTh *funArTh = new FunArTh();
				copySimplifiedArgumentsTo(funArTh);
				return funArTh;
			}
	};

	class FunArCth : public Function
	{
		public:
			virtual ~FunArCth() { }
			virtual double getValue();
			virtual Computable *getDerivative();
			virtual int getArity() { return 1; }
			virtual const char *getMask() { return "Arcth(@)"; }
			
			virtual Computable *clone() {
				FunArCth *funArCth = new FunArCth();
				copyArgumentsTo(funArCth);
				return funArCth;
			}

			virtual Computable *getSimplified() {
				FunArCth *funArCth = new FunArCth();
				copySimplifiedArgumentsTo(funArCth);
				return funArCth;
			}
	};

	//////////////////////////////////////////////////////

	class Error : public Position
	{
		public:
			Error(Position pos) : Position(pos) { }
			virtual const char *message() = 0;
			virtual const char *information() = 0;
	};

	class ESyntax : public Error
	{
		public:
			ESyntax(Position pos) : Error(pos) { }
			virtual const char *message() { return "Синтаксическая ошибка"; }
			virtual const char *information() { return "Эта лексема не ожидается в данной позиции"; }
	};

	class EInvalidChar : public Error
	{
		public:
			EInvalidChar(Position pos) : Error(pos) { }
			virtual const char *message() { return "Некорректный символ"; }
			virtual const char *information() { return "Этот символ нельзя использовать"; }
	};

	class EEmptyExpression : public Error
	{
		public:
			EEmptyExpression() : Error(Position(0, 1)) { }
			virtual const char *message() { return "Пустое выражение"; }
			virtual const char *information() { return ""; }
	};

	class EUndefinedSymbol : public Error
	{
		public:
			EUndefinedSymbol(Position pos, const char *name) : Error(pos) {
				this->name = strdup(name);
				this->infoText = NULL;
			}

			virtual ~EUndefinedSymbol() {
				free(name);
				if (infoText) free(infoText);
			}

			virtual const char *message() { return "Символ не определен"; }

			virtual const char *information() {
				if (infoText) free(infoText);
				int len = snprintf(NULL, 0, "\"%s\"", name);
				infoText = (char*) malloc(sizeof(char) * len + 1);
				sprintf(infoText, "\"%s\"", name);
				return infoText;
			}
	
		protected:
			char *name;
			char *infoText;
	};

	class EParenthesesMismatchOpen : public Error
	{
		public:
			EParenthesesMismatchOpen(Position pos) : Error(pos) { }
			virtual const char *message() { return "Несоответствие круглых скобок"; }
			virtual const char *information() { return "Не хватает открывающей"; }
	};

	class EParenthesesMismatchClose : public Error
	{
		public:
			EParenthesesMismatchClose(Position pos) : Error(pos) { }
			virtual const char *message() { return "Несоответствие круглых скобок"; }
			virtual const char *information() { return "Не хватает закрывающей"; }
	};

	class ENotEnoughArgs : public Error
	{
		public:
			ENotEnoughArgs(Position pos, int expected, int got) : Error(pos) {
				this->expected = expected;
				this->got = got;
				infoText = NULL;
			}
			virtual ~ENotEnoughArgs() {
				if (infoText) free(infoText);
			}

			virtual const char *message() { return "Не хватает аргументов для оператора или функции"; }

			virtual const char *information() {
				if (infoText) free(infoText);
				int len = snprintf(NULL, 0, "Ожидается: %d, получено: %d", expected, got);
				infoText = (char*) malloc(sizeof(char) * len + 1);
				sprintf(infoText, "Ожидается: %d, получено: %d", expected, got); 
				return infoText;
			}

		protected:
			int expected;
			int got;
			char *infoText;
	};

	//////////////////////////////////////////////////////

	class Taylor : public Scope
	{
		public:
			Taylor();
			~Taylor();

			Error *parse(const char *expr, Computable *&to);
			Computable *getSeries(Computable *source, int degree,
					Link *x, Symbol *arg, Symbol *a);

			/** Создает список независимых разложений степеней 0..degree
			 */
			void getSeriesList(std::list<Computable*> &tsList, Computable *source, int degree,
					Link *x, Symbol *arg, Symbol *a);

			static void setAccuracy(int acc);
			static int getAccuracy();
			static void setDisplayAccuracy(int acc);
			static inline int getDisplayAccuracy();
			static void setExpoPres(bool ep);
			static inline bool getExpoPres();
			static inline bool iszero(double x);
			static inline bool isone(double x);

		private:
			static int accuracy;
			static double epsilon;
			static int displayAccuracy;
			static bool expoPres;

			Computable *getObject(const char *expr, int argCount = -1);
			Error *checkSyntax(std::list<Computable*> &infix);
	};

} // namespace Taylor


