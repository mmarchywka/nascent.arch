#ifndef MJM_CURSING_LIST_H__
#define MJM_CURSING_LIST_H__
#include "mjm_globals.h"
#include <vector>
#include <list>
#include <sstream>
#include <string>


/*
// probably be better off with own list... 
template <class Ty> class mjm_cursing_list : private std::list<Ty>
{
typedef mjm_cursing_list Myt;
typedef typename std::list<Ty> Super;
typedef typename Super::iterator Cursor;
typedef typename Super::iterator iterator;
typedef size_t IdxTy;
typedef Ty ValueTy;
iterator begin()  { return Super::begin(); }
iterator end()  { return Super::end(); }
iterator rbegin()  { return Super::rbegin(); }
iterator rend() { return Super::rend(); }


public:
mjm_cursing_list(): m_cursor(begin()),m_size(0) {}
ValueTy  & next() 
{
if (m_cursor!=end()) ++m_cursor();
if (m_cursor==end()) m_cursor=begin();
return **m_cursor; 
}
ValueTy  & prev() 
{
if (m_cursor!=rend()) --m_cursor();
if (m_cursor==rend()) m_cursor=begin();
return **m_cursor; 
}
ValueTy& front() { return Super::front(); }
const ValueTy & front() const { return Super::front(); }
ValueTy &  back() { return Super::back(); }
const ValueTy & back() const { return Super::back(); }

void push_front(const Ty & ty) { Super::push_front(ty); }
void pop_back() { Super::pop_back(); }

Myt & add(const Ty & ty)
{
if (m_cursor==end()) {push_back(ty); return *this;}
// insert Ty BEFORE iterator, not really good?
insert(m_cursor,ty); 
++m_size;
return *this;
}

Myt & remove()
{
erase(m_cursor);
--m_size;
return *this;

}
// this will be slow anyway lol. 
template <class Tx> void find(Tx (Ty::*ptr)() , const Tx & val,IdxTy flags)
{
const Cursor start=m_cursor;
if (size()==0) return; 
do {
if ( (&next())->ptr()>val) break;
} while (start!=m_cursor);

}


IdxTy size() const { return m_size; } 
// we should make our own size tracker but then need to 
// keep Super private. 
Cursor m_cursor;
IdxTy m_size;


}; // mjm_cursing_list
*/

template <class Ty> class mjm_cursing_list2
{
typedef unsigned int IdxTy;
template <class Tx> class entry
{
public:
typedef entry<Tx> * Eptr;
entry():next(0),prior(0){}
entry(const Tx & ty):value(ty),next(0),prior(0){}
Tx value;
Eptr next;
Eptr prior;

};

typedef entry<Ty> * Eptr;
Eptr first;
Eptr last;
Eptr cursor;
IdxTy m_size;
public:
typedef Ty ValueTy;
	mjm_cursing_list2():first(0),last(0),cursor(0),m_size(0) {}
	// ideally we just delete our block... 
	// this is also the ideal test case for situations where
	// malloc and delete patterns are easy to predict
	~mjm_cursing_list2() {while (m_size) pop_back(); }
	void clear() {while (m_size) pop_back(); }

	IdxTy size() const { return m_size; } 
	ValueTy& front() { return (*first).value; }
	const ValueTy & front() const { return (*first).value; }
	ValueTy &  back() { return (*last).value; }
	const ValueTy & back() const { return (*last).value; }

	ValueTy& frontprior() { return (*(*first).prior).value; }
	const ValueTy & frontprior() const { return (*(*first).prior).value; }
	ValueTy &  backnext() { return (*(*last).next).value; }
	const ValueTy & backnext() const { return (*(*last).next).value; }
	// the idea that containers are or have iterators may be confusing
	// but here it is likely the container will require exactly one or
	// a defined group of iterators. 
	ValueTy& current() { return (*cursor).value; }
	const ValueTy & current() const { return (*cursor).value; }

	void push_front(const Ty & ty) 
	{ 
		Eptr oldfront=first;
		// this needs to find a spot in block not call OS.
		first = new  entry<Ty>(ty);
		first->prior=oldfront;
		if (oldfront==0) last=first;
		else oldfront->next=first;
		++m_size;	
	}
	void pop_back() 
	{	if (m_size==0) return; // don't mess this up. 
		Eptr newlast=(*last).next;
		if (newlast!=0) newlast->prior=0; 
		else first=0;
		// ultimately we will manage our own blocks,not call OS 
		// for all this crap 
		delete last;   
		last=newlast; 
		--m_size;
	}

// this will be slow anyway lol. 
template <class Tx> void find(Tx (Ty::*ptr)() , const Tx & val,IdxTy flags)
{
if (size()==0) return; 
//const Cursor start=m_cursor;
Eptr x=cursor;
do {
//if ( (&next())->ptr()>val) break;
if ( (x->next)->ptr()>val) break;
//} while (start!=m_cursor);
	x=x->next;
} while (x!=last);

}

}; // mjm_cursing_list2


#endif

