/*
	File:			Array.h

	Function:		Defines an array type that manages its own storage space,
					and can be used as a stack or a list.
					
	Author(s):		Andrew Willmott

	Copyright:		Copyright (c) 1995-1996, 1998, Andrew Willmott
 */

#ifndef __Array__
#define __Array__

#include <iostream>
#include "Basics.h"

using namespace std;

#define TMPLArray	template<class T>
#define TArray		Array<T>

const Int kFirstAllocation = 16; // Default number of items to initially 
								 // allocate to the array
/*
	Note: We grow the array in size exponentially. That is, every time we
	need to increase the size of the array, we double its size rather than
	increasing it by a fixed amount. For appending n items to an empty list,
	this gives us O(n) copies, as opposed to the O(n^2) copies required if
	we increment by a fixed size each time. 
	
	It would be useful to have an array data structure that utilises a tree of
	fixed-size arrays, thus trading off access time (const vs. o(logbn))
	against eliminating copies. 
*/


TMPLArray class Array
{
public:
					Array();
					Array(Int size, Int alloc = kFirstAllocation);		
					Array(const TArray &array);
				   ~Array();
	
//	Array operators
	
	inline T		&operator [] (Int i);		// indexing operator
	inline const T	&operator [] (Int i) const; // indexing operator
	inline Int		NumItems() const;			// Number of items in the array
	
	TArray			&operator = (const TArray &array);	// Assignment!
	
//	Useful for stack implementations

	inline T		&Top();						// Return top of stack
	inline Void		Pop();						// Delete top of stack
	inline Void		Push(const T &t);			// Push item onto stack  
	
// List Operations

	inline Void		Append(const T &t);			// Append single item to array
	inline T		&Last();					// Return last item in array
	Void			Clear();					// Delete all items

	Void			PreAllocate(Int numItems);	// Preallocate space for array
	Void			SetSize(Int newSize);		// Set array size directly.
	Void			Add(Int n = 1);				// Add n items to the array
	Void			Shrink(Int n = 1);			// shrink the array by n items
	Void 			Insert(Int i, Int n = 1);	// Insert n items at i
	Void			Delete(Int i, Int n = 1);	// Delete n items at i
	Void			ShrinkWrap();				// Ensure allocated space =
												// space being used
	Void			SwapWith(TArray &a);		// swaps this array with a
	Void			Replace(TArray &a);			// replace this array with a
												// & clear a.

	const T			&Item(Int i) const
					{ return(SELF[i]); };
	T				&Item(Int i)
					{ return(SELF[i]); };
	
// Low level access

	inline T		*Ref() const;				// Return pointer to array
	inline T		*Detach();					// As above, but the array 
												// no longer owns the data.
	Void			MemMap(const Char* filename);
	Void 			WriteFile(const Char *filename);
	Void 			ReadFile(const Char *filename);

//	Private...

protected:
	T				*item;		// pointer to array
	UInt32 			items;		// items in the array
	UInt32			allocated;	// number of items we have space allocated for. 
	
	Void 			Grow();
};	

TMPLArray ostream &operator << (ostream &s, TArray &array);
TMPLArray istream &operator >> (istream &s, TArray &array);


// --- Inlines ----------------------------------------------------------------


TMPLArray inline TArray::Array() : items(0), item(0), allocated(0)
{
}

TMPLArray inline Int TArray::NumItems() const
{
	return(items);
}

TMPLArray inline T &TArray::operator [] (Int i)
{
	CheckRange(i, 0, items, "(Array::[]) index out of range");

	return(item[i]);
}

TMPLArray inline const T &TArray::operator [] (Int i) const
{
	CheckRange(i, 0, items, "(Array::[]) index out of range");

	return(item[i]);
}

TMPLArray inline T &TArray::Top()
{
	return(item[items - 1]);
}

TMPLArray inline T &TArray::Last()
{
	return(item[items - 1]);
}

TMPLArray inline Void TArray::Push(const T &t)
{
	if (items >= allocated)
		Grow();
	
	item[items++] = t;
}

TMPLArray inline Void TArray::Append(const T &t)
{
	if (items >= allocated)
		Grow();
	
	item[items++] = t;
}

TMPLArray inline Void TArray::Pop()
{	
	items--;
}

TMPLArray inline Void TArray::Clear()
{	
	items = 0;
	allocated = 0;
	delete[] item;
	item = 0;
}

TMPLArray inline T *TArray::Ref() const
{
	return(item);
}

TMPLArray inline T *TArray::Detach()
{
	T* result = item;

	items = 0;
	allocated = 0;
	item = 0;

	return(result);
}


#ifndef __ArrayTmpl__
#include "Array.cxx"
#endif

#endif
