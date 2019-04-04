/*
	File:			Array.cxx

	Function:		Template definitions for Array.h

	Author(s):		Andrew Willmott

	Copyright:		Copyright (c) 1995-1996, Andrew Willmott

	Notes:			

*/

#define __ArrayTmpl__
#include "Array.h"
#include <ctype.h>

TMPLArray TArray::Array(Int size, Int alloc) : items(size),
					allocated(alloc)
{
	Assert(size > 0, "(Array) Initial array size must be positive!");
	if (allocated < size)
		allocated = size;
	
	item = new T[allocated];
	Assert(item != 0, "(Array) Out of memory");
}

TMPLArray TArray::Array(const TArray &array) : items(array.items), 
	allocated(array.allocated)
{
	Int i;
	
	item = new T[allocated];
	Assert(item != 0, "(Array) Out of memory");
	
	for (i = 0; i < array.items; i++)
		item[i] = array.item[i];
}

TMPLArray TArray::~Array()
{
	delete[] item;
}

TMPLArray TArray &TArray::operator = (const TArray &array)
{
	Int i;

	if (allocated < array.allocated)
	{
		delete[] item;
		allocated = array.allocated;	
		item = new T[allocated];	
		Assert(item != 0, "(Array) Out of memory");
	}
			
	for (i = 0; i < array.items; i++)
		item[i] = array.item[i];

	items = array.items;
	
	return(SELF);
}

TMPLArray ostream &operator << (ostream &s, TArray &array)
{	
	Int		i;
	Char	sepChar;

	s << '[';
	if (array.NumItems() >= 16)
		sepChar = '\n';
	else
		sepChar = ' ';
	
	if (array.NumItems() > 0)
	{
		s << array[0];

		for (i = 1; i < array.NumItems(); i++)
			s << sepChar << array[i];
	}
	
	s << ']';

	return(s);
}

TMPLArray istream &operator >> (istream &s, TArray &array)
{
    Char	c;
	
	//	Expected format: [a b c d ...]
	
    while (isspace(s.peek()))			// 	chomp white space
		s.get(c);
		
    if (s.peek() == '[')						
    {
    	s.get(c);
    	array.Clear();
    	
	    while (isspace(s.peek()))		// 	chomp white space
			s.get(c);
    	
		while (s.peek() != ']')
		{			
			array.Add(1);
			s >> array.Top();			//	read an item
    	
			if (!s)
			{
				_Warning("Couldn't read array component");
				return(s);
			}
	
		    while (isspace(s.peek()))	// 	chomp white space
				s.get(c);
		}			
		s.get(c);
	}
    else
	{
	    s.clear(ios::failbit);
	    _Warning("Error: Expected '[' while reading array");
	    return(s);
	}
	
    return(s);
}

TMPLArray Void TArray::PreAllocate(Int newSize)
{
	Int	i;
	T	*newArray;
	
	if (newSize > allocated)
	{
		if (allocated == 0)
			allocated = kFirstAllocation;
		else
			allocated *= 2;	
		
		while (newSize > allocated)
			allocated *= 2;	
		
		newArray = new T[allocated];
		Assert(newArray != 0, "(Array::PreAllocate) Out of memory");
	
		for (i = 0; i < items; i++)
			newArray[i] = item[i];	
		
		delete[] item;
		item = newArray;
	}
}

TMPLArray Void TArray::SetSize(Int newSize)
{
	PreAllocate(newSize);
	items = newSize;
}

TMPLArray Void TArray::Add(Int n)
{
	SetSize(items + n);
}

TMPLArray Void TArray::Shrink(Int n)
//	take away n items.
{
	items -= n;
}

TMPLArray Void TArray::Insert(Int i, Int n)
//	Make space at position i for n items.
{
	Assert(i >= 0 && i <= items, "(Array:InsertSpace) Illegal index");

	Int j;
	
	Add(n);
	
	for (j = items - 1; j >= i + n; j--)
		item[j] = (item - n)[j];
}

TMPLArray Void TArray::Delete(Int i, Int n)
//	Delete n items at position i.
{
	Assert(i >= 0 && i <= items, "(Array:InsertSpace) Illegal index");

	Int j;
	
	items -= n;
		
	for (j = i; j < items; j++)
		item[j] = (item + n)[j];
}

TMPLArray Void TArray::ShrinkWrap()
//	Shrink allocated space to be only the current size of array
// There is no realloc version of new in C++, so this involves another copy.
{
	Int	i;
	T	*newArray;
	
	allocated = items;	
	
	newArray = new T[allocated];
	Assert(newArray != 0, "(Array::ShrinkWrap) Out of memory");

	for (i = 0; i < items; i++)
		newArray[i] = item[i];	
	
	delete[] item;
	item = newArray;
}

TMPLArray Void TArray::Grow()
//	Allocate more space for the array. Used internally prior to an items++.
{
	Int	i;
	T	*newArray;
	
	if (allocated == 0)
		allocated = kFirstAllocation;
	else
		allocated *= 2;	
	
	newArray = new T[allocated];
	Assert(newArray != 0, "(Array::Grow) Out of memory");

	for (i = 0; i < items; i++)
		newArray[i] = item[i];	
	
	delete[] item;
	item = newArray;
}

TMPLArray Void TArray::SwapWith(TArray &a)
{
	Int a1, b1;
	
	Swap(a1, b1);
	
	Swap(items, a.items);
	Swap(allocated, a.allocated);
	Swap(item, a.item);
}

TMPLArray Void TArray::Replace(TArray &a)
{
	delete[] item;
	item = a.item;
	items = a.items;
	allocated = a.allocated;

	a.item = 0;
}

#ifndef CL_NO_MEM_MAP

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

TMPLArray Void TArray::MemMap(const Char *filename)
{
	Int		fd;

	Clear();
	fd = open(filename, O_RDWR);
	
	if (fd != -1)
	{
		Int			fsize;
		struct stat	fstats;
		
		fstat(fd, &fstats);
		fsize = fstats.st_size;

		items = fsize / sizeof(T);
		Assert(items * sizeof(T) == fsize, "(Array::MemMap) bad file size");

		item = (T*) mmap(0,  fsize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
		allocated = -1;

		if (item == MAP_FAILED)
			_Error("MemMap failed");
		close(fd);
	}
}

#endif

#include <stdio.h>

TMPLArray Void TArray::WriteFile(const Char *filename)
{
	FILE	*file;

	file = fopen(filename, "wb");
	if (file)
	{
		fwrite(item, items, sizeof(T), file);
		fclose(file);
	}
}

TMPLArray Void TArray::ReadFile(const Char *filename)
{
	FILE	*file = fopen(filename, "rb");

	Clear();

	if (file)
	{
		Int		fsize;

		fseek(file, 0, SEEK_END);
		fsize = ftell(file);
		rewind(file);
		items = fsize / sizeof(T);
		Assert(items * sizeof(T) == fsize, "(Array::ReadFile) bad file size");
		item = new T[items];
		allocated = items;
		fread(item, items, sizeof(T), file);
		fclose(file);
	}
}
