/*
	nameTable.h

	Class for a string-indexed table.

	Brett Allen
	May, 2002
*/

#ifndef NAME_TABLE_H
#define NAME_TABLE_H

#include <vector>
#include <iostream>
using namespace std;

class NameTable {
protected:
	// should be protected, but there's a bug in VCC:
public:
	vector<char*> names;
	vector<void*> items;

public:
	char fieldName[256];
	bool shallowSave;
	int type;

	NameTable(char *name = NULL);
	NameTable(NameTable &nt);
	virtual ~NameTable();

	void add(const char *name, void *item);
	int lookupName(const char *name);
	char *getName(int i);
	void *&getItem(int i);
	void *&getItem(const char *name);

	int size() {
		return (int)names.size();
	}

	virtual void clear();
};

class NameTableInt : public NameTable {
public:
	NameTableInt(char *name = NULL) : NameTable(name) { }
	NameTableInt(NameTableInt &nt) : NameTable(nt) { }

	void addInt(const char *name, int item);
	int &getInt(int i);
	int &getInt(const char *name);

//	friend ostream& operator <<(ostream& os, const NameTableInt& v);
//	friend istream& operator >>( istream& is, NameTableInt& v );
};

class NameTableString : public NameTable {
public:
	NameTableString(char *name = NULL) : NameTable(name) { }
	NameTableString(NameTableString &nt);
	~NameTableString();

	void addString(const char *name, const char *item);
	char *getString(int i);
	char *getString(const char *name);

	static char *dupString(const char* s);
};

template <class X>
class NameTableX : public NameTable {
public:
	NameTableX(char *name = NULL) : NameTable(name) { }
	NameTableX(NameTableInt &nt) : NameTable(nt) { }
	virtual ~NameTableX() {
		int i;
		for (i=0; i < items.size(); i++) {
			delete ((X*)(items[i]));
		}
	}

	void addX(const char *name, X item);
	X &getX(int i);
	X &getX(const char *name);
	
	/*
#ifdef WIN32
	// VCC is non-standard
	friend ostream& operator <<(ostream& os, const NameTableX<X>& v);
	friend istream& operator >>( istream& is, NameTableX<X>& v );
#else
	friend ostream& operator << <>(ostream& os, const NameTableX<X>& v);
	friend istream& operator >> <>( istream& is, NameTableX<X>& v );
#endif
	*/
};


template <class T>
class NameTableT : public NameTable {
public:
	NameTableT(char *name = NULL) : NameTable(name) { }
	NameTableT(NameTableT<T> &nt) : NameTable(nt) { }

	void addT(const char *name, T* item) {
		add(name, (void*)item);
	}

	T *&getT(int i) {
		return (T*&)getItem(i);
	}

	T *&getT(const char *name) {
		return (T*&)getItem(name);
	}
};

ostream& operator <<( ostream& os, const NameTableInt& v );
istream& operator >>( istream& is, NameTableInt& v );

template <class X>
void NameTableX<X>::addX(const char *name, X item) {
	add(name, new X(item));
}

template <class X>
X &NameTableX<X>::getX(int i) {
	return *((X*)getItem(i));
}

template <class X>
X &NameTableX<X>::getX(const char *name) {
	return *((X*)getItem(name));
}

template<class X>
ostream& operator <<( ostream& os, const NameTableX<X>& v ) {
	os << v.fieldName << endl;

	int i;
	for (i=0; i < v.names.size(); i++) {
//		os << v.names[i] << " " << *((X*)v.items[i]) << endl;
	}

	return os;
}

template<class X>
istream& operator >>( istream& is, NameTableX<X>& v ) {
	char temp[256], *copy;

	// a little hack: if the field name has been set already,
	//  we assume we don't need to read it in.
	if (strlen(v.fieldName) < 2) {
		is >> temp;
		strcpy(v.fieldName, temp);
	}

	while (true) {
		is >> temp;
		if (strcmp(temp, "end") == 0)
			break;
		if (!is.good()) {
			cerr << "parser warning: couldn't find end of " << v.fieldName << endl;
			break;
		}
		copy = new char[strlen(temp)+1];
		strcpy(copy, temp);
		v.names.push_back(copy);

		X curT;
		is >> curT;
		v.items.push_back(new X(curT));
	}

	return is;
}

#endif
