/*
	nameTable.cpp

	Class for a string-indexed table.

	Brett Allen
	May, 2002
*/

#include "nameTable.h"


NameTable::NameTable(char *name) {
	if (name == NULL)
		name = "?";
	strcpy(fieldName, name);
	shallowSave = false;
	type = 0;
}

NameTable::NameTable(NameTable &nt) {
	strcpy(fieldName, nt.fieldName);
	int i;
	for (i=0; i < nt.names.size(); i++) {
		names.push_back(strdup(nt.names[i]));
		items.push_back(nt.items[i]);
	}
}

NameTable::~NameTable() {
	int i;
	for (i=0; i < names.size(); i++) {
		free(names[i]);
	}
}

void NameTable::add(const char *name, void *item) {
	if (name == NULL)
		name = "";
	names.push_back(strdup(name));
	items.push_back(item);
}

int NameTable::lookupName(const char *name) {
	int i;

	// inefficient, but what the hey
	for (i=0; i < names.size(); i++) {
		if (strcmp(name, names[i]) == 0)
			return i;
	}
	return -1;
}

char *NameTable::getName(int i) {
	return names[i];
}

void *&NameTable::getItem(int i) {
	return items[i];
}

void *&NameTable::getItem(const char *name) {
	static void *n = NULL;

	int i = lookupName(name);
	if (i >= 0)
		return items[i];
	else
		return n;
}

void NameTable::clear() {
	int i;
	for (i=0; i < names.size(); i++)
		free(names[i]);
	names.clear();
	items.clear();
}


void NameTableInt::addInt(const char *name, int item) {
	add(name, (void*)item);
}

int &NameTableInt::getInt(int i) {
	return (int&)getItem(i);
}

int &NameTableInt::getInt(const char *name) {
	return (int&)getItem(name);
}

ostream& operator <<( ostream& os, const NameTableInt& v ) {
	os << v.fieldName << endl;

	int i;
	for (i=0; i < v.names.size(); i++) {
		os << v.names[i] << " " << ((int)v.items[i]) << endl;
	}

	return os;
}

istream& operator >>( istream& is, NameTableInt& v ) {
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

		int curT;
		is >> curT;
		v.items.push_back((void*)curT);
	}

	return is;
}



NameTableString::NameTableString(NameTableString &nt) : NameTable(nt) {
	int i;
	for (i=0; i < items.size(); i++) {
		items[i] = dupString((char*)items[i]);
	}
}

NameTableString::~NameTableString() {
	int i;
	for (i=0; i < items.size(); i++)
		if (items[i])
			delete []items[i];
}

void NameTableString::addString(const char *name, const char *item) {
	add(name, dupString(item));
}

char *NameTableString::getString(int i) {
	return (char*)getItem(i);
}

char *NameTableString::getString(const char *name) {
	return (char*)getItem(name);
}

char *NameTableString::dupString(const char* s) {
	char *ret = new char[strlen(s)+1];
	strcpy(ret, s);
	return ret;
}


/*
template<class T>
ostream& operator <<( ostream& os, const NameTable<T>& v ) {
	os << v.fieldName << endl;

	int i;
	for (i=0; i < v.names.size(); i++) {
		os << v.names[i] << " " << v.items[i] << endl;
	}

	return os;
}

template<class T>
istream& operator >>( istream& is, NameTable<T>& v ) {
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

		T curT;
		is >> curT;
		v.items.push_back(curT);
	}

	return is;
}
*/