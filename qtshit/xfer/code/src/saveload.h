/*
	saveload.h

	These classes enable convenient saving and loading of 
	data structures.

	Brett Allen
	May, 2002
*/

#ifndef SAVELOAD_H
#define SAVELOAD_H

#include "nameTable.h"

const int SL_VOID          = -1;
const int SL_INT           = -2;
const int SL_NAMED_INT     = -3;
const int SL_ARRAY         = -4;
const int SL_VEC3D         = -5;
const int SL_CHAR256       = -6;
const int SL_CHAR256_ARRAY = -7;
const int SL_VECTOR_INT    = -8;
const int SL_QUATNORM      = -9;
const int SL_DOUBLE        = -10;
const int SL_NAMETABLE_T   = -11;
const int SL_VECTOR_T      = -12;
const int SL_DOUBLE_ARRAY  = -13;
const int SL_INSTANCE      = -14;
const int SL_BOOL          = -15;

typedef char Char256[256];

class SLProperty {
public:
	int type;
	void *globalPtr;
	int offset;
	int offset2;

	SLProperty();
	void *getPtr(void *obj);
	void *getPtr2(void *obj);
};

class SaveLoad;

class SLInterface {
public:
	char *className;

	virtual void alteredData() { }
};

class SLClassInfo {
public:
	int size;
	NameTableT <SLProperty> members;
	SLInterface *(*newInstance)(int arrayLen);
/*
	SLClassInfo() { }
	SLClassInfo(SLClassInfo &ci) : size(ci.size), members(ci.members) {
		newInstance = ci.newInstance;
	}*/
};

typedef NameTableT<SLInterface> SLInterfaceNT;

class SaveLoad {
private:
	SLClassInfo *loadClassInfo(istream &is);

	void loadProperty(istream &is, void *object, SLProperty *prop);
	void saveProperty(ostream &os, void *object, SLProperty *prop);

public:
	NameTableT<SLClassInfo> classes;

	SLClassInfo *getClassInfo(SLInterface *object);

	void save(ostream &os, SLInterface *object);
	void *load(istream &is);

	void saveNameTable(ostream &os, SLInterfaceNT *nt);
	void loadNameTable(istream &is, SLInterfaceNT *nt);

	void loadInstance(istream &is, SLInterface *object, SLClassInfo *classInfo = NULL);
	void saveInstance(ostream &os, SLInterface *object, SLClassInfo *classInfo = NULL);
};

extern SaveLoad saveLoad;

#endif
