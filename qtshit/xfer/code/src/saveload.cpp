/*
	saveload.cpp

	These classes enable convenient saving and loading of 
	data structures.

	Brett Allen
	May, 2002
*/

#include "saveload.h"
#include "vec.h"
#include "quatnorm.h"

SLProperty::SLProperty() : globalPtr(NULL), offset(0), offset2(0) {
}

void *SLProperty::getPtr(void *obj) {
	if (globalPtr)
		return (char*)globalPtr + offset;
	return (char*)obj + offset;
}

void *SLProperty::getPtr2(void *obj) {
	if (globalPtr)
		return (char*)globalPtr + offset2;
	return (char*)obj + offset2;
}


SLClassInfo *SaveLoad::getClassInfo(SLInterface *object) {
	if (object == NULL)
		return NULL;

	int index = classes.lookupName(object->className);
	if (index < 0) {
		cerr << "unknown class: '" << object->className << "'" << endl;
		return NULL;
	}

	return classes.getT(index);	
}

void SaveLoad::save(ostream &os, SLInterface *object) {
	SLClassInfo *info = getClassInfo(object);

	if (info == NULL) {
		os << "NULL" << endl;
		return;
	}

	os << object->className << endl;
	saveInstance(os, object, info);
}

void SaveLoad::loadInstance(istream &is, SLInterface *object, SLClassInfo *classInfo) {
	char name[255];
	int index;

	if (classInfo == NULL)
		classInfo = getClassInfo(object);

	if (object == NULL || classInfo == NULL)
		return;

	while (is.good()) {
		is >> name;

		if (strcmp(name, "end") == 0)
			break;

		index = classInfo->members.lookupName(name);
		if (index < 0) {
			cerr << "unknown property: '" << name << "'" << endl;
			continue;
		}

		loadProperty(is, object, classInfo->members.getT(index));
	}
}

void SaveLoad::saveInstance(ostream &os, SLInterface *object, SLClassInfo *classInfo) {
	int i;

	if (classInfo == NULL)
		classInfo = getClassInfo(object);

	if (object == NULL || classInfo == NULL) {
		os << "NULL" << endl;
		return;
	}

	for (i=0; i < classInfo->members.size(); i++) {
		os << classInfo->members.getName(i) << " ";
		saveProperty(os, object, classInfo->members.getT(i));
	}
	os << "end" << endl;
}

void SaveLoad::loadProperty(istream &is, void *object, SLProperty *prop) {
	int i, max;
	Vec3d *v;
	QuatNorm *q;
	void *items;
	SLClassInfo *info;
	vector<int> *vint;
	vector<SLInterface *> *vt;
	Char256 str;
	double d;

	switch (prop->type) {
	case SL_VOID:
		*((void**)prop->getPtr(object)) = load(is);
		break;

	case SL_INSTANCE:
		loadInstance(is, (SLInterface*)prop->getPtr(object));
		break;

	case SL_INT:
		int intVal;
		is >> intVal;
		*((int*)prop->getPtr(object)) = intVal;
		break;

	case SL_DOUBLE:
		is >> d;
		*((double*)prop->getPtr(object)) = d;
		break;

	case SL_ARRAY:
		is >> max;
		*((int*)prop->getPtr2(object)) = max;

		info = loadClassInfo(is);
		if (!is)
			break;

		items = info->newInstance(max);
		for (i=0; i < max; i++) {
			loadInstance(is, (SLInterface*)((char*)items + i * info->size), info);
		}
		*(void**)prop->getPtr(object) = items;

		break;

	case SL_VEC3D:
		v = ((Vec3d*)prop->getPtr(object));
		is >> (*v)[0] >> (*v)[1] >> (*v)[2];
		break;

	case SL_CHAR256:
		is >> ((char*)prop->getPtr(object));
		break;

	case SL_CHAR256_ARRAY:
		is >> max;
		*((int*)prop->getPtr2(object)) = max;

		items = new Char256[max];
		for (i=0; i < max; i++) {
			is >> ((char*)items + 256*i);
		}
		*(void**)prop->getPtr(object) = items;
		break;

	case SL_DOUBLE_ARRAY:
		is >> max;
		*((int*)prop->getPtr2(object)) = max;

		items = new double[max];
		for (i=0; i < max; i++) {
			is >> *((double*)items + i);
		}
		*(void**)prop->getPtr(object) = items;
		break;

	case SL_VECTOR_INT:
		vint = (vector<int>*)prop->getPtr(object);
		is >> max;
		vint->resize(max);
		for (i=0; i < max; i++)
			is >> (*vint)[i];
		break;

	case SL_QUATNORM:
		q = (QuatNorm*)prop->getPtr(object);
		is >> q->x >> q->y >> q->z >> q->w;
		break;

	case SL_NAMETABLE_T:
		is >> max;
		for (i=0; i < max; i++) {
			is >> str;
			((NameTable*)prop->getPtr(object))->add(str, load(is));
		}
		break;

	case SL_VECTOR_T:
		vt = (vector<SLInterface*>*)prop->getPtr(object);
		is >> max;
		vt->clear();
		for (i=0; i < max; i++)
			vt->push_back((SLInterface*)load(is));
		break;

	case SL_BOOL:
		is >> str;
		if (str[0] == '1' || str[0] == 'T' || str[0] == 't')
			*((bool*)prop->getPtr(object)) = true;
		else
			*((bool*)prop->getPtr(object)) = false;
		break;

	default:
		if (prop->type < 0 || prop->type >= classes.size()) {
			cerr << "bad property type: " << prop->type << endl;
			return;
		}
		loadInstance(is, (SLInterface*)prop->getPtr(object), classes.getT(prop->type));
	}
}

void SaveLoad::saveProperty(ostream &os, void *object, SLProperty *prop) {
	int i, max;
	SLClassInfo *info;
	char *type;
	Vec3d *v;
	QuatNorm *q;
	char *str;
	vector<int> *vint;
	vector<SLInterface *> *vt;
	NameTable *nt;

	switch (prop->type) {
	case SL_VOID:
		save(os, *(SLInterface**)prop->getPtr(object));
		break;

	case SL_INSTANCE:
		saveInstance(os, (SLInterface*)prop->getPtr(object));
		break;

	case SL_INT:
		os << *((int*)prop->getPtr(object)) << endl;
		break;

	case SL_DOUBLE:
		os << *((double*)prop->getPtr(object)) << endl;
		break;

	case SL_ARRAY:
		max = *((int*)prop->getPtr2(object));
		type = (*(SLInterface**)prop->getPtr(object))->className;
		info = classes.getT(type);
		if (info == NULL) {
			cerr << "unknown class: " << type << endl;
			return;
		}
		os << max << " " << type << endl;
		os << type << endl;
		for (i=0; i < max; i++) {
			saveInstance(os, (SLInterface*)((*(char**)prop->getPtr(object)) + i*info->size), info);
		}
		break;

	case SL_VEC3D:
		v = (Vec3d*)prop->getPtr(object);
		os << (*v)[0] << " " << (*v)[1] << " " << (*v)[2] << endl;
		break;

	case SL_CHAR256:
		os << ((char*)prop->getPtr(object)) << endl;
		break;

	case SL_DOUBLE_ARRAY:
		max = *((int*)prop->getPtr2(object));

		os << max << endl;
		for (i=0; i < max; i++) {
			os << (*((double**)prop->getPtr(object)))[i] << endl;
		}
		break;

	case SL_CHAR256_ARRAY:
		max = *((int*)prop->getPtr2(object));

		str = *((char**)prop->getPtr(object));
		os << max << endl;
		for (i=0; i < max; i++) {
			os << (str + 256*i) << endl;
		}
		break;

	case SL_VECTOR_INT:
		vint = (vector<int>*)prop->getPtr(object);
		max = vint->size();
		os << max << endl;
		for (i=0; i < max; i++)
			os << (*vint)[i] << " ";
		os << endl;
		break;

	case SL_QUATNORM:
		q = (QuatNorm*)prop->getPtr(object);
		os << q->x << " " << q->y << " " << q->z << " " << q->w << endl;
		break;

	case SL_NAMETABLE_T:
		nt = (NameTable*)prop->getPtr(object);
		max = nt->size();
		os << max << endl;
		for (i=0; i < max; i++) {
			os << nt->getName(i) << endl;
			if (nt->shallowSave || nt->getItem(i) == NULL)
				os << "NULL" << endl;
			else
				save(os, (SLInterface*)nt->getItem(i));
		}
		break;

	case SL_VECTOR_T:
		vt = (vector<SLInterface*>*)prop->getPtr(object);
		max = vt->size();
		os << max << endl;
		for (i=0; i < max; i++)
			save(os, (*vt)[i]);
		break;

	case SL_BOOL:
		if (*((bool*)prop->getPtr(object)))
			os << "1" << endl;
		else
			os << "0" << endl;
		break;

	default:
		if (prop->type < 0 || prop->type >= classes.size()) {
			cerr << "bad property type: " << prop->type << endl;
			return;
		}
		saveInstance(os, (SLInterface*)prop->getPtr(object), classes.getT(prop->type));
	}
}

SLClassInfo *SaveLoad::loadClassInfo(istream &is) {
	char name[255];
	is >> name;

	if (strcmp(name, "NULL") == 0)
		return NULL;

	int index = classes.lookupName(name);
	if (index < 0) {
		cerr << "unknown class: '" << name << "'" << endl;
		return NULL;
	}

	return classes.getT(index);
}

void *SaveLoad::load(istream &is) {
	SLClassInfo *info = loadClassInfo(is);
	if (!info)
		return NULL;

	SLInterface *ret = info->newInstance(1);
	loadInstance(is, ret, info);

	return ret;
}

void SaveLoad::saveNameTable(ostream &os, SLInterfaceNT *nt) {
	int i;

	for (i=0; i < nt->size(); i++) {
		os << nt->getName(i) << endl;
		save(os, nt->getT(i));
	}
	os << "end" << endl;
}

void SaveLoad::loadNameTable(istream &is, SLInterfaceNT *nt) {
	char name[255];
	SLInterface *item;

	while (is.good()) {
		is >> name;
		if (strcmp(name, "end") == 0 || !is.good())
			break;
	
		item = (SLInterface*)load(is);
		nt->add(name, item);
	}
}
