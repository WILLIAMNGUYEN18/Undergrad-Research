/*
	doppel_config.h

	A class for storing the current configuration of the doppel
	tool.

	Brett Allen
	May, 2002
*/

#ifndef DOPPEL_CONFIG_H
#define DOPPEL_CONFIG_H

#include "doppel2.h"
#include "saveload.h"

class DoppelConfig : public SLInterface {
public:
	DoppelConfig();

	void copyToGlobal();
	void copyFromGlobal();

	void save(char *fname);
	bool load(char *fname);

	static void registerClass(SaveLoad *sl);
	static SLInterface *newInstance(int);
};

extern DoppelConfig config;

#endif