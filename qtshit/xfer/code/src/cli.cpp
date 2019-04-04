#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

#include "cli.h"
#include "ba.h"
#include "nameTable.h"

NameTable functionTable;
NameTableString varTable;


void registerFunction(cliFunc f, const char *name) {
	functionTable.add(name, (void*)f);
}

void processCommand(const char *command) {
	char cmdName[80];
	char curParam[80];
	int i;

	cout << "> " << command << endl;

	// check if this is a comment
	if (command[0] == '/')
		return;

	// extract the command name
	const char *params = extractString(command, cmdName, 80);

	// is this a special command?
	if (strcmp(cmdName, "repeat") == 0) {
		int rep = 0;
		params = extractInt(params, &rep);
		for (i=0; i < rep; i++) {
			cout << (i+1) << "/" << rep;
			processCommand(params);
		}
		return;
	}
	else if (strcmp(cmdName, "x") == 0) {
		params = extractString(params, curParam, 80);
		processFile(curParam);
		return;
	}
	else if (cmdName[0] == '!') {
		sprintf(curParam, "scripts/%s.txt", cmdName+1);
		processFile(curParam);
		return;
	}
	else if (strcmp(cmdName, "set") == 0) {
		char varName[80], value[80];

		params = extractString(params, varName, 80);
		params = extractString(params, value, 80);

		int oldIndex = varTable.lookupName(varName);
		if (oldIndex >= 0) {
			if (varTable.items[oldIndex])
				delete [](varTable.items[oldIndex]);
			varTable.items[oldIndex] = varTable.dupString(value);
		}
		else
			varTable.addString(varName, value);

		return;
	}

	// look up the command
	for (i=0; i < functionTable.size(); i++) {
		if (strncmp(cmdName, functionTable.getName(i), 80) == 0) {
			((cliFunc)(functionTable.getItem(i)))(params);
			return;
		}
	}

	cerr << "WARNING: unknown command '" << cmdName << "'" << endl;
}

void processFile(const char *fname) {
	ifstream in;
	if (!openIFStream(&in, fname, "script"))
		return;

	char str[1024];
	while (!in.eof()) {
		in.getline(str, 1023);
		if (strlen(str) > 0)
			processCommand(str);
	}
}

const char *extractString(const char *s, char *str, int maxLen) {
	int index = 0;
	int destIndex = 0;
	bool quoted = false;

	if (s[index] == '"') {
		quoted = true;
		index++;
	}

	for (; destIndex < maxLen-1; index++) {
		if (s[index] == 0 || s[index] == '\n' || s[index] == '\r') {
			if (str) str[destIndex] = 0;
			return s + index;
		}
		else if (quoted && s[index] == '"') {
			if (str) str[destIndex] = 0;
			index++;
			if (s[index] == ' ')
				index++;
			return s + index;
		}
		else if (!quoted && s[index] == ' ') {
			if (str) str[destIndex] = 0;
			index++;
			return s + index;
		}
		else if (s[index] == '{') {
			char varName[80];
			int i;
			for (i=0; i < 79; i++) {
				index++;
				if (s[index] == 0 || s[index] == '}')
					break;
				varName[i] = s[index];
			}
			if (s[index] == '0')
				index--;
			varName[i] = 0;

			char *s = varTable.getString(varName);
			if (s) {
//				cout << "|" << varName << "| -> |" << s << "|" << endl;
				for (i=0; i < strlen(s); i++) {
					if (str) str[destIndex++] = s[i];
					if (destIndex >= maxLen-1)
						break;
				}
			}
			else {
				cout << "WARNING: unknown variable, '" << varName << "'" << endl;
			}
		}
		else 
			if (str) str[destIndex++] = s[index];
	}

	if (str) str[destIndex] = 0;
	return s + index;
}

const char *extractDouble(const char *s, double *d) {
	char tempVal[80];
	const char *ret = extractString(s, tempVal, 80);
	
	if (d && tempVal[0] != 0)
		*d = atof(tempVal);

	return ret;
}

const char *extractInt(const char *s, int *i) {
	char tempVal[80];
	const char *ret = extractString(s, tempVal, 80);

	if (i && tempVal[0] != 0)
		*i = atoi(tempVal);

	return ret;
}

const char *extractBool(const char *s, bool *b) {
	char tempVal[80];
	const char *ret = extractString(s, tempVal, 80);

	if (b && tempVal[0] != 0)
		*b = (tempVal[0] != '0' && tempVal[0] != 'f' && tempVal[0] != 'F');

	return ret;
}

const char *cliLookupVar(const char *s) {
	return varTable.getString(s);
}