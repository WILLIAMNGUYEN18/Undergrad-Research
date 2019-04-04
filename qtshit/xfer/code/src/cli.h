#ifndef CLI_H
#define CLI_H

typedef void (*cliFunc)(const char*);

void registerFunction(cliFunc f, const char *name);
void processCommand(const char *command);
void processFile(const char *fname);

const char *extractString(const char *s, char *str, int maxLen = 1024);
const char *extractDouble(const char *s, double *d);
const char *extractInt(const char *s, int *i);
const char *extractBool(const char *s, bool *b);

const char *cliLookupVar(const char *s);

#endif