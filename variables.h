

void InitVars();
void ClearVars();
int lookupvar(char *name);
int LookupFilter(char *name, filter *F);
int LookupFilterSet(char *name, filterset *FS);
int LookupImage(char *name, image *I);
void AddFilter(char *name, filter F);
void AddFilterSet(char *name, filterset FS);
void AddImage(char *name, image I);
void RMVar(char *name);

void ListVars();
