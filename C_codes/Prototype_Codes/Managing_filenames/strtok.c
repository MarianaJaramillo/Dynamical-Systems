#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main ()
{
    char str[] ="Archivo.dat";
    const char s[5] = ".";
    char src[150];
    char dest[150];
    char *pt;

    strcpy(src, "_new.dat");

    pt = strtok (str, s);

    printf("%s\n", pt);

    strcpy(dest, pt);

    strcat(dest, src);

    printf("%s\n", dest);

    return 0;
}
