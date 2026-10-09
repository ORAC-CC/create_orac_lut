/* Harness around ORAC's unmodified common/nc_get_string_att.c: argv = file variable attribute */
#include <stdio.h>
void nc_get_string_att(const char *fname, const char *vname, const char *aname, char str[128]);
int main(int argc, char **argv) {
    char str[128] = {0};
    if (argc != 4) { fprintf(stderr, "usage: %s file variable attribute\n", argv[0]); return 2; }
    nc_get_string_att(argv[1], argv[2], argv[3], str);
    printf("%s:%s = \"%s\"\n", argv[2], argv[3], str);
    return 0;
}
