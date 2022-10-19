
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <poik.h>

int main(int argc, char **argv)
{
    Progress_Indicator_File = stderr;

    packed_fp2 E0 = {0};

    if (argc <= 1 || strcmp(argv[1], "--initial")) {
        if (!Parse_Fp2(stdin, &E0)) {
            fputs("please provide curve on stdin", stderr);
            return -1;
        }
    }

    struct secret sec;
    if (!Generate(&E0, &sec)) {
        fputs("invalid curve\n", stderr);
        return 1;
    }

    struct proof *prf = malloc(sizeof(*prf));
    if (!prf) {
        fputs("could not allocate proof structure\n", stderr);
        return 2;
    }

    Prove(&sec, prf);

    if (!Dump(stdout, prf)) {
        fputs("could not write proof to output file\n", stderr);
        free(prf);
        return 3;
    }

    free(prf);
    return 0;
}

