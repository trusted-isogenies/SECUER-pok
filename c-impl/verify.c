
#include <stdlib.h>
#include <stdio.h>

#include <poik.h>

int main()
{
    Progress_Indicator_File = stderr;

    struct proof *prf = malloc(sizeof(*prf));
    if (!prf) {
        fputs("could not allocate proof structure\n", stderr);
        return 1;
    }

    if (!Parse(stdin, prf)) {
        fputs("please provide proof on stdin\n", stderr);
        free(prf);
        return -1;
    }

    if (!Verify(prf)) {
        fputs("invalid proof\n", stderr);
        free(prf);
        return 2;
    }

    if (!Dump_Fp2(stdout, &prf->E0) \
     || !Dump_Fp2(stdout, &prf->E1)) {
        fputs("could not write curve to output file\n", stderr);
        free(prf);
        return 3;
    }

    free(prf);
    return 0;
}

