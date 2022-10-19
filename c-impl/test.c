#define _POSIX_C_SOURCE 199309

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <poik.h>

////////////////////////////////////////////////////////////////

#include <stdint.h>

#include <unistd.h>
#include <time.h>

static double _timer(clockid_t which)
{
    struct timespec ts;
    clock_gettime(which, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}
#define wall_timer() _timer(CLOCK_REALTIME)
#define core_timer() _timer(CLOCK_PROCESS_CPUTIME_ID)

FILE *Bench_Progress_Indicator_File = NULL;
static signed int BENCH_LOOPS = 1;

double gen_wall = 0;
double gen_core = 0;
double prf_wall = 0;
double prf_core = 0;
double ver_wall = 0;
double ver_core = 0;


////////////////////////////////////////////////////////////////

static void render_progress(signed int idx)
{
    if (!Bench_Progress_Indicator_File)
        return;

    if (idx == -1) {
        fputs("\r\x1b[K", Bench_Progress_Indicator_File);
        fflush(Bench_Progress_Indicator_File);
        return;
    }

    static size_t last = -1;
    unsigned const done = 64 * idx / BENCH_LOOPS;
    if (done == last)
        return;

    char buf[65] = {0};
    memset(buf, '.', 64);
    memset(buf, '=', done);
    buf[done] = '>';
    fprintf(Bench_Progress_Indicator_File, "\r\x1b[K%s", buf);
    fflush(Bench_Progress_Indicator_File);
    last = done;
}


////////////////////////////////////////////////////////////////

static void measure(char const *name)
{
    static double t0_wall, t0_core;

    if (!name) {    // pre
        t0_wall = wall_timer();
        t0_core = core_timer();
    }
    else {          // post
        double const wall = wall_timer() - t0_wall;
        double const core = core_timer() - t0_core;

        if (*name == 'G') {
            gen_wall += wall;
            gen_core += core;
        } else if (*name == 'P') {
            prf_wall += wall;
            prf_core += core;
        } else if (*name == 'V') {
            ver_wall += wall;
            ver_core += core;
        }
    }
}
#define pre() measure(NULL)
#define post(S) measure(S)

////////////////////////////////////////////////////////////////

char const *const msg[2] = {" \x1b[1;31mNO\x1b[0m", "\x1b[1;32mYES\x1b[0m"};

bool test_poik(packed_fp2 const *startE)
{
    Progress_Indicator_File = stderr;

    fprintf(stderr, "\n\x1b[35m      ╔═══════════════════════════════════════╗\n      ║ ");
    fprintf(stderr, "\x1b[0;1;4;36mProof of Isogeny Knowledge\x1b[0m — \x1b[1;33m%u bits\x1b[0m", POIK_BITS);
    fprintf(stderr, "\x1b[35m ║\n      ╚═══════════════════════════════════════╝\x1b[0m\n");

fputc('\n', stderr);
fprintf(stderr, "%11s: \x1b[34m%9lu\x1b[0m bytes\n", "secret", sizeof(struct secret));
fprintf(stderr, "%11s: \x1b[34m%9lu\x1b[0m bytes\n", "proof", sizeof(struct proof));

    bool ok = true, res;

    struct secret sec;

fputc('\n', stderr);
    ok &= res = Generate(startE, &sec);


fputc('\n', stderr);
    fprintf(stderr, "       ║ could generate secret?            %s\n", msg[res]);
    if (!ok) return ok;


    struct proof *prf = malloc(sizeof(*prf));
    if (!(ok &= !!prf)) {
        perror("could not allocate proof structure");
        return ok;
    }

    Prove(&sec, prf);

    ok &= res = Verify(prf);

    fprintf(stderr, "       ║ proof verifies?                   %s\n", msg[res]);

    345[(char*)prf] ^= 1;
    ok &= !(res = Verify(prf));
    fprintf(stderr, "       ║ proof fails after bitflip?        %s\n", msg[!res]);
fputc('\n', stderr);
fputc('\n', stderr);

    return ok;
}


bool bench_poik(packed_fp2 const *startE)
{
    Progress_Indicator_File = NULL;
    Bench_Progress_Indicator_File = stderr;

    bool ok = true, res;

    struct secret sec;
pre();
    ok &= res = Generate(startE, &sec);
post("Generate()");

    if (!ok) return ok;

    struct proof *prf = malloc(sizeof(*prf));
    if (!(ok &= !!prf)) {
        perror("could not allocate proof structure");
        return ok;
    }

pre();
    Prove(&sec, prf);
post("Prove()");

pre();
    ok &= res = Verify(prf);
post("Verify()");

    return ok;
}

////////////////////////////////////////////////////////////////

packed_fp2 const E0 = {0};   // y^2 = x^3 + x

int main()
{

    bool ok = test_poik(&E0);
    if (!ok) return 1;

    fprintf(stderr, "Benchmarking across %d iterations...\n", BENCH_LOOPS);
    // fflush(stderr);

    for (signed int i = 0; i < BENCH_LOOPS; ++i) {
        render_progress(i);
        bench_poik(&E0);
    }
    render_progress(-1);
    
    fputc('\n', stderr);
    fprintf(stderr, "%11s: \x1b[34m%9.4lf\x1b[0m wall seconds / \x1b[34m%9.4lf\x1b[0m core seconds / ratio: \x1b[34m%5.2lf\x1b[0m\n", "Generate()", gen_wall/BENCH_LOOPS, gen_core/BENCH_LOOPS, gen_core / gen_wall);
    fprintf(stderr, "%11s: \x1b[34m%9.4lf\x1b[0m wall seconds / \x1b[34m%9.4lf\x1b[0m core seconds / ratio: \x1b[34m%5.2lf\x1b[0m\n", "Prove()", prf_wall/BENCH_LOOPS, prf_core/BENCH_LOOPS, prf_core / prf_wall);
    fprintf(stderr, "%11s: \x1b[34m%9.4lf\x1b[0m wall seconds / \x1b[34m%9.4lf\x1b[0m core seconds / ratio: \x1b[34m%5.2lf\x1b[0m\n", "Verify()", ver_wall/BENCH_LOOPS, ver_core/BENCH_LOOPS, ver_core / ver_wall);

    return 0;
}
