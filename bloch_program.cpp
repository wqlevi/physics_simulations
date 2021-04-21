#include "bloch_pTx.c"
#include <stdio.h>
#include <stdlib.h>
#include "data.cpp"

#define free2(a,b) free(a);free(b);
#define free3(a,b,c) free(a);free(b);free(c);

int main(int argc, char* argv[])
{

    const char *infile, *outfile;
    if (argc<2)
        infile = "/tmp/infile.dat";
    else
        infile = argv[1];
    if (argc < 3)
        outfile = "/tmp/outfile.dat";
    else
        outfile = argv[2];
    struct data mdata; // apply data.cpp to format binary files
    if(!mdata.read(infile))
        return -1;
    printf("md: %d, ntime: %d, nTx: %d, npos: %d, nf: %d prop:%d\n", mdata.mode, mdata.ntime, mdata.nTx, mdata.npos, mdata.nfreq, mdata.properties);

    /* This is where the magic happens */
    int ret = blochsimfz(mdata.b1r,mdata.b1i,mdata.xgrad,mdata.ygrad,mdata.zgrad,mdata.tsteps,mdata.ntime,mdata.t1,mdata.t2, \
                         mdata.dfreq,mdata.nfreq,mdata.dxpos,mdata.dypos,mdata.dzpos,mdata.npos,mdata.mx,mdata.my,mdata.mz, \
                         mdata.mode,mdata.b0,mdata.nTx,mdata.sr,mdata.si);
    if (ret != 0)
        return ret;

    /* Write resulting magnetization to file */
    mdata.write(outfile);
}
