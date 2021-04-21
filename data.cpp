#include <stdio.h>
#include <stdlib.h>
#define tfread(a, b, c, d)                                       \
    if (fread(a, b, c, d) != (size_t)c)                          \
    {                                                            \
        printf("fread failed in line %d. Exiting.\n", __LINE__); \
        return (false);                                              \
    };
#define tfwrite(a, b, c, d)                                       \
    if (fwrite(a, b, c, d) != (size_t)c)                          \
    {                                                            \
        printf("fwrite failed in line %d. Exiting.\n", __LINE__); \
        return (false);                                              \
    };

struct data                                 // define struct data
{
    bool read(const char *infile);
    bool write(const char *outfile);
    ~data();                                // ?
    data& operator=(const data& a);         // ?

    double *xgrad = NULL, *ygrad = NULL, *zgrad = NULL;
    double *tsteps = NULL;
    int nfreq, npos, ntime, nTx, properties;
    double t1, t2;
    double *dfreq = NULL;
    double *dxpos = NULL, *dypos = NULL, *dzpos = NULL;
    double *mx = NULL, *my = NULL, *mz = NULL;
    double *b0 = NULL;
    double *b1r = NULL, *b1i = NULL;
    int mode;
    double *sr = NULL, *si = NULL;
    double *target = NULL;
    int iter;
    double *targetmag = NULL;
};

data::~data(){                  // release memory of variables
    free(xgrad);
    free(ygrad);
    free(zgrad);
    free(tsteps);
    free(dfreq);
    free(dxpos);
    free(dypos);
    free(dzpos);
    free(mx);
    free(my);
    free(mz);
    free(b0);
    free(b1r);
    free(b1i);
    free(sr);
    free(si);
    free(targetmag);
}

bool data::read(const char *infile){
        FILE *fh;                                 // define file type pointer
        fh = fopen(infile, "rb");
        if (fh == NULL)
        {
            printf("ERROR:  input file couldn't be opened. \n");
            return false;
        }
        iter = 0;
        tfread(&(mode), (size_t)4, 1, fh);        // read first 4 byte of a element from fh and store to var &mode
        tfread(&(ntime), (size_t)4, 1, fh);
        tfread(&(nTx), (size_t)4, 1, fh);
        tfread(&(npos), (size_t)4, 1, fh);
        tfread(&(nfreq), (size_t)4, 1, fh);
        tfread(&(properties), (size_t)4, 1, fh);
        fseek(fh, (size_t)8 * 4, SEEK_SET);       // set cursor at position with (size_t)8*4 offset from SEEK_SET(beginning of file)

        tfread(&(t1), (size_t)8, 1, fh);
        tfread(&(t2), (size_t)8, 1, fh);

        // read b1
        free(b1r);free(b1i);
        b1r = (double *)malloc(ntime * nTx * sizeof(double));
        b1i = (double *)malloc(ntime * nTx * sizeof(double));
        tfread(b1r, (size_t)8, ntime * nTx, fh);
        tfread(b1i, (size_t)8, ntime * nTx, fh);

        // read g
        //double *gx,*gy,*gz;
        free(xgrad);free(ygrad);free(zgrad);
        xgrad = (double *)malloc(ntime * sizeof(double));
        ygrad = (double *)malloc(ntime * sizeof(double));
        zgrad = (double *)malloc(ntime * sizeof(double));
        tfread(xgrad, (size_t)8, ntime, fh);
        tfread(ygrad, (size_t)8, ntime, fh);
        tfread(zgrad, (size_t)8, ntime, fh);
        // read b0
        free(b0);
        b0 = (double *)malloc(npos * sizeof(double));
        tfread(b0, (size_t)8, npos, fh);
        // read tp
        free(tsteps);
        tsteps = (double *)malloc(ntime * sizeof(double));
        tfread(tsteps, (size_t)8, ntime, fh);
        // read nf
        free(dfreq);
        dfreq = (double *)malloc(nfreq * sizeof(double));
        tfread(dfreq, (size_t)8, nfreq, fh);
        // read pos
        free(dxpos);
        free(dypos);
        free(dzpos);
        dxpos = (double *)malloc(npos * sizeof(double));
        dypos = (double *)malloc(npos * sizeof(double));
        dzpos = (double *)malloc(npos * sizeof(double));
        tfread(dxpos, (size_t)8, npos, fh);
        tfread(dypos, (size_t)8, npos, fh);
        tfread(dzpos, (size_t)8, npos, fh);

        // read mx,my,mz
        unsigned int ntout;
        if (mode & 2)
            ntout = ntime; /* Include time points.	*/
        else
            ntout = 1;

        unsigned int ntnfnpos = ntout * nfreq * npos;
        //double *mx,*my,*mz;
        free(mx);free(my);free(mz);
        mx = (double *)malloc(ntnfnpos * sizeof(double));
        my = (double *)malloc(ntnfnpos * sizeof(double));
        mz = (double *)malloc(ntnfnpos * sizeof(double));
        tfread(mx, (size_t)8, ntnfnpos, fh);
        tfread(my, (size_t)8, ntnfnpos, fh);
        tfread(mz, (size_t)8, ntnfnpos, fh);

        // read Tx coil sensitivities
        free(sr);free(si);
        sr = (double *)malloc(npos * nTx * sizeof(double));
        si = (double *)malloc(npos * nTx * sizeof(double));
        tfread(sr, (size_t)8, nTx * npos, fh);
        tfread(si, (size_t)8, nTx * npos, fh);
        bool has_target = (properties >> 0) & 0x1;  // from here on:[06:03 am 03/April/2021]
        target = (double*) malloc(sizeof(double)*npos);
        if (has_target){
            tfread(target, (size_t)8, npos, fh);
        } else {
            for (int i=0;i<npos;i++)
                target[i] = 7;
        }

        free(targetmag);
        targetmag = (double*) malloc(sizeof(double)*npos);
        for (int i=0;i<npos;i++) {
            targetmag[i] = sin(target[i]*M_PI/180);
        }

        fclose(fh);
        return true;

    }

    bool data::write(const char *outfile){
        FILE *fh;
        fh = fopen(outfile, "wb");
        if (fh == NULL)
        {
            printf("ERROR:  output file couldn't be opened. \n");
            return false;
        }
        tfwrite(&(mode), (size_t)4, 1, fh);    // pointer to array to be written
        tfwrite(&(ntime), (size_t)4, 1, fh);
        tfwrite(&(nTx), (size_t)4, 1, fh);
        tfwrite(&(npos), (size_t)4, 1, fh);
        tfwrite(&(nfreq), (size_t)4, 1, fh);
        tfwrite(&(properties), (size_t)4, 1, fh);
        // fill up 2x any number, we use mode for simplicity.
        // the empty bytes
        tfwrite(&(mode), (size_t)4, 1, fh);
        tfwrite(&(mode), (size_t)4, 1, fh);
        // write relax times
        tfwrite(&(t1), (size_t)8, 1, fh);
        tfwrite(&(t2), (size_t)8, 1, fh);

        // write b1
        tfwrite(b1r, (size_t)8, ntime * nTx, fh);
        tfwrite(b1i, (size_t)8, ntime * nTx, fh);

        // write g
        //double *gx,*gy,*gz;
        tfwrite(xgrad, (size_t)8, ntime, fh);
        tfwrite(ygrad, (size_t)8, ntime, fh);
        tfwrite(zgrad, (size_t)8, ntime, fh);
        // write b0
        tfwrite(b0, (size_t)8, npos, fh);
        // write tp
        tfwrite(tsteps, (size_t)8, ntime, fh);
        // write nf
        tfwrite(dfreq, (size_t)8, nfreq, fh);
        // write pos
        tfwrite(dxpos, (size_t)8, npos, fh);
        tfwrite(dypos, (size_t)8, npos, fh);
        tfwrite(dzpos, (size_t)8, npos, fh);

        // read mx,my,mz
        unsigned int ntout;
        if (mode & 2) //anything but mode=1 nor 0
            ntout = ntime; /* Include time points.	*/
        else
            ntout = 1;

        unsigned int ntnfnpos = ntout * nfreq * npos;
        //double *mx,*my,*mz;
        tfwrite(mx, (size_t)8, ntnfnpos, fh);
        tfwrite(my, (size_t)8, ntnfnpos, fh);
        tfwrite(mz, (size_t)8, ntnfnpos, fh);

        // read Tx coil sensitivities
        tfwrite(sr, (size_t)8, nTx * npos, fh);
        tfwrite(si, (size_t)8, nTx * npos, fh);

bool has_target = (properties >> 0) & 0x1;
        if (has_target)
            tfwrite(target, (size_t)8, npos, fh);


        fclose(fh);
        return true;
    }

    data& data::operator=(const data& a)   // OK, but there is a cost
{
    if (this == &a) return *this;
    nfreq = a.nfreq;
    npos = a.npos;
    ntime = a.ntime;
    nTx = a.nTx;
    properties = a.properties;
    t1 = a.t1;
    t2 = a.t2;
    mode = a.mode;
    iter = a.iter;
    free(xgrad);xgrad = (double*)malloc(sizeof(double) * ntime);
    free(ygrad);ygrad = (double*)malloc(sizeof(double) * ntime);
    free(zgrad);zgrad = (double*)malloc(sizeof(double) * ntime);
    free(tsteps);tsteps = (double*)malloc(sizeof(double) * ntime);
    for(int i=0;i<ntime;i++){
        xgrad[i] = a.xgrad[i];
        ygrad[i] = a.ygrad[i];
        zgrad[i] = a.zgrad[i];
        tsteps[i] = a.tsteps[i];
    }
    free(b1r);b1r = (double*)malloc(sizeof(double) * ntime*nTx);
    free(b1i);b1i = (double*)malloc(sizeof(double) * ntime*nTx);
    for(int i=0;i<ntime*nTx;i++){
        b1r[i] = a.b1r[i];
        b1i[i] = a.b1i[i];
    }
    free(b0);b0 = (double*)malloc(sizeof(double) * npos);
    free(dxpos);dxpos = (double*)malloc(sizeof(double) * npos);
    free(dypos);dypos = (double*)malloc(sizeof(double) * npos);
    free(dzpos);dzpos = (double*)malloc(sizeof(double) * npos);
    free(target);target = (double*)malloc(sizeof(double) * npos);
    free(targetmag);
    if (a.targetmag != NULL)
        targetmag = (double*)malloc(sizeof(double) * npos);
    else
        targetmag = NULL;
    for(int i=0;i<npos;i++){
        b0[i] = a.b0[i];
        dxpos[i] = a.dxpos[i];
        dypos[i] = a.dypos[i];
        dzpos[i] = a.dzpos[i];
        target[i] = a.target[i];
        if (a.targetmag != NULL)
            targetmag[i] = a.targetmag[i];
    }
    free(sr);sr = (double*)malloc(sizeof(double) * npos*nTx);
    free(si);si = (double*)malloc(sizeof(double) * npos*nTx);
    for(int i=0;i<npos*nTx;i++){
        sr[i] = a.sr[i];
        si[i] = a.si[i];
    }
    free(dfreq);dfreq = (double*)malloc(sizeof(double) * nfreq);
    for(int i=0;i<nfreq;i++){
        dfreq[i] = a.dfreq[i];
    }

    unsigned int ntout;
    if (mode & 2)
        ntout = ntime; /* Include time points.	*/
    else
        ntout = 1;
    unsigned int ntnfnpos = ntout * nfreq * npos;
    free(mx);mx = (double*)malloc(sizeof(double) * ntnfnpos);
    free(my);my = (double*)malloc(sizeof(double) * ntnfnpos);
    free(mz);mz = (double*)malloc(sizeof(double) * ntnfnpos);
    for(int i=0;i<nfreq;i++){
        mx[i] = a.mx[i];
        my[i] = a.my[i];
        mz[i] = a.mz[i];
    }

   return *this;
}
