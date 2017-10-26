#include "Python.h"
#include "numpy/arrayobject.h"
#include <fcntl.h>
#include <math.h>
#include <omp.h>
#include <sys/param.h>

#define VERSION "0.1"
#define FREE_ARG char*
#define NR_END 1

#define fillincrement 0.01
#define oneoversqrt2 0.707106781187

double pexp;
double **topo,**flow,**flow1,**flow2,**flow3,**flow4,**flow5,**flow6,**flow7,**flow8;
int lattice_size_x,lattice_size_y,*iup,*idown,*jup,*jdown;

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

double *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        double *v;

        v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
        return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh) {
        free((FREE_ARG) (v+nl-NR_END));
}

void free_vector(double *v, long nl, long nh) {
        free((FREE_ARG) (v+nl-NR_END));
}

double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	m += NR_END;
	m -= nrl;

	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	return m;
}

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 100000

void indexx(n,arr,indx)
double arr[];
int indx[],n;
{
        int i,indxt,ir=n,itemp,j,k,l=1;
        int *istack;
        int jstack=0;
        double a;

        istack=ivector(1,NSTACK);
        for (j=1;j<=n;j++) indx[j]=j;
        for (;;) {
                if (ir-l < M) {
                        for (j=l+1;j<=ir;j++) {
                                indxt=indx[j];
                                a=arr[indxt];
                                for (i=j-1;i>=1;i--) {
                                        if (arr[indx[i]] <= a) break;
                                        indx[i+1]=indx[i];
                                }
                                indx[i+1]=indxt;
                        }
                        if (jstack == 0) break;
                        ir=istack[jstack--];
                        l=istack[jstack--];
                } else {
                        k=(l+ir) >> 1;
                        SWAP(indx[k],indx[l+1]);
                        if (arr[indx[l+1]] > arr[indx[ir]]) {
                                SWAP(indx[l+1],indx[ir])
                        }
                        if (arr[indx[l]] > arr[indx[ir]]) {
                                SWAP(indx[l],indx[ir])
                        }
                        if (arr[indx[l+1]] > arr[indx[l]]) {
                                SWAP(indx[l+1],indx[l])
                        }
                        i=l+1;
                        j=ir;
                        indxt=indx[l];
                        a=arr[indxt];
                        for (;;) {
                                do i++; while (arr[indx[i]] < a);
                                do j--; while (arr[indx[j]] > a);
                                if (j < i) break;
                                SWAP(indx[i],indx[j])
                        }
                        indx[l]=indx[j];
                        indx[j]=indxt;
                        jstack += 2;
                        if (ir-i+1 >= j-l) {
                                istack[jstack]=ir;
                                istack[jstack-1]=i;
                                ir=j-1;
                        } else {
                                istack[jstack]=j-1;
                                istack[jstack-1]=l;
                                l=i;
                        }
                }
        }
        free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP

void setupgridneighbors()
{    int i,j;

     idown=ivector(1,lattice_size_x);
     iup=ivector(1,lattice_size_x);
     jup=ivector(1,lattice_size_y);
     jdown=ivector(1,lattice_size_y);
     for (i=1;i<=lattice_size_x;i++)
      {idown[i]=i-1;
       iup[i]=i+1;}
     idown[1]=1;
     iup[lattice_size_x]=lattice_size_x;
     for (j=1;j<=lattice_size_y;j++)
      {jdown[j]=j-1;
       jup[j]=j+1;}
     jdown[1]=1;
     jup[lattice_size_y]=lattice_size_y;
}

void fillinpitsandflats(i,j)
int i,j;
{    double min;

     min=topo[i][j];
     if (topo[iup[i]][j]<min) min=topo[iup[i]][j];
     if (topo[idown[i]][j]<min) min=topo[idown[i]][j];
     if (topo[i][jup[j]]<min) min=topo[i][jup[j]];
     if (topo[i][jdown[j]]<min) min=topo[i][jdown[j]];
     if (topo[iup[i]][jup[j]]<min) min=topo[iup[i]][jup[j]];
     if (topo[idown[i]][jup[j]]<min) min=topo[idown[i]][jup[j]];
     if (topo[idown[i]][jdown[j]]<min) min=topo[idown[i]][jdown[j]];
     if (topo[iup[i]][jdown[j]]<min) min=topo[iup[i]][jdown[j]];
     if ((topo[i][j]<=min)&&(i>1)&&(j>1)&&(i<lattice_size_x)&&(j<lattice_size_y))
      {topo[i][j]=min+fillincrement;
       fillinpitsandflats(i,j);
       fillinpitsandflats(iup[i],j);
       fillinpitsandflats(idown[i],j);
       fillinpitsandflats(i,jup[j]);
       fillinpitsandflats(i,jdown[j]);
       fillinpitsandflats(iup[i],jup[j]);
       fillinpitsandflats(idown[i],jup[j]);
       fillinpitsandflats(idown[i],jdown[j]);
       fillinpitsandflats(iup[i],jdown[j]);}
}

void mfdflowroute(i,j)
int i,j;
{    double tot;
 
     tot=0;
     if (topo[i][j]>topo[iup[i]][j]) 
      tot+=pow(topo[i][j]-topo[iup[i]][j],pexp);
     if (topo[i][j]>topo[idown[i]][j]) 
      tot+=pow(topo[i][j]-topo[idown[i]][j],pexp);
     if (topo[i][j]>topo[i][jup[j]]) 
      tot+=pow(topo[i][j]-topo[i][jup[j]],pexp);
     if (topo[i][j]>topo[i][jdown[j]]) 
      tot+=pow(topo[i][j]-topo[i][jdown[j]],pexp);
     if (topo[i][j]>topo[iup[i]][jup[j]]) 
      tot+=pow((topo[i][j]-topo[iup[i]][jup[j]])*oneoversqrt2,pexp);
     if (topo[i][j]>topo[iup[i]][jdown[j]]) 
      tot+=pow((topo[i][j]-topo[iup[i]][jdown[j]])*oneoversqrt2,pexp);
     if (topo[i][j]>topo[idown[i]][jup[j]]) 
      tot+=pow((topo[i][j]-topo[idown[i]][jup[j]])*oneoversqrt2,pexp);
     if (topo[i][j]>topo[idown[i]][jdown[j]]) 
      tot+=pow((topo[i][j]-topo[idown[i]][jdown[j]])*oneoversqrt2,pexp);
     if (topo[i][j]>topo[iup[i]][j]) 
      flow1[i][j]=pow(topo[i][j]-topo[iup[i]][j],pexp)/tot; 
       else flow1[i][j]=0;
     if (topo[i][j]>topo[idown[i]][j]) 
      flow2[i][j]=pow(topo[i][j]-topo[idown[i]][j],pexp)/tot; 
       else flow2[i][j]=0;
     if (topo[i][j]>topo[i][jup[j]]) 
      flow3[i][j]=pow(topo[i][j]-topo[i][jup[j]],pexp)/tot; 
       else flow3[i][j]=0;
     if (topo[i][j]>topo[i][jdown[j]]) 
      flow4[i][j]=pow(topo[i][j]-topo[i][jdown[j]],pexp)/tot; 
       else flow4[i][j]=0;
     if (topo[i][j]>topo[iup[i]][jup[j]]) 
      flow5[i][j]=pow((topo[i][j]-topo[iup[i]][jup[j]])*oneoversqrt2,pexp)/tot;
       else flow5[i][j]=0;
     if (topo[i][j]>topo[iup[i]][jdown[j]]) 
      flow6[i][j]=pow((topo[i][j]-topo[iup[i]][jdown[j]])*oneoversqrt2,pexp)/tot;
       else flow6[i][j]=0;
     if (topo[i][j]>topo[idown[i]][jup[j]]) 
      flow7[i][j]=pow((topo[i][j]-topo[idown[i]][jup[j]])*oneoversqrt2,pexp)/tot;
       else flow7[i][j]=0;
     if (topo[i][j]>topo[idown[i]][jdown[j]]) 
      flow8[i][j]=pow((topo[i][j]-topo[idown[i]][jdown[j]])*oneoversqrt2,pexp)/tot;
       else flow8[i][j]=0;
     flow[iup[i]][j]+=flow[i][j]*flow1[i][j];
     flow[idown[i]][j]+=flow[i][j]*flow2[i][j];
     flow[i][jup[j]]+=flow[i][j]*flow3[i][j];
     flow[i][jdown[j]]+=flow[i][j]*flow4[i][j];
     flow[iup[i]][jup[j]]+=flow[i][j]*flow5[i][j];
     flow[iup[i]][jdown[j]]+=flow[i][j]*flow6[i][j];
     flow[idown[i]][jup[j]]+=flow[i][j]*flow7[i][j];
     flow[idown[i]][jdown[j]]+=flow[i][j]*flow8[i][j];
}

void
SCA(double *sca, const double *z,
    const double w) {
    int i, j, k, t;
    int *topovecind;
    double *topovec, delta;
    
    delta = w * w;
    setupgridneighbors();
    topo=matrix(1,lattice_size_x,1,lattice_size_y);
    topovec=vector(1,lattice_size_x*lattice_size_y);
    topovecind=ivector(1,lattice_size_x*lattice_size_y);
    flow=matrix(1,lattice_size_x,1,lattice_size_y);
    flow1=matrix(1,lattice_size_x,1,lattice_size_y);
    flow2=matrix(1,lattice_size_x,1,lattice_size_y);
    flow3=matrix(1,lattice_size_x,1,lattice_size_y);
    flow4=matrix(1,lattice_size_x,1,lattice_size_y);
    flow5=matrix(1,lattice_size_x,1,lattice_size_y);
    flow6=matrix(1,lattice_size_x,1,lattice_size_y);
    flow7=matrix(1,lattice_size_x,1,lattice_size_y);
    flow8=matrix(1,lattice_size_x,1,lattice_size_y);
    
    k = 0;
    for(j = 1; j <= lattice_size_y; j++) {
        for(i = 1; i <= lattice_size_x; i++) {
            topo[i][j] = z[k++];
            flow[i][j] = delta;
        }
    }
    for(j = 1; j <= lattice_size_y; j++)
        for(i = 1; i <= lattice_size_x; i++)
            fillinpitsandflats(i,j);
    for(j = 1; j <= lattice_size_y; j++)
        for(i = 1; i <= lattice_size_x; i++)
            topovec[(j-1)*lattice_size_x+i] = topo[i][j];
    
    indexx(lattice_size_x * lattice_size_y, topovec, topovecind);
    t = lattice_size_x * lattice_size_y + 1;
    while(t > 1 ) {
        t--;
        i = (topovecind[t]) % lattice_size_x;
        if(i == 0)
            i = lattice_size_x;
        j = (topovecind[t]) / lattice_size_x+1;
        if(i == lattice_size_x)
            j--;
        mfdflowroute(i,j);
    }
    k = 0;
    for(j = 1; j <= lattice_size_y; j++)
        for(i = 1;i <= lattice_size_x; i++)
            sca[k++] = flow[i][j] / w;
}

static PyObject *
mfdrouting_SCA(PyObject *self, PyObject* args) {
    PyObject *zarg;
    PyArrayObject *z, *sca;
    double w;

    pexp = 1.1;

    // parse input
    if(!PyArg_ParseTuple(args, "Od|d", &zarg, &w, &pexp))
        return NULL;
    z = (PyArrayObject *) PyArray_ContiguousFromObject(zarg, PyArray_DOUBLE, 2, 2);
    sca = (PyArrayObject *) PyArray_ZEROS(2, z->dimensions, PyArray_DOUBLE, 0);
    if(!z || !sca)
        return NULL;
    
    lattice_size_y = z->dimensions[0];
    lattice_size_x = z->dimensions[1];

    // get specific catchment area estimates
    SCA((double *)sca->data,
        (double *)z->data, w);
    Py_DECREF(z);
    return PyArray_Return(sca);
}

static PyMethodDef mfdrouting_Methods[] = {
    {"SCA", mfdrouting_SCA, METH_VARARGS, "SCA(z, cellwidth, exponent = 1.1)\nReturns the specific catchment area of the topography (DEM) in z.\n\nParameters\n----------\nz : array_like\n\tA 2-D array containing the elevations of the topography (DEM).\ncellwidth : float\n\tThe cell width of the grid in meters.\nexponent : float, optional\n\tThe exponent determining the flow divergence in relation to the slope.\n\tSee http://dx.doi.org/10.1016/0098-3004(91)90048-i.\n\nReturns\n-------\na : ndarray\n\tThe specific catchment area matrix on the same grid as z (DEM).\n"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef ModDef = {
    PyModuleDef_HEAD_INIT,
    "mfdrouting",
    NULL,
    -1,
    mfdrouting_Methods
};

PyMODINIT_FUNC
PyInit_mfdrouting(void) {
    PyObject *mod;

    mod = PyModule_Create(&ModDef);
    PyModule_AddStringConstant(mod, "__author__", "Aljoscha Rheinwalt <aljoscha.rheinwalt@uni-potsdam.de>");
    PyModule_AddStringConstant(mod, "__version__", VERSION);
    import_array();

    return mod;
}

int
main(int argc, char **argv) {
    wchar_t pname[255];

    PyImport_AppendInittab("mfdrouting", PyInit_mfdrouting);
    mbstowcs(pname, argv[0], strlen(argv[0])+1);
    Py_SetProgramName(pname);
    Py_Initialize();
    PyImport_ImportModule("mfdrouting");
    PyMem_RawFree(argv[0]);
    return 0;
}
