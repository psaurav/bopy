#include <cstdio>
#include <iostream>
#include <vector>
#include <Python.h>
#include <numpy/arrayobject.h>
#include "cgen3fit.h"
#include "clinfit.h"

using namespace std;

static PyObject* fitAPhiDC(PyObject* self, PyObject* args) {
    //printf ("Inside getAPhiDC\n"); 
    //int n;
    PyObject* dict;
    PyArrayObject* data;
    //if (!PyArg_ParseTuple(args, "O!ddIddI", &PyArray_Type, &data, &del_t, &freq, &nOfData, &absErr, &relErr, &maxIter)) {
    if (!PyArg_ParseTuple(args, "O!O", &PyArray_Type, &data, &dict)) {
        return NULL;
    }
    double del_t = PyFloat_AsDouble(PyDict_GetItemString(dict, "delt"));
    double freq = PyFloat_AsDouble(PyDict_GetItemString(dict, "freq"));
    unsigned int nOfData = PyInt_AsLong(PyDict_GetItemString(dict, "nOfData"));
    double absErr = PyFloat_AsDouble(PyDict_GetItemString(dict, "absErr"));
    double relErr = PyFloat_AsDouble(PyDict_GetItemString(dict, "relErr"));
    unsigned int maxIter = PyInt_AsLong(PyDict_GetItemString(dict, "maxIter"));

    // 
    // Perform sanity checks 
    // 
    if (PyArray_NDIM(data) != 3) {
        printf("ERROR: data array not 3D");
        return NULL;
    }
    if (PyArray_DIM(data,0) < nOfData) {
        printf("ERROR: z-size of data is less than nOfData");
        return NULL;
    }
    unsigned int nn = PyArray_DIM(data,0);  // number of frames (outer most loop)
    unsigned int nr = PyArray_DIM(data,1);  // number of rows
    unsigned int nc = PyArray_DIM(data,2);  // number of columns (inner most loop)

    cerr << "Number of frames: " << nn << endl;
    cerr << "Number of rows: " << nr << endl;
    cerr << "Number of cols: " << nc << endl;

    //unsigned int *dim = (unsigned int*) PyArray_DIMS(data);
    std::vector<unsigned int> dim(3);
    dim[0] = nc;
    dim[1] = nr;
    dim[2] = nn;
    //npy_intp npy_dim2d[2] = {nr, nc};
    
    //PyArrayObject* A = (PyArrayObject *) PyArray_SimpleNew(2, npy_dim2d, PyArray_DOUBLE);
    //PyArrayObject* phi = (PyArrayObject *) PyArray_SimpleNew(2, npy_dim2d, PyArray_DOUBLE);
    //PyArrayObject* DC = (PyArrayObject *) PyArray_SimpleNew(2, npy_dim2d, PyArray_DOUBLE);
    //PyArrayObject* iters = (PyArrayObject *) PyArray_SimpleNew(2, npy_dim2d, PyArray_UINT);
    //PyArrayObject* status = (PyArrayObject *) PyArray_SimpleNew(2, npy_dim2d, PyArray_INT);
    //PyArrayObject* delA = (PyArrayObject *) PyArray_SimpleNew(2, npy_dim2d, PyArray_DOUBLE);
    //PyArrayObject* delPhi = (PyArrayObject *) PyArray_SimpleNew(2, npy_dim2d, PyArray_DOUBLE);
    //PyArrayObject* delDC = (PyArrayObject *) PyArray_SimpleNew(2, npy_dim2d, PyArray_DOUBLE);
    //PyArrayObject* initA = (PyArrayObject *) PyArray_SimpleNew(2, npy_dim2d, PyArray_DOUBLE);
    //PyArrayObject* initPhi = (PyArrayObject *) PyArray_SimpleNew(2, npy_dim2d, PyArray_DOUBLE);
    //PyArrayObject* initDC = (PyArrayObject *) PyArray_SimpleNew(2, npy_dim2d, PyArray_DOUBLE);
    //PyArrayObject* chi2ByDOF = (PyArrayObject *) PyArray_SimpleNew(2, npy_dim2d, PyArray_DOUBLE);

    double *dataPtr = (double*) PyArray_DATA(data); 
    double *APtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "A")); 
    double *phiPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "phi")); 
    double *DCPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "DC")); 
    unsigned long int* itersPtr = (unsigned long int*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "iters"));
    long int* statusPtr = (long int*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "status"));
    double* delDCPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "delDC"));
    double* delPhiPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "delPhi"));
    double* delAPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "delA"));
    //double* initAPtr = (double*) PyArray_DATA(initA);
    //double* initPhiPtr = (double*) PyArray_DATA(initPhi);
    //double* initDCPtr = (double*) PyArray_DATA(initDC);
    double* chi2ByDOFPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "chi2ByDOF"));


    //
    // Debug
    //
    //cerr << "First element: " << dataPtr[0] << endl;
    //long int index;
    //index = 400+nc*(0+nr*0);
    //cerr << "First [0,0,400]: " << dataPtr[index] << endl;
    //index = 400+nc*(2+nr*0);
    //cerr << "First [0,2,400]: " << dataPtr[index] << endl;
    
    int err = getAPhiDC(dataPtr, 
                &dim[0], 
                del_t, 
                freq, 
                nOfData, 
                absErr, 
                relErr, 
                maxIter, 
                APtr, 
                phiPtr, 
                DCPtr, 
                itersPtr, 
                statusPtr,
                delDCPtr,
                delAPtr,
                delPhiPtr,
                chi2ByDOFPtr);
                //initAPtr,
                //initPhiPtr,
                //initDCPtr);

    return Py_BuildValue("i", 0);
}

static PyObject* fitAPhiDC2(PyObject* self, PyObject* args) {
    PyObject* dict;
    PyArrayObject* data;
    if (!PyArg_ParseTuple(args, "O!O", &PyArray_Type, &data, &dict)) {
        return NULL;
    }
    double del_t = PyFloat_AsDouble(PyDict_GetItemString(dict, "delt"));
    double freq = PyFloat_AsDouble(PyDict_GetItemString(dict, "freq"));
    double phi_init = PyFloat_AsDouble(PyDict_GetItemString(dict, "phi_init"));
    unsigned int nOfData = PyInt_AsLong(PyDict_GetItemString(dict, "nOfData"));
    double absErr = PyFloat_AsDouble(PyDict_GetItemString(dict, "absErr"));
    double relErr = PyFloat_AsDouble(PyDict_GetItemString(dict, "relErr"));
    unsigned int maxIter = PyInt_AsLong(PyDict_GetItemString(dict, "maxIter"));

    // 
    // Perform sanity checks 
    // 
    if (PyArray_NDIM(data) != 3) {
        printf("ERROR: data array not 3D");
        return NULL;
    }
    if (PyArray_DIM(data,0) < nOfData) {
        printf("ERROR: z-size of data is less than nOfData");
        return NULL;
    }
    unsigned int nn = PyArray_DIM(data,0);  // number of frames (outer most loop)
    unsigned int nr = PyArray_DIM(data,1);  // number of rows
    unsigned int nc = PyArray_DIM(data,2);  // number of columns (inner most loop)

    //cerr << "Number of frames: " << nn << endl;
    //cerr << "Number of rows: " << nr << endl;
    //cerr << "Number of cols: " << nc << endl;

    //unsigned int *dim = (unsigned int*) PyArray_DIMS(data);
    std::vector<unsigned int> dim(3);
    dim[0] = nc;
    dim[1] = nr;
    dim[2] = nn;

    double *dataPtr = (double*) PyArray_DATA(data); 
    double *APtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "A")); 
    double *phiPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "phi")); 
    double *DCPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "DC")); 
    unsigned long int* itersPtr = (unsigned long int*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "iters"));
    long int* statusPtr = (long int*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "status"));
    double* delDCPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "delDC"));
    double* delPhiPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "delPhi"));
    double* delAPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "delA"));
    double* chi2ByDOFPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "chi2ByDOF"));
    double* pcts = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "pcts"));
    double* rse = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "rse"));
    double* nrse = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "nrse"));

    int err = getAPhiDC2(dataPtr, 
                &dim[0], 
                del_t, 
                freq, 
                phi_init,
                nOfData, 
                absErr, 
                relErr, 
                maxIter, 
                APtr, 
                phiPtr, 
                DCPtr, 
                itersPtr, 
                statusPtr,
                delDCPtr,
                delAPtr,
                delPhiPtr,
                chi2ByDOFPtr,
                pcts,
                rse,
                nrse);

    return Py_BuildValue("i", 0);
}

static PyObject* fitLinAPhiDC(PyObject* self, PyObject* args) {
    PyObject* dict;
    PyArrayObject* data;
    if (!PyArg_ParseTuple(args, "O!O", &PyArray_Type, &data, &dict)) {
        return NULL;
    }
    unsigned int nOfData = PyInt_AsLong(PyDict_GetItemString(dict, "nOfData"));

    // 
    // Perform sanity checks 
    // 
    if (PyArray_NDIM(data) != 3) {
        printf("ERROR: data array not 3D");
        return NULL;
    }
    if (PyArray_DIM(data,0) < nOfData) {
        printf("ERROR: z-size of data is less than nOfData");
        return NULL;
    }
    unsigned int nn = PyArray_DIM(data,0);  // number of frames (outer most loop)
    unsigned int nr = PyArray_DIM(data,1);  // number of rows
    unsigned int nc = PyArray_DIM(data,2);  // number of columns (inner most loop)

    //unsigned int *dim = (unsigned int*) PyArray_DIMS(data);
    std::vector<unsigned int> dim(3);
    dim[0] = nc;
    dim[1] = nr;
    dim[2] = nn;

    double *dataPtr = (double*) PyArray_DATA(data); 
    double *pinvPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "pinv")); 
    double *workPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "work")); 
    double *APtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "A")); 
    double *phiPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "phi")); 
    double *DCPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "DC")); 
    double* delDCPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "delDC"));
    double* delPhiPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "delPhi"));
    double* delAPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "delA"));
    double* chi2ByDOFPtr = (double*) PyArray_DATA((PyArrayObject*) PyDict_GetItemString(dict, "chi2ByDOF"));

    int err = getLinAPhiDC(pinvPtr,
            dataPtr,
            &dim[0],
            nOfData,
            workPtr,
            APtr,
            phiPtr,
            DCPtr,
            delDCPtr,
            delAPtr,
            delPhiPtr,
            chi2ByDOFPtr);

    return Py_BuildValue("i", 0);
}


static PyObject* fitPixelAPhiDC2(PyObject* self, PyObject* args) {
    PyObject* dict;
    PyArrayObject* data;
    if (!PyArg_ParseTuple(args, "O!O", &PyArray_Type, &data, &dict)) {
        return NULL;
    }
    
    std::cerr << "just before row: " << std::endl;
    unsigned long int row = (unsigned long int) PyInt_AsLong(PyDict_GetItemString(dict, "row"));
    std::cerr << "row: " << row << std::endl;
    unsigned long int col = (unsigned long int) PyInt_AsLong(PyDict_GetItemString(dict, "col"));
    double del_t = PyFloat_AsDouble(PyDict_GetItemString(dict, "delt"));
    double freq = PyFloat_AsDouble(PyDict_GetItemString(dict, "freq"));
    double phi_init = PyFloat_AsDouble(PyDict_GetItemString(dict, "phi_init"));
    unsigned int nOfData = PyInt_AsLong(PyDict_GetItemString(dict, "nOfData"));
    double absErr = PyFloat_AsDouble(PyDict_GetItemString(dict, "absErr"));
    double relErr = PyFloat_AsDouble(PyDict_GetItemString(dict, "relErr"));
    unsigned int maxIter = PyInt_AsLong(PyDict_GetItemString(dict, "maxIter"));


    // 
    // Perform sanity checks 
    // 
    if (PyArray_NDIM(data) != 3) {
        printf("ERROR: data array not 3D");
        return NULL;
    }
    if (PyArray_DIM(data,0) < nOfData) {
        printf("ERROR: z-size of data is less than nOfData");
        return NULL;
    }
    unsigned int nn = PyArray_DIM(data,0);  // number of frames (outer most loop)
    unsigned int nr = PyArray_DIM(data,1);  // number of rows
    unsigned int nc = PyArray_DIM(data,2);  // number of columns (inner most loop)

    cerr << "Number of frames: " << nn << endl;
    cerr << "Number of rows: " << nr << endl;
    cerr << "Number of cols: " << nc << endl;

    //unsigned int *dim = (unsigned int*) PyArray_DIMS(data);
    std::vector<unsigned int> dim(3);
    dim[0] = nc;
    dim[1] = nr;
    dim[2] = nn;

    double *dataPtr = (double*) PyArray_DATA(data); 

    double A, phi, DC, delA, delPhi, delDC, chi2ByDOF, pcts, rse, nrse, ssr;
    unsigned long int iters;
    long int status;

    int err = getPixelAPhiDC2(dataPtr, 
                &dim[0], 
                row,
                col,
                del_t, 
                freq, 
                phi_init,
                nOfData, 
                absErr, 
                relErr, 
                maxIter, 
                A, 
                phi, 
                DC, 
                iters, 
                status,
                delDC,
                delA,
                delPhi,
                chi2ByDOF,
                pcts,
                rse,
                nrse,
                ssr);

    PyDict_SetItemString(dict, "A", PyFloat_FromDouble(A)); 
    PyDict_SetItemString(dict, "phi", PyFloat_FromDouble(phi)); 
    PyDict_SetItemString(dict, "DC", PyFloat_FromDouble(DC)); 
    PyDict_SetItemString(dict, "iters", PyInt_FromLong(iters));
    PyDict_SetItemString(dict, "status", PyInt_FromLong(status));
    PyDict_SetItemString(dict, "delDC", PyFloat_FromDouble(delDC));
    PyDict_SetItemString(dict, "delPhi", PyFloat_FromDouble(delPhi));
    PyDict_SetItemString(dict, "delA", PyFloat_FromDouble(delA));
    PyDict_SetItemString(dict, "chi2ByDOF", PyFloat_FromDouble(chi2ByDOF));
    PyDict_SetItemString(dict, "pcts", PyFloat_FromDouble(pcts));
    PyDict_SetItemString(dict, "rse", PyFloat_FromDouble(rse));
    PyDict_SetItemString(dict, "nrse", PyFloat_FromDouble(nrse));
    PyDict_SetItemString(dict, "ssr", PyFloat_FromDouble(ssr));

    return Py_BuildValue("i", 0);
}

static PyObject* fitPixelAPhiDCTr(PyObject* self, PyObject* args) {
    PyObject* dict;
    PyArrayObject* data;
    if (!PyArg_ParseTuple(args, "O!O", &PyArray_Type, &data, &dict)) {
        return NULL;
    }
    
    std::cerr << "just before row: " << std::endl;
    unsigned long int row = (unsigned long int) PyInt_AsLong(PyDict_GetItemString(dict, "row"));
    std::cerr << "row: " << row << std::endl;
    unsigned long int col = (unsigned long int) PyInt_AsLong(PyDict_GetItemString(dict, "col"));
    double del_t = PyFloat_AsDouble(PyDict_GetItemString(dict, "delt"));
    double freq = PyFloat_AsDouble(PyDict_GetItemString(dict, "freq"));
    double phi_init = PyFloat_AsDouble(PyDict_GetItemString(dict, "phi_init"));
    unsigned int nOfData = PyInt_AsLong(PyDict_GetItemString(dict, "nOfData"));
    double absErr = PyFloat_AsDouble(PyDict_GetItemString(dict, "absErr"));
    double relErr = PyFloat_AsDouble(PyDict_GetItemString(dict, "relErr"));
    unsigned int maxIter = PyInt_AsLong(PyDict_GetItemString(dict, "maxIter"));


    // 
    // Perform sanity checks 
    // 
    if (PyArray_NDIM(data) != 3) {
        printf("ERROR: data array not 3D");
        return NULL;
    }
    if (PyArray_DIM(data,0) < nOfData) {
        printf("ERROR: z-size of data is less than nOfData");
        return NULL;
    }
    unsigned int nn = PyArray_DIM(data,0);  // number of frames (outer most loop)
    unsigned int nr = PyArray_DIM(data,1);  // number of rows
    unsigned int nc = PyArray_DIM(data,2);  // number of columns (inner most loop)

    cerr << "Number of frames: " << nn << endl;
    cerr << "Number of rows: " << nr << endl;
    cerr << "Number of cols: " << nc << endl;

    //unsigned int *dim = (unsigned int*) PyArray_DIMS(data);
    std::vector<unsigned int> dim(3);
    dim[0] = nc;
    dim[1] = nr;
    dim[2] = nn;

    double *dataPtr = (double*) PyArray_DATA(data); 

    double A, phi, DC, alpha, delA, delPhi, delDC, delAlpha, chi2ByDOF, pcts, rse, nrse, ssr;
    unsigned long int iters;
    long int status;

    int err = getPixelAPhiDCTr(dataPtr, 
                &dim[0], 
                row,
                col,
                del_t, 
                freq, 
                phi_init,
                nOfData, 
                absErr, 
                relErr, 
                maxIter, 
                A, 
                phi, 
                DC, 
                alpha,
                iters, 
                status,
                delDC,
                delA,
                delPhi,
                delAlpha,
                chi2ByDOF,
                pcts,
                rse,
                nrse,
                ssr);

    PyDict_SetItemString(dict, "A", PyFloat_FromDouble(A)); 
    PyDict_SetItemString(dict, "phi", PyFloat_FromDouble(phi)); 
    PyDict_SetItemString(dict, "DC", PyFloat_FromDouble(DC)); 
    PyDict_SetItemString(dict, "alpha", PyFloat_FromDouble(alpha)); 
    PyDict_SetItemString(dict, "iters", PyInt_FromLong(iters));
    PyDict_SetItemString(dict, "status", PyInt_FromLong(status));
    PyDict_SetItemString(dict, "delDC", PyFloat_FromDouble(delDC));
    PyDict_SetItemString(dict, "delPhi", PyFloat_FromDouble(delPhi));
    PyDict_SetItemString(dict, "delA", PyFloat_FromDouble(delA));
    PyDict_SetItemString(dict, "delAlpha", PyFloat_FromDouble(delAlpha));
    PyDict_SetItemString(dict, "chi2ByDOF", PyFloat_FromDouble(chi2ByDOF));
    PyDict_SetItemString(dict, "pcts", PyFloat_FromDouble(pcts));
    PyDict_SetItemString(dict, "rse", PyFloat_FromDouble(rse));
    PyDict_SetItemString(dict, "nrse", PyFloat_FromDouble(nrse));
    PyDict_SetItemString(dict, "ssr", PyFloat_FromDouble(ssr));

    return Py_BuildValue("i", 0);
}

static PyMethodDef gen3fitMethods[] =
{
    {"fitAPhiDC", fitAPhiDC, METH_VARARGS, "Get Amplitude, Phase and DC for Gen3 Data."},
    {"fitAPhiDC2", fitAPhiDC2, METH_VARARGS, "Get Amplitude, Phase and DC for Gen3 Data (v2)."},
    {"fitLinAPhiDC", fitLinAPhiDC, METH_VARARGS, "Get Amplitude, Phase and DC for Gen3 Data (linear fit)."},
    {"fitPixelAPhiDC2", fitPixelAPhiDC2, METH_VARARGS, "Get Amplitude, Phase and DC for single pixel Gen3 Data (v2)."},
    {"fitPixelAPhiDCTr", fitPixelAPhiDCTr, METH_VARARGS, "Get Amplitude, Phase and DC for single pixel Gen3 Data with trend (v2)."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initgen3fit(void)
{
    (void) Py_InitModule("gen3fit", gen3fitMethods);
    import_array();
}

