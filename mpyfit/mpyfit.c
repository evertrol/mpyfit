/*
 * Copyright (c) SRON - Netherlands Institute for Space Research (2013).
 * All Rights Reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above
 *    copyright notice, this list of conditions and the following
 *    disclaimer in the documentation and/or other materials provided
 *    with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <Python.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <numpy/arrayobject.h>
#include "mpfit.h"

#if PY_MAJOR_VERSION >= 3
#define PY3K
#endif

#ifndef NPY_ARRAY_IN_ARRAY
#define NPY_ARRAY_IN_ARRAY (NPY_CONTIGUOUS|NPY_ALIGNED)
#endif

struct user_data {
	PyObject *function;  // User supplied Python function
	PyArrayObject *parameters;  // Function parameters to fit for
	PyTupleObject *arguments;  // Extra function arguments
};
typedef struct user_data user_data;


/**
 *  Call the Python function 'function' with 'parameters' as the first
 *  argument and arguments as further arguments (unless arguments is
 *  NULL).
 *
 *  Return the output values as a double array of y and length ny.
 *
 */
int
call_function(PyObject *function, PyArrayObject *parameters,
              PyTupleObject *arguments, size_t *ny, double *y)
{
	PyObject *result = NULL;
	if (PyTuple_Size((PyObject *)arguments) > 0) {
		result = PyObject_CallFunctionObjArgs(
			function, parameters, arguments, NULL);
	} else {
		result = PyObject_CallFunctionObjArgs(
			function, parameters, NULL);
	}
	if (!result) {
		return -1;
	}
	PyArrayObject *output = (PyArrayObject *)PyArray_FROM_OTF(
		result, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	if (!output) {
		PyErr_Format(PyExc_TypeError,
		             "return value was not an array");
		return 2;
	}
	if (PyArray_NDIM(output) == 0) {
		*ny = 1;
	} else if (PyArray_NDIM(output) > 1) {
		PyErr_Format(PyExc_TypeError,
		             "return value array is one dimensional");
		Py_XDECREF(result);
		Py_XDECREF(output);
		return 3;
	} else {
		*ny = (size_t)PyArray_DIM(output, 0);
	}
	// Copy data into y if y is not NULL
	if (y) {
		double *tmpy = PyArray_DATA(output);
		memcpy(y, tmpy, *ny * sizeof(*tmpy));
	}
	Py_XDECREF(result);
	Py_XDECREF(output);

	return 0;
}


int
wrapper_user_func(int nx, int np, double *x, double *fvec, double **dvec,
                  void *private_data)
/**
 *  This is the function that is supplied to mpfit to minimize (i.e.,
 *  an mp_func type of function).
 *
 *  This function will call the Python function with addition data
 *  that the user supplied. Thus, mpfit is not directly calling the
 *  Python function (since that would require mpfit to create Python
 *  objects. Instead, the Python objects that serve as input arguments
 *  to the user supplied function.
 *
 *  The Python function returns an array of values that are passed
 *  into fvec.
 *
 */
{
	user_data *data = (user_data *)private_data;
	size_t ny;

	npy_intp dims[1];
	dims[0] = np;
	PyObject* array = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE,
	                                            (void *)x);
	if (!array) {
		return -11;
	}
	//Py_INCREF(array);
	int result = call_function(data->function, (PyArrayObject *)array,
	                           data->arguments, &ny, fvec);
	Py_XDECREF(array);
	if (result) {
		return result;
	}
	if (ny != nx) {
		PyErr_Format(PyExc_TypeError,
		             "return value array has incorrect size");
		return -12;
	}

	return 0;
}


static int
parse_limits(PyObject *infodict, mp_par *parinfo)
{
	PyObject *tuple = PyDict_GetItemString(infodict, "limits");
	PyObject *value = NULL;

	if (!tuple) {
		return 0;
	}

	if (!PyTuple_Check(tuple)) {
		PyErr_Format(PyExc_TypeError,
		             "'limits' is not a tuple");
		return 1;
	}
	if (PyTuple_Size(tuple) != 2) {
		PyErr_Format(PyExc_TypeError,
		             "'limits' does not have length 2");
		return 1;
	}
	value = PyTuple_GetItem(tuple, 0);
	if (!value) {
		PyErr_Format(PyExc_TypeError,
		             "can't parse 'limits'");
		return 1;
	}
	if (value != Py_None) {
		PyObject *floatval = PyNumber_Float(value);
		if (!floatval) {
			PyErr_Format(PyExc_TypeError,
			             "lower limit is not a number");
				}
		PyErr_Restore(NULL, NULL, NULL);
		parinfo->limits[0] = PyFloat_AsDouble(floatval);
		if (PyErr_Occurred()) {
			PyErr_Format(PyExc_TypeError,
			             "failed to parse lower limit");
		}
		parinfo->limited[0] = 1;
	}
	value = PyTuple_GetItem(tuple, 1);
	if (!value) {
		PyErr_Format(PyExc_TypeError,
		             "can't parse 'limits'");
		return 1;
	}
	if (value != Py_None) {
		PyObject *floatval = PyNumber_Float(value);
		if (!floatval) {
			PyErr_Format(PyExc_TypeError,
			             "upper limit is not a number");
				}
		PyErr_Restore(NULL, NULL, NULL);
		parinfo->limits[1] = PyFloat_AsDouble(floatval);
		if (PyErr_Occurred()) {
			PyErr_Format(PyExc_TypeError,
			             "failed to parse upper limit");
		}
		parinfo->limited[1] = 1;
	}

	return 0;
}


static PyObject*
mpfitwrap(PyObject *self, PyObject *args, PyObject *keywords)
{
	PyObject *arg1 = NULL;
	PyObject *parameters = NULL;
	PyTupleObject *user_args = NULL;
	PyListObject *parinfo = NULL;
	PyArrayObject *params = NULL;
	PyObject *fittedpars = NULL;
	PyArrayObject *parerrors = NULL;
	PyArrayObject *residuals = NULL;
	PyArrayObject *covariances = NULL;
	static char *keywordlist[] = {
		"function", "p", "args", "parinfo",
		"ftol", "xtol", "gtol", "epsfcn", "stepfactor",
		"covtol", "maxiter", "maxfev", NULL
	};
	mp_config mpconf = { 0 };
	mp_par *mpparinfo = NULL;


	if (!PyArg_ParseTupleAndKeywords(args, keywords,
#ifdef PY3K
	                                 "OO!|$O!O!ddddddii",
#else
	                                 "OO!|O!O!ddddddii",
#endif
	                                 keywordlist,
	                                 &arg1,
	                                 &PyArray_Type, &parameters,
	                                 &PyTuple_Type, &user_args,
	                                 &PyList_Type, &parinfo,
	                                 &mpconf.ftol,
	                                 &mpconf.xtol,
	                                 &mpconf.gtol,
	                                 &mpconf.epsfcn,
	                                 &mpconf.stepfactor,
	                                 &mpconf.covtol,
	                                 &mpconf.maxiter,
	                                 &mpconf.maxfev
		    )) {
		goto cleanup;
	}
	if (!PyCallable_Check(arg1)) {
		PyErr_Format(PyExc_TypeError, "not a callable function");
		goto cleanup;
	}
	PyObject *function = arg1;

	params = (PyArrayObject *)PyArray_FROM_OTF(
		parameters, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	if (!params) {
		PyErr_Format(PyExc_TypeError,
		             "failed to convert parameters array");
		goto cleanup;
	}
	if (PyArray_NDIM(params) != 1) {
		PyErr_Format(PyExc_TypeError,
		             "parameter array is not one dimensional");
		goto cleanup;
	}
	size_t np = (size_t)PyArray_DIM(params, 0);
	npy_intp dims[2];
	dims[0] = np;
	if (parinfo) {
		// Obtain the parameter information
		Py_ssize_t nparinfo = PyList_Size((PyObject *)parinfo);
		if (np != nparinfo) {
			PyErr_Format(PyExc_TypeError,
			             "the parameters and the parameter info "
			             "list do not have equal lengths");
			goto cleanup;
		}
		mpparinfo = malloc(np * sizeof(*mpparinfo));
		if (!mpparinfo) {
			PyErr_Format(PyExc_TypeError,
			             "memory allocation error");
			goto cleanup;
		}
		for (Py_ssize_t i = 0; i < nparinfo; i++) {
			PyObject *infodict;
			PyObject *value;
			// Set defaults
			mpparinfo[i].fixed = 0;
			mpparinfo[i].limited[0] = 0;
			mpparinfo[i].limited[1] = 0;
			mpparinfo[i].limits[0] = 0;
			mpparinfo[i].limits[1] = 0;
			mpparinfo[i].parname = NULL;
			mpparinfo[i].step = 0;
			mpparinfo[i].relstep = 0;
			mpparinfo[i].side = 0;
			mpparinfo[i].deriv_debug = 0;
			mpparinfo[i].deriv_reltol = 0;
			mpparinfo[i].deriv_abstol = 0;

			infodict = PyList_GetItem((PyObject *)parinfo, i);
			if (!PyDict_Check(infodict)) {
				PyErr_Format(PyExc_TypeError,
				             "parameter info items should be "
				             "dicts");
				goto cleanup;
			}

			value = PyDict_GetItemString(infodict, "fixed");
			if (value) {
				if (!PyBool_Check(value)) {
					PyErr_Format(PyExc_TypeError,
					             "'fixed' is not a bool");
					goto cleanup;
				}
				if (value == Py_True) {
					mpparinfo[i].fixed = 1;
				}
			}

			if (parse_limits(infodict, &mpparinfo[i])) {
				goto cleanup;
			}
			value = PyDict_GetItemString(infodict, "step");
			if (value) {
				PyObject *floatval = PyNumber_Float(value);
				if (!floatval) {
					PyErr_Format(PyExc_TypeError,
					             "'step' is not a float");
					goto cleanup;
				}
				PyErr_Restore(NULL, NULL, NULL);
				mpparinfo[i].step = PyFloat_AsDouble(floatval);
				if (PyErr_Occurred()) {
					PyErr_Format(PyExc_TypeError,
					             "failed to parse step");
					goto cleanup;
				}
			}
			value = PyDict_GetItemString(infodict, "relstep");
			if (value) {
				PyObject *floatval = PyNumber_Float(value);
				if (!floatval) {
					PyErr_Format(PyExc_TypeError,
					             "'relstep' is not a "
					             "float");
					goto cleanup;
				}
				PyErr_Restore(NULL, NULL, NULL);
				mpparinfo[i].relstep = PyFloat_AsDouble(
					floatval);
				if (PyErr_Occurred()) {
					PyErr_Format(PyExc_TypeError,
					             "failed to parse relstep");
					goto cleanup;
				}
			}
			value = PyDict_GetItemString(infodict, "deriv_debug");
			if (value) {
				if (!PyBool_Check(value)) {
					PyErr_Format(PyExc_TypeError,
					             "deriv_debug is not a "
					             "bool");
					goto cleanup;
				}
				if (value == Py_True) {
					mpparinfo[i].deriv_debug = 1;
				}
			}
			value = PyDict_GetItemString(infodict, "reltol");
			if (value) {
				if (!PyFloat_Check(value)) {
					PyErr_Format(PyExc_TypeError,
					             "'reltol' is not a "
					             "float");
					goto cleanup;
				}
				PyErr_Restore(NULL, NULL, NULL);
				mpparinfo[i].deriv_reltol = PyFloat_AsDouble(
					value);
				if (PyErr_Occurred()) {
					PyErr_Format(PyExc_TypeError,
					             "failed to parse reltol");
					goto cleanup;
				}
			}
			value = PyDict_GetItemString(infodict, "abstol");
			if (value) {
				if (!PyFloat_Check(value)) {
					PyErr_Format(PyExc_TypeError,
					             "'abstol' is not a "
					             "float");
					goto cleanup;
				}
				PyErr_Restore(NULL, NULL, NULL);
				mpparinfo[i].deriv_abstol = PyFloat_AsDouble(
					value);
				if (PyErr_Occurred()) {
					PyErr_Format(PyExc_TypeError,
					             "failed to parse abstol");
					goto cleanup;
				}
			}
		}
	}

	// We use these 3 steps, because PyArray_SimpleNewFromData
	// (which was originally used), doesn't own the data pointer
	// (pdata), and we can't free it once we returned from the
	// function.
	// See eg http://mail.python.org/pipermail/cplusplus-sig/2006-September/011051.html
	// Create a new array object (this allocates the memory)
	fittedpars = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
	// Get its data pointer
	double *pdata = PyArray_DATA((PyArrayObject *)fittedpars);
	// Copy the original data into the data pointer; pdata will be
	// altered by mpfit, and thus fittedpars will be set correctly.
	memcpy(pdata, PyArray_DATA(params), np * sizeof(*pdata));

	// Test that we can safely call the user function
	// This also obtains the size of the result array
	size_t ny;
	int status = call_function(function, params, user_args, &ny, NULL);
	if (status) {
		if (status < 0) {
			PyErr_Format(PyExc_TypeError,
			             "error calling function");
		}
		goto cleanup;
	}

	// Create an array to hold the parameter errors
	dims[0] = np;
	parerrors = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
	if (!parerrors) {
		goto cleanup;
	}
	// Create an array to hold the covariances
	dims[1] = np;
	covariances = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
	if (!covariances) {
		goto cleanup;
	}
	// Create an array to hold the residuals
	dims[0] = ny;
	residuals = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
	if (!residuals) {
		goto cleanup;
	}

	user_data private_data;
	private_data.function = function;
	private_data.parameters = params;
	private_data.arguments = user_args;

	mp_result result;
	memset(&result, 0, sizeof(result));
	result.xerror = PyArray_DATA((PyArrayObject *)parerrors);
	result.resid = PyArray_DATA((PyArrayObject *)residuals);
	result.covar = PyArray_DATA((PyArrayObject *)covariances);
	status = mpfit(wrapper_user_func, ny, np, pdata, mpparinfo, &mpconf,
	               (void *)&private_data, &result);
	if (status < 0) {
		if (status < -1) {
			PyErr_Format(
				PyExc_RuntimeError, "mpfit function error %d",
				status);
		}
		goto cleanup;
	}

cleanup:
	Py_XDECREF(params);
	free(mpparinfo);

	if (PyErr_Occurred()) {
		Py_XDECREF(fittedpars);
		Py_XDECREF(residuals);
		Py_XDECREF(parerrors);
		Py_XDECREF(covariances);
		return NULL;
	}

	return Py_BuildValue(
		"(N,{s:d,s:d,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:N,s:N,s:N})",
		fittedpars,
		"bestnorm", result.bestnorm,
		"orignorm", result.orignorm,
		"niter", result.niter,
		"nfev", result.nfev,
		"status", result.status,
		"npar", result.npar,
		"nfree", result.nfree,
		"npegged", result.npegged,
		"nfunc", result.nfunc,
		"residuals", residuals,
		"parerrors", parerrors,
		"covariances", covariances);
}


static PyMethodDef
mpfit_functions[] = {
	{ "mpfit", (PyCFunction)mpfitwrap, METH_VARARGS|METH_KEYWORDS,
	  "Wrapper around C mpfit"},
	{NULL, NULL, 0, NULL},
};


#ifdef PY3K
static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,
	"mpfit",
	"A wrapper around C mpfit.",
        -1,
	mpfit_functions,
        NULL,
        NULL,
        NULL,
	NULL,
};
#endif


PyMODINIT_FUNC
#ifdef PY3K
PyInit_mpfit(void)
#else
initmpfit(void)
#endif
{
	PyObject *module = NULL;

#ifdef PY3K
	module = PyModule_Create(&moduledef);
#else
	module = Py_InitModule3("mpfit", mpfit_functions,
	                        "Wrapper around C mpfit");
#endif
	if (!module) {
#ifdef PY3K
		return NULL;
#else
		return ;
#endif
	}
	import_array();

#ifdef PY3K
	return module;
#endif
}
