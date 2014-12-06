#include <Python.h>

static PyObject* ADMMModule(PyObject* self, PyObject* args);

static PyMethodef module_methods[] = {
  { "ADMMInterface", (PyCFunction)ADMMModule, METH_VARARGS, "Function to interface with the C++ module" }, {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initInterface() {
  Py_InitModule3( ADMMInterface, module_methods, "Interface for C++");
}
