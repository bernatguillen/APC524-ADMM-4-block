#include <Python.h>

/*Copying example from manual for now*/
static PyObject *ADMMModule(PyObject *self, PyObject *args){
  const char *command;
  int sts;
  if (!PyArg_ParseTuple(args, "s", &command)){
    return NULL;
  }
  
}
