#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>

namespace khnum {
void Test() {
    Py_Initialize();
    PyRun_SimpleString("from time import time,ctime\n"
                       "print('Today is', ctime(time()))\n");

}
} //namespace khnum