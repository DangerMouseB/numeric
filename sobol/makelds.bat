@echo off

del release\lds.*
del lds.dll

cl /nologo /c lds.cpp /Fo"release\\" /EHsc

link /nologo /OUT:release\lds.dll /def:ldsDLL.def /libpath:release lds.obj /DLL

move release\lds.dll

