@echo off

del release\distributions.*
del distributions.dll

cl /nologo /c distributions.cpp /Fo"release\\" /EHsc

link /nologo /OUT:release\distributions.dll /def:distributionsDLL.def /libpath:release distributions.obj /DLL

move release\distributions.dll

