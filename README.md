# EQeq Python Binding (pybind11) — README

## 概述
这是对原 **EQeq** 电荷平衡算法的改造版本，参考文献：[An Extended Charge Equilibration Method](https://doi.org/10.1021/jz3008485)  
使用 **pybind11** 封装为 Python 扩展模块 `eqeq`。  
模块提供一个便捷函数 `run()`，完成计算，并返回 **{label: charge}** 的 Python 字典。  
**chargecenters.dat** 与 **ionizationdata.dat** 已预编码进程序中。  
## 使用方法
下载 **Release** 里的 `.so` 文件，放到 `.py` 文件同目录下，然后在 Python 中引用即可。
```
import eqeq
charge = eqeq.run(xxx.cif)
print(charge)
```
或者自行编译
在根目录中执行
```
mkdir build
cd build
cmake ..
make -j4
```
## 可选参数
在调用 `run` 时可以自行添加参数，除 `cif` 路径之外的其他参数已在程序中预设，使用如下方法自定义，具体可阅读源文件。  
```
eqeq.run("mystructure.cif", precision=3, method="Ewald", lambda=1.2)
```
## Overview
This is a modified version of the original EQeq charge equilibration algorithm. Reference: [An Extended Charge Equilibration Method](https://doi.org/10.1021/jz3008485).  
The code is wrapped with **pybind11** as a Python extension module named `eqeq`.
The module provides a convenient function `run()` that performs the calculation and returns a Python dictionary {label: charge}.  
The files **chargecenters.dat** and **ionizationdata.dat** have been pre-encoded into the program.

## Usage
Download the `.so` file from the Release, place it in the same directory as your `.py` file, and then import it in Python:
```
import eqeq
charge = eqeq.run("xxx.cif")
print(charge)
```
Or build it yourself. From the project root directory run:
```
mkdir build
cd build
cmake ..
make -j4
```
## Optional Parameters
You can pass optional parameters to run. Besides the cif path, other parameters have sensible defaults in the program; you can override them as needed. For example:
```
eqeq.run("mystructure.cif", precision=3, method="Ewald", lambda=1.2)
```
See the source code for the full list of configurable parameters and their meanings.
