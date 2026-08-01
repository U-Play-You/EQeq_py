# EQeq Python Binding (pybind11) — README

## 概述
这是对原 **EQeq** 电荷平衡算法的改造版本，参考文献：[An Extended Charge Equilibration Method](https://doi.org/10.1021/jz3008485)  
使用 **pybind11** 封装为 Python 扩展模块 `eqeq`。  
模块提供一个便捷函数 `run()`，完成计算，并返回 **{label: charge}** 的 Python 字典。  
**chargecenters.dat** 与 **ionizationdata.dat** 已预编码进程序中。  
## 使用方法

本项目构建的扩展模块指定用于 **CPython 3.13**。请先将本仓库 Fork 到自己的 GitHub 账户，然后通过 GitHub Actions 编译：

1. 打开 Fork 后的仓库，进入 **Actions** 页面；如果 GitHub 提示工作流尚未启用，请先点击 **I understand my workflows, go ahead and enable them**。
2. 在左侧选择 **Build for CentOS 7**，点击 **Run workflow**，选择要构建的分支后再次点击 **Run workflow**。
3. 等待任务完成，打开该次运行记录，在 **Artifacts** 中下载 `eqeq-python313-centos7-x86_64`。
4. 解压 Artifact，将其中的 `eqeq*.so` 放到 Python 脚本所在目录。

该 Artifact 适用于 **CPython 3.13、Linux x86_64**，并基于 `manylinux2014` 构建。它不能直接用于其他 Python 版本、Windows 或 macOS。

然后在 Python 3.13 中导入并调用：

```python
import eqeq

charge = eqeq.run("xxx.cif")
print(charge)
```

## 可选参数

在调用 `run` 时可以自行添加参数，除 `cif` 路径之外的其他参数已在程序中预设，使用如下方法自定义，具体可阅读源文件。  

```python
eqeq.run("mystructure.cif", precision=3, method="Ewald", **{"lambda": 1.2})
```

## Overview

This is a modified version of the original EQeq charge equilibration algorithm. Reference: [An Extended Charge Equilibration Method](https://doi.org/10.1021/jz3008485).  
The code is wrapped with **pybind11** as a Python extension module named `eqeq`.  
The module provides a convenient function `run()` that performs the calculation and returns a Python dictionary {label: charge}.  
The files **chargecenters.dat** and **ionizationdata.dat** have been pre-encoded into the program.

## Usage

The extension produced by this project targets **CPython 3.13**. Fork this repository to your GitHub account, then build it with GitHub Actions:

1. Open the forked repository and go to **Actions**. If GitHub says workflows are disabled, click **I understand my workflows, go ahead and enable them** first.
2. Select **Build for CentOS 7**, click **Run workflow**, choose the branch to build, and click **Run workflow** again.
3. When the run finishes, open it and download `eqeq-python313-centos7-x86_64` from **Artifacts**.
4. Extract the Artifact and place the included `eqeq*.so` file in the same directory as your Python script.

The Artifact targets **CPython 3.13 on Linux x86_64** and is built against `manylinux2014`. It cannot be used directly with another Python version, Windows, or macOS.

Import and call it with Python 3.13:

```python
import eqeq

charge = eqeq.run("xxx.cif")
print(charge)
```

## Optional Parameters

You can pass optional parameters to run. Besides the cif path, other parameters have sensible defaults in the program; you can override them as needed. For example:

```python
eqeq.run("mystructure.cif", precision=3, method="Ewald", **{"lambda": 1.2})
```

See the source code for the full list of configurable parameters and their meanings.
