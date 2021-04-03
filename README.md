# Q4 Program

## 作者：Charlie Xia，Max Zhang，9527567

### 主要功能：

使用C++实现等参单元四节点的有限元程序的计算。
#### 2021/04/03更新 取消Windows10支持，优化性能
#### 2021/2/15更新 可以在Linux，Windows10，Max OS运行

##### 使用教程：Ubuntu：
 这里是高性能版本，使用英特尔编译器以及mkl并行加速eigen库，英特尔的编译器安装不再赘述。

安装eigen3库：

```bash
sudo apt install libeigen3-dev
```

克隆本项目到本地:

```bash
git clone https://github.com/charlie-C929/Q4_Program_Charlie.git
```

创建build文件夹，

```bash
mkdir build
cd build
CXX=icpc cmake .. 
make
```

进入 bin 文件夹：

```bash
cd .. # 这行语句是为了从刚才build文件夹返回上一级目录
cd bin
./demo 你的输入文件路径 想要的输出文件路径
#输入文件路径是必须的，输出文件并非必须，如果不写输出文件路径，结果将会直接打印到屏幕
#example
./demo mesh.k myresult.dat
./demo mesh.k
```

如果eigen库提示在使用空矩阵，请检查你的输入文件格式。

windows下建议使用Linux子系统。

~~##### 使用教程：Windows10：~~

~~使用Cmake GUI进行即可，也可用命令行，注意要在 build 文件夹下，文件结构不会乱。如果提示找不到eigen库，则查看src中的CMakeLists.txt中的提示。在Visual studio中运行Debug时要注意newRead.cpp中的文件路径，按自己的实际情况修改。~~


##### 使用教程：Mac OS：

参考Ubuntu
