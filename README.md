# 2-Dimensional-FFT without 1D FFT

This work etends the concept of FFT to 2 dimensions without applying FFT on rows and columns successively. This repository containes the c++ implementation of the 2D FFT.

**Usage**
```
g++ fft.cpp -L/usr/local/Cellar/opencv/3.3.0_3/lib -I/usr/local/Cellar/opencv/3.3.0_3/include/opencv -I/usr/local/Cellar/opencv/3.3.0_3/include -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs
```

Opencv include path can be found using

```
pkg-config opencv --cflags
```

and library path can be found using

```
pkg-config opencv --libs
```


