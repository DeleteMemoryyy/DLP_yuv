# YUV Converter

---

- 运行环境：

  - 64位Linux发行版，本机环境为Ubuntu 16.04 64-bit
  - Intel x86架构CPU（需要浮点指令集支持）

- 可执行文件： convert

- 可选参数
  - -ti：测试模式和输入文件
    - -t1 f1：在图片f1上执行测试1（单张图片alpha混合），f1为图片文件名，输入图片要求为YUV420P格式
    - -t2 f1 f2：在图片f1和f2上执行测试2（两张图片进行叠加），f1、f2为图片文件名，输入图片要求为YUV420P格式
  - -w width：手动输入图片宽度，默认为1920
  - -h height：手动输入图片高度，默认为1080:
  - -o IS：开启浮点指令集优化
    - -o mmx 或 -o MMX
    - -o sse2 或 -o SSE2
    - -o avx 或 -o AVX

- 输入实例：

  - ./convert -t1 demo/dem2.yuv -o SSE2
  - ./convert -t2 demo/dem1.yuv demo/dem2.yuv -o AVX