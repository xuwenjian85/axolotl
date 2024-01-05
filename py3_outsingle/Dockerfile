# 基础镜像
FROM python:3.11

# 设定工作目录
WORKDIR /app

# 复制应用程序代码到容器
COPY requirements.txt /app

# 安装Python软件包
RUN pip install -i https://pypi.tuna.tsinghua.edu.cn/simple -r requirements.txt

# 定义容器启动时执行的命令
# CMD ["python", "app.py"]
