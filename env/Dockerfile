FROM ubuntu
RUN apt-get update
RUN apt-get upgrade -y 
RUN apt-get install -y wget
RUN apt-get install -y zip
RUN apt-get install -y python3
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN apt-get install -y python3-pip
RUN pip3 install --break-system-packages h5py
RUN pip3 install --break-system-packages numpy
RUN pip3 install --break-system-packages pandas
RUN pip3 install --break-system-packages openpyxl
RUN pip3 install --break-system-packages matplotlib
RUN pip3 install --break-system-packages seaborn
RUN pip3 install --break-system-packages logomaker
RUN pip3 install --break-system-packages requests
RUN pip3 install --break-system-packages scipy
RUN pip3 install --break-system-packages scikit-learn
RUN pip3 install --break-system-packages tables
