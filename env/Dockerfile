FROM debian:bookworm-20250317
RUN apt-get update && \
    apt-get install -y \
        wget \
        zip \
        python3 \
        python3-pip \
        hmmer
RUN ln -s /usr/bin/python3 /usr/bin/python
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
RUN pip3 install --break-system-packages biopython
