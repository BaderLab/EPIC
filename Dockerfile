FROM jupyter/datascience-notebook:cfddc5a3163f
MAINTAINER Florian Goebels <florian.goebels@gmail.com>

RUN conda install -n python2 --quiet --yes \
    'rpy2=2.8*' \
    'r-base=3.3.2' \
    'r-irkernel=0.7*' \
    'r-plyr=1.8*' \
    'r-devtools=1.12*' \
    'r-tidyverse=1.0*' \
    'r-shiny=0.14*' \
    'r-rmarkdown=1.2*' \
    'r-forecast=7.3*' \
    'r-rsqlite=1.1*' \
    'r-reshape2=1.4*' \
    'r-nycflights13=0.2*' \
    'r-caret=6.0*' \
    'r-rcurl=1.95*' \
    'r-crayon=1.3*' \
    'r-kohonen=3.0.2*' \
    'r-randomforest=4.6*' && conda clean -tipsy

RUN /opt/conda/envs/python2/bin/Rscript -e "install.packages(\"https://cran.r-project.org/src/contrib/Archive/wccsom/wccsom_1.2.11.tar.gz\", repos=NULL, type=\"source\")"

USER root

RUN apt-get update
RUN apt-get install --fix-missing
RUN apt-get install -y default-jre


RUN pip2 install hide_code
RUN pip2 install fileupload
RUN pip2 install jupyter_contrib_nbextensions

RUN /opt/conda/envs/python2/bin/jupyter nbextension install --py fileupload
RUN /opt/conda/envs/python2/bin/jupyter nbextension enable --py fileupload

RUN /opt/conda/envs/python2/bin/jupyter nbextension install --py hide_code
RUN /opt/conda/envs/python2/bin/jupyter nbextension enable --py hide_code
RUN /opt/conda/envs/python2/bin/jupyter serverextension enable --py hide_code
RUN /opt/conda/envs/python2/bin/jupyter nbextension enable jupyter_contrib_nbextensions


USER $NB_USER
USER $NB_USER
USER $NB_USER
USER $NB_USER
USER $NB_USER
USER $NB_USER
USER $NB_USER
USER $NB_USER
USER $NB_USER
USER $NB_USER
USER $NB_USER
USER $NB_USER

RUN git clone https://github.com/BaderLab/EPIC $HOME/work/EPIC
RUN ln -s $HOME/work/EPIC/src/EPIC.ipynb $HOME/work/EPIC.ipynb



