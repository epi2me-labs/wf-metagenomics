ARG BASEIMAGE=git.oxfordnanolabs.local:4567/epi2melabs/workflow-containers/base-workflow-image:v0.0.1
FROM $BASEIMAGE

# Minimal install for example purposes
COPY environment.yaml $HOME/environment.yaml 
RUN \
    . $CONDA_DIR/etc/profile.d/mamba.sh \
    && micromamba activate \
    && micromamba install --help \
    && micromamba install --file $HOME/environment.yaml \
    && fix-permissions $CONDA_DIR \
    && fix-permissions $HOME

COPY environment_centrifuge.yaml $HOME/environment_centrifuge.yaml
ENV CONDA_CENTRIFUGE_DIR=/home/$WF_USER/conda_centrifuge

RUN \
    ./bin/micromamba shell init -s bash -p $CONDA_CENTRIFUGE_DIR \
    && source ~/.bashrc \
    && chown $WF_USER:$WF_GID $CONDA_CENTRIFUGE_DIR \
    && fix-permissions $CONDA_CENTRIFUGE_DIR \
    && fix-permissions $HOME

RUN \
    . $CONDA_CENTRIFUGE_DIR/etc/profile.d/mamba.sh \
    && micromamba activate $CONDA_CENTRIFUGE_DIR \
    && micromamba install --help \
    && micromamba install --file $HOME/environment_centrifuge.yaml \
    && fix-permissions $CONDA_CENTRIFUGE_DIR \
    && fix-permissions $HOME

USER $WF_UID
WORKDIR $HOME
