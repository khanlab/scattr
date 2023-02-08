# Stage: mrtrix
FROM python:3.9-slim-bullseye AS mrtrix
ARG MRTRIX_VER=3.0.4
RUN mkdir -p /opt \
    && apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
       git \
       g++ \
       libeigen3-dev \
       zlib1g-dev \
       libqt5opengl5-dev \
       libqt5svg5-dev \
       libgl1-mesa-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && git clone https://github.com/MRtrix3/mrtrix3.git /opt/mrtrix3 \
    && cd /opt/mrtrix3 \
    && git fetch --tags \
    && git checkout ${MRTRIX_VER} \
    && ./configure \ 
    && ./build \
    && apt-get purge -y -q g++ \
    && apt-get --purge -y -qq autoremove 
ENV PATH=/opt/mrtrix3/bin:$PATH

# Stage: freesurfer
FROM mrtrix AS freesurfer
ARG FS_VER=7.2.0
RUN mkdir -p /opt \
    && apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
       tcsh \
       wget \
       curl \
       libncurses5 \
       libxt6 \
    && wget  https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/${FS_VER}/freesurfer-linux-centos7_x86_64-${FS_VER}.tar.gz -O fs.tar.gz \
    && tar --no-same-owner -xzvf fs.tar.gz \
    && mv freesurfer /usr/local \
    && rm fs.tar.gz \
    && curl 'https://surfer.nmr.mgh.harvard.edu/fswiki/MatlabRuntime?action=AttachFile&do=get&target=runtime2014bLinux.tar.gz' -o 'fs_runtime2014b.tar.gz' \
    && tar --no-same-owner -xzvf fs_runtime2014b.tar.gz -C /usr/local/freesurfer \
    && rm fs_runtime2014b.tar.gz \
    && chmod -R 775 /usr/local/freesurfer/MCRv84 \
    && apt-get purge -y -q wget curl
# setup fs env
ENV OS=Linux \ 
    PATH=/usr/local/freesurfer/bin:/usr/local/freesurfer/fsfast/bin:/usr/local/freesurfer/tktools:/usr/local/freesurfer/mni/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH \
    FREESURFER_HOME=/usr/local/freesurfer \
    FREESURFER=/usr/local/freesurfer \
    SUBJECTS_DIR=/usr/local/freesurfer/subjects \
    LOCAL_DIR=/usr/local/freesurfer/local \
    FSFAST_HOME=/usr/local/freesurfer/fsfast \
    FMRI_ANALYSIS_DIR=/usr/local/freesurfer/fsfast \ 
    FUNCTIONALS_DIR=/usr/local/freesurfer/sessions \
    # set default fs options
    FS_OVERRIDE=0 \
    FIX_VERTEX_AREA="" \
    FSF_OUTPUT_FORMAT=nii.gz \
    # mni env requirements
    MINC_BIN_DIR=/usr/local/freesurfer/mni/bin \
    MINC_LIB_DIR=/usr/local/freesurfer/mni/lib \
    MNI_DIR=/usr/local/freesurfer/mni \
    MNI_DATAPATH=/usr/local/freesurfer/mni/data \
    MNI_PERL5LIB=/usr/local/freesurfer/mni/share/perl5 \
    PERL5LIB=/usr/local/freesurfer/mni/share/perl5

# Stage: ANTs
FROM freesurfer as ants
ARG ANTS_VER=2.4.3
RUN mkdir -p /opt \
    && apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
       wget \
       unzip \
    && wget https://github.com/ANTsX/ANTs/releases/download/v${ANTS_VER}/ants-${ANTS_VER}-centos7-X64-gcc.zip -O ants.zip \
    && unzip -qq ants.zip -d /opt \
    && rm ants.zip \
    && apt-get purge -y -q wget unzip 
# setup ants env
ENV ANTSPATH=/opt/ants-${ANTS_VER}/bin/ \
    PATH=/opt/ants-${ANTS_VER}/:/opt/ants-${ANTS_VER}/bin:$PATH

# Stage: build
# NOTE: g++ and libdatrie are required for poetry install
FROM ants AS build
COPY ./poetry.lock ./pyproject.toml /
RUN mkdir -p /opt \
    && apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
        g++=4:10.2.1-1 \
        libdatrie1=0.2.13-1 \
        parallel=20161222-1.1 \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && pip install --prefer-binary --no-cache-dir \
        poetry==1.2.2 \
    && poetry config virtualenvs.create false \
    && poetry install --only main \
    && apt-get purge -y -q g++ \
    && apt-get --purge -y -qq autoremove

# Stage: runtime
FROM build AS runtime
COPY ./scattr /opt/scattr
# ENTRYPOINT ["/opt/scattr/run.py"]
