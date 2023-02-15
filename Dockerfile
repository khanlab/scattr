# Stage: requirements
# Notes: g++ (snakebids), tcsh (freesurfer), parellel (scattr)
FROM python:3.9-slim-bullseye AS requirements 
RUN mkdir -p /opt \
    && apt-get update -qq \ 
    && apt-get install -y -q --no-install-recommends \
       libeigen3-dev \
       zlib1g-dev \
       libqt5opengl5-dev \
       libqt5svg5-dev \
       libgl1-mesa-dev \
       libncurses5 \
       libxt6 \
       parallel \
       rsync \
       tcsh \
       curl \
       git \
       g++ \
       unzip \
       wget \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Stage: mrtrix
FROM requirements AS mrtrix
ARG MRTRIX_VER=3.0.4
RUN git clone https://github.com/MRtrix3/mrtrix3.git /opt/mrtrix3 \
    && cd /opt/mrtrix3 \
    && git fetch --tags \
    && git checkout ${MRTRIX_VER} \
    && ./configure \ 
    && ./build

# Stage: freesurfer
FROM requirements AS freesurfer
ARG FS_VER=7.2.0
RUN wget https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/${FS_VER}/freesurfer-linux-centos7_x86_64-${FS_VER}.tar.gz -O fs.tar.gz \
    && tar --no-same-owner -xzvf fs.tar.gz \
    && mv freesurfer /usr/local \
    && rm fs.tar.gz \
    && curl 'https://surfer.nmr.mgh.harvard.edu/fswiki/MatlabRuntime?action=AttachFile&do=get&target=runtime2014bLinux.tar.gz' -o 'fs_runtime2014b.tar.gz' \
    && tar --no-same-owner -xzvf fs_runtime2014b.tar.gz -C /usr/local/freesurfer \
    && rm fs_runtime2014b.tar.gz \
    && chmod -R 775 /usr/local/freesurfer/MCRv84 

# Stage: ANTs
# Copy over antsApplyTransforms and antsRegistrationSyNQuick.sh to "minify"
FROM requirements as ants
ARG ANTS_VER=2.4.3
RUN wget https://github.com/ANTsX/ANTs/releases/download/v${ANTS_VER}/ants-${ANTS_VER}-centos7-X64-gcc.zip -O ants.zip \
    && unzip -qq ants.zip -d /opt \
    && mv /opt/ants-${ANTS_VER} /opt/ants \
    && rm ants.zip

# Stage: build
FROM requirements AS build
COPY . /opt/scattr/
RUN cd /opt/scattr \
    && pip install --prefer-binary --no-cache-dir \
        poetry==1.2.2 \
    && poetry build -f wheel

# Stage: runtime
# NOTE: use fsl stage for library dependencies
FROM requirements AS runtime
COPY --from=mrtrix /opt/mrtrix3 /opt/mrtrix3
COPY --from=freesurfer /usr/local/freesurfer /usr/local/freesurfer
COPY --from=ants /opt/ants/bin/antsApplyTransforms /opt/ants/bin/antsRegistration /opt/ants/bin/antsRegistrationSyNQuick.sh /opt/ants/bin/
COPY --from=build /opt/scattr/dist/*.whl /opt/scattr/
RUN WHEEL=`ls /opt/scattr | grep whl` \
    && pip install /opt/scattr/$WHEEL \
    && rm -r /opt/scattr \
    && apt-get purge -y -q curl g++ unzip wget \ 
    && apt-get --purge -y -qq autoremove
# Setup environments
ENV OS=Linux \
    # Path: mrtrix > freesurfer > ants
    PATH=/opt/mrtrix3/bin:/usr/local/freesurfer/bin:/usr/local/freesurfer/fsfast/bin:/usr/local/freesurfer/tktools:/usr/local/freesurfer/mni/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/ants/:/opt/ants/bin:$PATH \
    # freesurfer specific
    FREESURFER_HOME=/usr/local/freesurfer \
    FREESURFER=/usr/local/freesurfer \
    SUBJECTS_DIR=/usr/local/freesurfer/subjects \
    LOCAL_DIR=/usr/local/freesurfer/local \
    FSFAST_HOME=/usr/local/freesurfer/fsfast \
    FMRI_ANALYSIS_DIR=/usr/local/freesurfer/fsfast \ 
    FUNCTIONALS_DIR=/usr/local/freesurfer/sessions \
    FS_OVERRIDE=0 \
    FIX_VERTEX_AREA="" \
    FSF_OUTPUT_FORMAT=nii.gz \
    MINC_BIN_DIR=/usr/local/freesurfer/mni/bin \
    MINC_LIB_DIR=/usr/local/freesurfer/mni/lib \
    MNI_DIR=/usr/local/freesurfer/mni \
    MNI_DATAPATH=/usr/local/freesurfer/mni/data \
    MNI_PERL5LIB=/usr/local/freesurfer/mni/share/perl5 \
    PERL5LIB=/usr/local/freesurfer/mni/share/perl5 \
    # ants specific
    ANTSPATH=/opt/ants/bin/ \
ENTRYPOINT ["scattr"]
