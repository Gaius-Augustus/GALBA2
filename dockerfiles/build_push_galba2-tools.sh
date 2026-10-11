#!/bin/bash
# Build, check and push the galba2-tools image.
# Run on a machine with Docker and `docker login` done (sudo if your user is
# not in the docker group):
#
#   sudo bash dockerfiles/build_push_galba2-tools.sh            # build + check + push
#   PUSH=0 bash dockerfiles/build_push_galba2-tools.sh          # build + check only
#   TAG=v0.5.1 bash dockerfiles/build_push_galba2-tools.sh      # other version tag
#
# Why this rebuild: compleasm 0.2.7 crashed in Downloader.download_placement
# on the file_versions.tsv the BUSCO server serves since the odb12.2 release
# (compleasm issue #61); 0.2.10 parses it. The image also gets the miniprot
# PATH order from the Dockerfile (/usr/local/bin/miniprot 0.14 before
# compleasm's bundled 0.12).
#
# After the push, Snakemake keeps using the cached image of the same URI
# until the cache file is removed (its name is the md5 of the URI):
#
#   rm /home/katharina/git/GALBA2/.singularity_cache/840f53cf2aa7e7561f3f106d30090e52.simg
#   ssh brain rm /home/hoffk83/git/GALBA2/.singularity_cache/840f53cf2aa7e7561f3f106d30090e52.simg
#
# or point galba_tools_image at the version tag instead of :latest.

set -euo pipefail

IMAGE=katharinahoff/galba2-tools
TAG=${TAG:-v0.5.0}
PUSH=${PUSH:-1}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "[$(date)] Building $IMAGE:$TAG and $IMAGE:latest"
docker build -f "$SCRIPT_DIR/Dockerfile.galba2-tools" -t "$IMAGE:$TAG" -t "$IMAGE:latest" "$REPO_DIR"

echo "[$(date)] Checking the image"
docker run --rm "$IMAGE:$TAG" bash -c '
    set -e
    v=$(sed -n "s/^__version__ = \"\(.*\)\"/\1/p" /opt/compleasm_kit/_version.py)
    echo "compleasm $v"
    [ "$v" = "0.2.10" ] || { echo "ERROR: expected compleasm 0.2.10"; exit 1; }
    echo "miniprot $(miniprot --version) at $(command -v miniprot)"
    [ "$(command -v miniprot)" = /usr/local/bin/miniprot ] || { echo "ERROR: compleasm miniprot shadows /usr/local/bin/miniprot"; exit 1; }
    minisplice 2>&1 | head -1 || true
    miniprothint.py --help > /dev/null && echo "miniprothint ok"
    python3 -c "import Bio, matplotlib, pandas; print(\"biopython\", Bio.__version__, \"matplotlib\", matplotlib.__version__, \"pandas\", pandas.__version__)"
'

if [ "$PUSH" = "1" ]; then
    echo "[$(date)] Pushing $IMAGE:$TAG and $IMAGE:latest"
    docker push "$IMAGE:$TAG"
    docker push "$IMAGE:latest"
    echo "[$(date)] Done. Remove the cached .simg files (see header) or switch configs to $IMAGE:$TAG."
else
    echo "[$(date)] Built without pushing (PUSH=0)."
fi
