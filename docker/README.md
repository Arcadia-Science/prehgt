# Example instructions for how to build a docker container from a Dockerfile

The image needed to be run with an explicit linux platform if built on M1 Mac.

```
docker build --progress=plain --no-cache --platform linux/x86_64 -t tidy-prehgt tidy-prehgt
docker login
docker images
docker tag tidy-prehgt arcadiascience/tidy-prehgt:2.0.0
docker push arcadiascience/tidy-prehgt:2.0.0
```
