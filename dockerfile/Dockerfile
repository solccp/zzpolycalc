FROM intel/oneapi-hpckit:latest AS builder
WORKDIR /app
COPY . .
RUN ./build.sh

FROM ubuntu:22.04
COPY --from=builder /app/ZZPolyCalc /app/ZZPolyCalc
ENTRYPOINT ["/app/ZZPolyCalc"]
