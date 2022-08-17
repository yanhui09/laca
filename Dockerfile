# Container image that runs your code
FROM continuumio/miniconda3:latest
COPY entrypoint.sh /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]