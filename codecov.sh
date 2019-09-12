#!/bin/bash

apt-get install -y git
R -e "install.packages('covr'); covr::codecov(token = '$CODECOV_TOKEN')"
