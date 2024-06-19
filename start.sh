#!/bin/bash
source .venv/bin/activate \
&& cd guidance_Linux \
&& python3 script/guidance_main.py $@