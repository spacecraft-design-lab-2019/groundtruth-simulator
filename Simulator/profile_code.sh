#!/bin/bash
# profile code

python3 -m cProfile -o detumble_profile.prof -s cumtime compare_detumble_algorithms.py
snakeviz detumble_profile.prof