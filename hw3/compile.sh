#!/bin/bash
`root-config --cxx --cflags` -o main main.c `root-config --libs`
# g++ main.c -o main -g -Wall -Wextra -Werror -I /home/asatk/Downloads/setups/root_v6.26.06.Linux-ubuntu22-x86_64-gcc11.2/root/include