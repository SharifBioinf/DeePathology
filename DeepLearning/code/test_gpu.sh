#!/usr/bin/env bash
nvidia-smi | awk '/GeForce/ {getline;print $13}' | sed 's/%//g'
