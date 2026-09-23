#!/bin/bash
echo ""

./launch_create_PLANETi.sh &
BACK_PID=$!
wait $BACK_PID

python create_Timings.py
