#!/bin/bash
echo "Into create-"

echo "Launch create"
./launch_create_PLANETi.sh &
BACK_PID=$!
wait $BACK_PID

echo " go to create timing"
python create_Timings.py
