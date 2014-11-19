echo 'Warming up...'
python mess.py -m testPWM.txt testPWM2.txt -s testSeqs.txt testSeqs2.txt testSeqs3.txt -p 0.05 -t 8 > /dev/null

echo 'Testing 8 threads...'
time python mess.py -m testPWM.txt testPWM2.txt -s testSeqs.txt testSeqs2.txt testSeqs3.txt -p 0.05 -t 8 > /dev/null
#echo 'Threaded (8 threads)'

echo 'Testing 4 threads...'
time python mess.py -m testPWM.txt testPWM2.txt -s testSeqs.txt testSeqs2.txt testSeqs3.txt -p 0.05 -t 4 > /dev/null

echo 'Testing unthreaded...'
time python mess.py -m testPWM.txt testPWM2.txt -s testSeqs.txt testSeqs2.txt testSeqs3.txt -p 0.05 > /dev/null
#echo 'Unthreaded'

