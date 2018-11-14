#TAVP-WRK Instruction Execution Log Analysis script

echo "TAVP-WRK Instruction Execution State:"

echo "DECODER:IN"
grep "DECODER:IN" ./qforce.*.log | wc

echo "RESOURCER:IN"
grep "RESOURCER:IN" ./qforce.*.log | wc

echo "RESOURCER:OK"
grep "RESOURCER:OK" ./qforce.*.log | wc

echo "RESOURCER:OD"
grep "RESOURCER:OD" ./qforce.*.log | wc

echo "RESOURCER:DEP"
grep "RESOURCER:DEP" ./qforce.*.log | wc

echo "RESOURCER:LIM"
grep "RESOURCER:LIM" ./qforce.*.log | wc

echo "COMMUNICATOR:FET"
grep "COMMUNICATOR:FET" ./qforce.*.log | wc

echo "COMMUNICATOR:FSNC"
grep "COMMUNICATOR:FSNC" ./qforce.*.log | wc

echo "DISPATCHER:ISS"
grep "DISPATCHER:ISS" ./qforce.*.log | wc

echo "DISPATCHER:CML"
grep "DISPATCHER:CML" ./qforce.*.log | wc

echo "COMMUNICATOR:UPL"
grep "COMMUNICATOR:UPL" ./qforce.*.log | wc

echo "COMMUNICATOR:USNC"
grep "COMMUNICATOR:USNC" ./qforce.*.log | wc

echo "RESOURCER:OUT"
grep "RESOURCER:OUT" ./qforce.*.log | wc

echo "RETIRER:OUT"
grep "RETIRER:OUT" ./qforce.*.log | wc

echo ""
grep "" ./qforce.*.log | wc

grep -i error ./qforce.*.log
grep -i fatal ./qforce.*.log
