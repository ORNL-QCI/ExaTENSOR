#TAVP Instruction Execution Analysis script

echo "Deferred due to data dependency:"
grep "Tensor instruction deferred due to data dependency" ./qforce.*.log | wc

echo "Deferred due to lack of resources:"
grep "Tensor instruction deferred due to lack of resources" ./qforce.*.log | wc

echo "Resources acquired:"
grep "Resources acquired for tensor instruction" ./qforce.*.log | wc

echo "Issued to Communicator directly:"
grep "Tensor instruction issued from main queue" ./qforce.*.log | wc

echo "Issued to Communicator after deferrence:"
grep "Tensor instruction issued from deferred queue" ./qforce.*.log | wc

echo "Input prefetch initiated:"
grep "Initiated input prefetch for tensor instruction" ./qforce.*.log | wc

echo "Input prefetch synced:"
grep "Synced input prefetch for tensor instruction" ./qforce.*.log | wc

echo "Issued for execution by Dispatcher:"
grep "Issued tensor instruction" ./qforce.*.log | wc

echo "Completed execution:"
grep "Completed tensor instruction" ./qforce.*.log | wc

echo "Output upload initiated:"
grep "Initiated output upload for tensor instruction" ./qforce.*.log | wc

echo "Output upload synced:"
grep "Synced output upload for tensor instruction" ./qforce.*.log | wc

echo "Resources released:"
grep "Resources released for tensor instruction" ./qforce.*.log | wc

echo "Retired:"
grep "Retired tensor instruction" ./qforce.*.log | wc
