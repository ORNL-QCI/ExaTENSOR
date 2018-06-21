#TAVP Instruction Execution Analysis script

echo "Resources acquired:"
grep "Resources acquired for tensor instruction" ./qforce.*.log | wc

echo "Issued to Communicator:"
grep "Tensor instruction issued from main queue" ./qforce.*.log | wc

echo "Issued to Communicator after deferrence:"
grep "Tensor instruction issued from deferred queue" ./qforce.*.log | wc

echo "Unissued from deferred queue due to persisting dependency:"
grep "Tensor instruction kept in deferred queue due to dependency:" ./qforce.*.log | wc

echo "Input prefetch initiated:"
grep "Initiated input prefetch for tensor instruction" ./qforce.*.log | wc

echo "Input prefetch synced:"
grep "Synced input prefetch for tensor instruction" ./qforce.*.log | wc

echo "Issued for execution:"
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
