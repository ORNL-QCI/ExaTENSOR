#TAVP-MNG Instruction Execution Analysis script

echo "TAVP-MNG Instruction Execution State:"

echo "Decoded (all):"
grep "(TAVP-MNG:uDecoder): Decoded a new tensor instruction:" ./qforce.*.log | wc

echo "Received by Locator (all):"
grep "An instruction is received for operand location:" ./qforce.*.log | wc

echo "Deferred by Locator (tensor only):"
grep "An instruction is deferred (operands not ready):" ./qforce.*.log | wc

echo "Decomposed by Decomposer (tensor only):"
grep "(TAVP-MNG:Decomposer): A new instruction was decomposed:" ./qforce.*.log | wc

echo "Received back by bottom-TAVP Decomposer (tensor only):"
grep "A previously decomposed instruction is received back (bottom):" ./qforce.*.log | wc

echo "Dispatched (all):"
grep "A new instruction is dispatched:" ./qforce.*.log | wc

echo "Received by Collector from Decomposer (all):"
grep "A parent instruction was received:" ./qforce.*.log | wc

echo "Matched by Collector (tensor only):"
grep "Matched a child instruction:" ./qforce.*.log | wc

echo "Retired (all):"
grep "(TAVP-MNG:Retirer): Retired instruction:" ./qforce.*.log | wc

grep -i error ./qforce.*.log
grep -i fatal ./qforce.*.log
