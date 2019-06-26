#!/usr/bin/env bats

@test "acceptance testing" {
	tempdir="$BATS_TEST_DIRNAME/results"
	sjaracne local -e ./tests/inputs/Tcell1170.exp -g ./tests/inputs/TcellTF.txt -n 5 -o $tempdir
	run diff $tempdir/consensus_network_ncol_.txt ./tests/answerkey/acceptance/cnn_5.txt
	[ "$status" -eq 0 ]
	rm -r $tempdir

}
