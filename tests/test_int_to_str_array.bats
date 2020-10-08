#!/usr/bin/env bats

setup() {
	ymlPath="$BATS_TEST_DIRNAME/inputs/int_to_int.yml"
	touch $ymlPath
}

teardown() {
	rm $ymlPath
}


@test "int to string array" {
	echo 'number: 3' > $ymlPath
	cwltool ./SJARACNe/cwl/int_to_str_array.cwl $ymlPath | tr -d "\n" | tr -s " " | tee results
	run grep '"TF_run_000.adj", "TF_run_001.adj", "TF_run_002.adj"' results
	[ "$status" -eq 0 ]
	rm results
}
