#!/usr/bin/env bats

setup() {
	ymlPath="$BATS_TEST_DIRNAME/inputs/int_to_int.yml"
	touch $ymlPath
}

teardown() {
	rm $ymlPath
}

@test "int=0 to int[]" {
	echo 'number: 0' >> $ymlPath

	cwltool ./SJARACNe/cwl/int_to_int_array.cwl $ymlPath | tr -d '\n' | tr -s " " | tee results0
	run grep '[ ]' results0
	[ "$status" -eq 0 ]
	rm results0
}

@test "int=3 to int[0, 1, 2]" {
	echo 'number: 3' >> $ymlPath

	cwltool ./SJARACNe/cwl/int_to_int_array.cwl $ymlPath | tr -d '\n' | tr -s " " | tee results3
	run grep '0, 1, 2' results3
	[ "$status" -eq 0 ]
	rm results3
}
