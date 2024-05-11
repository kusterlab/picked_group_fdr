# exit on first error
set -e

./tests/system_tests/test_digest.sh
./tests/system_tests/test_picked_group_fdr.sh
./tests/system_tests/test_entrapment_fdr.sh
./tests/system_tests/test_pipeline.sh
./tests/system_tests/test_quantification.sh
./tests/system_tests/test_fragpipe.sh
./tests/system_tests/test_sage.sh

echo ""
echo "All system tests passed"